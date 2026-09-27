#' Separate individual colonies, the bacterium and its halo
#'
#' Turns the pixel classes into labelled objects: one region per expected
#' fungal colony, the central bacterial colony and its halo.
#'
#' @details
#' **Fungal mask clean-up.** The fungus class is closed and opened with a
#' disk of radius `morph_mm` (bridging wrinkle shadows, removing hairline
#' specks), holes are filled (unclassified craters and highlights inside a
#' colony belong to the colony) and fragments smaller than
#' `min_area_mm2` are dropped.
#'
#' **Seeded separation of cultures.** For each expected inoculation point
#' (from the [plate_layout()]), the colony *core* is the fungal pixel with
#' the largest distance to the colony edge (maximum of the Euclidean
#' distance transform) within `search_frac` x dish radius of that point.
#' Cores closer than `min_core_mm` to the colony edge are rejected (the
#' inoculation failed, or only a satellite is present) and the colony is
#' reported as not detected. The cores then seed a marker-controlled
#' watershed on the negated distance transform, restricted to the fungal
#' mask: each fungal pixel is assigned to the colony from which it can be
#' reached through fungal tissue with the least "constriction". Colonies
#' that touch are split along their neck, and a colony can never claim
#' fungal tissue that is not connected to it.
#'
#' **Edge placement.** Colour features are lightly smoothed for robust
#' classification, which blurs colony outlines outwards by about one pixel.
#' Outlines are therefore re-placed at the half-contrast position (the
#' classical 50 % edge criterion for a blurred step edge) using the
#' unsmoothed colour difference from agar.
#'
#' **Satellites.** Fungal fragments not connected to any core (e.g. spread
#' inoculum, a secondary micro-colony next to the main one) are labelled
#' separately and excluded from colony statistics, but their number and area
#' are reported as a QC measure.
#'
#' **Bacterium and halo.** The bacterial colony is the bacteria-class
#' component closest to the dish centre (holes filled). The halo is the set
#' of halo-class pixels connected to (within `halo_link_mm` of) the
#' bacterium, excluding bacterial and fungal pixels.
#'
#' @param plate A plate processed by [classify_pixels()].
#' @param layout The [plate_layout()] used for classification.
#' @param min_area_mm2 Minimum area of a fungal fragment (mm^2).
#' @param min_core_mm Minimum distance from core to colony edge (mm).
#' @param morph_mm Radius of morphological clean-up (mm).
#' @param halo_link_mm Maximum gap between bacterium and halo (mm).
#' @param bridge_mm Opening radius (mm) that separates the bacterial colony
#'   from bacteria-like tissue connected to it by thin bridges.
#' @param edge_criterion Sub-pixel-consistent edge placement: within a
#'   3-pixel band inside each colony outline, pixels are kept only if their
#'   colour difference from agar is at least this fraction of the colony's
#'   interior median ("half-maximum" criterion, unbiased for a blurred step
#'   edge). Set 0 to disable.
#'
#' @return The plate with element `objects`: a list containing the colony
#'   label matrix (`labels`, values = row index of `colony_table`),
#'   `satellites`, `bacteria` and `halo` masks and `colony_table`.
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 6)
#' p <- read_plate(sim$image) |> detect_plate() |> model_background() |>
#'   classify_pixels(plate_layout()) |> segment_colonies(plate_layout())
#' p$objects$colony_table
#' @export
segment_colonies <- function(plate, layout = plate_layout(), min_area_mm2 = 0.5,
                             min_core_mm = 0.4, morph_mm = 0.08, halo_link_mm = 1,
                             edge_criterion = 0.5, bridge_mm = 0.5) {
  stopifnot(is_mycohalo_plate(plate), !is.null(plate$classes))
  d <- plate$dish
  mmpp <- d$mm_per_px
  map <- plate$classes$map
  m <- dish_masks(plate)
  rpx <- morph_mm / mmpp

  # ---- bacterium first: one compact colony nearest the dish centre.
  # Bacteria-class pixels elsewhere (e.g. light colony margins, faint halo
  # fringe) are re-assigned to the most probable remaining class.
  bact <- matrix(FALSE, nrow(map), ncol(map))
  if (layout$centre == "bacteria" && any(map == 2L)) {
    b <- map == 2L
    if (rpx >= 1) b <- close_disk(b, rpx)
    # cut thin bridges (e.g. along a halo into a light colony) so that only
    # the compact central colony is kept
    b <- open_disk(b, max(1, bridge_mm / mmpp))
    b <- remove_small(b, min_area_mm2 / mmpp^2)
    if (any(b)) {
      lb <- label_components(b)
      ids <- seq_len(max(lb))
      dc <- vapply(ids, function(k) {
        w <- which(lb == k)
        min(sqrt((m$x[w] - d$x)^2 + (m$y[w] - d$y)^2))
      }, 1)
      if (min(dc) <= 0.25 * d$r_outer) bact <- fill_holes(lb == ids[which.min(dc)])
    }
  }
  stray <- map == 2L & !bact
  if (any(stray)) map <- reassign_pixels(plate, map, stray, exclude = "bacteria")
  map[bact] <- 2L

  # ---- fungus
  fung <- map == 1L
  if (rpx >= 1) fung <- open_disk(close_disk(fung, rpx), rpx)
  # fill holes, but never across bacterial or halo tissue
  fung <- fill_holes(fung, forbid = map %in% c(2L, 3L)) & m$agar & !bact
  fung <- remove_small(fung, min_area_mm2 / mmpp^2)
  dist <- edt(fung)

  exp_pos <- layout_pixels(layout, d)
  search_r <- layout$search_frac * d$r_outer
  n <- nrow(exp_pos)
  markers <- matrix(0L, nrow(fung), ncol(fung))
  core <- data.frame(id = exp_pos$id, expected_x = exp_pos$x, expected_y = exp_pos$y,
                     core_x = NA_real_, core_y = NA_real_, core_depth_mm = NA_real_,
                     detected = FALSE, stringsAsFactors = FALSE)
  if (layout$centre == "fungus") {
    core <- rbind(core, data.frame(id = "centre", expected_x = d$x, expected_y = d$y,
                                   core_x = NA_real_, core_y = NA_real_, core_depth_mm = NA_real_,
                                   detected = FALSE))
  }
  # candidate colony cores: greedy non-maximum suppression on the distance
  # transform (each accepted core suppresses its own inscribed disk)
  cand <- find_cores(dist, fung, min_depth_px = min_core_mm / mmpp,
                     max_n = 4L * nrow(core) + 4L)
  if (nrow(cand)) {
    # match expected positions to cores: closest pairs first, each core used
    # once, cores much shallower than the deepest (satellites) ignored
    cand <- cand[cand$depth >= 0.3 * max(cand$depth), , drop = FALSE]
    dmat <- sqrt(outer(core$expected_x, cand$x, "-")^2 + outer(core$expected_y, cand$y, "-")^2)
    dmat[dmat > search_r] <- Inf
    while (any(is.finite(dmat))) {
      k <- which(dmat == min(dmat), arr.ind = TRUE)[1, ]
      i <- k[1]; j <- k[2]
      core$core_x[i] <- cand$x[j]; core$core_y[i] <- cand$y[j]
      core$core_depth_mm[i] <- cand$depth[j] * mmpp
      core$detected[i] <- TRUE
      dmat[i, ] <- Inf; dmat[, j] <- Inf
    }
  }
  # two expected positions must not share one core
  dup <- duplicated(paste(core$core_x, core$core_y)) & core$detected
  if (any(dup)) {
    cli::cli_warn("Colonies {.val {core$id[dup]}} share a core with another colony and were not detected.")
    core$detected[dup] <- FALSE
  }
  for (i in which(core$detected)) {
    # seed = disk of half the local depth around the core (robust to plateaus)
    rr <- max(1, 0.5 * core$core_depth_mm[i] / mmpp)
    sel <- (m$x - core$core_x[i])^2 + (m$y - core$core_y[i])^2 <= rr^2 & fung
    markers[sel] <- i
  }
  cost <- -gauss_blur(dist, 1)
  # (the centroid check after the watershed rejects "colonies" that are
  # really bands of mis-classified tissue passing through the search area)
  labels <- if (any(markers > 0)) watershed_cpp(cost, markers, fung) else markers
  for (i in which(core$detected)) {
    w <- which(labels == i)
    if (sqrt((mean(m$x[w]) - core$expected_x[i])^2 + (mean(m$y[w]) - core$expected_y[i])^2) > search_r) {
      labels[w] <- 0L
      core$detected[i] <- FALSE
    }
  }
  satellites <- fung & labels == 0L
  # unassigned "fungus" fragments hugging the bacterium are its darker
  # margin (colony edge / swarming front), not fungal satellites
  if (any(bact) && any(satellites)) {
    ls <- label_components(satellites)
    touch <- unique(ls[dilate_disk(bact, 2) & ls > 0])
    if (length(touch)) {
      add <- ls %in% touch
      bact <- fill_holes(bact | add)
      satellites <- satellites & !add
      fung <- fung & !add
    }
  }
  # fragments touching the edge of the analysed disk are dish-wall
  # reflections or meniscus artefacts, not satellite colonies
  if (any(satellites)) {
    ls <- label_components(satellites)
    rim <- m$rr > d$r_agar - 3 & m$agar
    edge_ids <- unique(ls[rim & ls > 0])
    if (length(edge_ids)) satellites <- satellites & !(ls %in% edge_ids)
  }
  if (edge_criterion > 0) labels <- refine_edges(plate, labels, edge_criterion)

  # ---- halo: halo-class pixels connected to the bacterium
  halo <- matrix(FALSE, nrow(map), ncol(map))
  if (any(bact) && any(map == 3L)) {
    h <- map == 3L & !fung
    if (rpx >= 1) h <- open_disk(close_disk(h, rpx), rpx)
    near_b <- edt(!bact) * mmpp <= halo_link_mm
    lh <- label_components(h | bact)
    keep <- unique(lh[bact | (near_b & h)])
    keep <- keep[keep > 0]
    halo <- (lh %in% keep) & h
    halo <- fill_holes(halo | bact) & !bact & !fung
  }
  bact <- bact & !(labels > 0)
  plate$classes$map_final <- map

  plate$objects <- list(labels = labels, satellites = satellites, bacteria = bact,
                        halo = halo, fungus = fung, colony_table = core)
  plate
}

#' Re-assign pixels to their most probable class, excluding some classes
#' @keywords internal
#' @noRd
reassign_pixels <- function(plate, map, which_px, exclude) {
  cl <- plate$classes
  keep_k <- which(!cl$post_names %in% exclude)
  if (!length(keep_k)) {
    map[which_px] <- 4L
    return(map)
  }
  pos <- match(which(which_px), cl$obj_idx)
  ok <- !is.na(pos)
  P <- cl$post[pos[ok], keep_k, drop = FALSE]
  code <- c(fungus = 1L, bacteria = 2L, halo = 3L)[cl$post_names[keep_k]]
  newlab <- code[max.col(P, ties.method = "first")]
  idx <- which(which_px)
  map[idx[ok]] <- newlab
  map[idx[!ok]] <- 4L
  map
}

#' Re-place colony outlines at half contrast (internal)
#'
#' Local version of the 50 % edge criterion: each pixel in a thin band inside
#' the outline is compared with the *local* colony colour (inside the band)
#' and the *local* surrounding colour (just outside the outline, which may
#' be agar, halo or a neighbouring culture), and is kept only if it is
#' closer to the colony. `frac` shifts the threshold (0.5 = half way).
#' @keywords internal
#' @noRd
refine_edges <- function(plate, labels, frac = 0.5, band = 3) {
  S <- plate$bg$smooth
  mmpp <- plate$dish$mm_per_px
  sig <- max(2, 0.25 / mmpp)
  for (i in setdiff(unique(as.vector(labels)), 0L)) {
    cm <- labels == i
    w <- which(cm, arr.ind = TRUE)
    pad <- ceiling(3 * sig + band + 2)
    r1 <- max(1, min(w[, 1]) - pad); r2 <- min(nrow(cm), max(w[, 1]) + pad)
    c1 <- max(1, min(w[, 2]) - pad); c2 <- min(ncol(cm), max(w[, 2]) + pad)
    sub <- cm[r1:r2, c1:c2, drop = FALSE]
    dep_in <- edt(sub)
    dep_out <- edt(!sub)
    inner <- sub & dep_in > band + 1
    outer <- !sub & dep_out > 1 & dep_out <= band + 3
    if (sum(inner) < 50 || sum(outer) < 20) next
    bandm <- sub & dep_in <= band
    d_in <- 0; d_out <- 0
    for (ch in c("L", "a", "b")) {
      x <- S[[ch]][r1:r2, c1:c2, drop = FALSE]
      mi <- gauss_blur(x, sig, w = inner * 1)
      mo <- gauss_blur(x, sig, w = outer * 1)
      d_in <- d_in + (x - mi)^2
      d_out <- d_out + (x - mo)^2
    }
    d_in <- sqrt(d_in); d_out <- sqrt(d_out)
    # keep if the pixel lies on the colony side of the frac point
    drop <- bandm & is.finite(d_in) & is.finite(d_out) &
      d_in / pmax(d_in + d_out, 1e-9) > frac
    keep <- sub & !drop
    keep <- open_disk(keep, 1)
    keep <- largest_component(keep)
    loss <- sub & !keep
    tmp <- labels[r1:r2, c1:c2]
    tmp[loss] <- 0L
    labels[r1:r2, c1:c2] <- tmp
  }
  labels
}

#' Re-assign pixels to their most probable class, excluding some classes
#' @keywords internal
#' @noRd
reassign_pixels <- function(plate, map, which_px, exclude) {
  cl <- plate$classes
  keep_k <- which(!cl$post_names %in% exclude)
  if (!length(keep_k)) {
    map[which_px] <- 4L
    return(map)
  }
  pos <- match(which(which_px), cl$obj_idx)
  ok <- !is.na(pos)
  P <- cl$post[pos[ok], keep_k, drop = FALSE]
  code <- c(fungus = 1L, bacteria = 2L, halo = 3L)[cl$post_names[keep_k]]
  newlab <- code[max.col(P, ties.method = "first")]
  idx <- which(which_px)
  map[idx[ok]] <- newlab
  map[idx[!ok]] <- 4L
  map
}


#' Candidate colony cores by non-maximum suppression of the EDT (internal)
#' @keywords internal
#' @noRd
find_cores <- function(dist, mask, min_depth_px, max_n = 20L) {
  idx <- which(mask & dist >= min_depth_px)
  if (!length(idx)) return(data.frame(x = numeric(), y = numeric(), depth = numeric()))
  H <- nrow(dist)
  xs <- (idx - 1) %/% H + 1; ys <- (idx - 1) %% H + 1
  dv <- dist[idx]
  o <- order(dv, decreasing = TRUE)
  xs <- xs[o]; ys <- ys[o]; dv <- dv[o]
  alive <- rep(TRUE, length(dv))
  out <- list()
  while (any(alive) && length(out) < max_n) {
    j <- which(alive)[1]
    out[[length(out) + 1]] <- c(xs[j], ys[j], dv[j])
    r <- max(dv[j], min_depth_px) * 1.5
    alive[alive] <- (xs[alive] - xs[j])^2 + (ys[alive] - ys[j])^2 > r^2
  }
  o <- do.call(rbind, out)
  data.frame(x = o[, 1], y = o[, 2], depth = o[, 3])
}
