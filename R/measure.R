#' Measure colonies: size, shape, melanization and interaction metrics
#'
#' @details
#' ## Melanization readouts
#' All colour statistics use unsmoothed, flat-field-corrected pixels of the
#' colony *interior*: a rim of `edge_exclude_mm` is removed (mixed
#' colony/agar pixels, "partial-volume" effect) and specular highlights
#' (\eqn{L^* > 97} or a saturated channel) are excluded.
#' \describe{
#'   \item{`L_mean`, `L_median`, `L_sd`, `L_q05`...`L_q95`}{CIE lightness
#'     \eqn{L^*} (0 black - 100 white). The primary, perceptually uniform
#'     measure: a difference of 1 unit is the same visual darkening
#'     anywhere on the scale.}
#'   \item{`MI_mean`, `MI_median`}{Melanization index \eqn{MI = 100 - L^*}.
#'     Higher = darker = more melanized. Linear in \eqn{L^*}, so MI
#'     differences are comparable across the whole range.}
#'   \item{`MI_lo`, `MI_hi`}{95 % spatial block-bootstrap confidence interval
#'     of `MI_mean` (see [colony_profiles()]).}
#'   \item{`gray_imagej`}{Mean of (R + G + B) / 3 on the 0-255 scale,
#'     identical to ImageJ/Fiji's default RGB -> 8-bit conversion, for
#'     comparison with earlier *Z. tritici* melanization studies that report
#'     mean grey values.}
#'   \item{`a_mean`, `b_mean`, `chroma_mean`, `hue_deg`}{Chromatic
#'     components; DHN-melanin is dark brown-black, so strongly melanized
#'     colonies drift towards low \eqn{L^*} with slightly positive
#'     \eqn{a^*, b^*} (brown), whereas "grey" colonies are nearly neutral.}
#'   \item{`dark_fraction`}{Fraction of interior pixels with
#'     \eqn{L^* <} `dark_L` (only if `dark_L` is given). Choose `dark_L`
#'     from control plates, e.g. the 10th percentile of control colonies.}
#' }
#'
#' ## Interaction readouts
#' The *interaction direction* of a colony is the vector from its centroid
#' to the centroid of the bacterial colony (confrontations), or to the dish
#' centre (controls, so that the same metrics give the null expectation).
#' \describe{
#'   \item{`MI_facing`, `MI_away`, `delta_MI_facing`}{MI of the colony half
#'     facing the bacterium (angle < 90 deg to the interaction direction)
#'     and of the opposite half, and their difference. A localized
#'     melanization response to the bacterium gives `delta_MI_facing > 0`.}
#'   \item{`R_toward_mm`, `R_away_mm`, `growth_inhibition_pct`}{Colony radius
#'     from the centroid towards / away from the bacterium (maximum radius
#'     in a +/- 20 deg sector), and the percent inhibition of radial growth
#'     PIRG = 100 (R_away - R_toward) / R_away, the standard dual-culture
#'     antagonism index.}
#'   \item{`gap_bacteria_mm`, `gap_halo_mm`}{Shortest edge-to-edge distance
#'     from the colony to the bacterial colony and to the halo (0 = contact).}
#' }
#'
#' @param plate A plate processed by [segment_colonies()].
#' @param edge_exclude_mm Width of the colony margin excluded from colour
#'   statistics (mm).
#' @param dark_L Optional \eqn{L^*} threshold for `dark_fraction`.
#' @param n_boot Bootstrap replicates for confidence intervals.
#' @param block_mm Side of the square spatial blocks for the bootstrap (mm).
#' @param seed RNG seed for the bootstrap.
#'
#' @return Data frame with one row per expected colony.
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 7)
#' res <- analyze_plate(sim$image, verbose = FALSE)
#' res$colonies[, c("colony_id", "area_mm2", "L_mean", "MI_mean")]
#' @export
measure_colonies <- function(plate, edge_exclude_mm = 0.1, dark_L = NULL,
                             n_boot = 200L, block_mm = 0.5, seed = 1L) {
  stopifnot(is_mycohalo_plate(plate), !is.null(plate$objects))
  ob <- plate$objects
  d <- plate$dish
  mmpp <- d$mm_per_px
  m <- dish_masks(plate)
  L <- plate$lab[, , 1]; A <- plate$lab[, , 2]; B <- plate$lab[, , 3]
  rgb <- plate$rgb
  gray <- (rgb[, , 1] + rgb[, , 2] + rgb[, , 3]) / 3 * 255
  sat <- rgb[, , 1] >= 254.5 / 255 | rgb[, , 2] >= 254.5 / 255 | rgb[, , 3] >= 254.5 / 255
  spec <- L > 97 | sat
  has_bact <- any(ob$bacteria)
  if (has_bact) {
    bw <- which(ob$bacteria)
    tx <- mean(m$x[bw]); ty <- mean(m$y[bw])
    dist_b <- edt(!ob$bacteria) * mmpp
    bh <- ob$bacteria | ob$halo
    dist_h <- edt(!bh) * mmpp
  } else {
    tx <- d$x; ty <- d$y
  }
  block_id <- (floor(m$x * mmpp / block_mm)) * 1e5 + floor(m$y * mmpp / block_mm)
  ct <- ob$colony_table

  rows <- lapply(seq_len(nrow(ct)), function(i) {
    base <- data.frame(colony_id = ct$id[i], detected = ct$detected[i],
                       stringsAsFactors = FALSE)
    if (!ct$detected[i]) return(base)
    cm <- ob$labels == i
    w <- which(cm)
    area_px <- length(w)
    cx <- mean(m$x[w]); cy <- mean(m$y[w])
    per <- perimeter_px(cm)
    hull <- convex_hull_area(m$x[w], m$y[w])
    depth <- edt_local(cm)
    interior <- cm & depth * mmpp > edge_exclude_mm & !spec
    wi <- which(interior)
    if (length(wi) < 30) wi <- which(cm & !spec)
    Lv <- L[wi]
    q <- stats::quantile(Lv, c(0.05, 0.25, 0.75, 0.95), names = FALSE)
    am <- mean(A[wi]); bm <- mean(B[wi])
    ch <- lab_chroma_hue(am, bm)
    ci <- block_boot_mean(100 - Lv, block_id[wi], n_boot, seed)

    # interaction geometry
    vx <- tx - cx; vy <- ty - cy
    nv <- sqrt(vx^2 + vy^2)
    ang <- atan2(m$y[wi] - cy, m$x[wi] - cx) - atan2(vy, vx)
    ang <- abs(atan2(sin(ang), cos(ang)))
    facing <- ang < pi / 2
    angb <- atan2(m$y[w] - cy, m$x[w] - cx) - atan2(vy, vx)
    angb <- atan2(sin(angb), cos(angb))
    rad <- sqrt((m$x[w] - cx)^2 + (m$y[w] - cy)^2) * mmpp
    cone <- 20 * pi / 180
    r_to <- if (any(abs(angb) < cone)) max(rad[abs(angb) < cone]) else NA_real_
    r_aw <- if (any(abs(angb) > pi - cone)) max(rad[abs(angb) > pi - cone]) else NA_real_

    out <- data.frame(
      base,
      centroid_x_px = cx, centroid_y_px = cy,
      centroid_x_mm = (cx - d$x) * mmpp, centroid_y_mm = (cy - d$y) * mmpp,
      offset_from_expected_mm = sqrt((cx - ct$expected_x[i])^2 + (cy - ct$expected_y[i])^2) * mmpp,
      area_mm2 = area_px * mmpp^2,
      eq_diameter_mm = 2 * sqrt(area_px / pi) * mmpp,
      perimeter_mm = per * mmpp,
      circularity = min(1, 4 * pi * area_px / per^2),
      solidity = area_px / hull,
      max_depth_mm = max(depth[w]) * mmpp,
      n_pixels = length(wi),
      L_mean = mean(Lv), L_median = stats::median(Lv), L_sd = stats::sd(Lv),
      L_q05 = q[1], L_q25 = q[2], L_q75 = q[3], L_q95 = q[4],
      MI_mean = 100 - mean(Lv), MI_median = 100 - stats::median(Lv),
      MI_lo = ci[1], MI_hi = ci[2],
      gray_imagej = mean(gray[wi]),
      a_mean = am, b_mean = bm, chroma_mean = mean(sqrt(A[wi]^2 + B[wi]^2)),
      hue_deg = ch$hue,
      dark_fraction = if (is.null(dark_L)) NA_real_ else mean(Lv < dark_L),
      specular_fraction = mean(spec[w]),
      interaction_target = if (has_bact) "bacteria" else "dish_centre",
      dist_to_target_mm = nv * mmpp,
      MI_facing = if (any(facing)) 100 - mean(Lv[facing]) else NA_real_,
      MI_away = if (any(!facing)) 100 - mean(Lv[!facing]) else NA_real_,
      R_toward_mm = r_to, R_away_mm = r_aw,
      stringsAsFactors = FALSE
    )
    out$delta_MI_facing <- out$MI_facing - out$MI_away
    out$growth_inhibition_pct <- 100 * (r_aw - r_to) / r_aw
    out$gap_bacteria_mm <- if (has_bact) max(0, min(dist_b[w]) - mmpp) else NA_real_
    out$gap_halo_mm <- if (has_bact) max(0, min(dist_h[w]) - mmpp) else NA_real_
    out
  })
  bind_rows_fill(rows)
}

#' Edge-to-centre melanization profiles
#'
#' Melanization of a *Z. tritici* colony is not uniform: the growing margin
#' consists of young, often lighter hyphae, while older central tissue is
#' darker. Profiles as a function of distance from the colony **edge**
#' (not from the centroid) follow this age gradient for irregular, lobed
#' colonies too.
#'
#' @details
#' The Euclidean distance transform of each colony gives every pixel its
#' distance to the nearest colony edge. Pixels are grouped into rings of
#' width `ring_width_mm` (absolute mode, comparable between colonies of
#' different size) or into `n_rings` rings of equal relative depth
#' (relative mode, comparable in developmental stage).
#'
#' **Confidence intervals.** Neighbouring pixels are strongly correlated
#' (same hypha, JPEG blocks, smoothing in the camera), so a naive pixel
#' bootstrap grossly overstates precision. MycoHalo uses a *spatial block
#' bootstrap* (Künsch 1989; Lahiri 2003): the colony is divided into square
#' blocks of `block_mm`, and whole blocks are resampled with replacement.
#' Choose `block_mm` at least as large as the wrinkle / texture scale
#' (default 0.5 mm).
#'
#' @param result A `mycohalo_result` from [analyze_plate()], or a plate
#'   processed by [segment_colonies()].
#' @param ring_width_mm Ring width in absolute mode (mm).
#' @param mode `"absolute"` or `"relative"`.
#' @param n_rings Number of rings in relative mode.
#' @param n_boot,block_mm,seed Block-bootstrap settings.
#' @return Data frame with one row per colony x ring.
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 8, edge_lightening = 8)
#' res <- analyze_plate(sim$image, verbose = FALSE)
#' head(res$profiles)
#' @export
colony_profiles <- function(result, ring_width_mm = 0.5, mode = c("absolute", "relative"),
                            n_rings = 10L, n_boot = 200L, block_mm = 0.5, seed = 1L) {
  mode <- match.arg(mode)
  plate <- if (inherits(result, "mycohalo_result")) result$plate else result
  ob <- plate$objects
  mmpp <- plate$dish$mm_per_px
  m <- dish_masks(plate)
  L <- plate$lab[, , 1]
  rgb <- plate$rgb
  spec <- L > 97 | rgb[, , 1] >= 254.5 / 255 | rgb[, , 2] >= 254.5 / 255 | rgb[, , 3] >= 254.5 / 255
  block_id <- (floor(m$x * mmpp / block_mm)) * 1e5 + floor(m$y * mmpp / block_mm)
  ct <- ob$colony_table
  rows <- lapply(which(ct$detected), function(i) {
    cm <- ob$labels == i
    dep <- edt_local(cm) * mmpp
    w <- which(cm & !spec)
    dv <- dep[w] - 0.5 * mmpp          # pixel-centre distance from edge
    dmax <- max(dv)
    if (mode == "absolute") {
      br <- seq(0, ceiling(dmax / ring_width_mm) * ring_width_mm, by = ring_width_mm)
    } else {
      br <- seq(0, dmax, length.out = n_rings + 1)
    }
    ring <- findInterval(dv, br, rightmost.closed = TRUE, all.inside = TRUE)
    do.call(rbind, lapply(sort(unique(ring)), function(k) {
      s <- ring == k
      Lv <- L[w[s]]
      ci <- block_boot_mean(100 - Lv, block_id[w[s]], n_boot, seed)
      data.frame(colony_id = ct$id[i], ring = k,
                 from_edge_mm = br[k], to_edge_mm = br[k + 1],
                 mid_mm = (br[k] + br[k + 1]) / 2,
                 rel_depth = ((br[k] + br[k + 1]) / 2) / dmax,
                 n_pixels = sum(s), n_blocks = length(unique(block_id[w[s]])),
                 L_mean = mean(Lv), L_sd = stats::sd(Lv),
                 MI_mean = 100 - mean(Lv), MI_lo = ci[1], MI_hi = ci[2],
                 stringsAsFactors = FALSE)
    }))
  })
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Spatial block bootstrap CI of a mean
#' @keywords internal
#' @noRd
block_boot_mean <- function(v, block, n_boot = 200L, seed = 1L, level = 0.95) {
  if (length(v) < 2 || n_boot < 10) return(c(NA_real_, NA_real_))
  f <- factor(block)
  sums <- as.vector(tapply(v, f, sum))
  cnts <- as.vector(tapply(v, f, length))
  nb <- length(sums)
  if (nb < 3) return(c(NA_real_, NA_real_))
  Wm <- with_seed(if (is.null(seed)) sample.int(1e6, 1) else seed,
                  stats::rmultinom(n_boot, nb, rep(1 / nb, nb)))
  bm <- colSums(Wm * sums) / colSums(Wm * cnts)
  a <- (1 - level) / 2
  stats::quantile(bm, c(a, 1 - a), names = FALSE)
}

#' Plate-level measurements: dish, background, bacterium, halo, QC
#' @keywords internal
#' @noRd
measure_plate <- function(plate) {
  ob <- plate$objects
  d <- plate$dish
  mmpp <- d$mm_per_px
  m <- dish_masks(plate)
  L <- plate$lab[, , 1]; B <- plate$lab[, , 3]
  mu <- plate$bg$agar_lab
  ba <- sum(ob$bacteria) * mmpp^2
  ha <- sum(ob$halo) * mmpp^2
  halo_r <- NA_real_; halo_r95 <- NA_real_; bact_r <- NA_real_
  halo_db <- NA_real_; halo_dL <- NA_real_; bact_L <- NA_real_
  if (ba > 0) {
    bw <- which(ob$bacteria)
    bx <- mean(m$x[bw]); by <- mean(m$y[bw])
    bact_r <- sqrt(ba / pi)
    bact_L <- mean(L[bw])
    if (ha > 0) {
      hw <- which(ob$halo)
      halo_r <- sqrt((ba + ha) / pi)
      halo_r95 <- stats::quantile(sqrt((m$x[hw] - bx)^2 + (m$y[hw] - by)^2), 0.95, names = FALSE) * mmpp
      halo_db <- mean(B[hw]) - mu[["b"]]
      halo_dL <- mean(L[hw]) - mu[["L"]]
    }
  }
  data.frame(
    dish_x_px = d$x, dish_y_px = d$y, dish_r_px = d$r_outer,
    mm_per_px = mmpp, dish_fit_rmse_px = d$fit_rmse_px, dish_arc_coverage = d$arc_coverage,
    agar_L = mu[["L"]], agar_a = mu[["a"]], agar_b = mu[["b"]],
    agar_sd_L = plate$bg$agar_sd[["L"]],
    illumination_range_pct = plate$bg$illumination_range_pct,
    n_colonies_detected = sum(ob$colony_table$detected),
    n_satellites = if (any(ob$satellites)) max(label_components(ob$satellites)) else 0L,
    satellite_area_mm2 = sum(ob$satellites) * mmpp^2,
    bacteria_area_mm2 = ba, bacteria_radius_mm = bact_r, bacteria_L = bact_L,
    halo_area_mm2 = ha, halo_outer_radius_mm = halo_r, halo_r95_mm = halo_r95,
    halo_delta_b = halo_db, halo_delta_L = halo_dL,
    other_area_mm2 = sum((if (is.null(plate$classes$map_final)) plate$classes$map else plate$classes$map_final) == 4L) * mmpp^2,
    overexposed_fraction = plate$clipped_fraction,
    classifier = plate$classes$source,
    stringsAsFactors = FALSE
  )
}

#' rbind data frames with differing columns
#' @keywords internal
#' @noRd
bind_rows_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, TRUE)]
  if (length(rows) == 0) return(data.frame())
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) {
    miss <- setdiff(cols, names(r))
    for (mm in miss) r[[mm]] <- NA
    r[cols]
  })
  do.call(rbind, rows)
}
