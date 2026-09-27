#' Detect the Petri dish and set the spatial scale
#'
#' Locates the dish as a circle and derives the pixel-to-millimetre scale
#' from its known physical diameter.
#'
#' @details
#' **Algorithm.**
#' \enumerate{
#'   \item The background (e.g. black paper) lightness is estimated from the
#'     four image corners, which lie outside a round dish, as the median
#'     \eqn{L^*} with a robust (MAD) noise estimate \eqn{\sigma}.
#'   \item Pixels brighter than `median + max(6 sigma, 3)` (dark background)
#'     or darker than `median - max(6 sigma, 3)` (light background) form the
#'     foreground; it is cleaned by a morphological opening, the largest
#'     connected component is kept and its holes filled.
#'   \item Boundary pixels of that region, excluding those on the image
#'     border (dishes are often cropped at the sides), are fitted with a
#'     circle by RANSAC (Fischler & Bolles 1981) followed by an algebraic
#'     least-squares refinement (Kasa 1976) on the inliers. Partial arcs
#'     are therefore handled correctly.
#' }
#' The analysed agar disk has radius `agar_fraction * r_outer`, which
#' excludes the dish wall, its reflections and the agar meniscus.
#'
#' **Scale.** `mm_per_px = dish_diameter_mm / (2 * r_outer)`, where
#' `r_outer` is the radius of the edge that was *detected*. Depending on
#' lighting this is the outer rim of the dish or the inner wall / agar edge,
#' so check the dashed circle in [plot_qc()] once for your photo set-up and
#' measure *that* diameter with calipers (standard "90 mm" dishes are
#' ~88-90 mm outside and ~85-86 mm inside). Alternatively photograph a ruler
#' and pass `mm_per_px`.
#'
#' @param plate A `mycohalo_plate` from [read_plate()].
#' @param dish Optional manual circle `c(x, y, r)` in working-image pixels;
#'   skips detection.
#' @param dish_diameter_mm Physical outer diameter of the dish (mm).
#' @param agar_fraction Fraction of the outer radius analysed (0-1).
#' @param background `"auto"`, `"dark"` or `"light"` photographic background.
#' @param mm_per_px Optional scale (mm per *working-image* pixel), e.g. from
#'   a ruler photographed with the plates; overrides `dish_diameter_mm`.
#'   Note: if you measured the scale on the original file, multiply by
#'   `1 / plate$resize_factor`.
#' @param seed RNG seed for RANSAC (results are deterministic by default).
#'
#' @return The plate with element `dish`: a list with `x`, `y`, `r_outer`,
#'   `r_agar` (px), `mm_per_px`, `fit_rmse_px`, `arc_coverage` and `method`.
#'
#' @examples
#' sim <- simulate_plate(seed = 1)
#' plate <- detect_plate(read_plate(sim$image))
#' str(plate$dish)
#' sim$truth$dish
#' @export
detect_plate <- function(plate, dish = NULL, dish_diameter_mm = 90,
                         agar_fraction = 0.88, background = c("auto", "dark", "light"),
                         mm_per_px = NULL, seed = 1L) {
  stopifnot(is_mycohalo_plate(plate))
  background <- match.arg(background)
  if (!is.null(dish)) {
    stopifnot(length(dish) == 3)
    circ <- list(x = dish[1], y = dish[2], r = dish[3], rmse = NA_real_, coverage = NA_real_)
    method <- "manual"
  } else {
    circ <- fit_dish_circle(plate$lab[, , 1], background, seed)
    method <- "auto"
  }
  H <- dim(plate$rgb)[1]; W <- dim(plate$rgb)[2]
  if (method == "auto" && (circ$r < 0.15 * min(H, W) || circ$r > 0.8 * max(H, W))) {
    cli::cli_abort(c(
      "Dish detection failed (implausible radius {round(circ$r)} px).",
      "i" = "Photograph the dish on a uniform, contrasting background, or pass {.arg dish = c(x, y, r)}.",
      "i" = "Use {.fn plot_qc} to check the detected circle."
    ))
  }
  plate$dish <- list(
    x = circ$x, y = circ$y, r_outer = circ$r,
    r_agar = circ$r * agar_fraction,
    mm_per_px = if (is.null(mm_per_px)) dish_diameter_mm / (2 * circ$r) else mm_per_px,
    scale_source = if (is.null(mm_per_px)) "dish diameter" else "user",
    dish_diameter_mm = dish_diameter_mm,
    agar_fraction = agar_fraction,
    fit_rmse_px = circ$rmse,
    arc_coverage = circ$coverage,
    method = method
  )
  plate
}

#' @keywords internal
#' @noRd
fit_dish_circle <- function(L, background, seed) {
  H <- nrow(L); W <- ncol(L)
  Ls <- gauss_blur(L, 1.5)
  # Background = the darkest (or lightest) cluster of the image border.
  # Corners alone are not safe: dishes are often cropped by the frame, so
  # some corners show dish, not background.
  k <- max(3L, round(0.03 * min(H, W)))
  border <- c(Ls[1:k, ], Ls[(H - k + 1):H, ], Ls[, 1:k], Ls[, (W - k + 1):W])
  if (background == "auto") {
    centre <- Ls[round(H * 0.3):round(H * 0.7), round(W * 0.3):round(W * 0.7)]
    background <- if (stats::median(centre) >= stats::median(border)) "dark" else "light"
  }
  side <- if (background == "dark") border[border <= stats::quantile(border, 0.3)]
          else border[border >= stats::quantile(border, 0.7)]
  bg <- stats::median(side)
  sg <- stats::mad(side)
  delta <- max(6 * sg, 3)
  fg <- if (background == "dark") Ls > bg + delta else Ls < bg - delta
  r_open <- max(1, round(0.004 * min(H, W)))
  fg <- open_disk(fg, r_open)
  fg <- largest_component(fg)
  fg <- fill_holes(fg)
  b <- mask_boundary(fg)
  pts <- which(b, arr.ind = TRUE)
  # drop boundary pixels created by the image frame
  keep <- pts[, 1] > 3 & pts[, 1] < H - 2 & pts[, 2] > 3 & pts[, 2] < W - 2
  pts <- pts[keep, , drop = FALSE]
  if (nrow(pts) < 30) cli::cli_abort("Too few dish-edge pixels found for circle fitting.")
  x <- pts[, 2]; y <- pts[, 1]
  fit <- ransac_circle(x, y, tol = max(2, 0.004 * min(H, W)), seed = seed)
  ang <- atan2(y[fit$inliers] - fit$y, x[fit$inliers] - fit$x)
  coverage <- length(unique(cut(ang, breaks = seq(-pi, pi, length.out = 73)))) / 72
  r_ref <- refine_radius(L, fit, bg, background, ang)
  list(x = fit$x, y = fit$y, r = r_ref, r_mask = fit$r, rmse = fit$rmse, coverage = coverage)
}

#' Sub-pixel refinement of the dish radius
#'
#' The binary mask edge is biased outwards by smoothing and by the
#' low threshold. Radial intensity profiles are sampled (bilinear
#' interpolation, 0.25 px steps) across the rim at the inlier angles; the
#' median profile is located where it crosses half-way between the
#' background level and the rim peak (full width at half maximum edge
#' criterion).
#' @keywords internal
#' @noRd
refine_radius <- function(L, fit, bg, background, ang) {
  H <- nrow(L); W <- ncol(L)
  ang <- sort(unique(round(ang, 2)))
  if (length(ang) > 180) ang <- ang[round(seq(1, length(ang), length.out = 180))]
  rs <- seq(fit$r * 0.94, fit$r * 1.04, by = 0.25)
  prof <- sapply(ang, function(a) bilinear(L, fit$x + rs * cos(a), fit$y + rs * sin(a)))
  pm <- apply(prof, 1, stats::median, na.rm = TRUE)
  if (all(is.na(pm))) return(fit$r)
  if (background == "light") pm <- -pm
  bgl <- if (background == "light") -bg else bg
  pk <- which.max(pm)
  half <- (pm[pk] + bgl) / 2
  out <- pk:length(pm)
  j <- out[which(pm[out] < half)[1]]
  if (is.na(j) || j <= 1) return(fit$r)
  # linear interpolation between j-1 and j
  r <- rs[j - 1] + (half - pm[j - 1]) / (pm[j] - pm[j - 1]) * (rs[j] - rs[j - 1])
  if (!is.finite(r) || abs(r - fit$r) > 0.03 * fit$r) fit$r else r
}

#' Bilinear interpolation of a matrix at (x = column, y = row)
#' @keywords internal
#' @noRd
bilinear <- function(M, x, y) {
  H <- nrow(M); W <- ncol(M)
  ok <- x >= 1 & x <= W & y >= 1 & y <= H
  out <- rep(NA_real_, length(x))
  x <- x[ok]; y <- y[ok]
  x0 <- pmin(floor(x), W - 1); y0 <- pmin(floor(y), H - 1)
  fx <- x - x0; fy <- y - y0
  v <- (1 - fx) * (1 - fy) * M[cbind(y0, x0)] + fx * (1 - fy) * M[cbind(y0, x0 + 1)] +
    (1 - fx) * fy * M[cbind(y0 + 1, x0)] + fx * fy * M[cbind(y0 + 1, x0 + 1)]
  out[ok] <- v
  out
}

#' Circle through three points
#' @keywords internal
#' @noRd
circle3 <- function(x, y) {
  A <- cbind(2 * (x[2:3] - x[1]), 2 * (y[2:3] - y[1]))
  bb <- (x[2:3]^2 - x[1]^2) + (y[2:3]^2 - y[1]^2)
  if (abs(det(A)) < 1e-9) return(NULL)
  cxy <- solve(A, bb)
  c(cxy, sqrt((x[1] - cxy[1])^2 + (y[1] - cxy[2])^2))
}

#' Algebraic (Kasa) least-squares circle fit
#' @keywords internal
#' @noRd
kasa_fit <- function(x, y) {
  A <- cbind(x, y, 1)
  bb <- x^2 + y^2
  p <- qr.solve(A, bb)
  cx <- p[1] / 2; cy <- p[2] / 2
  c(cx, cy, sqrt(p[3] + cx^2 + cy^2))
}

#' RANSAC circle fit with least-squares refinement
#' @keywords internal
#' @noRd
ransac_circle <- function(x, y, tol, n_iter = 400L, seed = 1L) {
  n <- length(x)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  best <- NULL; best_n <- 0L
  for (i in seq_len(n_iter)) {
    s <- sample.int(n, 3)
    cc <- circle3(x[s], y[s])
    if (is.null(cc) || !is.finite(cc[3])) next
    d <- abs(sqrt((x - cc[1])^2 + (y - cc[2])^2) - cc[3])
    ni <- sum(d < tol)
    if (ni > best_n) { best_n <- ni; best <- cc }
  }
  if (is.null(best)) cli::cli_abort("Circle fitting failed.")
  inl <- abs(sqrt((x - best[1])^2 + (y - best[2])^2) - best[3]) < tol
  for (k in 1:3) {
    cc <- kasa_fit(x[inl], y[inl])
    res <- sqrt((x - cc[1])^2 + (y - cc[2])^2) - cc[3]
    inl <- abs(res) < tol
  }
  list(x = cc[1], y = cc[2], r = cc[3],
       rmse = sqrt(mean(res[inl]^2)), inliers = which(inl))
}

#' Pixel masks derived from the dish circle
#' @keywords internal
#' @noRd
dish_masks <- function(plate) {
  H <- dim(plate$rgb)[1]; W <- dim(plate$rgb)[2]
  d <- plate$dish
  yy <- matrix(rep(seq_len(H), W), H, W)
  xx <- matrix(rep(seq_len(W), each = H), H, W)
  rr <- sqrt((xx - d$x)^2 + (yy - d$y)^2)
  list(x = xx, y = yy, rr = rr, agar = rr <= d$r_agar)
}
