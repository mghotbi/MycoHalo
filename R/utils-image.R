# Internal image-processing helpers built on the C++ kernels in src/.
# Morphology uses the exact Euclidean distance transform, which is equivalent
# to erosion / dilation with a perfect Euclidean disk of any radius.

#' Euclidean distance transform
#' @param mask Logical matrix.
#' @return Distance (px) from each TRUE pixel to the nearest FALSE pixel.
#' @keywords internal
#' @noRd
edt <- function(mask) {
  storage.mode(mask) <- "logical"
  edt_cpp(mask)
}

#' Connected components
#' @keywords internal
#' @noRd
label_components <- function(mask, connectivity = 8L) {
  storage.mode(mask) <- "logical"
  label_cpp(mask, as.integer(connectivity))
}

#' Morphological erosion with a Euclidean disk of radius r (px)
#' @keywords internal
#' @noRd
erode_disk <- function(mask, r) {
  if (r <= 0) return(mask)
  d <- edt(mask)
  d > r
}

#' Morphological dilation with a Euclidean disk of radius r (px)
#' @keywords internal
#' @noRd
dilate_disk <- function(mask, r) {
  if (r <= 0) return(mask)
  d <- edt(!mask)
  d <= r
}

#' Opening (erosion then dilation): removes objects / spurs thinner than 2r
#' @keywords internal
#' @noRd
open_disk <- function(mask, r) dilate_disk(erode_disk(mask, r), r)

#' Closing (dilation then erosion): bridges gaps narrower than 2r
#' @keywords internal
#' @noRd
close_disk <- function(mask, r) erode_disk(dilate_disk(mask, r), r)

#' Fill holes: background regions not connected to the image border
#'
#' @param forbid Optional logical matrix; holes containing any forbidden
#'   pixel (e.g. bacterial tissue enclosed by a ring of fungus) stay open.
#' @keywords internal
#' @noRd
fill_holes <- function(mask, forbid = NULL) {
  bg <- label_components(!mask, 4L)
  border <- unique(c(bg[1, ], bg[nrow(bg), ], bg[, 1], bg[, ncol(bg)]))
  border <- border[border > 0]
  if (!is.null(forbid) && any(forbid & bg > 0)) {
    border <- unique(c(border, unique(bg[forbid & bg > 0])))
  }
  out <- mask
  out[bg > 0 & !(bg %in% border)] <- TRUE
  out
}

#' Remove connected components smaller than min_px
#' @keywords internal
#' @noRd
remove_small <- function(mask, min_px) {
  lab <- label_components(mask)
  if (max(lab) == 0) return(mask)
  sz <- tabulate(lab[lab > 0], nbins = max(lab))
  keep <- which(sz >= min_px)
  mask & (lab %in% keep)
}

#' Keep the largest connected component
#' @keywords internal
#' @noRd
largest_component <- function(mask) {
  lab <- label_components(mask)
  if (max(lab) == 0) return(mask)
  sz <- tabulate(lab[lab > 0], nbins = max(lab))
  lab == which.max(sz)
}

#' Local mean over a (2r+1)^2 window, optionally weighted (e.g. by a mask)
#' @keywords internal
#' @noRd
box_mean <- function(x, r, w = NULL) {
  r <- as.integer(max(0, round(r)))
  if (r == 0) return(x)
  if (is.null(w)) {
    return(box_blur_cpp(x, r, 1L))
  } else {
    w <- w * 1
    num <- box_sum_cpp(x * w, r)
    den <- box_sum_cpp(w, r)
  }
  out <- num / den
  out[den == 0] <- NA_real_
  out
}

#' Gaussian blur approximated by three successive box filters
#'
#' Three passes of a box filter of width w = 2r + 1 approximate a Gaussian
#' with variance 3 (w^2 - 1) / 12 (central-limit theorem; Wells 1986,
#' IEEE PAMI 8:234-239). With weights `w`, normalised convolution
#' blur(x * w) / blur(w) is returned, so masked pixels do not leak in.
#' @keywords internal
#' @noRd
gauss_blur <- function(x, sigma, w = NULL) {
  if (sigma <= 0.5) return(x)
  r <- max(1L, as.integer(round((sqrt(4 * sigma^2 + 1) - 1) / 2)))
  if (is.null(w)) return(box_blur_cpp(x, r, 3L))
  num <- x * w
  num[is.na(num)] <- 0
  num <- box_blur_cpp(num, r, 3L)
  den <- box_blur_cpp(w * 1, r, 3L)
  out <- num / den
  out[den < 1e-6] <- NA_real_
  out
}

#' Local standard deviation (texture) over a (2r+1)^2 window
#' @keywords internal
#' @noRd
local_sd <- function(x, r, w = NULL) {
  m1 <- box_mean(x, r, w)
  m2 <- box_mean(x^2, r, w)
  sqrt(pmax(m2 - m1^2, 0))
}

#' Otsu's threshold (Otsu 1979, IEEE Trans. Syst. Man Cybern. 9:62-66)
#'
#' Maximises the between-class variance of a 1-D histogram.
#' @keywords internal
#' @noRd
otsu_threshold <- function(x, nbins = 256L) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || diff(range(x)) == 0) return(mean(x))
  br <- seq(min(x), max(x), length.out = nbins + 1)
  h <- graphics::hist(x, breaks = br, plot = FALSE)
  p <- h$counts / sum(h$counts)
  mids <- h$mids
  w0 <- cumsum(p)
  mu <- cumsum(p * mids)
  muT <- mu[length(mu)]
  sb <- (muT * w0 - mu)^2 / (w0 * (1 - w0))
  sb[!is.finite(sb)] <- 0
  br[which.max(sb) + 1]
}

#' Boundary pixels of a mask (4-neighbourhood)
#' @keywords internal
#' @noRd
mask_boundary <- function(mask) {
  nr <- nrow(mask); nc <- ncol(mask)
  up    <- rbind(FALSE, mask[-nr, , drop = FALSE])
  down  <- rbind(mask[-1, , drop = FALSE], FALSE)
  left  <- cbind(FALSE, mask[, -nc, drop = FALSE])
  right <- cbind(mask[, -1, drop = FALSE], FALSE)
  mask & !(up & down & left & right)
}

#' Perimeter estimate of a binary region (px)
#'
#' Counts "crack" edges (pixel sides separating region from background) and
#' applies the Cauchy-Crofton isotropy correction pi/4: for boundaries with
#' random orientation the city-block crack length overestimates the
#' Euclidean length by 4/pi on average (e.g. Dorst & Smeulders 1987).
#' @keywords internal
#' @noRd
perimeter_px <- function(mask) {
  m <- rbind(FALSE, cbind(FALSE, mask, FALSE), FALSE)
  nr <- nrow(m); nc <- ncol(m)
  cracks <- sum(m[-1, ] != m[-nr, ]) + sum(m[, -1] != m[, -nc])
  cracks * pi / 4
}

#' Area of the convex hull of a set of pixels (px^2)
#' @keywords internal
#' @noRd
convex_hull_area <- function(x, y) {
  if (length(x) < 3) return(length(x))
  # use pixel corners so that hull encloses whole pixels
  xx <- c(x - 0.5, x + 0.5, x - 0.5, x + 0.5)
  yy <- c(y - 0.5, y - 0.5, y + 0.5, y + 0.5)
  h <- grDevices::chull(xx, yy)
  px <- xx[h]; py <- yy[h]
  abs(sum(px * c(py[-1], py[1]) - c(px[-1], px[1]) * py)) / 2
}

#' Resample a logical/integer mask with nearest neighbour
#' @keywords internal
#' @noRd
resize_nn <- function(m, nr, nc) {
  ri <- pmin(nrow(m), pmax(1L, round(seq(0.5, nrow(m) - 0.5, length.out = nr) + 0.5)))
  ci <- pmin(ncol(m), pmax(1L, round(seq(0.5, ncol(m) - 0.5, length.out = nc) + 0.5)))
  m[ri, ci, drop = FALSE]
}

#' EDT restricted to the bounding box of a mask (fast for small objects)
#' @keywords internal
#' @noRd
edt_local <- function(mask, pad = 2L) {
  w <- which(mask, arr.ind = TRUE)
  out <- matrix(0, nrow(mask), ncol(mask))
  if (!nrow(w)) return(out)
  r1 <- max(1, min(w[, 1]) - pad); r2 <- min(nrow(mask), max(w[, 1]) + pad)
  c1 <- max(1, min(w[, 2]) - pad); c2 <- min(ncol(mask), max(w[, 2]) + pad)
  sub <- mask[r1:r2, c1:c2, drop = FALSE]
  if (r1 == 1 || c1 == 1 || r2 == nrow(mask) || c2 == ncol(mask)) {
    sub <- rbind(FALSE, cbind(FALSE, sub, FALSE), FALSE)
    d <- edt(sub)[-c(1, nrow(sub)), -c(1, ncol(sub)), drop = FALSE]
  } else {
    d <- edt(sub)
  }
  out[r1:r2, c1:c2] <- d
  out
}
