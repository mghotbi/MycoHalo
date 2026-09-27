#' Convert sRGB values to CIE 1976 L*a*b*
#'
#' Converts gamma-encoded sRGB (IEC 61966-2-1) values in \[0, 1\] to CIELAB
#' under the D65 reference white (CIE 15:2004).
#'
#' @details
#' Steps: (1) inverse sRGB companding to linear light,
#' (2) linear RGB -> CIE XYZ with the sRGB/D65 primaries matrix,
#' (3) XYZ -> \eqn{L^*a^*b^*} with the CIE cube-root function and its linear
#' segment below \eqn{(6/29)^3}.
#'
#' \eqn{L^*} (0 = black, 100 = diffuse white) is the perceptual lightness
#' used by MycoHalo as the primary melanization readout; \eqn{a^*}
#' (green-red) and \eqn{b^*} (blue-yellow) carry the chromatic information
#' that separates achromatic melanized hyphae from warm-coloured bacteria,
#' pigments and halos.
#'
#' @param rgb Numeric matrix with three columns (R, G, B) in \[0, 1\], or a
#'   H x W x 3 array.
#'
#' @return Same shape as the input, with channels L, a, b.
#'
#' @examples
#' srgb_to_lab(rbind(c(1, 1, 1), c(0.5, 0.5, 0.5), c(0, 0, 0)))
#'
#' @export
srgb_to_lab <- function(rgb) {
  is_arr <- length(dim(rgb)) == 3L
  d <- dim(rgb)
  if (is_arr) rgb <- matrix(rgb, ncol = 3L)
  if (!is.matrix(rgb) || ncol(rgb) != 3L) {
    cli::cli_abort("{.arg rgb} must be an n x 3 matrix or a H x W x 3 array.")
  }
  lin <- srgb_linearize(rgb)
  xyz <- lin %*% t(.srgb_to_xyz)
  lab <- xyz_to_lab(xyz)
  if (is_arr) {
    lab <- array(lab, dim = d)
  } else {
    colnames(lab) <- c("L", "a", "b")
  }
  lab
}

#' Convert CIE L*a*b* to sRGB
#'
#' Inverse of [srgb_to_lab()]; out-of-gamut values are clipped to \[0, 1\].
#'
#' @param lab Numeric n x 3 matrix (L, a, b) or H x W x 3 array.
#' @return sRGB values in \[0, 1\], same shape as input.
#' @examples
#' lab_to_srgb(cbind(L = 50, a = 0, b = 0))
#' @export
lab_to_srgb <- function(lab) {
  is_arr <- length(dim(lab)) == 3L
  d <- dim(lab)
  if (is_arr) lab <- matrix(lab, ncol = 3L)
  fy <- (lab[, 1] + 16) / 116
  fx <- fy + lab[, 2] / 500
  fz <- fy - lab[, 3] / 200
  finv <- function(t) ifelse(t > 6 / 29, t^3, 3 * (6 / 29)^2 * (t - 4 / 29))
  xyz <- cbind(finv(fx) * .d65[1], finv(fy) * .d65[2], finv(fz) * .d65[3])
  lin <- xyz %*% t(solve(.srgb_to_xyz))
  rgb <- srgb_compand(pmin(pmax(lin, 0), 1))
  rgb <- matrix(rgb, ncol = 3L)
  if (is_arr) array(rgb, dim = d) else rgb
}

# sRGB (D65) primaries -> XYZ, IEC 61966-2-1
.srgb_to_xyz <- matrix(c(
  0.4124564, 0.3575761, 0.1804375,
  0.2126729, 0.7151522, 0.0721750,
  0.0193339, 0.1191920, 0.9503041
), nrow = 3, byrow = TRUE)

# D65 reference white (CIE 1931 2 degree observer), Y normalised to 1
.d65 <- c(0.95047, 1, 1.08883)

#' @keywords internal
#' @noRd
srgb_linearize <- function(v) {
  out <- ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  matrix(out, ncol = ncol(as.matrix(v)))
}

#' @keywords internal
#' @noRd
srgb_compand <- function(v) {
  ifelse(v <= 0.0031308, 12.92 * v, 1.055 * v^(1 / 2.4) - 0.055)
}

#' @keywords internal
#' @noRd
xyz_to_lab <- function(xyz) {
  eps <- (6 / 29)^3
  f <- function(t) ifelse(t > eps, t^(1 / 3), t / (3 * (6 / 29)^2) + 4 / 29)
  fx <- f(xyz[, 1] / .d65[1])
  fy <- f(xyz[, 2] / .d65[2])
  fz <- f(xyz[, 3] / .d65[3])
  cbind(L = 116 * fy - 16, a = 500 * (fx - fy), b = 200 * (fy - fz))
}

#' Chroma and hue angle from a* and b*
#'
#' \eqn{C^*_{ab} = \sqrt{a^{*2} + b^{*2}}}, \eqn{h_{ab} = atan2(b^*, a^*)} in
#' degrees \[0, 360). Melanized hyphae have low chroma; yellow/cream
#' bacterial growth and diffusible pigments have hue angles near 70-100 deg.
#'
#' @param a,b Numeric vectors.
#' @return A list with elements `chroma` and `hue` (degrees).
#' @examples
#' lab_chroma_hue(a = c(0, 5), b = c(0, 20))
#' @export
lab_chroma_hue <- function(a, b) {
  h <- atan2(b, a) * 180 / pi
  list(chroma = sqrt(a^2 + b^2), hue = ifelse(h < 0, h + 360, h))
}

#' Colour difference CIE76 (Euclidean distance in L*a*b*)
#' @param lab1,lab2 n x 3 matrices (or vectors of length 3, recycled).
#' @return Numeric vector of \eqn{\Delta E^*_{ab}}; ~2.3 is a just-noticeable
#'   difference (Sharma 2003).
#' @examples
#' delta_e76(c(50, 0, 0), c(52, 1, -1))
#' @export
delta_e76 <- function(lab1, lab2) {
  lab1 <- matrix(lab1, ncol = 3)
  lab2 <- matrix(lab2, ncol = 3)
  sqrt(rowSums((lab1 - lab2[rep_len(seq_len(nrow(lab2)), nrow(lab1)), , drop = FALSE])^2))
}
