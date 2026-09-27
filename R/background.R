#' Model the agar background and correct uneven illumination
#'
#' Estimates the colour of the agar, removes the spatial illumination
#' gradient (vignetting, oblique light) by flat-field correction, and flags
#' every dish pixel that differs significantly from agar.
#'
#' @details
#' **1. Agar colour by mode seeking.** The agar is the most frequent colour
#' inside the dish even on crowded confrontation plates, because agar pixels
#' are concentrated in colour space while colonies, bacteria and halos are
#' spread along gradients. The mode of the (smoothed) \eqn{L^*a^*b^*}
#' distribution of low-texture pixels is found with a coarse 3-D histogram
#' followed by mean-shift refinement (Comaniciu & Meer 2002). It works on
#' dark and on light agar.
#'
#' **2. Flat-field correction.** Illumination acts multiplicatively on
#' *linear* light. For each linear RGB channel c a quadratic surface
#' \eqn{\log I_c(x, y) = \beta_0 + \beta_1 x + \beta_2 y + \beta_3 x^2 +
#' \beta_4 xy + \beta_5 y^2} is fitted to agar pixels with iterative
#' 3-MAD trimming, and every pixel is multiplied by
#' \eqn{\exp(\hat\ell_c(x_0, y_0) - \hat\ell_c(x, y))}, i.e. normalised to
#' the illumination at the dish centre. Because the correction is estimated
#' from agar *around* the colonies and applied to the colonies, colony
#' lightness becomes independent of where on the dish the colony sits.
#'
#' **3. Object detection.** After correction the agar has a spatially
#' constant colour \eqn{\mu} with per-channel robust SD \eqn{\sigma} (MAD).
#' A pixel is *object* (non-agar) if its standardised distance
#' \eqn{D^2 = \sum_c ((x_c - \mu_c)/\sigma_c)^2} exceeds the 0.999 quantile
#' of \eqn{\chi^2_3} **and** its colour difference exceeds `min_delta_e`
#' (default 3, just above the ~2.3 just-noticeable difference), so that
#' neither sensor noise nor imperceptible tints are called objects.
#' In addition, a pixel is object when its local texture (SD of \eqn{L^*}
#' in a ~0.2 mm window, log scale) exceeds the agar texture by more than
#' `texture_z` robust standard deviations: strongly melanized, wrinkled
#' colonies can have almost the *mean* colour of a dark agar but are far
#' rougher. Cast shadows of raised colonies (darker than agar, agar-like
#' chromaticity, smooth) are removed from the objects: under directional
#' light they would otherwise be counted as melanized tissue on one side of
#' every colony.
#'
#' @param plate A `mycohalo_plate` with a detected dish.
#' @param flat_field Logical; apply the illumination correction.
#' @param degree Polynomial degree of the illumination surface (1 or 2).
#' @param min_delta_e Minimum CIE76 colour difference from agar for a pixel
#'   to be considered object.
#' @param texture_z Robust z-score of log local texture (relative to agar)
#'   above which a pixel is an object even if its colour matches the agar.
#' @param shadow_texture_z Pixels darker than agar, with agar-like
#'   chromaticity and texture z below this value are cast shadows of raised
#'   colonies and are excluded (set `-Inf` to disable).
#' @param smooth_mm Gaussian smoothing (mm) of colour used for object
#'   detection (measurements always use unsmoothed pixels).
#' @param feature_smooth_mm Gaussian smoothing (mm) of colour used as
#'   classification feature; averages over hyphal ridges and valleys so a
#'   wrinkled colony is described by its mean colour.
#' @param reference Optional exposure / white-balance reference applied
#'   after flat-field correction: a grey-card region `c(x, y, r)` in relative
#'   image coordinates (see [calibrate_color()]), or `"agar"` to scale the
#'   image so that the agar takes the colour `reference_lab` (only valid when
#'   all plates use the same medium batch and the agar is not tinted).
#' @param reference_lab Known CIELAB of the reference.
#'
#' @return The plate with corrected `rgb`/`lab` and an element `bg`.
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 3, vignetting = 0.3)
#' p <- detect_plate(read_plate(sim$image))
#' p <- model_background(p)
#' p$bg$agar_lab
#' p$bg$illumination_range_pct
#' @export
model_background <- function(plate, flat_field = TRUE, degree = 2L,
                              min_delta_e = 3, texture_z = 4, shadow_texture_z = 2.5,
                              smooth_mm = 0.08,
                              feature_smooth_mm = 0.25,
                              reference = NULL, reference_lab = c(50, 0, 0)) {
  stopifnot(is_mycohalo_plate(plate), !is.null(plate$dish))
  m <- dish_masks(plate)
  mmpp <- plate$dish$mm_per_px
  sig_px <- max(1, smooth_mm / mmpp)
  tex_r <- max(1L, round(0.1 / mmpp))
  H <- dim(plate$rgb)[1]; W <- dim(plate$rgb)[2]

  agar_idx <- which(m$agar)
  smooth_lab <- function(lab, sg = sig_px) {
    list(L = gauss_blur(lab[, , 1], sg), a = gauss_blur(lab[, , 2], sg),
         b = gauss_blur(lab[, , 3], sg))
  }
  S <- smooth_lab(plate$lab)
  tex <- local_sd(plate$lab[, , 1], tex_r)

  # ---- 1. agar colour: histogram mode + mean shift on low-texture pixels
  sub <- agar_idx
  if (length(sub) > 200000) sub <- sub[round(seq(1, length(sub), length.out = 200000))]
  X <- cbind(S$L[sub], S$a[sub], S$b[sub])
  low_tex <- tex[sub] <= stats::quantile(tex[sub], 0.5)
  mu <- lab_mode(X[low_tex, , drop = FALSE])
  agar_set <- agar_set_from(X, mu)

  # ---- 2. flat-field correction in linear RGB
  illum_range <- NA_real_
  gain_rng <- c(1, 1)
  if (flat_field) {
    xn <- (m$x - plate$dish$x) / plate$dish$r_outer
    yn <- (m$y - plate$dish$y) / plate$dish$r_outer
    design <- function(x, y) {
      if (degree >= 2) cbind(1, x, y, x^2, x * y, y^2) else cbind(1, x, y)
    }
    fit_idx <- sub[agar_set]
    if (length(fit_idx) > 60000) fit_idx <- fit_idx[round(seq(1, length(fit_idx), length.out = 60000))]
    Dfit <- design(xn[fit_idx], yn[fit_idx])
    lin <- srgb_linearize(matrix(plate$rgb, ncol = 3))
    Dall <- design(as.vector(xn), as.vector(yn))
    d0 <- design(0, 0)
    ysurf <- NULL
    for (ch in 1:3) {
      yv <- log(lin[fit_idx, ch] + 1e-4)
      keep <- rep(TRUE, length(yv))
      for (it in 1:4) {
        beta <- qr.coef(qr(Dfit[keep, , drop = FALSE]), yv[keep])
        beta[is.na(beta)] <- 0
        res <- yv - Dfit %*% beta
        s <- stats::mad(res[keep])
        keep <- abs(res) <= 3 * max(s, 1e-6)
      }
      surf <- as.vector(Dall %*% beta)
      # normalise to the illumination at the dish centre: a fixed geometric
      # reference, identical for every plate. The agar surrounds the centre,
      # so the surface is interpolated there even when a bacterium covers it.
      g <- exp(as.numeric(d0 %*% beta) - surf)
      g <- pmin(pmax(g, 0.5), 2)
      lin[, ch] <- lin[, ch] * g
      if (ch == 2) ysurf <- exp(surf)
      gain_rng <- range(c(gain_rng, g[m$agar]))
    }
    ys <- ysurf[as.vector(m$agar)]
    illum_range <- 100 * (max(ys) - min(ys)) / mean(ys)
    rgb <- srgb_compand(pmin(pmax(lin, 0), 1))
    plate$rgb <- array(rgb, dim = dim(plate$rgb))
    plate$lab <- srgb_to_lab(plate$rgb)
    S <- smooth_lab(plate$lab)
    X <- cbind(S$L[sub], S$a[sub], S$b[sub])
    mu <- lab_mode(X[low_tex, , drop = FALSE])
    agar_set <- agar_set_from(X, mu)
  }

  # ---- 2b. exposure / white-balance calibration (after flat-field, so the
  # reference is expressed at the illumination of the dish centre)
  if (!is.null(reference)) {
    if (identical(reference, "agar")) {
      ref_mask <- matrix(FALSE, H, W)
      ref_mask[sub[agar_set]] <- TRUE
      where <- "agar"
    } else {
      stopifnot(is.numeric(reference), length(reference) == 3)
      ref_mask <- (m$x - reference[1] * W)^2 + (m$y - reference[2] * H)^2 <=
        (reference[3] * min(H, W))^2
      where <- "reference patch"
      if (flat_field && any(ref_mask & m$rr > plate$dish$r_outer)) {
        cli::cli_inform(c("i" = "The reference patch lies outside the dish; its illumination is extrapolated from the flat-field model. Place the card close to the dish and evenly lit."))
      }
    }
    plate <- apply_reference_gain(plate, ref_mask, reference_lab, where)
    S <- smooth_lab(plate$lab)
    X <- cbind(S$L[sub], S$a[sub], S$b[sub])
    mu <- lab_mode(X[low_tex, , drop = FALSE])
    agar_set <- agar_set_from(X, mu)
  }

  # ---- 3. robust agar statistics and object mask
  for (it in 1:2) {
    mu <- apply(X[agar_set, , drop = FALSE], 2, stats::median)
    sdv <- pmax(apply(X[agar_set, , drop = FALSE], 2, stats::mad), 0.25)
    D2 <- rowSums(sweep(sweep(X, 2, mu), 2, sdv, "/")^2)
    dE <- sqrt(rowSums(sweep(X, 2, mu)^2))
    agar_set <- !(D2 > stats::qchisq(0.999, 3) & dE > min_delta_e)
  }
  D2all <- ((S$L - mu[1]) / sdv[1])^2 + ((S$a - mu[2]) / sdv[2])^2 + ((S$b - mu[3]) / sdv[3])^2
  dEall <- sqrt((S$L - mu[1])^2 + (S$a - mu[2])^2 + (S$b - mu[3])^2)
  colour_obj <- D2all > stats::qchisq(0.999, 3) & dEall > min_delta_e

  # texture: heavily melanized colonies can have almost the same *mean*
  # colour as dark agar, but their wrinkled surface has a local SD of L*
  # several times that of agar. Texture outliers are objects too.
  tex <- local_sd(plate$lab[, , 1], tex_r)
  ltex <- log(tex + 0.1)
  at <- ltex[sub][agar_set]
  tmed <- stats::median(at); tmad <- max(stats::mad(at), 0.05)
  tz <- gauss_blur((ltex - tmed) / tmad, max(1, tex_r))
  tex_obj <- tz > texture_z
  obj <- m$agar & (colour_obj | tex_obj)
  # cast shadows of raised colonies: agar under less light is darker than
  # agar, keeps the agar's chromaticity and stays smooth. Melanized tissue
  # can be darker than agar too, but it is rough (hyphal texture).
  dLs <- S$L - mu[1]
  dC <- sqrt((S$a - mu[2])^2 + (S$b - mu[3])^2)
  shadow <- obj & dLs < 0 & tz < shadow_texture_z & dC < pmax(3, 0.5 * abs(dLs))
  if (any(shadow)) shadow <- remove_small(shadow, max(5, round(0.05 / mmpp^2)))
  obj <- obj & !shadow
  obj <- remove_small(obj, max(5, round(0.05 / mmpp^2)))
  S_feat <- smooth_lab(plate$lab, max(1, feature_smooth_mm / mmpp))

  plate$bg <- list(
    agar_lab = stats::setNames(mu, c("L", "a", "b")),
    agar_sd = stats::setNames(sdv, c("L", "a", "b")),
    object = obj,
    delta_e = dEall,
    smooth = S,
    smooth_feat = S_feat,
    texture = tex,
    texture_z = tz,
    shadow = shadow,
    flat_field = flat_field,
    illumination_range_pct = illum_range,
    gain_range = gain_rng,
    min_delta_e = min_delta_e
  )
  plate
}

#' Mode of a 3-D colour distribution: coarse histogram + mean shift
#' @keywords internal
#' @noRd
lab_mode <- function(X, bin = 2, h = 4) {
  key <- paste(floor(X[, 1] / bin), floor(X[, 2] / bin), floor(X[, 3] / bin))
  tb <- table(key)
  top <- names(tb)[which.max(tb)]
  centre <- colMeans(X[key == top, , drop = FALSE])
  for (it in 1:20) {
    d <- sqrt(rowSums(sweep(X, 2, centre)^2))
    new <- colMeans(X[d < h, , drop = FALSE])
    if (sqrt(sum((new - centre)^2)) < 0.01) break
    centre <- new
  }
  centre
}

#' Initial agar pixel set around a colour centre
#' @keywords internal
#' @noRd
agar_set_from <- function(X, mu) {
  near <- sqrt(rowSums(sweep(X, 2, mu)^2)) < 8
  s <- pmax(apply(X[near, , drop = FALSE], 2, stats::mad), 0.25)
  D2 <- rowSums(sweep(sweep(X, 2, mu), 2, s, "/")^2)
  D2 < stats::qchisq(0.999, 3) | sqrt(rowSums(sweep(X, 2, mu)^2)) < 3
}
