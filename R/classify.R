#' Classify non-agar pixels into fungus, bacteria and halo
#'
#' This step solves the core problem of confrontation plates: grey fungal
#' colonies, the whitish bacterial colony and the yellowish halo overlap in
#' brightness but not in *colour and texture*.
#'
#' @details
#' **Features** (per object pixel, after flat-field correction and light
#' Gaussian smoothing). The colour difference to the agar of the *same*
#' plate, \eqn{(\Delta L^*, \Delta a^*, \Delta b^*)}, is split into its
#' magnitude \eqn{\log \Delta E} and its direction (unit vector). A
#' diffusible pigment at decreasing concentration moves the colour along a
#' fixed direction from the agar towards the pigment colour, so the faint
#' fringe of a halo shares the direction of its core, while grey fungal
#' tissue points along the achromatic \eqn{L^*} axis at every
#' melanization level. The fifth feature is texture, the log local standard
#' deviation of \eqn{L^*} in a 0.4 mm window (wrinkled hyphal colonies vs.
#' smooth bacterial growth). Referencing colour to the agar of the same
#' plate also cancels most white-balance differences between photographs.
#'
#' **Class model.** Each class has a bivariate Gaussian for
#' (\eqn{\log \Delta E}, texture) and a mean colour *direction* with an
#' isotropic spread. The direction likelihood of every pixel is widened by
#' its own noise-induced uncertainty \eqn{(\sigma_{noise} / \Delta E)^2}, so
#' faint pixels (colony margins, halo fringe) are not forced into a class by
#' noise (a heteroscedastic, measurement-error-aware mixture).
#'
#' The halo is a *diffusion zone*: pigment concentration, and hence colour
#' magnitude, decays continuously with distance from the bacterium. Its
#' magnitude likelihood is therefore uniform in \eqn{\log \Delta E} up to
#' its core level (with a Gaussian fall-off above), so the faint outer fringe
#' is as plausible a halo pixel as the saturated core.
#'
#' **Unsupervised model (default).** A mixture with one component
#' per class is fitted by EM (Dempster, Laird & Rubin 1977). The EM is
#' *initialised from the plate layout*: pixels near the expected fungal
#' inoculation points seed the fungus class, pixels at the centre seed the
#' bacterium, and pixels in the annulus between them seed the halo. The
#' components therefore keep their biological meaning, and the halo class
#' is dropped automatically when the bacterium produces no halo (its colour
#' direction must differ from that of the fungus by > 12 degrees).
#'
#' **Supervised model (optional).** Pass a classifier from
#' [train_classifier()] (Gaussian / quadratic discriminant model built from
#' pixels you annotated), which is recommended when a new
#' bacterium-fungus pair has unusual colours.
#'
#' **Spatial regularisation.** Class posteriors are smoothed before the
#' maximum-a-posteriori decision (a fast approximation of a Markov random
#' field prior). Each pixel gets a confidence
#' \eqn{c = 1 / (1 + d^2 / 0.01)} from its direction uncertainty
#' \eqn{d^2}; its final posterior is \eqn{c} x (posterior smoothed at
#' `posterior_smooth_mm`) + \eqn{(1 - c)} x (confidence-weighted posterior
#' of its neighbourhood at `context_mm`). Faint pixels such as the fading
#' fringe of a halo thus take the class of the confident tissue they belong
#' to, instead of being mistaken for pale fungal margin.
#' Pixels far from every class (Mahalanobis distance above the
#' \eqn{1 - 10^{-6}} quantile of \eqn{\chi^2_4}) are labelled `"other"`
#' (dust, bubbles, contaminants) and never measured.
#'
#' @param plate A plate processed by [model_background()].
#' @param layout A [plate_layout()].
#' @param classifier Optional `mycohalo_classifier` from [train_classifier()].
#' @param posterior_smooth_mm Smoothing of class posteriors for confident
#'   pixels (mm).
#' @param context_mm Neighbourhood (mm) from which low-confidence (faint)
#'   pixels borrow their class.
#' @param max_em_pixels Maximum number of pixels used to fit the mixture.
#'
#' @return The plate with element `classes`: an integer matrix
#'   (0 agar/outside, 1 fungus, 2 bacteria, 3 halo, 4 other) and the fitted
#'   model.
#' @examples
#' sim <- simulate_plate(seed = 4)
#' p <- read_plate(sim$image) |> detect_plate() |> model_background()
#' p <- classify_pixels(p, plate_layout())
#' table(p$classes$map)
#' @export
classify_pixels <- function(plate, layout = plate_layout(), classifier = NULL,
                            posterior_smooth_mm = 0.15, context_mm = 0.5,
                            max_em_pixels = 60000L) {
  stopifnot(is_mycohalo_plate(plate), !is.null(plate$bg))
  F <- pixel_features(plate)
  obj_idx <- which(plate$bg$object)
  H <- dim(plate$rgb)[1]; W <- dim(plate$rgb)[2]
  map <- matrix(0L, H, W)
  if (length(obj_idx) == 0) {
    plate$classes <- list(map = map, model = NULL)
    return(plate)
  }
  Xo <- F[obj_idx, , drop = FALSE]

  if (!is.null(classifier)) {
    model <- classifier$model
    model <- model[intersect(names(model), c("fungus", "bacteria", "halo"))]
    if (!"fungus" %in% names(model)) cli::cli_abort("Classifier must contain a {.val fungus} class.")
  } else {
    model <- fit_layout_mixture(plate, layout, F, obj_idx, max_em_pixels)
  }

  # posterior for all object pixels
  lp <- sapply(model, function(g) log(g$prior) + gauss_logdens(Xo, g))
  lp <- matrix(lp, ncol = length(model))
  mx <- apply(lp, 1, max)
  post <- exp(lp - mx)
  post <- post / rowSums(post)
  md <- sapply(model, function(g) gauss_logdens(Xo, g, parts = TRUE))
  md <- matrix(md, ncol = length(model))
  outlier <- apply(md, 1, min) > stats::qchisq(1 - 1e-6, 4)

  # confidence-weighted, two-scale spatial smoothing of posteriors.
  # conf is high where the colour direction is well determined (strong
  # contrast) and low for faint pixels (halo fringe, pale margins); faint
  # pixels therefore borrow their class from confident neighbours within
  # context_mm, confident pixels are only lightly smoothed.
  mmpp <- plate$dish$mm_per_px
  conf <- 1 / (1 + Xo[, "d2"] / 0.01)
  objm <- plate$bg$object * 1
  wimg <- matrix(0, H, W); wimg[obj_idx] <- conf
  sig_s <- max(1, posterior_smooth_mm / mmpp)
  sig_l <- max(2, context_mm / mmpp)
  P <- matrix(0, length(obj_idx), length(model))
  for (k in seq_along(model)) {
    img <- matrix(0, H, W)
    img[obj_idx] <- post[, k]
    ps <- gauss_blur(img, sig_s, w = objm)[obj_idx]
    pl <- gauss_blur(img, sig_l, w = wimg)[obj_idx]
    ps[is.na(ps)] <- post[is.na(ps), k]
    pl[is.na(pl)] <- ps[is.na(pl)]
    P[, k] <- conf * ps + (1 - conf) * pl
  }
  P[is.na(P)] <- 0
  cls <- max.col(P, ties.method = "first")
  code <- c(fungus = 1L, bacteria = 2L, halo = 3L)[names(model)]
  lab <- code[cls]
  lab[outlier] <- 4L
  map[obj_idx] <- lab

  plate$classes <- list(map = map, model = model, post = P, obj_idx = obj_idx,
                        post_names = names(model),
                        n_classes = length(model),
                        source = if (is.null(classifier)) "layout-initialised EM" else "supervised")
  plate
}

#' Per-pixel classification features (n_pixels x 5 matrix, all pixels)
#'
#' Colour is described relative to the agar as a difference vector
#' (dL, da, db) and decomposed into its magnitude (log Delta E) and its
#' direction (unit vector uL, ua, ub). A diffusible pigment (halo) at
#' varying concentration moves the colour along one direction from agar
#' towards the pigment colour, so its *direction* is concentration-
#' invariant: the faint outer fringe of a halo points the same way as its
#' saturated core, whereas grey fungal tissue points along the achromatic
#' L axis at any melanization level. Texture (log local SD of L*) is the
#' fifth feature. `d2` is the per-pixel variance of the direction caused by
#' sensor noise, (sigma_noise / Delta E)^2: faint pixels have uncertain
#' direction, and the classifier widens their direction likelihood
#' accordingly (heteroscedastic model), instead of trusting noise.
#' @keywords internal
#' @noRd
pixel_features <- function(plate) {
  S <- if (!is.null(plate$bg$smooth_feat)) plate$bg$smooth_feat else plate$bg$smooth
  mu <- plate$bg$agar_lab
  dL <- as.vector(S$L) - mu[1]
  da <- as.vector(S$a) - mu[2]
  db <- as.vector(S$b) - mu[3]
  E <- pmax(sqrt(dL^2 + da^2 + db^2), 0.5)
  sn <- sqrt(mean(plate$bg$agar_sd^2))       # noise of smoothed colour (per channel)
  cbind(logE = log(E), uL = dL / E, ua = da / E, ub = db / E,
        tex = log1p(as.vector(plate$bg$texture)),
        d2 = (sn / E)^2)
}


#' Class log-density: magnitude/texture Gaussian x noise-aware direction
#'
#' log p(x | k) = log N((logE, tex); mu_k, S_k)
#'              - ||u - m_k||^2 / (2 v) - log(2 pi v),  v = s2_k + d2_i
#' The direction term is an isotropic Gaussian on the tangent plane of the
#' unit sphere (2 degrees of freedom) whose variance adds the class spread
#' s2_k and the pixel's own noise-induced direction variance d2_i.
#' @keywords internal
#' @noRd
gauss_logdens <- function(X, g, parts = FALSE) {
  M <- X[, c("logE", "tex"), drop = FALSE]
  U <- X[, c("uL", "ua", "ub"), drop = FALSE]
  v <- g$s2 + X[, "d2"]
  qd <- rowSums(sweep(U, 2, g$dir)^2) / v
  if (isTRUE(g$gradient)) {
    # diffusion zone: any intensity up to the core level is equally likely
    # (uniform in log dE), Gaussian fall-off above it; texture Gaussian.
    sdE <- sqrt(g$S[1, 1]); sdT <- sqrt(g$S[2, 2])
    zE <- pmax(M[, 1] - g$mu[1], 0) / sdE
    zT <- (M[, 2] - g$mu[2]) / sdT
    qm <- zE^2 + zT^2
    if (parts) return(qm + qd)
    width <- max(g$mu[1] - g$logE_min, 0) + sdE * sqrt(pi / 2)
    return(-log(width) - 0.5 * zE^2 - 0.5 * (zT^2 + log(2 * pi * sdT^2)) -
             0.5 * qd - log(2 * pi * v))
  }
  qm <- stats::mahalanobis(M, g$mu, g$S)
  if (parts) return(qm + qd)
  -0.5 * (qm + g$logdet + 2 * log(2 * pi)) - 0.5 * qd - log(2 * pi * v)
}

#' Weighted fit of a class model
#' @keywords internal
#' @noRd
fit_gauss <- function(X, w = NULL, prior = 1, gradient = FALSE) {
  if (is.null(w)) w <- rep(1, nrow(X))
  sw <- sum(w)
  M <- X[, c("logE", "tex"), drop = FALSE]
  mu <- colSums(M * w) / sw
  S <- crossprod(sweep(M, 2, mu) * sqrt(w)) / sw
  if (gradient) {
    # core level = weighted 75th percentile of log dE; spread above it
    o <- order(M[, 1]); cw <- cumsum(w[o]) / sw
    mu[1] <- M[o, 1][which(cw >= 0.75)[1]]
    up <- M[, 1] > mu[1]
    S[1, 1] <- if (any(up)) sum(w[up] * (M[up, 1] - mu[1])^2) / sum(w[up]) else 0.09
    S[1, 2] <- S[2, 1] <- 0
  }
  diag(S) <- pmax(diag(S), c(0.3, 0.1)^2) + 1e-6
  U <- X[, c("uL", "ua", "ub"), drop = FALSE]
  dir <- colSums(U * w) / sw
  dir <- dir / sqrt(sum(dir^2))
  dev2 <- rowSums(sweep(U, 2, dir)^2)
  s2 <- max(sum(w * pmax(dev2 - 2 * X[, "d2"], 0)) / (2 * sw), 0.03^2)
  list(mu = mu, S = S, logdet = as.numeric(determinant(S, logarithm = TRUE)$modulus),
       dir = dir, s2 = s2, prior = prior, gradient = gradient,
       logE_min = log(3))
}

#' Angle (degrees) between the colour directions of two classes
#' @keywords internal
#' @noRd
dir_angle <- function(g1, g2) acos(min(1, max(-1, sum(g1$dir * g2$dir)))) * 180 / pi

#' Layout-initialised Gaussian mixture (EM)
#' @keywords internal
#' @noRd
fit_layout_mixture <- function(plate, layout, F, obj_idx, max_em_pixels, min_halo_angle = 12) {
  d <- plate$dish
  m <- dish_masks(plate)
  R <- d$r_outer
  exp_pos <- layout_pixels(layout, d)
  xo <- m$x[obj_idx]; yo <- m$y[obj_idx]
  rc <- sqrt((xo - d$x)^2 + (yo - d$y)^2) / R
  dnear <- rep(Inf, length(obj_idx))
  for (i in seq_len(nrow(exp_pos))) {
    dnear <- pmin(dnear, sqrt((xo - exp_pos$x[i])^2 + (yo - exp_pos$y[i])^2) / R)
  }
  X <- F[obj_idx, , drop = FALSE]
  tz <- plate$bg$texture_z[obj_idx]
  # Initial sets combine the layout (where) with texture (what): fungal
  # colonies are rough, a diffusion halo is smooth. This keeps the
  # initialisation right even when inoculation points are a few mm off.
  near_f <- dnear < layout$search_frac & rc > 0.15
  if (layout$centre == "bacteria") {
    ann <- rc > 0.12 & rc < 0.28
    tz_f <- if (any(near_f)) stats::quantile(tz[near_f], 0.5) else Inf
    tz_h <- if (any(ann)) stats::quantile(tz[ann], 0.5) else -Inf
    init <- list(fungus = which(near_f & tz >= tz_f & dnear < 0.25),
                 bacteria = which(rc < 0.07),
                 halo = which(ann & tz <= tz_h & !(near_f & tz >= tz_f)))
  } else {
    init <- list(fungus = which(dnear < 0.2))
    if (layout$centre == "fungus") init$fungus <- c(init$fungus, which(rc < 0.07))
  }
  min_n <- 200L
  init <- init[vapply(init, length, 1L) >= min_n]
  if (!"fungus" %in% names(init)) {
    cli::cli_warn("Few object pixels near the expected fungal positions; check {.fn plate_layout} and {.fn plot_qc}.")
    init$fungus <- seq_along(obj_idx)
  }
  if (length(init) == 1L) {
    # single class: fit it to all object pixels
    return(list(fungus = fit_gauss(X[sample_idx(nrow(X), max_em_pixels), , drop = FALSE])))
  }
  model <- lapply(names(init), function(k) fit_gauss(X[init[[k]], , drop = FALSE],
                                                     prior = 1 / length(init),
                                                     gradient = k == "halo"))
  names(model) <- names(init)
  # Semi-supervised EM: a subsample of all object pixels plus the seed
  # pixels, whose class is clamped. Clamping keeps each component tied to
  # its biological meaning (a free EM can let e.g. the halo component drift
  # onto pale colony margins on real plates).
  s <- sample_idx(nrow(X), max_em_pixels)
  seeds <- lapply(init, function(ix) ix[sample_idx(length(ix), 15000L)])
  Xs <- rbind(X[s, , drop = FALSE], X[unlist(seeds), , drop = FALSE])
  fixed <- c(rep(NA_integer_, length(s)),
             rep(seq_along(seeds), vapply(seeds, length, 1L)))
  model <- run_em(model, Xs, fixed = fixed)
  # drop classes that collapsed to (almost) nothing
  keep <- vapply(model, function(g) g$prior > 0.002, TRUE)
  keep["fungus"] <- TRUE
  # a halo is a chromatic zone: it must differ in colour *direction* from
  # the fungus; otherwise the "halo" component only modelled colony margins
  if ("halo" %in% names(model) && dir_angle(model$fungus, model$halo) < min_halo_angle) {
    keep["halo"] <- FALSE
  }
  if (!all(keep)) {
    model <- model[keep]
    pr <- vapply(model, `[[`, 1, "prior")
    for (k in seq_along(model)) model[[k]]$prior <- pr[k] / sum(pr)
    kept <- which(keep)
    fx <- match(fixed, kept)
    fx[is.na(fx) & !is.na(fixed)] <- NA_integer_
    model <- if (length(model) > 1) run_em(model, Xs, fixed = fx) else list(fungus = fit_gauss(Xs))
  }
  model
}

#' EM iterations for the class mixture
#' @keywords internal
#' @noRd
run_em <- function(model, Xs, max_iter = 30L, fixed = NULL) {
  ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    lp <- sapply(model, function(g) log(g$prior) + gauss_logdens(Xs, g))
    mx <- apply(lp, 1, max)
    rsp <- exp(lp - mx)
    tot <- rowSums(rsp)
    rsp <- rsp / tot
    if (!is.null(fixed)) {
      fr <- which(!is.na(fixed))
      rsp[fr, ] <- 0
      rsp[cbind(fr, fixed[fr])] <- 1
    }
    ll <- sum(mx + log(tot))
    for (k in seq_along(model)) {
      w <- rsp[, k]
      if (sum(w) < 50) next
      model[[k]] <- fit_gauss(Xs, w, prior = mean(w), gradient = isTRUE(model[[k]]$gradient))
    }
    if (abs(ll - ll_old) < 1e-4 * abs(ll)) break
    ll_old <- ll
  }
  model
}

#' Deterministic evenly spaced subsample of indices
#' @keywords internal
#' @noRd
sample_idx <- function(n, k) {
  if (n <= k) return(seq_len(n))
  unique(round(seq(1, n, length.out = k)))
}

#' Extract annotated training pixels
#'
#' Collects classification features of pixels inside circular regions you
#' marked on a processed plate (see [annotate_plate()] for interactive
#' annotation).
#'
#' @param plate A plate processed by [model_background()].
#' @param regions Data frame with columns `class` (`"fungus"`,
#'   `"bacteria"`, `"halo"`), `x`, `y`, `r` in working-image pixels.
#' @return Data frame of features with a `class` column.
#' @examples
#' sim <- simulate_plate(seed = 5)
#' p <- read_plate(sim$image) |> detect_plate() |> model_background()
#' reg <- data.frame(class = c("fungus", "bacteria"),
#'                   x = c(sim$truth$colonies$x[1], p$dish$x),
#'                   y = c(sim$truth$colonies$y[1], p$dish$y), r = 8)
#' head(extract_training_pixels(p, reg))
#' @export
extract_training_pixels <- function(plate, regions) {
  stopifnot(is_mycohalo_plate(plate), !is.null(plate$bg),
            all(c("class", "x", "y", "r") %in% names(regions)))
  F <- pixel_features(plate)
  m <- dish_masks(plate)
  out <- lapply(seq_len(nrow(regions)), function(i) {
    ins <- which((m$x - regions$x[i])^2 + (m$y - regions$y[i])^2 <= regions$r[i]^2)
    data.frame(class = regions$class[i], F[ins, , drop = FALSE], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Train a supervised pixel classifier
#'
#' Fits one multivariate Gaussian per class (quadratic discriminant
#' analysis) to annotated pixels from one or several plates.
#'
#' @param training Data frame from [extract_training_pixels()] (or a list
#'   of them, which are row-bound).
#' @param equal_priors Use equal class priors (recommended, because
#'   annotation effort rather than biology determines class frequencies).
#' @return An object of class `mycohalo_classifier`.
#' @examples
#' sim <- simulate_plate(seed = 5)
#' p <- read_plate(sim$image) |> detect_plate() |> model_background()
#' reg <- data.frame(class = c("fungus", "bacteria"),
#'                   x = c(sim$truth$colonies$x[1], p$dish$x),
#'                   y = c(sim$truth$colonies$y[1], p$dish$y), r = 8)
#' clf <- train_classifier(extract_training_pixels(p, reg))
#' clf
#' @export
train_classifier <- function(training, equal_priors = TRUE) {
  if (!is.data.frame(training)) training <- do.call(rbind, training)
  feats <- c("logE", "uL", "ua", "ub", "tex", "d2")
  cl <- unique(training$class)
  model <- lapply(cl, function(k) {
    X <- as.matrix(training[training$class == k, feats])
    fit_gauss(X, prior = if (equal_priors) 1 / length(cl) else nrow(X) / nrow(training),
              gradient = k == "halo")
  })
  names(model) <- cl
  structure(list(model = model, n = table(training$class)), class = "mycohalo_classifier")
}

#' @export
print.mycohalo_classifier <- function(x, ...) {
  cli::cli_text("{.strong <mycohalo_classifier>} Gaussian (QDA) pixel classifier")
  for (k in names(x$model)) {
    mu <- signif(x$model[[k]]$mu, 3)
    dd <- signif(x$model[[k]]$dir, 2)
    cli::cli_bullets(c("*" = "{k}: n = {x$n[[k]]}; median |dE| = {signif(exp(mu[1]), 3)}, colour direction (L,a,b) = ({dd[1]}, {dd[2]}, {dd[3]})"))
  }
  invisible(x)
}

#' Interactively annotate training regions
#'
#' Displays the plate and lets you click region centres for each class with
#' [graphics::locator()] (finish each class with Esc / right-click). Works
#' in RStudio's plot pane and in regular graphics devices.
#'
#' @param plate A plate processed by [model_background()].
#' @param classes Classes to annotate, in order.
#' @param r_mm Radius of each clicked region (mm).
#' @return Data frame of regions suitable for [extract_training_pixels()].
#' @examples
#' \dontrun{
#' reg <- annotate_plate(p)
#' clf <- train_classifier(extract_training_pixels(p, reg))
#' }
#' @export
annotate_plate <- function(plate, classes = c("fungus", "bacteria", "halo"), r_mm = 0.5) {
  if (!interactive()) cli::cli_abort("{.fn annotate_plate} needs an interactive session.")
  H <- dim(plate$rgb)[1]; W <- dim(plate$rgb)[2]
  r <- r_mm / plate$dish$mm_per_px
  out <- list()
  for (k in classes) {
    show_rgb(plate$rgb, main = sprintf("Click centres of '%s' regions; Esc when done", k))
    pts <- graphics::locator(type = "p", col = "red", pch = 3)
    if (!is.null(pts)) {
      out[[k]] <- data.frame(class = k, x = pts$x, y = pts$y, r = r)
    }
  }
  do.call(rbind, out)
}
