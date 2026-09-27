#' Simulate a plate photograph with known ground truth
#'
#' Generates a synthetic confrontation (or control) plate that reproduces
#' the difficulties of real photographs: grey, wrinkled fungal colonies with
#' irregular lobed outlines, a whitish central bacterium with a diffuse
#' yellow halo that overlaps the colonies in brightness, satellite
#' micro-colonies, a bright dish wall, vignetting, exposure errors and
#' sensor noise. Because every pixel's true class and colour are known, the
#' simulator is used to validate MycoHalo's accuracy (see the unit tests and
#' the vignette) and to test new settings before applying them to real data.
#'
#' @param width,height Image size (px).
#' @param dish_r_frac Dish outer radius as a fraction of the image height.
#' @param dish_diameter_mm Physical dish diameter (sets the true scale).
#' @param layout A [plate_layout()] defining colony positions.
#' @param colony_radius_mm Mean colony radius (mm), recycled per colony.
#' @param colony_L Mean \eqn{L^*} of each colony (the melanization truth;
#'   lower = darker), recycled.
#' @param colony_ab \eqn{a^*, b^*} of the colonies.
#' @param edge_lightening \eqn{L^*} units by which the outer 1.5 mm margin
#'   is lighter than the interior (age gradient).
#' @param facing_darkening \eqn{L^*} units by which the colony half facing
#'   the centre is darker (simulated localized melanization response).
#' @param texture_sd SD (\eqn{L^*}) of the wrinkle texture.
#' @param bacteria Include a central bacterium (default: from `layout`).
#' @param bacteria_radius_mm,halo_radius_mm Radii of bacterium and halo.
#' @param bacteria_lab,halo_lab CIELAB of bacterium and of the halo core.
#' @param agar_lab,background_lab CIELAB of agar and photographic background.
#' @param satellite Add a small satellite colony next to the first colony.
#' @param vignetting Relative light fall-off at the dish edge (0 = none).
#' @param exposure Multiplicative exposure error in linear light.
#' @param noise_sd Gaussian sensor noise (sRGB units, 0-1).
#' @param grey_card Add an 18 % grey card (L* = 50) in the lower-left corner.
#' @param seed Random seed.
#'
#' @return A list with `image` (H x W x 3 sRGB array) and `truth` (dish
#'   circle, `mm_per_px`, per-colony true area and mean \eqn{L^*} before
#'   illumination effects, label matrix, bacteria and halo masks, grey-card
#'   region).
#' @examples
#' sim <- simulate_plate(seed = 1)
#' show <- function(a) { plot(as.raster(a)) }
#' show(sim$image)
#' sim$truth$colonies
#' @export
simulate_plate <- function(width = 900, height = 1200, dish_r_frac = 0.47,
                           dish_diameter_mm = 90, layout = plate_layout(),
                           colony_radius_mm = 7, colony_L = c(52, 46, 40, 58),
                           colony_ab = c(0.5, -1), edge_lightening = 4,
                           facing_darkening = 0, texture_sd = 5,
                           bacteria = NULL, bacteria_radius_mm = 6, halo_radius_mm = 17,
                           bacteria_lab = c(86, 0, 9), halo_lab = c(42, -2, 22),
                           agar_lab = c(22, -1, -4), background_lab = c(8, 0, 1),
                           satellite = TRUE, vignetting = 0.25, exposure = 1,
                           noise_sd = 0.012, grey_card = FALSE, seed = 1) {
  set.seed(seed)
  if (is.null(bacteria)) bacteria <- layout$centre == "bacteria"
  H <- height; W <- width
  R <- dish_r_frac * H
  cx <- W / 2 + 0.5; cy <- H / 2 + 0.5
  mmpp <- dish_diameter_mm / (2 * R)
  xx <- matrix(rep(seq_len(W), each = H), H, W)
  yy <- matrix(rep(seq_len(H), W), H, W)
  rr <- sqrt((xx - cx)^2 + (yy - cy)^2)

  Lm <- matrix(background_lab[1], H, W)
  Am <- matrix(background_lab[2], H, W)
  Bm <- matrix(background_lab[3], H, W)
  set_px <- function(mask, lab) {
    Lm[mask] <<- lab[1]; Am[mask] <<- lab[2]; Bm[mask] <<- lab[3]
  }
  # dish wall (bright ring) and agar
  set_px(rr <= R, c(62, -1, -3))
  set_px(rr <= 0.965 * R, agar_lab)

  smooth_noise <- function(sc) {
    z <- matrix(stats::rnorm(H * W), H, W)
    z <- gauss_blur(z, sc)
    (z - mean(z)) / stats::sd(as.vector(z))
  }

  # halo (smooth radial blend towards halo colour) + bacterium
  bact <- matrix(FALSE, H, W); halo_true <- bact
  if (bacteria) {
    rb <- bacteria_radius_mm / mmpp
    rh <- halo_radius_mm / mmpp
    th <- atan2(yy - cy, xx - cx)
    wobble <- 1 + 0.06 * cos(4 * th + 0.3) + 0.04 * cos(7 * th + 1)
    w <- 1 / (1 + exp((rr / wobble - rh) / (0.04 * rh)))
    inside <- rr <= 0.965 * R
    Lm[inside] <- (1 - w[inside]) * Lm[inside] + w[inside] * halo_lab[1]
    Am[inside] <- (1 - w[inside]) * Am[inside] + w[inside] * halo_lab[2]
    Bm[inside] <- (1 - w[inside]) * Bm[inside] + w[inside] * halo_lab[3]
    halo_true <- w > 0.5 & inside
    rbe <- rb * (1 + 0.08 * cos(5 * th + 2) + 0.05 * cos(9 * th))
    bact <- rr <= rbe
    tb <- smooth_noise(2) * 3
    Lm[bact] <- bacteria_lab[1] + tb[bact] - 6 * (rr[bact] / rb)^2
    Am[bact] <- bacteria_lab[2]; Bm[bact] <- bacteria_lab[3] + 4 * (rr[bact] / rb)
    halo_true <- halo_true & !bact
  }

  # fungal colonies
  pos <- layout_pixels(layout, list(x = cx, y = cy, r_outer = R, mm_per_px = mmpp))
  n <- nrow(pos)
  crad <- rep_len(colony_radius_mm, n) / mmpp
  cL <- rep_len(colony_L, n)
  labels <- matrix(0L, H, W)
  tex <- smooth_noise(1.2) * 0.6 + smooth_noise(3) * 0.8
  tex <- tex / stats::sd(as.vector(tex)) * texture_sd
  trueL <- numeric(n); area <- numeric(n)
  for (i in seq_len(n)) {
    th <- atan2(yy - pos$y[i], xx - pos$x[i])
    ri <- sqrt((xx - pos$x[i])^2 + (yy - pos$y[i])^2)
    k <- 3:9
    amp <- stats::runif(length(k), 0.01, 0.05)
    ph <- stats::runif(length(k), 0, 2 * pi)
    edge <- crad[i] * (1 + Reduce(`+`, lapply(seq_along(k), function(j) amp[j] * cos(k[j] * th + ph[j]))))
    cm <- ri <= edge & rr <= 0.95 * R
    labels[cm] <- i
    depth_mm <- (edge - ri) * mmpp
    Lc <- cL[i] + tex + edge_lightening * pmax(0, 1 - depth_mm / 1.5)
    if (facing_darkening != 0) {
      toward <- atan2(cy - pos$y[i], cx - pos$x[i])
      dang <- abs(atan2(sin(th - toward), cos(th - toward)))
      Lc <- Lc - facing_darkening * (dang < pi / 2)
    }
    Lm[cm] <- Lc[cm]; Am[cm] <- colony_ab[1]; Bm[cm] <- colony_ab[2]
    trueL[i] <- mean(Lc[cm]); area[i] <- sum(cm) * mmpp^2
  }
  sat_mask <- matrix(FALSE, H, W)
  if (satellite && n > 0) {
    ang <- atan2(pos$y[1] - cy, pos$x[1] - cx) + 0.9
    sx <- pos$x[1] + (crad[1] + 3.2 / mmpp) * cos(ang)
    sy <- pos$y[1] + (crad[1] + 3.2 / mmpp) * sin(ang)
    sat_mask <- (xx - sx)^2 + (yy - sy)^2 <= (1.3 / mmpp)^2 & labels == 0L
    Lm[sat_mask] <- cL[1] + tex[sat_mask]; Am[sat_mask] <- colony_ab[1]; Bm[sat_mask] <- colony_ab[2]
  }
  if (bacteria) {
    halo_true <- halo_true & labels == 0L & !sat_mask
    bact <- bact & labels == 0L
  }
  card <- NULL
  if (grey_card) {
    card <- c(0.1, 0.92, 0.05)
    cm <- (xx - card[1] * W)^2 + (yy - card[2] * H)^2 <= (card[3] * min(H, W) * 1.3)^2
    set_px(cm, c(50, 0, 0))
  }

  lab <- array(c(Lm, Am, Bm), c(H, W, 3))
  rgb <- lab_to_srgb(lab)
  lin <- srgb_linearize(matrix(rgb, ncol = 3))
  vig <- (1 - vignetting * pmin(rr / R, 1.2)^2) * exposure
  lin <- lin * as.vector(vig)
  rgb <- srgb_compand(pmin(pmax(lin, 0), 1))
  rgb <- rgb + stats::rnorm(length(rgb), 0, noise_sd)
  rgb <- array(pmin(pmax(rgb, 0), 1), c(H, W, 3))

  list(image = rgb,
       truth = list(
         dish = c(x = cx, y = cy, r = R), mm_per_px = mmpp,
         colonies = data.frame(colony_id = pos$id, x = pos$x, y = pos$y,
                               area_mm2 = area, L_mean = trueL, MI_mean = 100 - trueL),
         labels = labels, bacteria = bact, halo = halo_true, satellite = sat_mask,
         grey_card = card, vignetting = vignetting, exposure = exposure))
}
