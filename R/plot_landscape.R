#' Melanization landscape: 3D shaded surface and 2D map of a plate
#'
#' Turns a plate into a "melanization landscape". In 3D, each colony rises
#' from the agar with a height and colour given by its local melanization
#' index (MI = 100 - L*), so darker tissue forms higher, browner relief. The
#' bacterium and its halo are painted on the floor. In 2D, the same MI field
#' is drawn over a faded photograph, with the halo and bacterium outlined.
#' An arrow on each colony points to the bacterium and is annotated with
#' the side-specific response `delta_MI_facing`.
#'
#' @details
#' The MI field is smoothed with a Gaussian of `smooth_mm` to show tissue-
#' scale structure (margins, sectors, facing-side darkening) rather than
#' single hyphal ridges. In 3D the colony edges are tapered over
#' `taper_mm`, so colonies look like rounded relief instead of cliffs. The
#' colour scale is a single-hue sepia ramp from light to dark, because MI is
#' a magnitude and melanin is brown-black. Heights and colours share the same
#' scale, set by the observed MI range (or `mi_range`).
#'
#' @param x A `mycohalo_result` computed with `keep_images = TRUE`.
#' @param type `"both"` (3D + 2D side by side), `"3d"` or `"2d"`.
#' @param mi_range Optional common MI range `c(min, max)` for colours and
#'   heights. Fix it when comparing several plates.
#' @param smooth_mm Smoothing of the MI field (mm).
#' @param relief_smooth_mm Smoothing (mm) of the MI field used for the 3D
#'   heights (colours use `smooth_mm`).
#' @param taper_mm Width of the rounded colony edge in 3D (mm).
#' @param grid Approximate number of grid cells along the longer side of
#'   the 3D surface (speed vs. detail).
#' @param theta,phi Viewing angles of the 3D view (see [graphics::persp()]).
#' @param expand Vertical exaggeration of the 3D relief.
#' @param main Title.
#' @return Invisibly, a list with the MI matrix used and the colour ramp.
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 42,
#'                      facing_darkening = 5, edge_lightening = 8)
#' res <- analyze_plate(sim$image, verbose = FALSE)
#' plot_melanization_map(res)
#' @export
plot_melanization_map <- function(x, type = c("both", "3d", "2d"), mi_range = NULL,
                                  smooth_mm = 0.25, relief_smooth_mm = 1.5,
                                  taper_mm = 2.5, grid = 180,
                                  theta = -30, phi = 32, expand = 0.3,
                                  main = NULL) {
  type <- match.arg(type)
  if (!inherits(x, "mycohalo_result") || is.null(x$plate)) {
    cli::cli_abort("{.arg x} must be a {.cls mycohalo_result} computed with {.code keep_images = TRUE}.")
  }
  p <- x$plate
  ob <- p$objects
  d <- p$dish
  mmpp <- d$mm_per_px
  H <- dim(p$rgb)[1]; W <- dim(p$rgb)[2]
  col_mask <- ob$labels > 0
  if (!any(col_mask)) cli::cli_abort("No colonies to draw.")

  # ---- MI field, smoothed inside colonies only (normalised convolution)
  MI <- 100 - p$lab[, , 1]
  MIs <- gauss_blur(MI, max(1, smooth_mm / mmpp), w = col_mask * 1)
  MIs[!col_mask] <- NA
  if (is.null(mi_range)) {
    mi_range <- stats::quantile(MIs[col_mask], c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    mi_range <- c(floor(mi_range[1]), ceiling(mi_range[2]))
  }
  ramp <- melanin_ramp(100)
  to_col <- function(v) {
    k <- round(1 + 99 * pmin(pmax((v - mi_range[1]) / diff(mi_range), 0), 1))
    ramp[k]
  }

  if (is.null(main)) main <- paste("Melanization landscape -", p$id)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (type == "both") {
    graphics::layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.09))
  }
  # crop to the analysed dish for both views
  r <- d$r_agar
  rows <- max(1, floor(d$y - r)):min(H, ceiling(d$y + r))
  cols <- max(1, floor(d$x - r)):min(W, ceiling(d$x + r))

  if (type %in% c("both", "3d")) {
    graphics::par(mar = c(0, 0, 2.5, 0))
    step <- max(1L, round(max(length(rows), length(cols)) / grid))
    rs <- rev(rows[seq(1, length(rows), by = step)])     # bottom-to-top so y (mm) increases
    cs <- cols[seq(1, length(cols), by = step)]
    dep <- edt_local(col_mask) * mmpp
    taper <- pmin(dep / taper_mm, 1)
    taper <- taper^2 * (3 - 2 * taper)                 # smoothstep: rounded edges
    # relief from a more strongly smoothed field (tissue scale, not ridges)
    MIh <- gauss_blur(MI, max(1, relief_smooth_mm / mmpp), w = col_mask * 1)
    h <- (MIh - mi_range[1] + 0.15 * diff(mi_range))
    h[!col_mask] <- 0
    h <- pmax(h, 0) * taper
    m <- dish_masks(p)
    inside <- m$rr <= d$r_agar
    z <- t(h[rs, cs])                                  # persp: x = columns, y = rows
    z[!t(inside[rs, cs])] <- NA
    xs <- (cs - d$x) * mmpp; ys <- -(rs - d$y) * mmpp  # mm, y up
    # facet colours (nx-1) x (ny-1) from the lower-left vertex of each facet
    fc <- function(M) t(M[rs, cs])[-length(cs), -length(rs), drop = FALSE]
    floor_col <- matrix("#E6E8EB", nrow(z) - 1, ncol(z) - 1)
    floor_col[fc(ob$halo)] <- "#E9CF7A"
    floor_col[fc(ob$bacteria)] <- "#FFF6DC"
    floor_col[fc(ob$satellites)] <- "#B9B2A8"
    mif <- fc(MIs)
    colony_f <- fc(col_mask)
    floor_col[colony_f] <- to_col(mif[colony_f])
    graphics::persp(xs, ys, z, theta = theta, phi = phi, expand = expand,
                    col = floor_col, border = NA, shade = 0.55, ltheta = -110, lphi = 40,
                    box = FALSE, scale = TRUE, r = 4)
    graphics::title(main, cex.main = 1.1)
  }

  if (type %in% c("both", "2d")) {
    graphics::par(mar = c(0.5, 0.5, 2.5, 0.5))
    base <- p$rgb[rows, cols, , drop = FALSE]
    g <- (base[, , 1] + base[, , 2] + base[, , 3]) / 3
    g <- 0.55 + 0.45 * (g - min(g)) / max(1e-9, diff(range(g)))   # faded photo
    img <- array(c(g, g, g), c(length(rows), length(cols), 3))
    cm <- col_mask[rows, cols]
    rgbc <- grDevices::col2rgb(to_col(MIs[rows, cols][cm])) / 255
    for (k in 1:3) { ch <- img[, , k]; ch[cm] <- rgbc[k, ]; img[, , k] <- ch }
    paint_outline <- function(img, mask, colv, thick = 1) {
      b <- mask_boundary(mask)
      if (thick > 1) b <- dilate_disk(b, thick - 1)
      v <- grDevices::col2rgb(colv) / 255
      for (k in 1:3) { ch <- img[, , k]; ch[b] <- v[k]; img[, , k] <- ch }
      img
    }
    th <- max(1, round(length(rows) / 450))
    if (any(ob$halo)) img <- paint_outline(img, (ob$halo | ob$bacteria)[rows, cols], "#C9971C", th)
    if (any(ob$bacteria)) img <- paint_outline(img, ob$bacteria[rows, cols], "#7A6A3A", th)
    img <- paint_outline(img, cm, "#2B2118", th)
    show_rgb(img, main = if (type == "2d") main else "Melanization index map")
    cl <- x$colonies
    ct <- ob$colony_table
    mm <- dish_masks(p)
    if (any(ob$bacteria)) {
      bw <- which(ob$bacteria)
      tx <- mean(mm$x[bw]); ty <- mean(mm$y[bw])
    } else {
      tx <- d$x; ty <- d$y
    }
    for (i in which(ct$detected)) {
      cx <- cl$centroid_x_px[i] - cols[1] + 1
      cy <- cl$centroid_y_px[i] - rows[1] + 1
      vx <- tx - cl$centroid_x_px[i]; vy <- ty - cl$centroid_y_px[i]
      nv <- sqrt(vx^2 + vy^2)
      L <- 0.6 * cl$eq_diameter_mm[i] / 2 / mmpp
      graphics::arrows(cx, cy, cx + vx / nv * L, cy + vy / nv * L, length = 0.08,
                       lwd = 2, col = "white")
      # label block placed outside the colony (above top colonies, below
      # bottom ones), using the colony's actual vertical extent
      wy <- which(ob$labels == i) %% H
      wy[wy == 0] <- H
      top <- min(wy) - rows[1] + 1; bot <- max(wy) - rows[1] + 1
      dl <- 0.032 * length(rows)
      if (cy < length(rows) / 2) { y1 <- top - 1.9 * dl; y2 <- top - 0.9 * dl } else { y1 <- bot + 0.9 * dl; y2 <- bot + 1.9 * dl }
      graphics::text(cx, y1, sprintf("%s   MI %.1f", ct$id[i], cl$MI_mean[i]),
                     cex = 0.8, col = "#1F1A17", font = 2)
      dv <- sprintf("%+.1f", cl$delta_MI_facing[i])
      graphics::text(cx, y2, bquote(Delta * MI[facing] == .(dv)),
                     cex = 0.8, col = "#1F1A17")
    }
  }

  if (type == "both") {
    graphics::par(mar = c(1.6, 18, 0.9, 18))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    graphics::rasterImage(grDevices::as.raster(matrix(ramp, nrow = 1)), 0, 0, 1, 1,
                          interpolate = TRUE)
    at <- pretty(mi_range, 5); at <- at[at >= mi_range[1] & at <= mi_range[2]]
    graphics::axis(1, at = (at - mi_range[1]) / diff(mi_range), labels = at,
                   tick = FALSE, line = -0.8, cex.axis = 0.8, col.axis = "#3A3530")
    graphics::mtext(expression("Melanization index (100 - " * L^"*" * ");   3D floor: halo = gold, bacterium = cream, agar = grey"),
                    side = 3, line = 0.1, cex = 0.7, col = "#3A3530")
  }
  invisible(list(MI = MIs, mi_range = mi_range, ramp = ramp))
}

#' Single-hue sepia ramp, light to dark (perceptually ordered in HCL)
#' @keywords internal
#' @noRd
melanin_ramp <- function(n) {
  t <- seq(0, 1, length.out = n)
  grDevices::hcl(h = 55 - 20 * t, c = 12 + 38 * sin(pi * pmin(t, 0.85) / 1.7),
                 l = 94 - 80 * t, fixup = TRUE)
}
