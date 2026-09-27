#' Quality-control figure
#'
#' Always look at the QC figure before trusting numbers. The left panel
#' shows the (flat-field corrected) photograph with the detected dish, the
#' analysed agar disk, expected inoculation points and search radii, colony
#' outlines, the bacterium and its halo; the right panel shows the pixel
#' classification.
#'
#' @param x A `mycohalo_result` (computed with `keep_images = TRUE`).
#' @param panels Which panels to draw.
#' @param colony_cols Colours for the colonies.
#' @return Invisibly `x`.
#' @examples
#' sim <- simulate_plate(seed = 11)
#' res <- analyze_plate(sim$image, verbose = FALSE)
#' plot_qc(res)
#' @export
plot_qc <- function(x, panels = c("overlay", "classes"),
                    colony_cols = c("#00E5FF", "#FF4081", "#76FF03", "#FFD740",
                                    "#E040FB", "#FF6E40", "#40C4FF", "#B2FF59")) {
  if (!inherits(x, "mycohalo_result") || is.null(x$plate)) {
    cli::cli_abort("{.arg x} must be a {.cls mycohalo_result} computed with {.code keep_images = TRUE}.")
  }
  p <- x$plate
  ob <- p$objects
  d <- p$dish
  H <- dim(p$rgb)[1]; W <- dim(p$rgb)[2]
  ncol_ <- length(panels)
  op <- graphics::par(mfrow = c(1, ncol_), mar = c(0.5, 0.5, 2.5, 0.5), bg = "white")
  on.exit(graphics::par(op), add = TRUE)
  ct <- ob$colony_table
  col_rgb <- grDevices::col2rgb(colony_cols) / 255
  thick <- max(1, round(min(H, W) / 700))

  paint <- function(img, mask, rgbv) {
    for (k in 1:3) {
      ch <- img[, , k]; ch[mask] <- rgbv[k]; img[, , k] <- ch
    }
    img
  }
  outline <- function(mask) {
    b <- mask_boundary(mask)
    if (thick > 1) b <- dilate_disk(b, thick - 1)
    b
  }

  if ("overlay" %in% panels) {
    img <- p$rgb
    for (i in which(ct$detected)) {
      img <- paint(img, outline(ob$labels == i), col_rgb[, (i - 1) %% ncol(col_rgb) + 1])
    }
    if (any(ob$satellites)) img <- paint(img, outline(ob$satellites), c(1, 0, 1))
    if (any(ob$halo)) img <- paint(img, outline(ob$halo | ob$bacteria), c(1, 0.6, 0))
    if (any(ob$bacteria)) img <- paint(img, outline(ob$bacteria), c(1, 1, 0.2))
    show_rgb(img, main = sprintf("%s  |  %.4f mm/px", p$id, d$mm_per_px))
    th <- seq(0, 2 * pi, length.out = 361)
    graphics::lines(d$x + d$r_outer * cos(th), d$y + d$r_outer * sin(th), col = "#00B8D4", lty = 2, lwd = 1.5)
    graphics::lines(d$x + d$r_agar * cos(th), d$y + d$r_agar * sin(th), col = "#00B8D4", lwd = 1)
    lay <- x$params$layout
    sr <- lay$search_frac * d$r_outer
    for (i in seq_len(nrow(ct))) {
      graphics::points(ct$expected_x[i], ct$expected_y[i], pch = 3, col = "white", cex = 1.2, lwd = 2)
      graphics::lines(ct$expected_x[i] + sr * cos(th), ct$expected_y[i] + sr * sin(th),
                      col = "grey70", lty = 3)
    }
    cl <- x$colonies
    for (i in seq_len(nrow(ct))) {
      if (!ct$detected[i]) {
        graphics::text(ct$expected_x[i], ct$expected_y[i], paste(ct$id[i], "\nnot detected"),
                       col = "red", font = 2, cex = 0.9)
        next
      }
      lbl <- sprintf("%s\nMI %.1f", ct$id[i], cl$MI_mean[i])
      graphics::text(cl$centroid_x_px[i], cl$centroid_y_px[i], lbl, col = "black", font = 2, cex = 0.95)
      graphics::text(cl$centroid_x_px[i] - 1.5, cl$centroid_y_px[i] - 1.5, lbl, col = "white", font = 2, cex = 0.95)
    }
    graphics::legend("bottomleft", bg = grDevices::adjustcolor("white", 0.8), box.col = NA, cex = 0.75,
                     legend = c("dish edge", "analysed agar", "expected position / search",
                                "bacterium", "halo", "satellite"),
                     col = c("#00B8D4", "#00B8D4", "grey70", "#FFFF33", "#FF9900", "magenta"),
                     lty = c(2, 1, 3, 1, 1, 1), lwd = 2)
  }
  if ("classes" %in% panels) {
    img <- array(0, c(H, W, 3))
    m <- dish_masks(p)
    cmap <- if (!is.null(p$classes$map_final)) p$classes$map_final else p$classes$map
    img <- paint(img, m$agar, c(0.25, 0.27, 0.30))
    img <- paint(img, cmap == 1L, c(0.55, 0.55, 0.55))
    img <- paint(img, cmap == 2L, c(0.98, 0.95, 0.85))
    img <- paint(img, cmap == 3L, c(0.85, 0.70, 0.15))
    img <- paint(img, cmap == 4L, c(0.9, 0.1, 0.1))
    if (!is.null(p$bg$shadow)) img <- paint(img, p$bg$shadow, c(0.10, 0.16, 0.45))
    for (i in which(ct$detected)) {
      img <- paint(img, ob$labels == i, col_rgb[, (i - 1) %% ncol(col_rgb) + 1] * 0.8)
    }
    show_rgb(img, main = sprintf("Pixel classes (%s)", p$classes$source))
    graphics::legend("bottomleft", bg = grDevices::adjustcolor("white", 0.8), box.col = NA, cex = 0.75,
                     legend = c("agar", "fungus (unassigned)", "bacteria", "halo", "other", "cast shadow"),
                     fill = c("#40454D", "#8C8C8C", "#FAF2D9", "#D9B326", "#E61A1A", "#1A2973"))
  }
  invisible(x)
}

#' Display an RGB array with image-row y coordinates (row 1 at top)
#' @keywords internal
#' @noRd
show_rgb <- function(rgb, main = "") {
  H <- dim(rgb)[1]; W <- dim(rgb)[2]
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, W + 0.5), ylim = c(H + 0.5, 0.5), asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(grDevices::as.raster(pmin(pmax(rgb, 0), 1)), 0.5, H + 0.5, W + 0.5, 0.5,
                        interpolate = FALSE)
  graphics::title(main)
}

#' Plot edge-to-centre melanization profiles
#'
#' @param x A `mycohalo_result`, or a profile data frame (e.g. `out$profiles`
#'   from [analyze_plates()]).
#' @param colour_by Column used for colour (default `colony_id`).
#' @param facet_by Optional column for facets (e.g. `plate_id`, treatment).
#' @return A ggplot object (requires \pkg{ggplot2}).
#' @examples
#' sim <- simulate_plate(seed = 12, edge_lightening = 10)
#' res <- analyze_plate(sim$image, verbose = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_profiles(res)
#' @export
plot_profiles <- function(x, colour_by = "colony_id", facet_by = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort("Install {.pkg ggplot2}.")
  df <- if (inherits(x, "mycohalo_result")) x$profiles else x
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mid_mm, y = .data$MI_mean,
                                        colour = .data[[colour_by]], fill = .data[[colour_by]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$MI_lo, ymax = .data$MI_hi),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(x = "Distance from colony edge (mm)",
                  y = expression("Melanization index (100 - " * L^"*" * ")"),
                  colour = colour_by, fill = colour_by) +
    ggplot2::theme_bw(base_size = 12)
  if (!is.null(facet_by)) g <- g + ggplot2::facet_wrap(facet_by)
  g
}
