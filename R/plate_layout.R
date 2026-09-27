#' Describe the inoculation layout of a plate
#'
#' MycoHalo uses the *known* inoculation geometry as a spatial prior. This
#' is what makes it possible to separate several cultures that share one
#' dish, including colonies that touch each other or grow into the halo of
#' the central bacterium.
#'
#' @details
#' Positions are given relative to the dish centre in units of the dish
#' outer radius, with x increasing to the right and **y increasing
#' downwards** (image convention). The default reproduces the standard
#' four-quadrant design: four fungal colonies on the diagonals at 45 % of the
#' radius (TL = top-left, TR, BR, BL) and, for confrontations, a bacterium
#' in the centre.
#'
#' Always photograph plates in the same orientation (e.g. a mark on the dish
#' rim at 12 o'clock) so that colony IDs are reproducible; otherwise use
#' `rotation_deg`.
#'
#' @param n_colonies Number of fungal inoculation points placed evenly on a
#'   circle (ignored if `positions` is given).
#' @param radius_frac Distance of the inoculation points from the dish
#'   centre, as a fraction of the outer dish radius.
#' @param radius_mm Alternatively, the distance in mm (overrides
#'   `radius_frac` once the dish scale is known).
#' @param start_deg Angle of the first colony, degrees clockwise from
#'   3 o'clock. The default -135 puts the first colony at top-left.
#' @param rotation_deg Extra clockwise rotation of the whole layout.
#' @param centre What is inoculated at the centre: `"bacteria"`, `"none"`
#'   or `"fungus"`.
#' @param ids Colony identifiers.
#' @param positions Optional data frame with columns `id`, `x_rel`, `y_rel`
#'   for arbitrary designs.
#' @param search_frac Radius (fraction of the dish radius) around each
#'   expected position within which the colony core is searched. Colonies
#'   are matched to expected positions closest-first, one colony per
#'   position, so search areas may overlap; inoculation points only need to
#'   be approximately where the layout says.
#'
#' @return An object of class `mycohalo_layout`.
#'
#' @examples
#' plate_layout()                      # 4 Zymoseptoria + central bacterium
#' plate_layout(centre = "none")       # control plate
#' plate_layout(n_colonies = 3, ids = c("A", "B", "C"))
#' @export
plate_layout <- function(n_colonies = 4L, radius_frac = 0.45, radius_mm = NULL,
                         start_deg = -135, rotation_deg = 0,
                         centre = c("bacteria", "none", "fungus"),
                         ids = NULL, positions = NULL, search_frac = 0.4) {
  centre <- match.arg(centre)
  if (is.null(positions)) {
    ang <- (start_deg + rotation_deg + (seq_len(n_colonies) - 1) * 360 / n_colonies) * pi / 180
    if (is.null(ids)) {
      ids <- if (n_colonies == 4L && start_deg == -135) c("TL", "TR", "BR", "BL")
             else sprintf("C%d", seq_len(n_colonies))
    }
    positions <- data.frame(id = ids, x_rel = radius_frac * cos(ang),
                            y_rel = radius_frac * sin(ang), stringsAsFactors = FALSE)
  } else {
    stopifnot(all(c("id", "x_rel", "y_rel") %in% names(positions)))
    positions <- as.data.frame(positions)
    if (rotation_deg != 0) {
      a <- rotation_deg * pi / 180
      xr <- positions$x_rel * cos(a) - positions$y_rel * sin(a)
      yr <- positions$x_rel * sin(a) + positions$y_rel * cos(a)
      positions$x_rel <- xr; positions$y_rel <- yr
    }
  }
  structure(list(positions = positions, centre = centre, radius_mm = radius_mm,
                 search_frac = search_frac),
            class = "mycohalo_layout")
}

#' @export
print.mycohalo_layout <- function(x, ...) {
  cli::cli_text("{.strong <mycohalo_layout>} {nrow(x$positions)} fungal colonies; centre: {x$centre}")
  print(x$positions, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Expected colony positions in working-image pixels
#' @keywords internal
#' @noRd
layout_pixels <- function(layout, dish) {
  pos <- layout$positions
  scale <- rep(1, nrow(pos))
  if (!is.null(layout$radius_mm)) {
    rf <- sqrt(pos$x_rel^2 + pos$y_rel^2)
    target <- (layout$radius_mm / dish$mm_per_px) / dish$r_outer
    scale <- ifelse(rf > 0, target / rf, 1)
  }
  data.frame(id = pos$id,
             x = dish$x + pos$x_rel * scale * dish$r_outer,
             y = dish$y + pos$y_rel * scale * dish$r_outer,
             stringsAsFactors = FALSE)
}
