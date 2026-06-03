#' Print a MycoHalo plate
#'
#' @param x A mycohalo_plate object.
#' @param ... Unused.
#'
#' @return Invisibly returns x.
#'
#' @export
print.mycohalo_plate <- function(x, ...) {

  cat("\n<MycoHalo Plate>\n\n")

  cat(
    "Dimensions:",
    x$metadata$width,
    "x",
    x$metadata$height,
    "\n"
  )

  cat(
    "Plate detected:",
    !is.null(x$plate),
    "\n"
  )

  cat(
    "Colonies detected:",
    !is.null(x$colonies),
    "\n"
  )

  cat(
    "Melanization computed:",
    !is.null(x$melanization),
    "\n"
  )

  invisible(x)
}
