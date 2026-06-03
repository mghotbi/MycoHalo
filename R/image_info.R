#' Extract image metadata
#'
#' Retrieve metadata from a plate image.
#'
#' @param image A magick-image object.
#'
#' @return
#' A tibble containing image metadata.
#'
#' @examples
#' \dontrun{
#'
#' img <- read_plate(
#'   "inst/extdata/control_01.jpg"
#' )
#'
#' image_info(img)
#'
#' }
#'
#' @export
#'
#' @importFrom magick image_info
#' @importFrom cli cli_abort
image_info <- function(image) {

  if (!inherits(image, "magick-image")) {
    cli::cli_abort(
      c(
        "{.arg image} must be a magick-image object.",
        "x" = "Received class: {.cls {class(image)}}."
      )
    )
  }

  magick::image_info(image)
}
