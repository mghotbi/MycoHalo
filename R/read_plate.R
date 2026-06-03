#' Read a plate image
#'
#' Import a fungal colony plate image and create a
#' MycoHalo plate object.
#'
#' @param path Character scalar.
#' Path to an image file.
#'
#' @return
#' A `mycohalo_plate` object.
#'
#' @examples
#' \dontrun{
#' plate <- read_plate(
#'   "inst/extdata/control/control_01.jpg"
#' )
#' }
#'
#' @export
#'
#' @importFrom magick image_read
#' @importFrom cli cli_abort
read_plate <- function(path) {

  if (!is.character(path)) {
    cli::cli_abort(
      c(
        "{.arg path} must be a character vector.",
        "x" = "Received class: {.cls {class(path)}}."
      )
    )
  }

  if (length(path) != 1) {
    cli::cli_abort(
      c(
        "{.arg path} must contain exactly one file path.",
        "x" = "Length is {length(path)}."
      )
    )
  }

  if (!file.exists(path)) {
    cli::cli_abort(
      c(
        "Image file does not exist.",
        "x" = "Cannot find file {.file {path}}."
      )
    )
  }

  image <- magick::image_read(path)

  create_plate(image)
}
