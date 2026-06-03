#' Convert image to grayscale
#'
#' Convert an EBImage image to grayscale.
#'
#' @param plate A mycohalo_plate object.
#'
#' @return
#' Updated mycohalo_plate object.
#'
#' @export
#'
#' @importFrom cli cli_abort
to_grayscale <- function(plate) {

  if (!inherits(plate, "mycohalo_plate")) {
    cli::cli_abort(
      "{.arg plate} must be a mycohalo_plate object."
    )
  }

  plate$gray <- EBImage::channel(
    plate$ebimage,
    mode = "gray"
  )

  plate
}
