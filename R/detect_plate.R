#' Detect Petri dish
#'
#' Detect the Petri dish region from the image background.
#'
#' @param plate A mycohalo_plate object.
#'
#' @return
#' Updated mycohalo_plate object.
#'
#' @examples
#' \dontrun{
#' plate <- read_plate("plate.jpg")
#' plate <- detect_plate(plate)
#' }
#'
#' @export
#'
#' @importFrom cli cli_abort
detect_plate <- function(plate) {

  if (!inherits(plate, "mycohalo_plate")) {
    cli::cli_abort(
      "{.arg plate} must be a mycohalo_plate object."
    )
  }

  gray <- EBImage::channel(
    plate$ebimage,
    mode = "gray"
  )

  plate$gray <- gray

  plate
}
