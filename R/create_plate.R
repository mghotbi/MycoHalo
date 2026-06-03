#' Create a MycoHalo plate object
#'
#' @param image A magick-image object.
#'
#' @return A mycohalo_plate object.
#'
#' @export
create_plate <- function(image) {

  if (!inherits(image, "magick-image")) {
    cli::cli_abort(
      "{.arg image} must be a magick-image object."
    )
  }

  plate <- list(
    image = image,
    metadata = magick::image_info(image),
    ebimage = image_to_ebimage(image),
    plate = NULL,
    colonies = NULL,
    halo = NULL,
    melanization = NULL,
    radial_profiles = NULL
  )

  class(plate) <- "mycohalo_plate"

  plate
}
