#' Read image as EBImage object
#'
#' @keywords internal
#' @noRd
image_to_ebimage <- function(image) {

  tmp <- tempfile(fileext = ".png")

  magick::image_write(
    image,
    path = tmp,
    format = "png"
  )

  EBImage::readImage(tmp)
}
