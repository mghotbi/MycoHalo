#' Convert image to grayscale matrix
#'
#' Convert a magick image into a numeric grayscale matrix.
#'
#' @param image A magick-image object.
#'
#' @return
#' A numeric matrix with values between 0 and 255.
#'
#' @keywords internal
#' @noRd
image_to_matrix <- function(image) {

  if (!inherits(image, "magick-image")) {
    cli::cli_abort(
      c(
        "{.arg image} must be a magick-image object.",
        "x" = "Received class: {.cls {class(image)}}."
      )
    )
  }

  img_data <- magick::image_data(
    image,
    channels = "gray"
  )

  matrix(
    as.integer(img_data),
    nrow = dim(img_data)[2],
    ncol = dim(img_data)[3]
  )
}
