#' Detect Petri dish
#'
#' Detect the Petri dish region from a fungal plate image.
#'
#' The function:
#'
#' \enumerate{
#'   \item Converts the image to grayscale.
#'   \item Applies intensity thresholding.
#'   \item Labels connected components.
#'   \item Selects the largest connected component.
#'   \item Stores the resulting plate mask and summary statistics.
#' }
#'
#' @param plate A `mycohalo_plate` object.
#'
#' @param threshold Numeric threshold between 0 and 1.
#'
#' @return
#' Updated `mycohalo_plate` object.
#'
#' @examples
#' \dontrun{
#'
#' plate <- read_plate(
#'   "plate.jpg"
#' )
#'
#' plate <- detect_plate(
#'   plate,
#'   threshold = 0.25
#' )
#'
#' }
#'
#' @export
#'
#' @importFrom cli cli_abort
detect_plate <- function(
    plate,
    threshold = 0.25
) {

  if (!inherits(
    plate,
    "mycohalo_plate"
  )) {

    cli::cli_abort(
      "{.arg plate} must be a mycohalo_plate object."
    )

  }

  gray <- EBImage::channel(
    plate$ebimage,
    mode = "gray"
  )

  mask <- gray > threshold

  labels <- EBImage::bwlabel(mask)

  features <- EBImage::computeFeatures.shape(
    labels
  )

  areas <- features[, "s.area"]

  largest_component <- which.max(
    areas
  )

  plate_mask <- labels == largest_component

  plate$gray <- gray

  plate$plate <- list(
    mask = plate_mask,
    area_px = max(areas, na.rm = TRUE),
    threshold = threshold
  )

  plate
}
