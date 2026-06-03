#' Diagnose grayscale intensity distribution
#'
#' Visualize the grayscale intensity distribution of an image and
#' optionally overlay threshold values.
#'
#' @param gray An EBImage grayscale image.
#'
#' @param thresholds Numeric vector of thresholds to display.
#'
#' @param bins Number of histogram bins.
#'
#' @return Invisibly returns the histogram object.
#'
#' @examples
#' \dontrun{
#' gray <- EBImage::channel(
#'   plate$ebimage,
#'   mode = "gray"
#' )
#'
#' diagnose_histogram(
#'   gray,
#'   thresholds = c(
#'     0.15,
#'     0.20,
#'     0.25
#'   )
#' )
#' }
#'
#' @export
diagnose_histogram <- function(
    gray,
    thresholds = NULL,
    bins = 100
) {

  if (!inherits(gray, "Image")) {
    cli::cli_abort(
      "{.arg gray} must be an EBImage Image object."
    )
  }

  h <- graphics::hist(
    as.vector(gray),
    breaks = bins,
    col = "grey80",
    border = "white",
    main = "Pixel intensity distribution",
    xlab = "Grayscale intensity",
    ylab = "Frequency"
  )

  if (!is.null(thresholds)) {

    graphics::abline(
      v = thresholds,
      col = "red",
      lwd = 2,
      lty = 2
    )

    graphics::text(
      x = thresholds,
      y = max(h$counts) * 0.95,
      labels = thresholds,
      col = "red",
      pos = 4
    )
  }

  invisible(h)
}
