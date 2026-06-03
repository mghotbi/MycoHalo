#' Diagnose Threshold Selection
#'
#' Visualize the effect of multiple threshold values on a grayscale image.
#'
#' This function is intended as a diagnostic and exploratory tool during
#' image-processing workflow development. It displays a series of binary
#' masks generated from a grayscale image using different threshold values,
#' allowing rapid assessment of threshold sensitivity.
#'
#' Threshold diagnostics are particularly useful when developing colony
#' segmentation, plate detection, halo detection, and melanization
#' quantification workflows.
#'
#' @param gray An EBImage grayscale image.
#'   Typically produced by:
#'
#'   \preformatted{
#'   gray <- EBImage::channel(
#'     plate$ebimage,
#'     mode = "gray"
#'   )
#'   }
#'
#' @param thresholds Numeric vector of threshold values.
#'   Values should generally fall between 0 and 1.
#'
#' @param nrow Number of rows in the diagnostic panel layout.
#'
#' @param ncol Number of columns in the diagnostic panel layout.
#'
#' @return
#' Invisibly returns a list containing:
#'
#' \describe{
#'   \item{thresholds}{
#'     Threshold values evaluated.
#'   }
#'
#'   \item{masks}{
#'     List of binary masks corresponding to each threshold.
#'   }
#' }
#'
#' @details
#' For each threshold value:
#'
#' \enumerate{
#'   \item Pixels greater than the threshold are classified as foreground.
#'   \item Pixels less than or equal to the threshold are classified as background.
#'   \item A binary mask is generated.
#'   \item The mask is displayed for visual inspection.
#' }
#'
#' This function is intended for exploratory image analysis and algorithm
#' development. It is not typically used in production workflows.
#'
#' @seealso
#' \code{\link{detect_plate}}
#'
#' \code{\link{detect_colonies}}
#'
#' @examples
#' \dontrun{
#'
#' plate <- read_plate(
#'   "inst/extdata/confrontation/2026_05_26_20.JPG"
#' )
#'
#' gray <- EBImage::channel(
#'   plate$ebimage,
#'   mode = "gray"
#' )
#'
#' diagnose_threshold(
#'   gray = gray,
#'   thresholds = c(
#'     0.10,
#'     0.15,
#'     0.20,
#'     0.25
#'   )
#' )
#'
#' }
#'
#' @export
diagnose_threshold <- function(
    gray,
    thresholds = c(
      0.10,
      0.15,
      0.20,
      0.25
    ),
    nrow = 2,
    ncol = 2
) {

  if (!inherits(gray, "Image")) {
    cli::cli_abort(
      c(
        "{.arg gray} must be an EBImage Image object.",
        "x" = "Received class: {.cls {class(gray)}}."
      )
    )
  }

  if (!is.numeric(thresholds)) {
    cli::cli_abort(
      "{.arg thresholds} must be numeric."
    )
  }

  old_par <- graphics::par(
    no.readonly = TRUE
  )

  on.exit(
    graphics::par(old_par),
    add = TRUE
  )

  graphics::par(
    mfrow = c(
      nrow,
      ncol
    )
  )

  masks <- vector(
    mode = "list",
    length = length(thresholds)
  )

  for (i in seq_along(thresholds)) {

    thr <- thresholds[i]

    masks[[i]] <- gray > thr

    graphics::image(
      masks[[i]],
      main = paste0(
        "Threshold = ",
        thr
      ),
      axes = FALSE,
      col = c(
        "black",
        "white"
      )
    )
  }

  invisible(
    list(
      thresholds = thresholds,
      masks = masks
    )
  )
}
