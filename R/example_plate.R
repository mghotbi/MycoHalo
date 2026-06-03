#' Load example plate image
#'
#' Load an example image distributed with MycoHalo.
#'
#' @param type Character scalar.
#' Either `"control"` or `"confrontation"`.
#'
#' @return
#' A `mycohalo_plate` object.
#'
#' @examples
#' \dontrun{
#' plate <- example_plate("control")
#' }
#'
#' @export
example_plate <- function(type = c(
  "control",
  "confrontation"
)) {

  type <- match.arg(type)

  path <- switch(
    type,
    control = system.file(
      "extdata",
      "control",
      "2026_05_26_37.JPG",
      package = "MycoHalo"
    ),
    confrontation = system.file(
      "extdata",
      "confrontation",
      "2026_05_26_20.JPG",
      package = "MycoHalo"
    )
  )

  read_plate(path)
}
