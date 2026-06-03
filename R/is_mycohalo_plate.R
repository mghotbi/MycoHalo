#' Check whether an object is a MycoHalo plate
#'
#' @param x An object.
#'
#' @return
#' A logical scalar.
#'
#' @examples
#' is_mycohalo_plate(mtcars)
#'
#' @export
is_mycohalo_plate <- function(x) {

  inherits(x, "mycohalo_plate")
}
