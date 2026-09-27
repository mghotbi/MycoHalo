#' MycoHalo: colour-calibrated quantification of fungal melanization in
#' confrontation assays
#'
#' MycoHalo measures melanization of *Zymoseptoria tritici* (or any other
#' fungus) colonies grown at fixed inoculation positions on Petri dishes,
#' alone (controls) or confronted with a bacterium inoculated in the centre.
#'
#' @section Why colour and not grayscale:
#' On a photographed plate the grey fungal colonies, the cream-white
#' bacterial colony and its yellowish diffusion halo overlap strongly in
#' grayscale brightness, so no single grey threshold can separate them.
#' They are, however, well separated in the perceptually uniform CIELAB
#' colour space: melanized hyphae are *achromatic* (chroma \eqn{C^*_{ab}}
#' close to 0) and textured, whereas bacterial growth and the halo are
#' *warm* (positive \eqn{b^*}, hue angle ~60-100 degrees) and smooth.
#' MycoHalo therefore
#' \enumerate{
#'   \item converts sRGB to CIE 1976 \eqn{L^*a^*b^*} (D65),
#'   \item models the agar background and its illumination gradient
#'     (flat-field) with a robust polynomial surface,
#'   \item classifies every non-agar pixel as fungus, bacterium or halo with a
#'     Gaussian mixture initialised from the known plate layout, and
#'   \item separates individual colonies with a marker-controlled watershed
#'     seeded at the expected inoculation positions.
#' }
#'
#' @section Main entry points:
#' \itemize{
#'   \item [analyze_plate()] — full pipeline for one image.
#'   \item [analyze_plates()] — batch processing of a folder / file vector.
#'   \item [plot_qc()] — quality-control overlay (always inspect!).
#'   \item [plot_melanization_map()] — 3D shaded melanization landscape and
#'     2D melanization map with the side-specific response.
#'   \item [plate_layout()] — describe where colonies were inoculated.
#'   \item [simulate_plate()] — synthetic plates with known ground truth,
#'     used for validation and unit tests.
#'   \item [compare_melanization()] — plate-aware statistics.
#' }
#'
#' @keywords internal
#' @useDynLib MycoHalo, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
