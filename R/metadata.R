#' Create a plate metadata sheet for a folder of photos
#'
#' Writes a CSV with one row per image, ready to fill in (in Excel, Numbers
#' or R), which [analyze_plates()] then merges into every output table.
#'
#' @details
#' Columns:
#' \describe{
#'   \item{file}{Image file name (do not edit).}
#'   \item{plate_id}{Identifier used in all outputs (file name by default).}
#'   \item{treatment}{`"control"` (Zymoseptoria alone, no bacterium) or the
#'     name of the confrontation, e.g. `"bacteria"`. This column decides the
#'     layout: `"control"` plates are analysed without a central bacterium.}
#'   \item{bacterium, zymo_strain, medium, day, replicate, photo_date, notes}{
#'     Free design variables; add or remove any columns you like.}
#' }
#' `treatment` can be pre-filled from the file names with `patterns`, e.g.
#' `patterns = c(control = "ctrl|control", bacteria = "bact|conf")`.
#'
#' @param dir Folder with plate images.
#' @param file Output CSV path (default `plate_metadata.csv` in `dir`).
#' @param treatment Default treatment for all plates (e.g. `"bacteria"`), or
#'   `NA` to leave empty.
#' @param patterns Optional named character vector of regular expressions:
#'   names are treatment labels, values are matched against file names.
#' @param extra Named list of additional constant columns, e.g.
#'   `list(bacterium = "isolate_7", zymo_strain = "IPO323", day = 7)`.
#' @param overwrite Overwrite an existing file.
#' @return The metadata data frame (invisibly), also written to `file`.
#' @examples
#' d <- tempfile(); dir.create(d)
#' for (f in c("ctrl_01.png", "ctrl_02.png", "bact_01.png")) {
#'   png::writePNG(array(0.5, c(10, 10, 3)), file.path(d, f))
#' }
#' meta <- make_metadata(d, patterns = c(control = "ctrl", bacteria = "bact"),
#'                       extra = list(bacterium = "isolate_7", day = 7))
#' meta
#' @export
make_metadata <- function(dir, file = file.path(dir, "plate_metadata.csv"),
                          treatment = NA_character_, patterns = NULL,
                          extra = list(), overwrite = FALSE) {
  if (!dir.exists(dir)) cli::cli_abort("Folder {.file {dir}} does not exist.")
  files <- list.files(dir, pattern = "\\.(jpe?g|png|tiff?|heic)$", ignore.case = TRUE)
  if (!length(files)) cli::cli_abort("No images (jpg, png, tif, heic) in {.file {dir}}.")
  trt <- rep(treatment, length(files))
  if (!is.null(patterns)) {
    for (k in names(patterns)) {
      hit <- grepl(patterns[[k]], files, ignore.case = TRUE)
      trt[hit] <- k
    }
  }
  meta <- data.frame(file = files,
                     plate_id = tools::file_path_sans_ext(files),
                     treatment = trt,
                     bacterium = NA_character_, zymo_strain = NA_character_,
                     medium = NA_character_, day = NA_integer_, replicate = NA_integer_,
                     photo_date = NA_character_, notes = NA_character_,
                     stringsAsFactors = FALSE)
  for (k in names(extra)) meta[[k]] <- extra[[k]]
  if (file.exists(file) && !overwrite) {
    cli::cli_abort(c("{.file {file}} already exists.", "i" = "Use {.code overwrite = TRUE} or another {.arg file}."))
  }
  utils::write.csv(meta, file, row.names = FALSE, na = "")
  cli::cli_alert_success("Wrote {nrow(meta)} plate{?s} to {.file {file}}.")
  if (any(is.na(meta$treatment) | meta$treatment == "")) {
    cli::cli_alert_info("Fill in the {.field treatment} column ({.val control} or e.g. {.val bacteria}) before running {.fn analyze_plates}.")
  }
  invisible(meta)
}

#' Layout chooser for a metadata row
#'
#' Convenience function for [analyze_plates()]: plates whose `treatment` is
#' `"control"` (case-insensitive) or empty get `plate_layout(centre = "none")`,
#' all others the confrontation layout with a central bacterium.
#'
#' @param row One metadata row (data frame) or `NULL`.
#' @param ... Passed to [plate_layout()] (e.g. `radius_mm = 20`).
#' @return A `mycohalo_layout`.
#' @examples
#' layout_by_treatment(data.frame(treatment = "control"))
#' @export
layout_by_treatment <- function(row, ...) {
  ctrl <- is.null(row) || is.null(row$treatment) || is.na(row$treatment[1]) ||
    tolower(row$treatment[1]) %in% c("control", "ctrl", "")
  plate_layout(centre = if (ctrl) "none" else "bacteria", ...)
}

#' Remove position-specific lighting bias from side-specific metrics
#'
#' Under directional light every raised colony is brighter on the side that
#' faces the lamp, whatever the treatment. On a fixed copy stand this bias
#' depends on the colony's position in the frame (TL, TR, BR, BL), so it is
#' estimated from control plates, per position, and subtracted from all
#' plates.
#'
#' @param colonies Colony table from [analyze_plates()].
#' @param control Value of `group` that marks control plates.
#' @param group Column holding the treatment.
#' @param vars Side-specific columns to correct.
#' @return `colonies` with extra columns `<var>_corrected`.
#' @examples
#' d <- data.frame(colony_id = rep(c("TL", "TR"), 4),
#'                 treatment = rep(c("control", "bacteria"), each = 4),
#'                 delta_MI_facing = c(4, -8, 5, -7, 9, -2, 10, -3))
#' correct_facing_bias(d)
#' @export
correct_facing_bias <- function(colonies, control = "control", group = "treatment",
                                vars = c("delta_MI_facing", "growth_inhibition_pct")) {
  vars <- intersect(vars, names(colonies))
  ctrl <- colonies[[group]] == control
  if (!any(ctrl, na.rm = TRUE)) cli::cli_abort("No rows with {group} == {.val {control}}.")
  for (v in vars) {
    bias <- tapply(colonies[[v]][ctrl], colonies$colony_id[ctrl], mean, na.rm = TRUE)
    colonies[[paste0(v, "_corrected")]] <- colonies[[v]] - unname(bias[colonies$colony_id])
  }
  colonies
}
