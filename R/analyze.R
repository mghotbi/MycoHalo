#' Analyse one plate photograph
#'
#' Runs the complete MycoHalo pipeline: import, dish detection and scale,
#' optional colour calibration, flat-field correction, colour/texture pixel
#' classification, seeded colony separation, measurement and QC.
#'
#' @param path Image file or H x W x 3 sRGB array.
#' @param layout A [plate_layout()]. Use `plate_layout(centre = "none")` for
#'   control plates without bacteria.
#' @param id Plate identifier (default: file name).
#' @param dish_diameter_mm Outer diameter of the dish (mm); sets the scale.
#' @param agar_fraction Fraction of the dish radius that is analysed.
#' @param dish Optional manual dish circle `c(x, y, r)` in working pixels.
#' @param mm_per_px Optional known scale (mm per working-image pixel);
#'   overrides `dish_diameter_mm`.
#' @param reference Optional exposure / white-balance reference: a grey-card
#'   region `c(x, y, r)` in relative image coordinates, or `"agar"`; see
#'   [model_background()].
#' @param reference_lab Known CIELAB of the reference.
#' @param max_dim Working resolution (longest side, px).
#' @param flat_field Apply illumination correction (recommended).
#' @param classifier Optional supervised classifier ([train_classifier()]).
#' @param min_delta_e Minimum colour difference from agar for objects.
#' @param min_area_mm2,min_core_mm,morph_mm,edge_criterion Segmentation settings, see
#'   [segment_colonies()].
#' @param edge_exclude_mm,dark_L,n_boot,block_mm Measurement settings, see
#'   [measure_colonies()].
#' @param ring_width_mm,profile_mode,n_rings Profile settings, see
#'   [colony_profiles()].
#' @param qc_file Optional path of a PNG quality-control figure.
#' @param keep_images Keep image arrays in the result (needed for
#'   [plot_qc()] later; set `FALSE` in large batches to save memory).
#' @param verbose Print progress.
#'
#' @return An object of class `mycohalo_result`, a list with
#'   \describe{
#'     \item{colonies}{One row per expected colony ([measure_colonies()]).}
#'     \item{profiles}{Edge-to-centre profiles ([colony_profiles()]).}
#'     \item{plate_summary}{One row of plate-level measures and QC.}
#'     \item{plate}{The processed plate (images and masks) if
#'       `keep_images = TRUE`.}
#'     \item{params, warnings}{Settings used and QC warnings.}
#'   }
#'
#' @examples
#' sim <- simulate_plate(seed = 10)
#' res <- analyze_plate(sim$image, id = "sim10", verbose = FALSE)
#' res
#' res$colonies[, c("colony_id", "MI_mean", "delta_MI_facing")]
#' \dontrun{
#' res <- analyze_plate("plates/2026_05_26_20.JPG", layout = plate_layout())
#' plot_qc(res)
#' ctrl <- analyze_plate("plates/2026_05_26_37.JPG",
#'                       layout = plate_layout(centre = "none"))
#' }
#' @export
analyze_plate <- function(path, layout = plate_layout(), id = NULL,
                          dish_diameter_mm = 90, agar_fraction = 0.88, dish = NULL,
                          mm_per_px = NULL,
                          reference = NULL, reference_lab = c(50, 0, 0),
                          max_dim = 2000, flat_field = TRUE, classifier = NULL,
                          min_delta_e = 3, min_area_mm2 = 0.5, min_core_mm = 0.4,
                          morph_mm = 0.08, edge_criterion = 0.5, edge_exclude_mm = 0.1, dark_L = NULL,
                          n_boot = 200L, block_mm = 0.5, ring_width_mm = 0.5,
                          profile_mode = c("absolute", "relative"), n_rings = 10L,
                          qc_file = NULL, keep_images = TRUE, verbose = TRUE) {
  profile_mode <- match.arg(profile_mode)
  stopifnot(inherits(layout, "mycohalo_layout"))
  params <- as.list(environment())
  params$path <- if (is.character(path)) path else "<array>"
  params$layout <- layout
  step <- function(msg) if (verbose) cli::cli_alert_info(msg)
  warn <- character()

  step("Reading image")
  plate <- read_plate(path, max_dim = max_dim, id = id)
  step("Detecting dish")
  plate <- detect_plate(plate, dish = dish, dish_diameter_mm = dish_diameter_mm,
                        agar_fraction = agar_fraction, mm_per_px = mm_per_px)
  step(sprintf("Dish r = %.0f px, %.4f mm/px", plate$dish$r_outer, plate$dish$mm_per_px))
  step("Modelling agar background and illumination")
  plate <- model_background(plate, flat_field = flat_field, min_delta_e = min_delta_e,
                            reference = reference, reference_lab = reference_lab)
  step("Classifying pixels (fungus / bacteria / halo)")
  plate <- classify_pixels(plate, layout, classifier = classifier)
  step("Separating colonies")
  plate <- segment_colonies(plate, layout, min_area_mm2 = min_area_mm2,
                            min_core_mm = min_core_mm, morph_mm = morph_mm,
                            edge_criterion = edge_criterion)
  step("Measuring")
  col <- measure_colonies(plate, edge_exclude_mm = edge_exclude_mm, dark_L = dark_L,
                          n_boot = n_boot, block_mm = block_mm)
  prof <- colony_profiles(plate, ring_width_mm = ring_width_mm, mode = profile_mode,
                          n_rings = n_rings, n_boot = n_boot, block_mm = block_mm)
  ps <- measure_plate(plate)

  # ---- QC warnings
  if (any(!col$detected)) warn <- c(warn, sprintf("Colony not detected: %s", paste(col$colony_id[!col$detected], collapse = ", ")))
  if (!is.na(plate$dish$arc_coverage) && plate$dish$arc_coverage < 0.4)
    warn <- c(warn, "Dish edge visible on < 40 % of its circumference; check the scale.")
  if (isTRUE(ps$illumination_range_pct > 30))
    warn <- c(warn, sprintf("Strong illumination gradient (%.0f %%); improve lighting.", ps$illumination_range_pct))
  if (any(col$specular_fraction > 0.02, na.rm = TRUE))
    warn <- c(warn, "Colonies contain > 2 % specular / over-exposed pixels; reduce exposure or use diffuse light.")
  if (any(col$offset_from_expected_mm > 8, na.rm = TRUE))
    warn <- c(warn, "A colony centroid is > 8 mm from its expected position; check layout / orientation.")
  if (layout$centre == "bacteria" && ps$bacteria_area_mm2 == 0)
    warn <- c(warn, "No bacterial colony found at the centre.")
  for (w in warn) cli::cli_warn(w)

  col <- cbind(plate_id = plate$id, col, stringsAsFactors = FALSE)
  if (nrow(prof)) prof <- cbind(plate_id = plate$id, prof, stringsAsFactors = FALSE)
  ps <- cbind(plate_id = plate$id, n_warnings = length(warn),
              warnings = paste(warn, collapse = " | "), ps, stringsAsFactors = FALSE)

  res <- structure(list(colonies = col, profiles = prof, plate_summary = ps,
                        plate = plate, params = params, warnings = warn,
                        version = as.character(utils::packageVersion("MycoHalo"))),
                   class = "mycohalo_result")
  if (!is.null(qc_file)) {
    grDevices::png(qc_file, width = 2400, height = 1250, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_qc(res)
  }
  if (!keep_images) res$plate <- NULL
  res
}

#' @export
print.mycohalo_result <- function(x, ...) {
  cli::cli_text("{.strong <mycohalo_result>} plate {.val {x$plate_summary$plate_id}}")
  ps <- x$plate_summary
  cli::cli_bullets(c(
    "*" = "Scale: {signif(ps$mm_per_px, 4)} mm/px; agar L* = {round(ps$agar_L, 1)}",
    "*" = "Colonies detected: {ps$n_colonies_detected} / {nrow(x$colonies)}; satellites: {ps$n_satellites}",
    "*" = if (!is.na(ps$bacteria_radius_mm)) "Bacterium r = {round(ps$bacteria_radius_mm, 2)} mm; halo outer r = {round(ps$halo_outer_radius_mm, 2)} mm" else "No bacterium"
  ))
  cols <- intersect(c("colony_id", "area_mm2", "L_mean", "MI_mean", "MI_lo", "MI_hi",
                      "delta_MI_facing", "growth_inhibition_pct"), names(x$colonies))
  print(format(x$colonies[, cols], digits = 3), row.names = FALSE)
  if (length(x$warnings)) cli::cli_alert_warning("{length(x$warnings)} QC warning{?s}: {x$warnings}")
  invisible(x)
}

#' Analyse many plates
#'
#' Batch wrapper around [analyze_plate()] that never stops on a single bad
#' image: failures are recorded and reported. A QC figure is written for
#' every plate — inspect them before using the numbers.
#'
#' @param files Character vector of image paths, or a directory (all
#'   JPEG/PNG/TIFF files are used).
#' @param metadata Optional data frame with a `file` column (base name) or
#'   `plate_id` column plus any design variables (treatment, strain,
#'   bacterium, day, replicate ...). Merged into all output tables.
#' @param layout A single [plate_layout()], or a function
#'   `function(metadata_row)` returning the layout for each plate (e.g.
#'   `centre = "none"` for controls).
#' @param qc_dir Directory for QC PNGs (created if needed); `NULL` to skip.
#' @param ... Further arguments passed to [analyze_plate()].
#' @param verbose Print progress.
#'
#' @return A list of three data frames (`colonies`, `profiles`, `plates`)
#'   and a `failed` data frame.
#' @examples
#' \dontrun{
#' meta <- data.frame(file = c("2026_05_26_20.JPG", "2026_05_26_37.JPG"),
#'                    treatment = c("bacteria", "control"))
#' lay <- function(row) plate_layout(centre = if (row$treatment == "control") "none" else "bacteria")
#' out <- analyze_plates("plates/", metadata = meta, layout = lay, qc_dir = "qc/")
#' write.csv(out$colonies, "colonies.csv", row.names = FALSE)
#' }
#' @export
analyze_plates <- function(files, metadata = NULL, layout = plate_layout(),
                           qc_dir = "mycohalo_qc", ..., verbose = TRUE) {
  if (length(files) == 1 && dir.exists(files)) {
    files <- list.files(files, pattern = "\\.(jpe?g|png|tiff?|heic)$",
                        ignore.case = TRUE, full.names = TRUE)
  }
  if (!length(files)) cli::cli_abort("No image files found.")
  if (!is.null(qc_dir)) dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  out_c <- list(); out_p <- list(); out_s <- list(); failed <- list()
  for (k in seq_along(files)) {
    f <- files[k]
    pid <- tools::file_path_sans_ext(basename(f))
    mrow <- NULL
    if (!is.null(metadata)) {
      if ("file" %in% names(metadata)) mrow <- metadata[basename(metadata$file) == basename(f), , drop = FALSE]
      else if ("plate_id" %in% names(metadata)) mrow <- metadata[metadata$plate_id == pid, , drop = FALSE]
      if (!is.null(mrow) && nrow(mrow) == 0) mrow <- NULL
    }
    lay <- if (is.function(layout)) layout(mrow) else layout
    if (verbose) cli::cli_alert("[{k}/{length(files)}] {basename(f)}")
    qc <- if (is.null(qc_dir)) NULL else file.path(qc_dir, paste0(pid, "_qc.png"))
    r <- tryCatch(analyze_plate(f, layout = lay, id = pid, qc_file = qc,
                                keep_images = FALSE, verbose = FALSE, ...),
                  error = function(e) e)
    if (inherits(r, "error")) {
      failed[[length(failed) + 1]] <- data.frame(file = f, error = conditionMessage(r))
      if (verbose) cli::cli_alert_danger(conditionMessage(r))
      next
    }
    addm <- function(df) {
      if (is.null(mrow) || !nrow(df)) return(df)
      extra <- mrow[1, setdiff(names(mrow), names(df)), drop = FALSE]
      cbind(df, extra[rep(1, nrow(df)), , drop = FALSE], row.names = NULL)
    }
    out_c[[pid]] <- addm(r$colonies)
    out_p[[pid]] <- addm(r$profiles)
    out_s[[pid]] <- addm(cbind(r$plate_summary, file = f))
  }
  list(colonies = bind_rows_fill(out_c), profiles = bind_rows_fill(out_p),
       plates = bind_rows_fill(out_s),
       failed = if (length(failed)) do.call(rbind, failed) else data.frame(file = character(), error = character()))
}
