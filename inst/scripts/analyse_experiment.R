# ---------------------------------------------------------------------------
# MycoHalo - lab template for analysing one confrontation experiment
#
# 1. Put all plate photos of the experiment in one folder (photo_dir).
# 2. Run STEP 1 once: it writes plate_metadata.csv into that folder.
# 3. Open plate_metadata.csv (Excel / Numbers / R), fill in the columns
#    (at least 'treatment': "control" or e.g. "bacteria"), save as CSV.
# 4. Run STEP 2.
# Get this file:  file.copy(system.file("scripts", "analyse_experiment.R", package = "MycoHalo"), ".")
# ---------------------------------------------------------------------------
library(MycoHalo)

## ---- SETTINGS -------------------------------------------------------------
photo_dir        <- "photos"            # folder with JPG/PNG/TIFF plate photos
output_dir       <- "mycohalo_results"
dish_diameter_mm <- 90                  # diameter of the edge detected (check QC once!)
grey_card        <- NULL                # e.g. c(0.08, 0.92, 0.03) if a grey card is in every photo
colony_distance_mm <- NULL              # e.g. 20 if colonies are inoculated 20 mm from the centre
## -------------------------------------------------------------------------

meta_file <- file.path(photo_dir, "plate_metadata.csv")

## ---- STEP 1: create the metadata sheet (only once) -------------------------
if (!file.exists(meta_file)) {
  make_metadata(photo_dir,
                treatment = "bacteria",                          # default for all plates
                patterns  = c(control = "ctrl|control"),         # file names containing ctrl -> control
                extra     = list(bacterium = "isolate_X", zymo_strain = "strain_Y", day = 7))
  stop("Now fill in ", meta_file, " (treatment = 'control' for plates without bacterium) and run again.")
}

## ---- STEP 2: analyse all plates ---------------------------------------------
dir.create(output_dir, showWarnings = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

out <- analyze_plates(file.path(photo_dir, meta$file), metadata = meta,
                      layout = function(row) layout_by_treatment(row, radius_mm = colony_distance_mm),
                      qc_dir = file.path(output_dir, "qc"),
                      dish_diameter_mm = dish_diameter_mm, reference = grey_card)

# remove the position-specific lighting bias of the side-specific metrics
if (any(out$colonies$treatment == "control")) out$colonies <- correct_facing_bias(out$colonies)

write.csv(out$colonies, file.path(output_dir, "colonies.csv"), row.names = FALSE)
write.csv(out$profiles, file.path(output_dir, "profiles.csv"), row.names = FALSE)
write.csv(out$plates,   file.path(output_dir, "plates.csv"),   row.names = FALSE)
if (nrow(out$failed)) write.csv(out$failed, file.path(output_dir, "failed.csv"), row.names = FALSE)

## ---- STEP 3: plate-aware statistics ---------------------------------------
if (length(unique(out$colonies$treatment)) > 1 && "control" %in% out$colonies$treatment) {
  for (resp in intersect(c("MI_mean", "delta_MI_facing_corrected", "growth_inhibition_pct_corrected", "area_mm2"),
                         names(out$colonies))) {
    cat("\n=====", resp, "=====\n")
    r <- compare_melanization(out$colonies, response = resp, group = "treatment", reference = "control")
    cat(r$method, "\n"); print(r$coefficients)
  }
}

## ---- STEP 4: figures ------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(out$profiles)) {
  g <- plot_profiles(out$profiles, colour_by = "treatment", facet_by = "colony_id")
  ggplot2::ggsave(file.path(output_dir, "profiles.pdf"), g, width = 9, height = 6)
}
dir.create(file.path(output_dir, "landscapes"), showWarnings = FALSE)
for (k in seq_len(nrow(meta))) {
  res <- analyze_plate(file.path(photo_dir, meta$file[k]), layout = layout_by_treatment(meta[k, ]),
                       dish_diameter_mm = dish_diameter_mm, reference = grey_card, verbose = FALSE)
  grDevices::png(file.path(output_dir, "landscapes", paste0(meta$plate_id[k], ".png")),
                 2400, 1300, res = 170)
  plot_melanization_map(res, mi_range = c(30, 90), main = paste(meta$plate_id[k], "-", meta$treatment[k]))
  grDevices::dev.off()
}
cat("\nDone. Inspect every QC image in", file.path(output_dir, "qc"), "before using the numbers.\n")
