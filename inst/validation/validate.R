# Validation of MycoHalo against simulated plates with known ground truth.
# Re-run after any change to the algorithms:  Rscript inst/validation/validate.R
library(MycoHalo)
run <- function(name, sim_args, an_args = list()) {
  sim <- do.call(simulate_plate, sim_args)
  t0 <- Sys.time()
  res <- suppressWarnings(do.call(analyze_plate, c(list(sim$image, verbose = FALSE, id = name), an_args)))
  el <- as.numeric(Sys.time() - t0, units = "secs")
  tr <- sim$truth$colonies
  co <- res$colonies
  m <- merge(tr, co, by = "colony_id", suffixes = c("_true", ""))
  data.frame(scenario = name, colony = m$colony_id, detected = m$detected,
             area_err_pct = round(100 * (m$area_mm2 - m$area_mm2_true) / m$area_mm2_true, 2),
             MI_err = round(m$MI_mean - m$MI_mean_true, 2),
             dMI_facing = round(m$delta_MI_facing, 2),
             scale_err_pct = round(100 * (res$plate_summary$mm_per_px - sim$truth$mm_per_px) / sim$truth$mm_per_px, 3),
             halo_r = round(res$plate_summary$halo_outer_radius_mm, 2),
             sec = round(el, 1), warn = res$plate_summary$n_warnings)
}
out <- rbind(
  run("default", list(seed = 1)),
  run("no_vignette", list(seed = 2, vignetting = 0)),
  run("strong_vignette", list(seed = 3, vignetting = 0.45)),
  run("control_plate", list(seed = 4, layout = plate_layout(centre = "none")), list(layout = plate_layout(centre = "none"))),
  run("facing_dark5", list(seed = 5, facing_darkening = 5)),
  run("touching", list(seed = 6, colony_radius_mm = 11, halo_radius_mm = 14)),
  run("noisy", list(seed = 7, noise_sd = 0.03, texture_sd = 8)),
  run("exposure0.7_card", list(seed = 8, exposure = 0.7, grey_card = TRUE), list(reference = c(0.1, 0.92, 0.05))),
  run("exposure0.7_nocard", list(seed = 8, exposure = 0.7, grey_card = TRUE)),
  run("very_dark_colonies", list(seed = 9, colony_L = c(30, 27, 25, 33))),
  run("missing_colony", list(seed = 10, colony_radius_mm = c(7, 7, 0.01, 7))),
  run("big_halo_over_colonies", list(seed = 11, halo_radius_mm = 26)),
  run("no_halo", list(seed = 12, halo_lab = c(22, -1, -4)))
)
print(out, row.names = FALSE)
write.csv(out, "simulation_validation.csv", row.names = FALSE)
