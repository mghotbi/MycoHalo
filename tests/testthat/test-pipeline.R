# Small simulated plates keep the tests fast (~ 3 s each)
small <- function(...) simulate_plate(width = 600, height = 800, ...)

test_that("dish detection recovers centre, radius and scale", {
  sim <- small(seed = 1)
  p <- detect_plate(read_plate(sim$image))
  expect_lt(abs(p$dish$x - sim$truth$dish[["x"]]), 1.5)
  expect_lt(abs(p$dish$y - sim$truth$dish[["y"]]), 1.5)
  expect_lt(abs(p$dish$mm_per_px / sim$truth$mm_per_px - 1), 0.005)
})

test_that("confrontation plate: colonies, area and melanization are accurate", {
  sim <- small(seed = 2)
  res <- analyze_plate(sim$image, verbose = FALSE)
  m <- merge(sim$truth$colonies, res$colonies, by = "colony_id", suffixes = c("_true", ""))
  expect_true(all(m$detected))
  expect_true(all(abs(m$area_mm2 / m$area_mm2_true - 1) < 0.04))
  expect_true(all(abs(m$MI_mean - m$MI_mean_true) < 1))
  expect_gt(res$plate_summary$bacteria_area_mm2, 0)
  expect_gt(res$plate_summary$halo_area_mm2, 0)
  expect_equal(res$plate_summary$n_satellites, 1L)
  # melanization ranking of the four colonies is recovered exactly
  expect_equal(order(m$MI_mean), order(m$MI_mean_true))
})

test_that("control plate without bacteria", {
  lay <- plate_layout(centre = "none")
  sim <- small(seed = 3, layout = lay)
  res <- analyze_plate(sim$image, layout = lay, verbose = FALSE)
  expect_true(all(res$colonies$detected))
  expect_equal(res$plate_summary$bacteria_area_mm2, 0)
  expect_true(all(res$colonies$interaction_target == "dish_centre"))
})

test_that("a missing colony is reported, not invented", {
  sim <- small(seed = 4, colony_radius_mm = c(7, 7, 0.01, 7))
  res <- suppressWarnings(analyze_plate(sim$image, verbose = FALSE))
  expect_false(res$colonies$detected[res$colonies$colony_id == "BR"])
  expect_true(any(grepl("not detected", res$warnings)))
})

test_that("side-specific melanization response is detected", {
  sim <- small(seed = 5, facing_darkening = 6)
  res <- analyze_plate(sim$image, verbose = FALSE)
  expect_true(all(res$colonies$delta_MI_facing > 3))
})

test_that("grey-card calibration removes exposure error", {
  sim <- small(seed = 6, exposure = 0.7, grey_card = TRUE)
  raw <- analyze_plate(sim$image, verbose = FALSE)
  cal <- suppressMessages(analyze_plate(sim$image, reference = sim$truth$grey_card, verbose = FALSE))
  err_raw <- mean(abs(raw$colonies$MI_mean - sim$truth$colonies$MI_mean))
  err_cal <- mean(abs(cal$colonies$MI_mean - sim$truth$colonies$MI_mean))
  expect_gt(err_raw, 4)
  expect_lt(err_cal, 1)
})

test_that("profiles and QC plot are produced", {
  sim <- small(seed = 7, edge_lightening = 10)
  res <- analyze_plate(sim$image, verbose = FALSE)
  pr <- res$profiles
  expect_true(all(c("colony_id", "mid_mm", "MI_mean", "MI_lo", "MI_hi") %in% names(pr)))
  # the lighter margin gives a lower MI in the outermost ring
  first <- pr[pr$ring == 1, ]
  inner <- pr[pr$ring == 4, ]
  expect_true(all(first$MI_mean < inner$MI_mean))
  f <- tempfile(fileext = ".png")
  grDevices::png(f, 1200, 600)
  expect_silent(plot_qc(res))
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("supervised classifier path runs", {
  sim <- small(seed = 8)
  p <- read_plate(sim$image) |> detect_plate() |> model_background()
  tr <- sim$truth
  reg <- data.frame(class = c("fungus", "fungus", "bacteria", "halo"),
                    x = c(tr$colonies$x[1:2], p$dish$x, p$dish$x + 150),
                    y = c(tr$colonies$y[1:2], p$dish$y, p$dish$y), r = 10)
  clf <- train_classifier(extract_training_pixels(p, reg))
  expect_s3_class(clf, "mycohalo_classifier")
  res <- analyze_plate(sim$image, classifier = clf, verbose = FALSE)
  expect_true(all(res$colonies$detected))
})

test_that("batch analysis with metadata and per-plate layouts", {
  dir <- withr::local_tempdir()
  s1 <- small(seed = 9)
  s2 <- small(seed = 10, layout = plate_layout(centre = "none"))
  png::writePNG(s1$image, file.path(dir, "plate_A.png"))
  png::writePNG(s2$image, file.path(dir, "plate_B.png"))
  meta <- data.frame(file = c("plate_A.png", "plate_B.png"),
                     treatment = c("bacteria", "control"))
  lay <- function(row) plate_layout(centre = if (row$treatment == "control") "none" else "bacteria")
  out <- analyze_plates(dir, metadata = meta, layout = lay,
                        qc_dir = file.path(dir, "qc"), verbose = FALSE)
  expect_equal(nrow(out$failed), 0)
  expect_equal(nrow(out$colonies), 8)
  expect_true(all(c("treatment", "plate_id") %in% names(out$colonies)))
  expect_equal(length(list.files(file.path(dir, "qc"))), 2)
})

test_that("melanization landscape renders", {
  sim <- small(seed = 11, facing_darkening = 4)
  res <- analyze_plate(sim$image, verbose = FALSE)
  f <- tempfile(fileext = ".png")
  grDevices::png(f, 1600, 900)
  out <- plot_melanization_map(res, grid = 80)
  grDevices::dev.off()
  expect_true(file.exists(f))
  expect_length(out$mi_range, 2)
})
