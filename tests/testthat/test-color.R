test_that("sRGB -> CIELAB reference values", {
  lab <- srgb_to_lab(rbind(c(1, 1, 1), c(0, 0, 0), c(0.5, 0.5, 0.5), c(1, 0, 0)))
  expect_equal(unname(lab[1, ]), c(100, 0, 0), tolerance = 1e-3)
  expect_equal(unname(lab[2, ]), c(0, 0, 0), tolerance = 1e-6)
  expect_equal(unname(lab[3, 1]), 53.389, tolerance = 1e-3)
  # sRGB red: L* 53.24, a* 80.09, b* 67.20 (Lindbloom reference)
  expect_equal(unname(lab[4, ]), c(53.24, 80.09, 67.20), tolerance = 1e-3)
})

test_that("Lab <-> sRGB round trip", {
  set.seed(1)
  rgb <- matrix(runif(300, 0.02, 0.98), ncol = 3)
  expect_equal(lab_to_srgb(srgb_to_lab(rgb)), rgb, tolerance = 1e-6)
})

test_that("chroma, hue and delta E", {
  ch <- lab_chroma_hue(a = 0, b = 10)
  expect_equal(ch$chroma, 10)
  expect_equal(ch$hue, 90)
  expect_equal(delta_e76(c(50, 0, 0), c(53, 4, 0)), 5)
})
