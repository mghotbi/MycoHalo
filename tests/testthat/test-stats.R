test_that("compare_melanization respects plates", {
  set.seed(1)
  d <- data.frame(plate_id = rep(sprintf("p%02d", 1:12), each = 4),
                  treatment = rep(c("control", "bacteria"), each = 24))
  d$MI_mean <- 50 + 4 * (d$treatment == "bacteria") + rep(rnorm(12, 0, 1), each = 4) + rnorm(48, 0, 0.5)
  r1 <- compare_melanization(d, reference = "control", method = "plate_means")
  expect_equal(r1$coefficients$term, "bacteria")
  expect_true(r1$coefficients$lower < 4 && r1$coefficients$upper > 4)
  skip_if_not_installed("lme4")
  r2 <- compare_melanization(d, reference = "control", method = "mixed")
  est <- r2$coefficients$estimate[r2$coefficients$term == "bacteria"]
  expect_equal(est, r1$coefficients$estimate, tolerance = 1e-6)
  expect_gt(r2$icc, 0.3)
})
