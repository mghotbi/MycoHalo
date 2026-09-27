test_that("exact EDT matches brute force", {
  set.seed(1)
  m <- matrix(runif(30 * 40) > 0.15, 30, 40)
  d <- MycoHalo:::edt(m)
  bg <- which(!m, arr.ind = TRUE)
  fg <- which(m, arr.ind = TRUE)
  bf <- apply(fg, 1, function(p) sqrt(min((bg[, 1] - p[1])^2 + (bg[, 2] - p[2])^2)))
  expect_equal(d[m], bf, tolerance = 1e-10)
  expect_true(all(d[!m] == 0))
})

test_that("connected components respect connectivity", {
  m <- matrix(FALSE, 5, 5)
  m[1, 1] <- m[2, 2] <- TRUE
  expect_equal(max(MycoHalo:::label_components(m, 8)), 1)
  expect_equal(max(MycoHalo:::label_components(m, 4)), 2)
})

test_that("box blur preserves constants and means", {
  x <- matrix(3, 20, 30)
  expect_equal(MycoHalo:::box_blur_cpp(x, 3L, 3L), x)
  set.seed(2)
  y <- matrix(rnorm(400), 20)
  b <- MycoHalo:::box_blur_cpp(y, 1L, 1L)
  expect_equal(b[10, 10], mean(y[9:11, 9:11]))
})

test_that("watershed splits two touching disks along the neck", {
  xx <- matrix(rep(1:80, each = 50), 50, 80)
  yy <- matrix(rep(1:50, 80), 50, 80)
  m <- (xx - 25)^2 + (yy - 25)^2 <= 15^2 | (xx - 52)^2 + (yy - 25)^2 <= 15^2
  mk <- matrix(0L, 50, 80); mk[25, 25] <- 1L; mk[25, 52] <- 2L
  lab <- MycoHalo:::watershed_cpp(-MycoHalo:::edt(m), mk, m)
  expect_true(all(lab[m] > 0))
  expect_true(all(lab[25, 10:36] == 1L))
  expect_true(all(lab[25, 41:66] == 2L))
})

test_that("fill_holes fills enclosed holes but not forbidden ones", {
  m <- matrix(TRUE, 9, 9); m[c(1, 9), ] <- FALSE; m[, c(1, 9)] <- FALSE
  m[5, 5] <- FALSE
  expect_true(MycoHalo:::fill_holes(m)[5, 5])
  f <- matrix(FALSE, 9, 9); f[5, 5] <- TRUE
  expect_false(MycoHalo:::fill_holes(m, forbid = f)[5, 5])
})

test_that("Otsu separates a bimodal sample", {
  set.seed(3)
  x <- c(rnorm(500, 0, 1), rnorm(500, 10, 1))
  expect_true(abs(MycoHalo:::otsu_threshold(x) - 5) < 1.5)
})
