testthat::test_that(
  "read_plate imports image",
  {

    test_path <- testthat::test_path(
      "testdata",
      "control_01.jpg"
    )

    testthat::expect_s3_class(
      read_plate(test_path),
      "magick-image"
    )

  }
)

