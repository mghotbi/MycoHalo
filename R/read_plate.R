#' Read a plate photograph
#'
#' Imports a plate image (JPEG, PNG, TIFF, HEIC ... anything ImageMagick
#' reads), applies the EXIF orientation, optionally down-samples it to a
#' working resolution and converts it to CIELAB.
#'
#' @details
#' **Working resolution.** Segmentation and colour statistics are computed on
#' an image whose longest side is at most `max_dim` pixels. For a 90 mm dish
#' photographed so that it fills the frame, `max_dim = 2000` gives roughly
#' 0.05 mm per pixel, i.e. a 1 cm colony is described by ~30 000 pixels —
#' far more than needed for sub-percent precision of mean colour and ~1 %
#' precision of area. Down-sampling uses an anti-aliasing (area-averaging)
#' filter, so it *reduces* JPEG block noise without biasing mean colour.
#' Use `max_dim = Inf` to keep the native resolution.
#'
#' **Clipping.** Pixels with any channel at 255 carry no colour
#' information (the sensor was saturated). Their fraction is stored and
#' [analyze_plate()] warns if colonies contain clipped pixels. Photograph
#' with manual exposure so that the brightest bacterial growth stays below
#' ~ 245.
#'
#' @param path Path to an image file, **or** a numeric H x W x 3 array with
#'   sRGB values in \[0, 1\] (e.g. from [simulate_plate()]).
#' @param max_dim Maximum size (px) of the longest image side used for
#'   analysis. Default 2000.
#' @param id Optional plate identifier; defaults to the file name without
#'   extension.
#'
#' @return An object of class `mycohalo_plate`: a list with
#'   \describe{
#'     \item{rgb}{H x W x 3 array of sRGB values in \[0, 1\].}
#'     \item{lab}{H x W x 3 array of CIELAB values.}
#'     \item{id, path}{Identifier and source path.}
#'     \item{original_dim}{Width and height of the file (px).}
#'     \item{resize_factor}{Working px per original px.}
#'     \item{clipped_fraction}{Fraction of pixels with a channel at 255.}
#'   }
#'
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 1)
#' plate <- read_plate(sim$image, id = "sim1")
#' plate
#'
#' @seealso [analyze_plate()], [calibrate_color()]
#' @export
read_plate <- function(path, max_dim = 2000, id = NULL) {
  if (is.array(path) && length(dim(path)) == 3L) {
    rgb <- path
    if (dim(rgb)[3] > 3) rgb <- rgb[, , 1:3]
    if (max(rgb, na.rm = TRUE) > 1) rgb <- rgb / 255
    src <- NA_character_
    orig <- c(width = dim(rgb)[2], height = dim(rgb)[1])
    fac <- 1
    if (max(dim(rgb)[1:2]) > max_dim) {
      fac <- max_dim / max(dim(rgb)[1:2])
      rgb <- resize_array(rgb, fac)
    }
    if (is.null(id)) id <- "plate"
  } else {
    if (!is.character(path) || length(path) != 1L) {
      cli::cli_abort("{.arg path} must be a single file path or an H x W x 3 array.")
    }
    if (!file.exists(path)) {
      cli::cli_abort(c("Image file does not exist.", "x" = "Cannot find {.file {path}}."))
    }
    img <- magick::image_read(path)
    img <- magick::image_orient(img)
    info <- magick::image_info(img)
    orig <- c(width = info$width, height = info$height)
    fac <- 1
    if (max(orig) > max_dim) {
      fac <- max_dim / max(orig)
      geom <- magick::geometry_size_pixels(
        width = round(orig[["width"]] * fac),
        height = round(orig[["height"]] * fac),
        preserve_aspect = FALSE
      )
      img <- magick::image_resize(img, geom, filter = "Triangle")
    }
    img <- magick::image_convert(img, colorspace = "sRGB", type = "TrueColor")
    arr <- as.integer(magick::image_data(img, channels = "rgb"))
    rgb <- array(as.numeric(arr) / 255, dim = dim(arr))
    src <- normalizePath(path)
    if (is.null(id)) id <- tools::file_path_sans_ext(basename(path))
  }
  storage.mode(rgb) <- "double"
  # highlight clipping (over-exposure); dark clipping of the background is harmless
  clipped <- mean(rgb[, , 1] >= 254.5 / 255 | rgb[, , 2] >= 254.5 / 255 | rgb[, , 3] >= 254.5 / 255)
  new_plate(rgb = rgb, id = id, path = src, original_dim = orig,
            resize_factor = fac, clipped_fraction = clipped)
}

#' @keywords internal
#' @noRd
new_plate <- function(rgb, id, path, original_dim, resize_factor, clipped_fraction) {
  lab <- srgb_to_lab(rgb)
  structure(
    list(
      rgb = rgb,
      lab = lab,
      id = id,
      path = path,
      original_dim = original_dim,
      resize_factor = resize_factor,
      clipped_fraction = clipped_fraction,
      calibration = NULL
    ),
    class = "mycohalo_plate"
  )
}

#' Area-averaging resize of an H x W x 3 array (internal, used for arrays)
#' @keywords internal
#' @noRd
resize_array <- function(rgb, fac) {
  H <- dim(rgb)[1]; W <- dim(rgb)[2]
  nh <- max(1L, round(H * fac)); nw <- max(1L, round(W * fac))
  ri <- findInterval(seq_len(H) - 0.5, seq(0, H, length.out = nh + 1), rightmost.closed = TRUE)
  ci <- findInterval(seq_len(W) - 0.5, seq(0, W, length.out = nw + 1), rightmost.closed = TRUE)
  out <- array(0, c(nh, nw, 3))
  for (k in 1:3) {
    m <- rgb[, , k]
    rs <- rowsum(m, ri)                   # nh x W
    cs <- t(rowsum(t(rs), ci))            # nh x nw
    cnt <- outer(tabulate(ri, nh), tabulate(ci, nw))
    out[, , k] <- cs / cnt
  }
  out
}

#' Test for a MycoHalo plate object
#' @param x Any object.
#' @return `TRUE` if `x` inherits from `mycohalo_plate`.
#' @examples
#' is_mycohalo_plate(1)
#' @export
is_mycohalo_plate <- function(x) inherits(x, "mycohalo_plate")

#' @export
print.mycohalo_plate <- function(x, ...) {
  cli::cli_text("{.strong <mycohalo_plate>} {.val {x$id}}")
  cli::cli_bullets(c(
    "*" = "Working size: {dim(x$rgb)[2]} x {dim(x$rgb)[1]} px (resize factor {signif(x$resize_factor, 3)})",
    "*" = "Source: {.file {x$path}}",
    "*" = "Over-exposed pixels: {signif(100 * x$clipped_fraction, 2)} %",
    "*" = "Colour calibration: {if (is.null(x$calibration)) 'none' else x$calibration$method}"
  ))
  invisible(x)
}

#' Colour / exposure calibration against a reference patch
#'
#' Rescales the image in *linear* RGB so that a neutral reference patch
#' photographed in every image (grey card, white standard, colour-checker
#' neutral patch) takes its known CIELAB value. This is a diagonal
#' (von Kries-type) correction: it removes differences in exposure and white
#' balance between photographs, which is essential for comparing absolute
#' \eqn{L^*} (melanization) across plates, days and cameras.
#'
#' @details
#' A photographic 18 % grey card has \eqn{L^* \approx 50}; the neutral 5
#' patch of an X-Rite ColorChecker Classic has \eqn{L^* \approx 51}. The patch
#' must be evenly lit and not clipped.
#'
#' @param plate A `mycohalo_plate`.
#' @param region Reference patch as `c(x, y, r)` in **relative** image
#'   coordinates (fractions of width, height and of the shorter side), e.g.
#'   `c(0.08, 0.92, 0.03)` for a card in the bottom-left corner.
#' @param ref_lab Known CIELAB of the patch (default neutral `c(50, 0, 0)`).
#' @return The calibrated `mycohalo_plate` (with `$calibration` gains).
#' @examples
#' sim <- simulate_plate(width = 300, height = 400, seed = 2, exposure = 0.8, grey_card = TRUE)
#' plate <- read_plate(sim$image)
#' plate <- calibrate_color(plate, region = sim$truth$grey_card)
#' @export
calibrate_color <- function(plate, region, ref_lab = c(50, 0, 0)) {
  stopifnot(is_mycohalo_plate(plate), length(region) == 3)
  H <- dim(plate$rgb)[1]; W <- dim(plate$rgb)[2]
  yy <- row(plate$rgb[, , 1]); xx <- col(plate$rgb[, , 1])
  inside <- (xx - region[1] * W)^2 + (yy - region[2] * H)^2 <= (region[3] * min(H, W))^2
  apply_reference_gain(plate, inside, ref_lab, "reference patch")
}

#' Diagonal gain correction in linear RGB (internal)
#' @keywords internal
#' @noRd
apply_reference_gain <- function(plate, mask, ref_lab, method) {
  if (sum(mask) < 20) cli::cli_abort("Reference region contains < 20 pixels.")
  lin <- srgb_linearize(matrix(plate$rgb, ncol = 3))
  meas <- apply(lin[as.vector(mask), , drop = FALSE], 2, stats::median)
  target <- srgb_linearize(lab_to_srgb(matrix(ref_lab, ncol = 3)))
  gain <- as.vector(target) / meas
  lin <- sweep(lin, 2, gain, `*`)
  if (mean(lin > 1) > 1e-4) {
    cli::cli_warn("{signif(100 * mean(lin > 1), 2)} % of channel values exceed 1 after calibration and were clipped.")
  }
  rgb <- srgb_compand(pmin(pmax(lin, 0), 1))
  plate$rgb <- array(rgb, dim = dim(plate$rgb))
  plate$lab <- srgb_to_lab(plate$rgb)
  plate$calibration <- list(method = method, gain = gain, ref_lab = ref_lab)
  plate
}
