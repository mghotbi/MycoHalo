## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 9, fig.height = 4.7,
                      dpi = 60, dev = "jpeg", out.width = "100%")
library(MycoHalo)

## ----quick--------------------------------------------------------------------
# A simulated confrontation plate (replace with your file path)
sim <- simulate_plate(seed = 42, facing_darkening = 4)
res <- analyze_plate(sim$image, id = "demo", verbose = FALSE)
res

## ----qc, fig.cap = "Quality-control figure. Left: detected dish (dashed), analysed agar disk (solid), expected inoculation points (+) with their search radius (dotted), colony outlines, bacterium (yellow), halo (orange), satellite (magenta). Right: pixel classes."----
plot_qc(res)

## ----real, eval = FALSE-------------------------------------------------------
#  # confrontation plate: 4 Zymoseptoria + central bacterium (default layout)
#  res <- analyze_plate("plates/2026_05_26_20.JPG",
#                       reference = c(0.08, 0.92, 0.03))   # grey card, bottom-left
#  # control plate: same layout, no bacterium
#  ctrl <- analyze_plate("plates/2026_05_26_37.JPG",
#                        layout = plate_layout(centre = "none"))

## ----real-plates, fig.width = 10, fig.height = 5.6----------------------------
treat <- analyze_plate(example_plate("confrontation"), verbose = FALSE)
ctrl  <- analyze_plate(example_plate("control"),
                       layout = plate_layout(centre = "none"), verbose = FALSE)
rbind(cbind(plate = "bacteria", treat$colonies[, c("colony_id", "area_mm2", "MI_mean")]),
      cbind(plate = "control",  ctrl$colonies[,  c("colony_id", "area_mm2", "MI_mean")]))
plot_qc(treat)
plot_melanization_map(treat, mi_range = c(35, 85), main = "Zymoseptoria + bacterium")
plot_melanization_map(ctrl,  mi_range = c(35, 85), main = "Zymoseptoria control")

## ----metadata, eval = FALSE---------------------------------------------------
#  # 1. create the sheet (file names containing "ctrl" become controls)
#  make_metadata("photos",
#                treatment = "bacteria",
#                patterns  = c(control = "ctrl|control"),
#                extra     = list(bacterium = "isolate_X", zymo_strain = "strain_Y", day = 7))
#  # -> photos/plate_metadata.csv:
#  #    file, plate_id, treatment, bacterium, zymo_strain, medium, day, replicate, photo_date, notes
#  
#  # 2. fix / complete it in Excel or in R, e.g.
#  meta <- read.csv("photos/plate_metadata.csv")
#  meta$treatment[meta$file == "2026_05_26_37.JPG"] <- "control"
#  write.csv(meta, "photos/plate_metadata.csv", row.names = FALSE, na = "")
#  
#  # 3. analyse everything; the layout follows the treatment column
#  out <- analyze_plates("photos", metadata = meta, layout = layout_by_treatment,
#                        qc_dir = "qc", reference = c(0.08, 0.92, 0.03))
#  out$colonies <- correct_facing_bias(out$colonies)   # lighting bias from controls
#  
#  write.csv(out$colonies, "colonies.csv", row.names = FALSE)
#  compare_melanization(out$colonies, group = "treatment", reference = "control")

## ----layouts------------------------------------------------------------------
plate_layout()                                   # default: TL, TR, BR, BL + bacterium
plate_layout(centre = "none")                    # controls
plate_layout(radius_mm = 20)                     # inoculated 20 mm from the centre
plate_layout(n_colonies = 3, ids = c("A", "B", "C"), start_deg = -90)

## ----supervised, eval = FALSE-------------------------------------------------
#  p   <- read_plate("plate.jpg") |> detect_plate() |> model_background()
#  reg <- annotate_plate(p)                       # click fungus, bacteria, halo regions
#  clf <- train_classifier(extract_training_pixels(p, reg))
#  res <- analyze_plate("plate.jpg", classifier = clf)

## ----profiles, fig.height = 3.5-----------------------------------------------
sim2 <- simulate_plate(seed = 7, edge_lightening = 10)
res2 <- analyze_plate(sim2$image, verbose = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE)) plot_profiles(res2)

## ----landscape, fig.width = 10, fig.height = 5.6, fig.cap = "Melanization landscape of a simulated confrontation plate in which each colony is 5 L* units darker on the side facing the bacterium."----
sim3 <- simulate_plate(width = 900, height = 1200, seed = 42,
                       facing_darkening = 5, edge_lightening = 8)
res3 <- analyze_plate(sim3$image, id = "simulated plate", verbose = FALSE)
plot_melanization_map(res3)

## ----stats--------------------------------------------------------------------
set.seed(1)
d <- data.frame(plate_id = rep(sprintf("p%02d", 1:12), each = 4),
                treatment = rep(c("control", "bacteria"), each = 24))
d$MI_mean <- 55 + 3 * (d$treatment == "bacteria") +
  rep(rnorm(12, 0, 1.5), each = 4) + rnorm(48, 0, 1)
cmp <- compare_melanization(d, response = "MI_mean", group = "treatment",
                            reference = "control")
cmp$method
cmp$coefficients
cmp$icc      # share of variance between plates

## ----validation---------------------------------------------------------------
v <- read.csv(system.file("validation", "simulation_validation.csv", package = "MycoHalo"))
agg <- do.call(rbind, lapply(split(v, v$scenario), function(s) data.frame(
  scenario = s$scenario[1],
  detected = paste0(sum(s$detected), "/", nrow(s)),
  max_abs_area_err_pct = max(abs(s$area_err_pct), na.rm = TRUE),
  max_abs_MI_err = max(abs(s$MI_err), na.rm = TRUE),
  mean_delta_MI_facing = round(mean(s$dMI_facing, na.rm = TRUE), 2))))
agg[order(agg$max_abs_MI_err), ]

