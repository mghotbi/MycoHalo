# MycoHalo 0.1.0

Complete rewrite. The earlier grayscale-threshold prototype could not
separate grey *Zymoseptoria* colonies from the bacterial colony and its halo
and did not separate several cultures in one dish.

* Colour-calibrated analysis in CIELAB with flat-field illumination
  correction and optional grey-card calibration (`calibrate_color()`,
  `analyze_plate(reference = )`).
* Texture-aware object detection, so that strongly melanized colonies on
  dark agar are found.
* Fungus / bacterium / halo classification by a layout-initialised,
  noise-aware colour-direction mixture (`classify_pixels()`), or by a
  supervised classifier (`train_classifier()`, `annotate_plate()`).
* Separation of cultures sharing one dish: plate layouts (`plate_layout()`),
  distance-transform colony cores, seeded watershed, satellite detection and
  half-contrast edge placement (`segment_colonies()`).
* Measurements: size, shape, melanization index with spatial block-bootstrap
  CIs, ImageJ-compatible grey value, side-specific melanization facing the
  bacterium, growth inhibition (PIRG), gaps to bacterium and halo, halo
  geometry and colour (`measure_colonies()`), edge-to-centre profiles
  (`colony_profiles()`).
* `plot_melanization_map()`: 3D shaded melanization landscape and 2D MI map.
* `make_metadata()`, `layout_by_treatment()` and `correct_facing_bias()` for
  experiment-level workflows; cast-shadow detection; local half-contrast
  edge placement; semi-supervised (seed-clamped) EM; calibrated on real
  Zymoseptoria control and confrontation photographs.
* Batch processing with metadata and per-plate layouts (`analyze_plates()`),
  QC figures (`plot_qc()`), profile plots (`plot_profiles()`).
* Plate-aware statistics with mixed models (`compare_melanization()`).
* Plate simulator with ground truth (`simulate_plate()`) and a validation
  suite (`inst/validation/`).
* No Bioconductor dependency: image kernels (exact Euclidean distance
  transform, connected components, watershed, separable filters) are
  implemented in C++ (Rcpp).
