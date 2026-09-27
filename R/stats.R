#' Compare melanization between treatments, respecting plate structure
#'
#' The four colonies on one dish share medium, bacterium, handling and
#' photograph; they are *not* independent replicates. Treating them as such
#' (pseudoreplication, Hurlbert 1984) inflates the evidence. This function
#' therefore fits a linear mixed model with a random intercept per plate
#' (\pkg{lme4}), or, if \pkg{lme4} is not installed or the design has too
#' few plates, analyses plate means with Welch's t-test / one-way ANOVA.
#'
#' @param colonies Colony table (e.g. `out$colonies` from [analyze_plates()])
#'   containing the response, the group variable and `plate_id`.
#' @param response Response column (default `"MI_mean"`).
#' @param group Treatment column.
#' @param plate Plate identifier column.
#' @param reference Optional reference level of `group` (e.g. "control").
#' @param method `"auto"`, `"mixed"` or `"plate_means"`.
#' @return A list with the fitted model, a coefficient table with 95 % CIs,
#'   the intra-plate correlation (ICC) and per-group descriptive statistics.
#' @examples
#' set.seed(1)
#' d <- data.frame(plate_id = rep(sprintf("p%02d", 1:10), each = 4),
#'                 treatment = rep(c("control", "bacteria"), each = 20))
#' d$MI_mean <- 50 + 3 * (d$treatment == "bacteria") +
#'   rep(rnorm(10, 0, 1.5), each = 4) + rnorm(40, 0, 1)
#' compare_melanization(d, group = "treatment", reference = "control")
#' @export
compare_melanization <- function(colonies, response = "MI_mean", group = "treatment",
                                 plate = "plate_id", reference = NULL,
                                 method = c("auto", "mixed", "plate_means")) {
  method <- match.arg(method)
  df <- colonies[!is.na(colonies[[response]]), , drop = FALSE]
  df$.y <- df[[response]]
  df$.g <- factor(df[[group]])
  if (!is.null(reference)) df$.g <- stats::relevel(df$.g, reference)
  df$.p <- factor(df[[plate]])
  n_plates <- tapply(df$.p, df$.g, function(z) length(unique(z)))
  desc <- do.call(rbind, lapply(split(df, df$.g), function(s) {
    pm <- tapply(s$.y, s$.p, mean)
    data.frame(group = s$.g[1], n_plates = length(pm), n_colonies = nrow(s),
               mean = mean(s$.y), sd_between_plates = stats::sd(pm),
               sd_within_plates = sqrt(mean(tapply(s$.y, s$.p, stats::var), na.rm = TRUE)))
  }))
  use_mixed <- method == "mixed" ||
    (method == "auto" && requireNamespace("lme4", quietly = TRUE) && min(n_plates) >= 3)
  if (use_mixed) {
    if (!requireNamespace("lme4", quietly = TRUE)) cli::cli_abort("Install {.pkg lme4}.")
    fit <- lme4::lmer(.y ~ .g + (1 | .p), data = df, REML = TRUE)
    co <- summary(fit)$coefficients
    ci <- suppressMessages(stats::confint(fit, method = "Wald", parm = "beta_"))
    vc <- as.data.frame(lme4::VarCorr(fit))
    icc <- vc$vcov[vc$grp == ".p"] / sum(vc$vcov)
    tab <- data.frame(term = sub("^\\.g", "", rownames(co)), estimate = co[, 1],
                      se = co[, 2], t = co[, 3], lower = ci[, 1], upper = ci[, 2],
                      row.names = NULL)
    tab$df_approx <- nlevels(df$.p) - nlevels(df$.g)
    tab$p_value <- 2 * stats::pt(-abs(tab$t), tab$df_approx)
    list(method = "linear mixed model: response ~ group + (1 | plate)",
         coefficients = tab, icc = icc, descriptives = desc, model = fit,
         note = "p-values use a conservative between-plate df approximation (n_plates - n_groups).")
  } else {
    pm <- stats::aggregate(.y ~ .p + .g, data = df, FUN = mean)
    if (nlevels(df$.g) == 2) {
      tt <- stats::t.test(.y ~ .g, data = pm)
      # t.test reports level1 - level2; we report level2 - reference (level1)
      tab <- data.frame(term = levels(df$.g)[2],
                        estimate = unname(tt$estimate[2] - tt$estimate[1]),
                        lower = -tt$conf.int[2], upper = -tt$conf.int[1],
                        t = -unname(tt$statistic), df = unname(tt$parameter),
                        p_value = tt$p.value, row.names = NULL)
    } else {
      fit <- stats::lm(.y ~ .g, data = pm)
      tab <- as.data.frame(stats::anova(fit))
    }
    list(method = "Welch test / ANOVA on plate means", coefficients = tab,
         icc = NA_real_, descriptives = desc, model = NULL)
  }
}
