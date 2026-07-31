suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(splines)
})

options(survey.lonely.psu = "adjust")
root <- normalizePath(file.path(getwd(), "work"), winslash = "/", mustWork = TRUE)
derived_dir <- file.path(root, "data", "derived")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
periods <- list(
  pre = readRDS(file.path(derived_dir, "pre_eligible.rds")),
  post = readRDS(file.path(derived_dir, "post_eligible.rds"))
)
period_labels <- c(
  pre = "2017-March 2020",
  post = "August 2021-August 2023"
)

make_design <- function(dat) {
  svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~diet_weight,
            nest = TRUE, data = dat)
}

extract_term <- function(model, term, exponentiate = FALSE) {
  beta <- coef(model)[term]
  se <- sqrt(diag(vcov(model)))[term]
  df <- degf(model$survey.design)
  crit <- qt(0.975, df)
  out <- data.frame(
    estimate = unname(beta),
    se_log_or_beta = unname(se),
    conf_low = unname(beta - crit * se),
    conf_high = unname(beta + crit * se),
    p_value = 2 * pt(abs(beta / se), df, lower.tail = FALSE),
    design_df = df,
    n = nobs(model)
  )
  if (exponentiate) {
    out[c("estimate", "conf_low", "conf_high")] <-
      exp(out[c("estimate", "conf_low", "conf_high")])
  }
  out
}

weighted_quantiles <- function(design, formula, probs) {
  as.numeric(unlist(svyquantile(formula, design, probs, ci = FALSE, na.rm = TRUE)))
}

descriptive_rows <- list()
self_dkr_rows <- list()
sensitivity_rows <- list()
spline_rows <- list()
spline_test_rows <- list()

for (p in names(periods)) {
  dat <- periods[[p]]
  label <- period_labels[[p]]

  alcohol_ok <- !is.na(dat$alcohol_g_day) &
    dat$alcohol_g_day <= ifelse(dat$sex == "Men", 30, 20)
  dat$masld_274 <- ifelse(
    dat$steatosis_274 == 1 & dat$any_cmrf == 1 & alcohol_ok, 1,
    ifelse(dat$steatosis_274 == 0 | dat$any_cmrf == 0 |
             dat$metald_range == 1 | dat$high_alcohol == 1, 0, NA)
  )
  dat$masld_248 <- ifelse(
    dat$steatosis_248 == 1 & dat$any_cmrf == 1 & alcohol_ok, 1,
    ifelse(dat$steatosis_248 == 0 | dat$any_cmrf == 0 |
             dat$metald_range == 1 | dat$high_alcohol == 1, 0, NA)
  )
  design <- make_design(dat)

  # Compact period-level descriptive table.
  continuous <- c(
    age = "Age, years",
    bmi = "Body mass index, kg/m2",
    energy_kcal = "Energy intake, kcal/day",
    carb_energy_pct = "Carbohydrate, % energy",
    dkr_total = "Dietary ketogenic ratio",
    cap = "Controlled attenuation parameter, dB/m"
  )
  binary <- c(
    "I(sex == 'Women')" = "Women",
    "I(race == 'Non-Hispanic White')" = "Non-Hispanic White",
    "I(education == 'College or above')" = "College education or above",
    "I(smoking == 'Current')" = "Current smoking",
    "I(meets_activity_guideline == 1)" = "Meeting activity guideline",
    masld = "MASLD (CAP >=285 dB/m)"
  )
  for (v in names(continuous)) {
    est <- svymean(as.formula(paste0("~", v)), design, na.rm = TRUE)
    descriptive_rows[[paste(p, v)]] <- data.frame(
      period = label, characteristic = continuous[[v]], type = "Mean",
      estimate = as.numeric(coef(est)[1]), se = as.numeric(SE(est)[1]),
      observed_n = sum(!is.na(dat[[v]]))
    )
  }
  for (v in names(binary)) {
    est <- svymean(as.formula(paste0("~as.numeric(", v, ")")),
                   design, na.rm = TRUE)
    descriptive_rows[[paste(p, v)]] <- data.frame(
      period = label, characteristic = binary[[v]], type = "Percent",
      estimate = 100 * as.numeric(coef(est)[1]),
      se = 100 * as.numeric(SE(est)[1]),
      observed_n = sum(!is.na(model.frame(as.formula(paste0("~", v)), dat)[[1]]))
    )
  }

  # Distribution of DKR specifically among self-reported low-carbohydrate eaters.
  low_design <- subset(design, self_lowcarb == 1)
  low_q <- weighted_quantiles(low_design, ~dkr_total, c(0.05, 0.25, 0.5, 0.75, 0.95))
  low_mean <- svymean(~dkr_total + carb_energy_pct, low_design, na.rm = TRUE)
  self_dkr_rows[[p]] <- data.frame(
    period = label,
    unweighted_n = sum(dat$self_lowcarb == 1),
    weighted_prevalence = as.numeric(coef(svymean(~self_lowcarb, design))),
    dkr_mean = as.numeric(coef(low_mean)["dkr_total"]),
    dkr_p05 = low_q[1], dkr_p25 = low_q[2], dkr_median = low_q[3],
    dkr_p75 = low_q[4], dkr_p95 = low_q[5],
    carb_energy_pct_mean = as.numeric(coef(low_mean)["carb_energy_pct"]),
    net_dkr_ge_1_n = sum(dat$self_lowcarb == 1 & dat$dkr_net >= 1),
    net_dkr_ge_1_5_n = sum(dat$self_lowcarb == 1 & dat$dkr_net >= 1.5)
  )

  # Standardize all continuous DKR variants with period-specific survey moments.
  for (v in c("dkr_total", "dkr_net", "dkr_mean_daily")) {
    mu <- as.numeric(coef(svymean(as.formula(paste0("~", v)), design, na.rm = TRUE)))
    sdv <- sqrt(as.numeric(svyvar(as.formula(paste0("~", v)), design, na.rm = TRUE)))
    dat[[paste0(v, "_z")]] <- (dat[[v]] - mu) / sdv
  }
  design <- make_design(dat)
  primary_covars <- paste(
    "age + sex + race + education + pir + smoking +",
    "meets_activity_guideline + alcohol_g_day + energy_kcal"
  )
  specs <- list(
    primary = list(outcome = "masld", exposure = "dkr_total_z",
                   extra = "", label = "Primary definition"),
    bmi_adjusted = list(outcome = "masld", exposure = "dkr_total_z",
                        extra = "+ bmi", label = "Additional adjustment for BMI"),
    waist_adjusted = list(outcome = "masld", exposure = "dkr_total_z",
                          extra = "+ waist", label = "Additional adjustment for waist circumference"),
    net_carbohydrate = list(outcome = "masld", exposure = "dkr_net_z",
                            extra = "", label = "Net-carbohydrate DKR"),
    mean_daily_ratio = list(outcome = "masld", exposure = "dkr_mean_daily_z",
                            extra = "", label = "Mean of daily DKR values"),
    cap274_masld = list(outcome = "masld_274", exposure = "dkr_total_z",
                        extra = "", label = "MASLD using CAP >=274 dB/m"),
    cap248_masld = list(outcome = "masld_248", exposure = "dkr_total_z",
                        extra = "", label = "MASLD using CAP >=248 dB/m"),
    cap_continuous = list(outcome = "cap", exposure = "dkr_total_z",
                          extra = "", label = "CAP as a continuous outcome")
  )
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    form <- as.formula(paste(
      sp$outcome, "~", sp$exposure, "+", primary_covars, sp$extra
    ))
    fam <- if (sp$outcome == "cap") gaussian() else quasibinomial()
    mod <- svyglm(form, design, family = fam, na.action = na.omit)
    tr <- extract_term(mod, sp$exposure, exponentiate = sp$outcome != "cap")
    sensitivity_rows[[paste(p, nm)]] <- cbind(
      data.frame(period = label, analysis = nm, description = sp$label,
                 outcome = sp$outcome), tr
    )
  }

  # Exclude MetALD and higher alcohol exposure from the comparison group.
  low_alcohol_design <- subset(design, alcohol_ok)
  low_alcohol_mod <- svyglm(
    as.formula(paste("masld ~ dkr_total_z +", primary_covars)),
    low_alcohol_design, family = quasibinomial(), na.action = na.omit
  )
  tr <- extract_term(low_alcohol_mod, "dkr_total_z", exponentiate = TRUE)
  sensitivity_rows[[paste(p, "exclude_higher_alcohol")]] <- cbind(
    data.frame(
      period = label, analysis = "exclude_higher_alcohol",
      description = "Excluding MetALD and higher alcohol exposure",
      outcome = "masld"
    ),
    tr
  )

  # Exclude participants with a known positive HBV surface antigen or HCV RNA.
  known_viral <- (
    (!is.na(dat$hbv_positive) & dat$hbv_positive == 1) |
      (!is.na(dat$hcv_positive) & dat$hcv_positive == 1)
  )
  viral_exclusion_design <- subset(design, !known_viral)
  viral_exclusion_mod <- svyglm(
    as.formula(paste("masld ~ dkr_total_z +", primary_covars)),
    viral_exclusion_design, family = quasibinomial(), na.action = na.omit
  )
  tr <- extract_term(viral_exclusion_mod, "dkr_total_z", exponentiate = TRUE)
  sensitivity_rows[[paste(p, "exclude_viral")]] <- cbind(
    data.frame(
      period = label, analysis = "exclude_viral",
      description = "Excluding known HBV or HCV positivity",
      outcome = "masld"
    ),
    tr
  )

  # Stabilized inverse-probability weighting for complete-case participation.
  main_vars <- c(
    "masld", "dkr_total_z", "age", "sex", "race", "education", "pir",
    "smoking", "meets_activity_guideline", "alcohol_g_day", "energy_kcal"
  )
  dat$complete_main <- complete.cases(dat[, main_vars])
  selection_design <- make_design(dat)
  selection_mod <- svyglm(
    complete_main ~ dkr_total_z + age + sex + race + energy_kcal,
    selection_design, family = quasibinomial()
  )
  p_complete <- pmin(pmax(as.numeric(predict(selection_mod, type = "response")), 0.05), 0.99)
  complete_prevalence <- as.numeric(
    coef(svymean(~as.numeric(complete_main), selection_design))
  )
  dat$ipw_weight <- ifelse(
    dat$complete_main,
    dat$diet_weight * complete_prevalence / p_complete,
    NA_real_
  )
  ipw_limits <- quantile(dat$ipw_weight[dat$complete_main], c(0.01, 0.99), na.rm = TRUE)
  dat$ipw_weight <- pmin(pmax(dat$ipw_weight, ipw_limits[1]), ipw_limits[2])
  ipw_design <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~ipw_weight,
    nest = TRUE, data = dat[dat$complete_main, ]
  )
  ipw_mod <- svyglm(
    as.formula(paste("masld ~ dkr_total_z +", primary_covars)),
    ipw_design, family = quasibinomial()
  )
  tr <- extract_term(ipw_mod, "dkr_total_z", exponentiate = TRUE)
  sensitivity_rows[[paste(p, "ipw_complete_case")]] <- cbind(
    data.frame(
      period = label, analysis = "ipw_complete_case",
      description = "Inverse-probability weighted complete-case analysis",
      outcome = "masld"
    ),
    tr
  )

  # Adjusted spline contrasts against the period-specific median DKR.
  knots <- weighted_quantiles(design, ~dkr_total, c(0.05, 0.35, 0.65, 0.95))
  basis <- ns(dat$dkr_total, knots = knots[c(2, 3)],
              Boundary.knots = knots[c(1, 4)])
  dat$s1 <- basis[, 1]
  dat$s2 <- basis[, 2]
  dat$s3 <- basis[, 3]
  spline_design <- make_design(dat)
  spline_mod <- svyglm(
    as.formula(paste("masld ~ s1 + s2 + s3 +", primary_covars)),
    spline_design, family = quasibinomial(), na.action = na.omit
  )
  # Rotate the spline basis so that the linear component is explicit; the two
  # remaining columns then provide a direct 2-df test for nonlinearity.
  residual_basis <- qr.resid(qr(cbind(1, dat$dkr_total)), basis)
  nonlinear_basis <- qr.Q(qr(residual_basis))[, 1:2, drop = FALSE]
  dat$nl1 <- nonlinear_basis[, 1]
  dat$nl2 <- nonlinear_basis[, 2]
  rotated_design <- make_design(dat)
  rotated_mod <- svyglm(
    as.formula(paste("masld ~ dkr_total + nl1 + nl2 +", primary_covars)),
    rotated_design, family = quasibinomial(), na.action = na.omit
  )
  wald_block <- function(model, terms) {
    b <- coef(model)[terms]
    v <- vcov(model)[terms, terms, drop = FALSE]
    chisq <- as.numeric(t(b) %*% solve(v, b))
    df1 <- length(terms)
    df2 <- degf(model$survey.design)
    c(statistic = chisq / df1, df1 = df1, df2 = df2,
      p_value = pf(chisq / df1, df1, df2, lower.tail = FALSE))
  }
  overall_test <- wald_block(rotated_mod, c("dkr_total", "nl1", "nl2"))
  nonlinear_test <- wald_block(rotated_mod, c("nl1", "nl2"))
  spline_test_rows[[p]] <- data.frame(
    period = label,
    p_overall = overall_test["p_value"],
    p_nonlinear = nonlinear_test["p_value"],
    overall_F = overall_test["statistic"],
    nonlinear_F = nonlinear_test["statistic"],
    design_df = overall_test["df2"]
  )
  grid <- seq(knots[1], knots[4], length.out = 121)
  bg <- ns(grid, knots = knots[c(2, 3)], Boundary.knots = knots[c(1, 4)])
  med <- weighted_quantiles(design, ~dkr_total, 0.5)
  bm <- ns(med, knots = knots[c(2, 3)], Boundary.knots = knots[c(1, 4)])
  contrast <- sweep(bg, 2, as.numeric(bm), "-")
  beta <- coef(spline_mod)[c("s1", "s2", "s3")]
  vv <- vcov(spline_mod)[c("s1", "s2", "s3"), c("s1", "s2", "s3")]
  log_or <- as.numeric(contrast %*% beta)
  se <- sqrt(rowSums((contrast %*% vv) * contrast))
  crit <- qt(0.975, degf(spline_mod$survey.design))
  spline_rows[[p]] <- data.frame(
    period = label, dkr = grid, odds_ratio = exp(log_or),
    conf_low = exp(log_or - crit * se), conf_high = exp(log_or + crit * se),
    reference_dkr = med
  )
}

descriptive <- bind_rows(descriptive_rows)
self_dkr <- bind_rows(self_dkr_rows)
sensitivity <- bind_rows(sensitivity_rows)
spline <- bind_rows(spline_rows)
spline_tests <- bind_rows(spline_test_rows)

# Cross-period heterogeneity, using independent period-specific estimates.
heterogeneity <- sensitivity |>
  filter(analysis %in% c("primary", "cap_continuous")) |>
  mutate(log_or_or_beta = ifelse(outcome == "cap", estimate, log(estimate))) |>
  select(period, analysis, outcome, log_or_or_beta, se_log_or_beta) |>
  tidyr::pivot_wider(
    names_from = period,
    values_from = c(log_or_or_beta, se_log_or_beta),
    names_sep = "__"
  )
pre_name <- period_labels[["pre"]]
post_name <- period_labels[["post"]]
heterogeneity$z_difference <- (
  heterogeneity[[paste0("log_or_or_beta__", pre_name)]] -
    heterogeneity[[paste0("log_or_or_beta__", post_name)]]
) / sqrt(
  heterogeneity[[paste0("se_log_or_beta__", pre_name)]]^2 +
    heterogeneity[[paste0("se_log_or_beta__", post_name)]]^2
)
heterogeneity$p_heterogeneity <- 2 * pnorm(
  abs(heterogeneity$z_difference), lower.tail = FALSE
)

write.csv(descriptive, file.path(table_dir, "table1_period_characteristics.csv"),
          row.names = FALSE)
write.csv(self_dkr, file.path(table_dir, "self_lowcarb_dkr_distribution.csv"),
          row.names = FALSE)
write.csv(sensitivity, file.path(table_dir, "sensitivity_associations.csv"),
          row.names = FALSE)
write.csv(spline, file.path(table_dir, "spline_predictions.csv"),
          row.names = FALSE)
write.csv(spline_tests, file.path(table_dir, "spline_tests.csv"),
          row.names = FALSE)
write.csv(heterogeneity, file.path(table_dir, "period_heterogeneity.csv"),
          row.names = FALSE)

p <- ggplot2::ggplot(
  spline,
  ggplot2::aes(x = dkr, y = odds_ratio, colour = period, fill = period)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf_low, ymax = conf_high),
    alpha = 0.12, colour = NA
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::labs(
    x = "Dietary ketogenic ratio",
    y = "Adjusted odds ratio for MASLD",
    colour = NULL, fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(figure_dir, "figure_spline.png"), p,
                 width = 7.2, height = 4.8, dpi = 400)

flow <- read.csv(file.path(table_dir, "sample_flow.csv"), check.names = FALSE)
flow_long <- bind_rows(lapply(seq_len(nrow(flow)), function(i) {
  data.frame(
    period = flow$period[i],
    stage = factor(
      c("All examined", "Age >=20 years", "Two reliable recalls",
        "Complete VCTE", "Eligible analysis sample"),
      levels = rev(c("All examined", "Age >=20 years", "Two reliable recalls",
                     "Complete VCTE", "Eligible analysis sample"))
    ),
    n = as.numeric(flow[i, c("total", "adult", "two_recall", "valid_vcte", "eligible")])
  )
}))
flow_plot <- ggplot2::ggplot(
  flow_long,
  ggplot2::aes(x = period, y = stage, group = period)
) +
  ggplot2::geom_line(colour = "#AAB7C4", linewidth = 0.7) +
  ggplot2::geom_label(
    ggplot2::aes(label = paste0(stage, "\n", format(n, big.mark = ","))),
    size = 3.0, linewidth = 0.25, fill = "white", colour = "#1F4D78"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    plot.margin = ggplot2::margin(10, 15, 10, 15),
    axis.text.x = ggplot2::element_text()
  )
ggplot2::ggsave(
  file.path(figure_dir, "figure_s1_selection_flow.png"),
  flow_plot, width = 7.2, height = 5.8, dpi = 400
)

message("Sensitivity analyses and manuscript tables complete.")
