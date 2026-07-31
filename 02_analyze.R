suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(ggplot2)
  library(splines)
  library(broom)
})

options(survey.lonely.psu = "adjust")

root <- normalizePath(file.path(getwd(), "work"), winslash = "/", mustWork = TRUE)
derived_dir <- file.path(root, "data", "derived")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

periods <- list(
  pre = readRDS(file.path(derived_dir, "pre_eligible.rds")),
  post = readRDS(file.path(derived_dir, "post_eligible.rds"))
)
period_labels <- c(
  pre = "2017-March 2020",
  post = "August 2021-August 2023"
)

make_design <- function(dat) {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~diet_weight,
    nest = TRUE,
    data = dat
  )
}

weighted_sd <- function(design, variable) {
  sqrt(as.numeric(svyvar(as.formula(paste0("~", variable)), design, na.rm = TRUE)))
}

fmt_ci <- function(est, low, high, digits = 2) {
  paste0(
    formatC(est, format = "f", digits = digits), " (",
    formatC(low, format = "f", digits = digits), ", ",
    formatC(high, format = "f", digits = digits), ")"
  )
}

weighted_auc <- function(score, outcome, weight) {
  ok <- complete.cases(score, outcome, weight) & weight > 0
  score <- score[ok]
  outcome <- outcome[ok]
  weight <- weight[ok]
  ord <- order(score)
  score <- score[ord]
  outcome <- outcome[ord]
  weight <- weight[ord]
  pos_w <- weight * (outcome == 1)
  neg_w <- weight * (outcome == 0)
  groups <- split(seq_along(score), score)
  neg_before <- 0
  concordant <- 0
  for (idx in groups) {
    gp <- sum(pos_w[idx])
    gn <- sum(neg_w[idx])
    concordant <- concordant + gp * neg_before + 0.5 * gp * gn
    neg_before <- neg_before + gn
  }
  concordant / (sum(pos_w) * sum(neg_w))
}

weighted_kappa <- function(a, b, w) {
  ok <- complete.cases(a, b, w) & w > 0
  a <- as.integer(a[ok])
  b <- as.integer(b[ok])
  w <- w[ok]
  tab <- matrix(0, 2, 2)
  for (i in 0:1) for (j in 0:1) {
    tab[i + 1, j + 1] <- sum(w[a == i & b == j]) / sum(w)
  }
  po <- sum(diag(tab))
  pe <- sum(rowSums(tab) * colSums(tab))
  c(kappa = (po - pe) / (1 - pe), agreement = po,
    pabak = 2 * po - 1,
    prevalence_index = abs(tab[1, 1] - tab[2, 2]),
    bias_index = abs(tab[2, 1] - tab[1, 2]),
    positive_agreement = 2 * tab[2, 2] / (2 * tab[2, 2] + tab[2, 1] + tab[1, 2]),
    negative_agreement = 2 * tab[1, 1] / (2 * tab[1, 1] + tab[2, 1] + tab[1, 2]))
}

weighted_correlation <- function(x, y, w) {
  ok <- complete.cases(x, y, w) & w > 0
  x <- x[ok]
  y <- y[ok]
  w <- w[ok] / sum(w[ok])
  mx <- sum(w * x)
  my <- sum(w * y)
  covariance <- sum(w * (x - mx) * (y - my))
  variance_x <- sum(w * (x - mx)^2)
  variance_y <- sum(w * (y - my)^2)
  c(day1_day2_correlation = covariance / sqrt(variance_x * variance_y))
}

replicate_interval <- function(rep_design, statistic, df) {
  fit <- withReplicates(rep_design, statistic)
  estimate <- as.numeric(coef(fit))
  se <- as.numeric(SE(fit))
  names(estimate) <- names(coef(fit))
  names(se) <- names(coef(fit))
  crit <- qt(0.975, df)
  list(
    estimate = estimate,
    se = se,
    low = estimate - crit * se,
    high = estimate + crit * se
  )
}

tidy_svy <- function(model, term, exponentiate = FALSE) {
  beta <- coef(model)
  if (!term %in% names(beta)) stop("Term not found: ", term)
  se <- sqrt(diag(vcov(model)))[term]
  df <- degf(model$survey.design)
  crit <- qt(0.975, df = df)
  row <- data.frame(
    term = term,
    estimate = unname(beta[term]),
    std.error = unname(se),
    statistic = unname(beta[term] / se),
    p.value = 2 * pt(abs(beta[term] / se), df = df, lower.tail = FALSE),
    conf.low = unname(beta[term] - crit * se),
    conf.high = unname(beta[term] + crit * se)
  )
  if (exponentiate) {
    row$estimate <- exp(row$estimate)
    row$conf.low <- exp(row$conf.low)
    row$conf.high <- exp(row$conf.high)
  }
  row
}

model_formula <- function(outcome, exposure, adjusted = TRUE) {
  base <- paste(outcome, "~", exposure, "+ age + sex + race")
  if (!adjusted) return(as.formula(base))
  as.formula(paste(
    base,
    "+ education + pir + smoking + meets_activity_guideline",
    "+ alcohol_g_day + energy_kcal"
  ))
}

distribution_rows <- list()
agreement_rows <- list()
association_rows <- list()
cutpoint_rows <- list()
performance_rows <- list()
analysis_counts <- list()
repeatability_rows <- list()

for (p in names(periods)) {
  dat <- periods[[p]]
  # The published 1.5 threshold is defined for net-carbohydrate DKR.
  dat$dkr_net_ge_1 <- as.integer(dat$dkr_net >= 1)
  dat$dkr_net_ge_1_5 <- as.integer(dat$dkr_net >= 1.5)
  dat$carb_lt_10 <- as.integer(dat$carb_energy_pct < 10)
  dat$carb_lt_20 <- as.integer(dat$carb_energy_pct < 20)
  dat$carb_lt_26 <- as.integer(dat$carb_energy_pct < 26)
  design <- make_design(dat)
  label <- period_labels[[p]]

  # Weighted exposure distribution.
  q_probs <- c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
  qs <- as.numeric(unlist(
    svyquantile(~dkr_total, design, quantiles = q_probs, ci = FALSE, na.rm = TRUE)
  ))
  mean_dkr <- svymean(~dkr_total, design, na.rm = TRUE)
  prevalence_vars <- c(
    "self_lowcarb" = "self_lowcarb",
    "dkr_net_ge_1" = "dkr_net_ge_1",
    "dkr_net_ge_1_5" = "dkr_net_ge_1_5",
    "carb_lt_10" = "carb_lt_10",
    "carb_lt_20" = "carb_lt_20",
    "carb_lt_26" = "carb_lt_26"
  )
  prev <- lapply(names(prevalence_vars), function(nm) {
    expr <- unname(prevalence_vars[nm])
    est <- svymean(as.formula(paste0("~", expr)), design, na.rm = TRUE)
    data.frame(
      measure = nm,
      estimate = as.numeric(coef(est)[1]),
      se = as.numeric(SE(est)[1])
    )
  }) |> bind_rows()
  distribution_rows[[p]] <- bind_rows(
    data.frame(
      period = label, measure = "dkr_mean",
      estimate = as.numeric(coef(mean_dkr)), se = as.numeric(SE(mean_dkr))
    ),
    data.frame(
      period = label,
      measure = paste0("dkr_p", q_probs * 100),
      estimate = qs, se = NA_real_
    ),
    transform(prev, period = label)
  )

  # Cross-proxy agreement.
  design <- update(
    design,
    dkr_net_ge_1 = as.integer(dkr_net >= 1),
    dkr_net_ge_1_5 = as.integer(dkr_net >= 1.5),
    carb_lt_20 = as.integer(carb_energy_pct < 20),
    carb_lt_26 = as.integer(carb_energy_pct < 26)
  )
  rep_design <- as.svrepdesign(design, type = "JKn", mse = TRUE)
  repeatability <- replicate_interval(
    rep_design,
    function(w, data) weighted_correlation(data$dkr_day1, data$dkr_day2, w),
    degf(design)
  )
  repeatability_rows[[p]] <- data.frame(
    period = label,
    measure = "Survey-weighted Pearson correlation between day-specific DKR values",
    n = sum(complete.cases(dat$dkr_day1, dat$dkr_day2)),
    estimate = repeatability$estimate["day1_day2_correlation"],
    se = repeatability$se["day1_day2_correlation"],
    conf_low = repeatability$low["day1_day2_correlation"],
    conf_high = repeatability$high["day1_day2_correlation"]
  )
  by_lowcarb <- svyby(
    ~dkr_total + carb_energy_pct,
    ~self_lowcarb,
    design,
    svymean,
    na.rm = TRUE,
    vartype = c("se", "ci")
  )
  write.csv(
    by_lowcarb,
    file.path(table_dir, paste0("proxy_means_", p, ".csv")),
    row.names = FALSE
  )

  for (proxy in c("dkr_net_ge_1", "dkr_net_ge_1_5", "carb_lt_20", "carb_lt_26")) {
    kap <- replicate_interval(
      rep_design,
      function(w, data) weighted_kappa(data$self_lowcarb, data[[proxy]], w),
      degf(design)
    )
    agreement_rows[[paste(p, proxy)]] <- data.frame(
      period = label,
      proxy = proxy,
      unweighted_self_lowcarb = sum(dat$self_lowcarb == 1, na.rm = TRUE),
      unweighted_proxy_positive = sum(dat[[proxy]] == 1, na.rm = TRUE),
      kappa = kap$estimate["kappa"],
      kappa_se = kap$se["kappa"],
      kappa_low = kap$low["kappa"],
      kappa_high = kap$high["kappa"],
      overall_agreement = kap$estimate["agreement"],
      overall_agreement_low = kap$low["agreement"],
      overall_agreement_high = kap$high["agreement"],
      pabak = kap$estimate["pabak"],
      pabak_low = kap$low["pabak"],
      pabak_high = kap$high["pabak"],
      prevalence_index = kap$estimate["prevalence_index"],
      prevalence_index_low = kap$low["prevalence_index"],
      prevalence_index_high = kap$high["prevalence_index"],
      bias_index = kap$estimate["bias_index"],
      bias_index_low = kap$low["bias_index"],
      bias_index_high = kap$high["bias_index"],
      positive_agreement = kap$estimate["positive_agreement"],
      positive_agreement_low = kap$low["positive_agreement"],
      positive_agreement_high = kap$high["positive_agreement"],
      negative_agreement = kap$estimate["negative_agreement"],
      auc = NA_real_, auc_se = NA_real_, auc_low = NA_real_, auc_high = NA_real_
    )
  }
  auc_fit <- replicate_interval(
    rep_design,
    function(w, data) weighted_auc(data$dkr_total, data$self_lowcarb, w),
    degf(design)
  )
  agreement_rows[[paste(p, "auc")]] <- data.frame(
    period = label,
    proxy = "DKR continuous predicting self-reported low-carbohydrate diet",
    unweighted_self_lowcarb = sum(dat$self_lowcarb == 1, na.rm = TRUE),
    unweighted_proxy_positive = NA_integer_,
    kappa = NA_real_, kappa_se = NA_real_, kappa_low = NA_real_, kappa_high = NA_real_,
    overall_agreement = NA_real_,
    overall_agreement_low = NA_real_, overall_agreement_high = NA_real_,
    pabak = NA_real_, pabak_low = NA_real_, pabak_high = NA_real_,
    prevalence_index = NA_real_, prevalence_index_low = NA_real_, prevalence_index_high = NA_real_,
    bias_index = NA_real_, bias_index_low = NA_real_, bias_index_high = NA_real_,
    positive_agreement = NA_real_,
    positive_agreement_low = NA_real_, positive_agreement_high = NA_real_,
    negative_agreement = NA_real_,
    auc = auc_fit$estimate[1], auc_se = auc_fit$se[1],
    auc_low = auc_fit$low[1], auc_high = auc_fit$high[1]
  )

  # Analysis counts before modeling.
  model_vars <- c(
    "masld", "dkr_total", "age", "sex", "race", "education", "pir",
    "smoking", "meets_activity_guideline", "alcohol_g_day", "energy_kcal",
    "diet_weight", "SDMVPSU", "SDMVSTRA"
  )
  complete_main <- complete.cases(dat[, model_vars])
  analysis_counts[[p]] <- data.frame(
    period = label,
    eligible_n = nrow(dat),
    masld_observed_n = sum(!is.na(dat$masld)),
    masld_cases = sum(dat$masld == 1, na.rm = TRUE),
    complete_main_n = sum(complete_main),
    complete_main_cases = sum(dat$masld[complete_main] == 1)
  )

  # Standardized DKR is computed using survey-weighted moments.
  dkr_mean <- as.numeric(coef(svymean(~dkr_total, design, na.rm = TRUE)))
  dkr_sd <- weighted_sd(design, "dkr_total")
  dat$dkr_z <- (dat$dkr_total - dkr_mean) / dkr_sd

  # Weighted quartiles and median split.
  qcuts <- as.numeric(unlist(
    svyquantile(~dkr_total, design, c(0.25, 0.5, 0.75), ci = FALSE, na.rm = TRUE)
  ))
  dat$dkr_q <- cut(
    dat$dkr_total,
    breaks = c(-Inf, qcuts, Inf),
    labels = c("Q1", "Q2", "Q3", "Q4"),
    include.lowest = TRUE
  )
  dat$dkr_median_high <- as.integer(dat$dkr_total >= qcuts[2])
  design <- make_design(dat)

  specifications <- list(
    continuous = "dkr_z",
    quartile = "dkr_q",
    median_split = "dkr_median_high"
  )
  for (outcome in c("masld", "steatosis_285", "cap")) {
    family <- if (outcome == "cap") gaussian() else quasibinomial()
    for (spec in names(specifications)) {
      exposure <- specifications[[spec]]
      f <- model_formula(outcome, exposure, adjusted = TRUE)
      mod <- svyglm(f, design = design, family = family, na.action = na.omit)
      terms <- if (spec == "quartile") paste0("dkr_q", c("Q2", "Q3", "Q4")) else exposure
      for (term in terms) {
        tr <- tidy_svy(mod, term, exponentiate = outcome != "cap")
        association_rows[[paste(p, outcome, spec, term)]] <- data.frame(
          period = label, outcome = outcome, specification = spec, term = term,
          estimate = tr$estimate, conf_low = tr$conf.low, conf_high = tr$conf.high,
          p_value = tr$p.value, n = nobs(mod)
        )
      }
      if (outcome == "masld") {
        pred <- predict(mod, type = "response")
        mf <- model.frame(mod)
        obs <- model.response(mf)
        w <- weights(mod, "prior")
        performance_rows[[paste(p, spec)]] <- data.frame(
          period = label, specification = spec, n = nobs(mod),
          weighted_brier = weighted.mean((obs - pred)^2, w, na.rm = TRUE)
        )
      }
    }
  }

  # Cutpoint scan: same model, cutpoints fixed by weighted percentiles.
  scan_probs <- seq(0.20, 0.80, by = 0.02)
  scan_cuts <- as.numeric(unlist(
    svyquantile(~dkr_total, design, scan_probs, ci = FALSE, na.rm = TRUE)
  ))
  for (i in seq_along(scan_probs)) {
    dat$scan_high <- as.integer(dat$dkr_total >= scan_cuts[i])
    scan_design <- make_design(dat)
    mod <- svyglm(
      model_formula("masld", "scan_high", adjusted = TRUE),
      design = scan_design,
      family = quasibinomial(),
      na.action = na.omit
    )
    tr <- tidy_svy(mod, "scan_high", exponentiate = TRUE)
    log_estimate <- as.numeric(coef(mod)["scan_high"])
    log_se <- sqrt(as.numeric(vcov(mod)["scan_high", "scan_high"]))
    cutpoint_rows[[paste(p, i)]] <- data.frame(
      period = label,
      percentile = scan_probs[i],
      cutpoint = scan_cuts[i],
      high_unweighted_n = sum(dat$scan_high == 1, na.rm = TRUE),
      estimate = tr$estimate, conf_low = tr$conf.low, conf_high = tr$conf.high,
      log_estimate = log_estimate, log_se = log_se,
      p_value = tr$p.value, n = nobs(mod)
    )
  }

  # Restricted cubic spline predictions for the main MASLD model.
  spline_knots <- as.numeric(unlist(
    svyquantile(~dkr_total, design, c(0.05, 0.35, 0.65, 0.95),
                ci = FALSE, na.rm = TRUE)
  ))
  dat$s1 <- ns(
    dat$dkr_total,
    knots = spline_knots[c(2, 3)],
    Boundary.knots = spline_knots[c(1, 4)]
  )[, 1]
  dat$s2 <- ns(
    dat$dkr_total,
    knots = spline_knots[c(2, 3)],
    Boundary.knots = spline_knots[c(1, 4)]
  )[, 2]
  dat$s3 <- ns(
    dat$dkr_total,
    knots = spline_knots[c(2, 3)],
    Boundary.knots = spline_knots[c(1, 4)]
  )[, 3]
  spline_design <- make_design(dat)
  spline_model <- svyglm(
    model_formula("masld", "s1 + s2 + s3", adjusted = TRUE),
    design = spline_design,
    family = quasibinomial(),
    na.action = na.omit
  )
  saveRDS(
    list(model = spline_model, knots = spline_knots, dkr_mean = dkr_mean, dkr_sd = dkr_sd),
    file.path(derived_dir, paste0("spline_model_", p, ".rds"))
  )
}

distribution <- bind_rows(distribution_rows)
agreement <- bind_rows(agreement_rows)
associations <- bind_rows(association_rows)
cutpoints <- bind_rows(cutpoint_rows)
cutpoints <- cutpoints |>
  group_by(period) |>
  mutate(
    p_holm = p.adjust(p_value, method = "holm"),
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_fdr = p.adjust(p_value, method = "BH")
  ) |>
  ungroup()
cutpoint_heterogeneity <- cutpoints |>
  select(period, percentile, log_estimate, log_se) |>
  tidyr::pivot_wider(
    names_from = period,
    values_from = c(log_estimate, log_se),
    names_sep = "__"
  )
pre_label <- period_labels[["pre"]]
post_label <- period_labels[["post"]]
cutpoint_heterogeneity <- cutpoint_heterogeneity |>
  mutate(
    z_period_difference = (
      .data[[paste0("log_estimate__", post_label)]] -
        .data[[paste0("log_estimate__", pre_label)]]
    ) / sqrt(
      .data[[paste0("log_se__", post_label)]]^2 +
        .data[[paste0("log_se__", pre_label)]]^2
    ),
    p_period_interaction = 2 * pnorm(abs(z_period_difference), lower.tail = FALSE),
    p_period_interaction_holm = p.adjust(p_period_interaction, method = "holm")
  ) |>
  select(percentile, z_period_difference, p_period_interaction,
         p_period_interaction_holm)
cutpoints <- cutpoints |>
  left_join(cutpoint_heterogeneity, by = "percentile")
performance <- bind_rows(performance_rows)
counts <- bind_rows(analysis_counts)
repeatability <- bind_rows(repeatability_rows)

write.csv(distribution, file.path(table_dir, "weighted_dkr_distribution.csv"), row.names = FALSE)
write.csv(agreement, file.path(table_dir, "proxy_agreement.csv"), row.names = FALSE)
write.csv(associations, file.path(table_dir, "main_associations.csv"), row.names = FALSE)
write.csv(cutpoints, file.path(table_dir, "cutpoint_scan.csv"), row.names = FALSE)
write.csv(performance, file.path(table_dir, "model_performance.csv"), row.names = FALSE)
write.csv(counts, file.path(table_dir, "analysis_counts.csv"), row.names = FALSE)
write.csv(repeatability, file.path(table_dir, "dkr_day_repeatability.csv"), row.names = FALSE)

# Figure 1: exposure distributions, scaled by period-specific weights.
plot_data <- bind_rows(lapply(names(periods), function(p) {
  dat <- periods[[p]]
  dat$plot_weight <- dat$diet_weight / sum(dat$diet_weight, na.rm = TRUE)
  dat$period_label <- period_labels[[p]]
  dat
}))
p1 <- ggplot(
  plot_data,
  aes(x = dkr_net, weight = plot_weight, colour = period_label, fill = period_label)
) +
  geom_density(alpha = 0.16, linewidth = 0.9, adjust = 1.1) +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "#7A0019") +
  coord_cartesian(xlim = c(0, 1.6)) +
  labs(
    x = "Net-carbohydrate dietary ketogenic ratio",
    y = "Survey-weighted density",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(figure_dir, "figure_dkr_distribution.png"), p1,
       width = 7.2, height = 4.6, dpi = 400)

# Figure 2: cutpoint scan.
p2 <- ggplot(cutpoints, aes(x = cutpoint, y = estimate, colour = period)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_ribbon(
    aes(ymin = conf_low, ymax = conf_high, fill = period),
    colour = NA, alpha = 0.12
  ) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  labs(
    x = "Candidate DKR cutpoint (weighted 20th-80th percentiles)",
    y = "Adjusted odds ratio for MASLD (high vs. low DKR)",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(figure_dir, "figure_cutpoint_scan.png"), p2,
       width = 7.2, height = 4.8, dpi = 400)

# Figure 3: main continuous associations by period and outcome.
forest <- associations |>
  filter(specification == "continuous") |>
  mutate(
    outcome_label = recode(
      outcome,
      masld = "MASLD",
      steatosis_285 = "Hepatic steatosis",
      cap = "CAP (dB/m)"
    )
  )
p3 <- ggplot(
  filter(forest, outcome != "cap"),
  aes(x = estimate, y = outcome_label, colour = period)
) +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),
                 height = 0.16, position = position_dodge(width = 0.45)) +
  geom_point(position = position_dodge(width = 0.45), size = 2.4) +
  scale_x_log10() +
  labs(x = "Adjusted odds ratio per 1-SD higher DKR", y = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(figure_dir, "figure_main_forest.png"), p3,
       width = 7.0, height = 3.8, dpi = 400)

message("Analysis complete. Main complete-case counts: ",
        paste(counts$period, counts$complete_main_n, collapse = "; "))
