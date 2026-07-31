# LSM >= 8 kPa sensitivity analysis (add-on to 02_analyze.R)
# Mirrors the design, standardization, covariates, and inference of the main analysis.
suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
})

options(survey.lonely.psu = "adjust")

root <- normalizePath(file.path(getwd(), "work"), winslash = "/", mustWork = TRUE)
derived_dir <- file.path(root, "data", "derived")
out_dir <- file.path(root, "results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

weighted_sd <- function(design, variable) {
  sqrt(as.numeric(svyvar(as.formula(paste0("~", variable)), design, na.rm = TRUE)))
}

tidy_svy <- function(model, term, exponentiate = FALSE) {
  beta <- coef(model)
  se <- sqrt(diag(vcov(model)))[term]
  df <- degf(model$survey.design)
  crit <- qt(0.975, df = df)
  row <- data.frame(
    term = term,
    estimate = unname(beta[term]),
    std.error = unname(se),
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

model_formula <- function(outcome, exposure) {
  as.formula(paste(
    outcome, "~", exposure,
    "+ age + sex + race + education + pir + smoking + meets_activity_guideline",
    "+ alcohol_g_day + energy_kcal"
  ))
}

rows <- list()
for (p in names(periods)) {
  dat <- periods[[p]]
  dat$lsm_ge8 <- as.integer(dat$lsm >= 8)
  design <- make_design(dat)
  dkr_mean <- as.numeric(coef(svymean(~dkr_total, design, na.rm = TRUE)))
  dkr_sd <- weighted_sd(design, "dkr_total")
  dat$dkr_z <- (dat$dkr_total - dkr_mean) / dkr_sd
  design <- make_design(dat)

  prev <- svymean(~lsm_ge8, design, na.rm = TRUE)
  mod <- svyglm(model_formula("lsm_ge8", "dkr_z"), design = design,
                family = quasibinomial(), na.action = na.omit)
  tr <- tidy_svy(mod, "dkr_z", exponentiate = TRUE)
  mf <- model.frame(mod)
  rows[[p]] <- data.frame(
    period = period_labels[[p]],
    n = nobs(mod),
    cases = sum(dat$lsm_ge8[complete.cases(
      dat[, c("lsm_ge8", "dkr_total", "age", "sex", "race", "education", "pir",
              "smoking", "meets_activity_guideline", "alcohol_g_day",
              "energy_kcal")])] == 1, na.rm = TRUE),
    weighted_prevalence = as.numeric(coef(prev)[1]),
    or = tr$estimate, conf_low = tr$conf.low, conf_high = tr$conf.high,
    log_estimate = as.numeric(coef(mod)["dkr_z"]),
    log_se = sqrt(as.numeric(vcov(mod)["dkr_z", "dkr_z"])),
    p_value = tr$p.value
  )
}

res <- bind_rows(rows)
# Between-period heterogeneity (inverse-variance difference of log ORs).
z <- (res$log_estimate[2] - res$log_estimate[1]) /
  sqrt(res$log_se[2]^2 + res$log_se[1]^2)
res$period_heterogeneity_p <- NA_real_
res$period_heterogeneity_p[1] <- 2 * pnorm(abs(z), lower.tail = FALSE)

write.csv(res, file.path(out_dir, "lsm_ge8_sensitivity.csv"), row.names = FALSE)
print(res, row.names = FALSE)
