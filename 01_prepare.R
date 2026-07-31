suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

options(survey.lonely.psu = "adjust")

root <- normalizePath(file.path(getwd(), "work"), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(root, "data", "raw")
derived_dir <- file.path(root, "data", "derived")
results_dir <- file.path(root, "results")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

read_component <- function(period, stem) {
  read_xpt(file.path(raw_dir, period, paste0(stem, ".xpt"))) |>
    zap_labels() |>
    mutate(SEQN = as.numeric(SEQN))
}

valid_num <- function(x, max_valid = Inf) {
  x <- as.numeric(x)
  x[is.na(x) | x < 0 | x > max_valid] <- NA_real_
  x
}

yes <- function(x) ifelse(is.na(x), NA, x == 1)

any_flag <- function(...) {
  m <- cbind(...)
  out <- apply(m, 1, function(z) {
    if (any(z == 1, na.rm = TRUE)) return(1)
    if (any(!is.na(z))) return(0)
    NA_real_
  })
  as.integer(out)
}

rowmean_valid <- function(df) {
  out <- rowMeans(df, na.rm = TRUE)
  out[rowSums(!is.na(df)) == 0] <- NA_real_
  out
}

freq_days_year <- function(x) {
  map <- c(
    `0` = 0, `1` = 365, `2` = 300, `3` = 3.5 * 52,
    `4` = 2 * 52, `5` = 52, `6` = 2.5 * 12, `7` = 12,
    `8` = 9, `9` = 4.5, `10` = 1.5
  )
  as.numeric(unname(map[as.character(x)]))
}

freq_week <- function(q, unit) {
  q <- valid_num(q, 1000)
  unit_chr <- toupper(trimws(as.character(unit)))
  mult <- case_when(
    unit_chr %in% c("D", "1") ~ 7,       # per day
    unit_chr %in% c("W", "2") ~ 1,       # per week
    unit_chr %in% c("M", "3") ~ 12 / 52, # per month
    unit_chr %in% c("Y", "4") ~ 1 / 52,  # per year
    TRUE ~ NA_real_
  )
  ifelse(q == 0, 0, q * mult)
}

mean_bp <- function(dat, prefix) {
  sys <- select(dat, starts_with(paste0(prefix, "SY")))
  dia <- select(dat, starts_with(paste0(prefix, "DI")))
  list(
    sbp = rowmean_valid(as.data.frame(lapply(sys, valid_num, max_valid = 300))),
    dbp = rowmean_valid(as.data.frame(lapply(dia, valid_num, max_valid = 200)))
  )
}

build_period <- function(period) {
  is_pre <- period == "pre"
  s <- if (is_pre) {
    c(DEMO="P_DEMO",DR1="P_DR1TOT",DR2="P_DR2TOT",LUX="P_LUX",
      BMX="P_BMX",BPX="P_BPXO",BPQ="P_BPQ",GHB="P_GHB",
      BIO="P_BIOPRO",HDL="P_HDL",DIQ="P_DIQ",ALQ="P_ALQ",
      SMQ="P_SMQ",PAQ="P_PAQ",HEPBD="P_HEPBD",HEPC="P_HEPC",HEQ="P_HEQ")
  } else {
    c(DEMO="DEMO_L",DR1="DR1TOT_L",DR2="DR2TOT_L",LUX="LUX_L",
      BMX="BMX_L",BPX="BPXO_L",BPQ="BPQ_L",GHB="GHB_L",
      BIO="BIOPRO_L",HDL="HDL_L",DIQ="DIQ_L",ALQ="ALQ_L",
      SMQ="SMQ_L",PAQ="PAQ_L",HEPBD="HEPBD_L",HEPC="HEPC_L",HEQ="HEQ_L")
  }
  x <- read_component(period, s["DEMO"])
  for (key in setdiff(names(s), "DEMO")) {
    z <- read_component(period, s[key])
    duplicate_names <- intersect(setdiff(names(z), "SEQN"), names(x))
    if (length(duplicate_names)) z <- select(z, -all_of(duplicate_names))
    x <- left_join(x, z, by = "SEQN")
  }

  weight_var <- if (is_pre) "WTDR2DPP" else "WTDR2D"
  x$diet_weight <- as.numeric(x[[weight_var]])

  # Core dietary exposure: ratio of the two-day mean macronutrient intakes.
  fat <- rowMeans(cbind(valid_num(x$DR1TTFAT, 1000), valid_num(x$DR2TTFAT, 1000)), na.rm = FALSE)
  protein <- rowMeans(cbind(valid_num(x$DR1TPROT, 1000), valid_num(x$DR2TPROT, 1000)), na.rm = FALSE)
  carbohydrate <- rowMeans(cbind(valid_num(x$DR1TCARB, 2000), valid_num(x$DR2TCARB, 2000)), na.rm = FALSE)
  fiber <- rowMeans(cbind(valid_num(x$DR1TFIBE, 200), valid_num(x$DR2TFIBE, 200)), na.rm = FALSE)
  energy <- rowMeans(cbind(valid_num(x$DR1TKCAL, 10000), valid_num(x$DR2TKCAL, 10000)), na.rm = FALSE)

  x$fat_g <- fat
  x$protein_g <- protein
  x$carbohydrate_g <- carbohydrate
  x$fiber_g <- fiber
  x$energy_kcal <- energy
  x$dkr_total <- (0.9 * fat + 0.46 * protein) /
    (0.1 * fat + 0.58 * protein + carbohydrate)
  net_carb <- pmax(carbohydrate - fiber, 0)
  x$dkr_net <- (0.9 * fat + 0.46 * protein) /
    (0.1 * fat + 0.58 * protein + net_carb)

  d1 <- (0.9 * x$DR1TTFAT + 0.46 * x$DR1TPROT) /
    (0.1 * x$DR1TTFAT + 0.58 * x$DR1TPROT + x$DR1TCARB)
  d2 <- (0.9 * x$DR2TTFAT + 0.46 * x$DR2TPROT) /
    (0.1 * x$DR2TTFAT + 0.58 * x$DR2TPROT + x$DR2TCARB)
  x$dkr_day1 <- d1
  x$dkr_day2 <- d2
  x$dkr_mean_daily <- rowMeans(cbind(d1, d2), na.rm = FALSE)
  x$carb_energy_pct <- 100 * carbohydrate * 4 / energy
  # NHANES checkbox items are blank when not selected; blanks are therefore
  # classified as 0 rather than missing.
  x$self_lowcarb <- as.integer(!is.na(x$DRQSDT9) & x$DRQSDT9 == 9)
  x$self_weightloss <- as.integer(!is.na(x$DRQSDT1) & x$DRQSDT1 == 1)
  x$self_any_special <- as.integer(!is.na(x$DRQSDIET) & x$DRQSDIET == 1)

  # Liver measurements and quality.
  x$cap <- valid_num(x$LUXCAPM, 500)
  x$lsm <- valid_num(x$LUXSMED, 100)
  x$valid_vcte <- x$LUAXSTAT == 1
  x$steatosis_285 <- as.integer(x$cap >= 285)
  x$steatosis_274 <- as.integer(x$cap >= 274)
  x$steatosis_248 <- as.integer(x$cap >= 248)
  x$cacld <- as.integer(x$lsm >= 10)

  # Demographic and socioeconomic covariates.
  x$age <- valid_num(x$RIDAGEYR, 120)
  x$sex <- factor(x$RIAGENDR, levels = c(1, 2), labels = c("Men", "Women"))
  x$race <- factor(
    x$RIDRETH3,
    levels = c(1, 2, 3, 4, 6, 7),
    labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
               "Non-Hispanic Black", "Non-Hispanic Asian", "Other/multiracial")
  )
  x$education <- factor(
    case_when(
      x$DMDEDUC2 %in% c(1, 2) ~ "Less than high school",
      x$DMDEDUC2 == 3 ~ "High school/GED",
      x$DMDEDUC2 %in% c(4, 5) ~ "College or above",
      TRUE ~ NA_character_
    ),
    levels = c("Less than high school", "High school/GED", "College or above")
  )
  x$pir <- valid_num(x$INDFMPIR, 5)

  # Cardiometabolic criteria.
  x$bmi <- valid_num(x$BMXBMI, 100)
  x$waist <- valid_num(x$BMXWAIST, 250)
  bp <- mean_bp(x, "BPXO")
  x$sbp <- bp$sbp
  x$dbp <- bp$dbp
  bp_med_var <- if (is_pre) "BPQ050A" else "BPQ150"
  lipid_med_var <- if (is_pre) "BPQ100D" else "BPQ101D"
  x$bp_med <- ifelse(
    x[[bp_med_var]] == 1, 1,
    ifelse(x[[bp_med_var]] == 2 | x$BPQ020 == 2, 0, NA)
  )
  x$lipid_med <- ifelse(
    x[[lipid_med_var]] == 1, 1,
    ifelse(x[[lipid_med_var]] == 2 | x$BPQ080 == 2, 0, NA)
  )
  x$hba1c <- valid_num(x$LBXGH, 30)
  x$triglycerides <- valid_num(x$LBXSTR, 3000)
  x$hdl <- valid_num(x$LBDHDD, 200)
  x$diabetes_dx <- ifelse(x$DIQ010 %in% c(1, 3), 1,
                          ifelse(x$DIQ010 == 2, 0, NA))
  x$diabetes_med <- ifelse(
    x$DIQ050 == 1, 1,
    ifelse(x$DIQ050 == 2 | x$DIQ010 %in% c(2, 3), 0, NA)
  )

  x$cmrf_adiposity <- any_flag(
    ifelse(is.na(x$bmi), NA, x$bmi >= ifelse(x$RIDRETH3 == 6, 23, 25)),
    ifelse(is.na(x$waist), NA, x$waist >= ifelse(x$RIAGENDR == 1, 94, 80))
  )
  x$cmrf_glycemia <- any_flag(
    ifelse(is.na(x$hba1c), NA, x$hba1c >= 5.7),
    x$diabetes_dx == 1,
    x$diabetes_med == 1
  )
  x$cmrf_bp <- any_flag(
    ifelse(is.na(x$sbp), NA, x$sbp >= 130),
    ifelse(is.na(x$dbp), NA, x$dbp >= 85),
    x$bp_med == 1
  )
  x$cmrf_tg <- any_flag(
    ifelse(is.na(x$triglycerides), NA, x$triglycerides >= 150),
    x$lipid_med == 1
  )
  x$cmrf_hdl <- any_flag(
    ifelse(is.na(x$hdl), NA, x$hdl <= ifelse(x$RIAGENDR == 1, 40, 50)),
    x$lipid_med == 1
  )
  cmrf_mat <- cbind(
    x$cmrf_adiposity, x$cmrf_glycemia, x$cmrf_bp, x$cmrf_tg, x$cmrf_hdl
  )
  x$cmrf_count <- rowSums(cmrf_mat, na.rm = TRUE)
  x$cmrf_observed <- rowSums(!is.na(cmrf_mat))
  x$any_cmrf <- ifelse(x$cmrf_count >= 1, 1,
                       ifelse(x$cmrf_observed == 5, 0, NA))

  # Alcohol: midpoint frequency mapping x drinks per drinking day x 14 g.
  drink_days <- freq_days_year(x$ALQ121)
  drinks_day <- valid_num(x$ALQ130, 15)
  x$alcohol_g_day <- ifelse(
    x$ALQ121 == 0 | x$ALQ111 == 2, 0,
    drink_days * drinks_day * 14 / 365
  )
  x$metald_range <- as.integer(
    (!is.na(x$alcohol_g_day)) &
      x$alcohol_g_day > ifelse(x$RIAGENDR == 1, 30, 20) &
      x$alcohol_g_day <= ifelse(x$RIAGENDR == 1, 60, 50)
  )
  x$high_alcohol <- as.integer(
    (!is.na(x$alcohol_g_day)) &
      x$alcohol_g_day > ifelse(x$RIAGENDR == 1, 60, 50)
  )
  x$masld <- ifelse(
    x$steatosis_285 == 1 & x$any_cmrf == 1 &
      !is.na(x$alcohol_g_day) &
      x$alcohol_g_day <= ifelse(x$RIAGENDR == 1, 30, 20),
    1,
    ifelse(
      x$steatosis_285 == 0 | x$any_cmrf == 0 |
        x$metald_range == 1 | x$high_alcohol == 1,
      0, NA
    )
  )

  # Smoking.
  x$smoking <- factor(
    case_when(
      x$SMQ020 == 2 ~ "Never",
      x$SMQ020 == 1 & x$SMQ040 == 3 ~ "Former",
      x$SMQ020 == 1 & x$SMQ040 %in% c(1, 2) ~ "Current",
      TRUE ~ NA_character_
    ),
    levels = c("Never", "Former", "Current")
  )

  # Harmonized leisure-time physical activity.
  if (is_pre) {
    mod <- ifelse(x$PAQ665 == 1, valid_num(x$PAQ670, 7) * valid_num(x$PAD675, 1440), 0)
    vig <- ifelse(x$PAQ650 == 1, valid_num(x$PAQ655, 7) * valid_num(x$PAD660, 1440), 0)
  } else {
    mod_freq <- freq_week(x$PAD790Q, x$PAD790U)
    vig_freq <- freq_week(x$PAD810Q, x$PAD810U)
    mod <- ifelse(mod_freq == 0, 0, mod_freq * valid_num(x$PAD800, 1440))
    vig <- ifelse(vig_freq == 0, 0, vig_freq * valid_num(x$PAD820, 1440))
  }
  x$ltpa_moderate_equiv_min <- mod + 2 * vig
  x$meets_activity_guideline <- as.integer(x$ltpa_moderate_equiv_min >= 150)

  # Viral hepatitis sensitivity flags.
  x$hbv_positive <- ifelse(x$LBDHBG == 1, 1,
                           ifelse(x$LBDHBG == 2, 0, NA))
  # Active HCV infection is defined by a positive RNA result. LBXHCR values
  # 2 and 3 denote RNA-negative and negative screening-antibody results.
  # LBXHCG contains genotype categories and is not a binary infection marker.
  x$hcv_positive <- ifelse(
    x$LBXHCR == 1, 1,
    ifelse(x$LBXHCR %in% c(2, 3), 0, NA)
  )

  x$period <- factor(
    if (is_pre) "2017-March 2020" else "August 2021-August 2023"
  )

  # Primary eligible sample. Implausible energy exclusions follow prior NHANES diet studies.
  x$energy_plausible <- (
    (x$RIAGENDR == 1 & x$energy_kcal >= 800 & x$energy_kcal <= 4200) |
      (x$RIAGENDR == 2 & x$energy_kcal >= 600 & x$energy_kcal <= 3500)
  )
  x$eligible_base <- (
    x$age >= 20 &
      x$DR1DRSTZ == 1 & x$DR2DRSTZ == 1 &
      x$valid_vcte &
      !is.na(x$dkr_total) & is.finite(x$dkr_total) &
      !is.na(x$diet_weight) & x$diet_weight > 0 &
      x$energy_plausible &
      (is.na(x$RIDEXPRG) | x$RIDEXPRG != 1)
  )
  x
}

pre <- build_period("pre")
post <- build_period("post")

saveRDS(pre, file.path(derived_dir, "pre_full.rds"))
saveRDS(post, file.path(derived_dir, "post_full.rds"))
saveRDS(filter(pre, eligible_base), file.path(derived_dir, "pre_eligible.rds"))
saveRDS(filter(post, eligible_base), file.path(derived_dir, "post_eligible.rds"))

flow <- bind_rows(
  data.frame(
    period = "2017-March 2020",
    total = nrow(pre),
    adult = sum(pre$age >= 20, na.rm = TRUE),
    two_recall = sum(pre$age >= 20 & pre$DR1DRSTZ == 1 & pre$DR2DRSTZ == 1, na.rm = TRUE),
    valid_vcte = sum(pre$age >= 20 & pre$DR1DRSTZ == 1 & pre$DR2DRSTZ == 1 & pre$valid_vcte, na.rm = TRUE),
    eligible = sum(pre$eligible_base, na.rm = TRUE)
  ),
  data.frame(
    period = "August 2021-August 2023",
    total = nrow(post),
    adult = sum(post$age >= 20, na.rm = TRUE),
    two_recall = sum(post$age >= 20 & post$DR1DRSTZ == 1 & post$DR2DRSTZ == 1, na.rm = TRUE),
    valid_vcte = sum(post$age >= 20 & post$DR1DRSTZ == 1 & post$DR2DRSTZ == 1 & post$valid_vcte, na.rm = TRUE),
    eligible = sum(post$eligible_base, na.rm = TRUE)
  )
)
write.csv(flow, file.path(results_dir, "tables", "sample_flow.csv"), row.names = FALSE)

variable_dictionary <- data.frame(
  construct = c(
    "Two-day dietary weight", "Total-carbohydrate DKR", "Net-carbohydrate DKR",
    "Self-reported low-carbohydrate diet", "Controlled attenuation parameter",
    "Complete elastography examination", "Age", "Sex", "Race and ethnicity",
    "Education", "Income-to-poverty ratio", "Body mass index", "Waist circumference",
    "Blood pressure", "Glycated hemoglobin", "Triglycerides", "HDL cholesterol",
    "Alcohol intake", "Smoking", "Leisure-time physical activity",
    "Hepatitis B surface antigen", "Hepatitis C RNA"
  ),
  variable = c(
    "WTDR2DPP / WTDR2D",
    "DR1TTFAT, DR2TTFAT, DR1TPROT, DR2TPROT, DR1TCARB, DR2TCARB",
    "Primary DKR variables plus DR1TFIBE and DR2TFIBE", "DRQSDT9", "LUXCAPM",
    "LUAXSTAT", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR",
    "BMXBMI", "BMXWAIST", "BPXO measurements and BPQ medication variables",
    "LBXGH", "LBXSTR", "LBDHDD", "ALQ111, ALQ121, ALQ130", "SMQ020, SMQ040",
    "PAQ/PAD leisure-time activity variables", "LBDHBG", "LBXHCR"
  ),
  operational_definition = c(
    "Official two-day dietary sample weight for the corresponding period",
    "(0.9 x fat + 0.46 x protein) / (0.1 x fat + 0.58 x protein + total carbohydrate)",
    "DKR calculated after subtracting dietary fiber from total carbohydrate",
    "Checkbox selected when DRQSDT9 = 9", "Median CAP, dB/m",
    "LUAXSTAT = 1", "Age at screening, years", "Men or women",
    "Six-category RIDRETH3 classification",
    "Less than high school, high school/GED, or college or above",
    "Family income-to-poverty ratio", "kg/m2", "cm",
    "Mean available oscillometric measurements plus medication use",
    "Percent", "mg/dL", "mg/dL",
    "Past-year frequency x drinks per drinking day x 14 g / 365",
    "Never, former, or current",
    "Moderate minutes plus twice vigorous minutes per week",
    "Positive surface antigen", "Positive HCV RNA"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  variable_dictionary,
  file.path(results_dir, "tables", "variable_dictionary.csv"),
  row.names = FALSE
)
message("Prepared eligible samples: pre=", sum(pre$eligible_base, na.rm=TRUE),
        "; post=", sum(post$eligible_base, na.rm=TRUE))
