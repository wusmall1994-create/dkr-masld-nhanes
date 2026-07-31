options(timeout = 300)

root <- normalizePath(file.path(getwd(), "work"), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(root, "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

files <- list(
  pre = c(
    DEMO = "P_DEMO", DR1 = "P_DR1TOT", DR2 = "P_DR2TOT", LUX = "P_LUX",
    BMX = "P_BMX", BPX = "P_BPXO", BPQ = "P_BPQ", GHB = "P_GHB",
    BIO = "P_BIOPRO", HDL = "P_HDL", DIQ = "P_DIQ", ALQ = "P_ALQ",
    SMQ = "P_SMQ", PAQ = "P_PAQ", HEPBD = "P_HEPBD",
    HEPC = "P_HEPC", HEQ = "P_HEQ"
  ),
  post = c(
    DEMO = "DEMO_L", DR1 = "DR1TOT_L", DR2 = "DR2TOT_L", LUX = "LUX_L",
    BMX = "BMX_L", BPX = "BPXO_L", BPQ = "BPQ_L", GHB = "GHB_L",
    BIO = "BIOPRO_L", HDL = "HDL_L", DIQ = "DIQ_L", ALQ = "ALQ_L",
    SMQ = "SMQ_L", PAQ = "PAQ_L", HEPBD = "HEPBD_L",
    HEPC = "HEPC_L", HEQ = "HEQ_L"
  )
)

bases <- c(
  pre = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles",
  post = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles"
)

manifest <- data.frame()
for (period in names(files)) {
  period_dir <- file.path(raw_dir, period)
  dir.create(period_dir, recursive = TRUE, showWarnings = FALSE)
  for (key in names(files[[period]])) {
    stem <- unname(files[[period]][key])
    dest <- file.path(period_dir, paste0(stem, ".xpt"))
    url <- paste0(bases[[period]], "/", stem, ".xpt")
    if (!file.exists(dest) || file.info(dest)$size < 1000) {
      message("Downloading ", url)
      download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
    }
    manifest <- rbind(
      manifest,
      data.frame(
        period = period, key = key, stem = stem, url = url,
        bytes = file.info(dest)$size,
        md5 = unname(tools::md5sum(dest)),
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(manifest, file.path(raw_dir, "download_manifest.csv"), row.names = FALSE)
message("Downloaded/verified ", nrow(manifest), " files.")
