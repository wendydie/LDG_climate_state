# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 00_data_preparation.R
# Last updated: 2026-08-01
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(palaeoverse)
library(dplyr)
library(stringr)
source("./R/options.R")
check_pbdb <- function(dat, step = "") {
  message(
    step,
    "\n  Occurrences : ", format(n_distinct(dat$occurrence_no), big.mark = ","),
    "\n  Collections : ", format(n_distinct(dat$collection_no), big.mark = ","),
    "\n  Genera      : ", format(n_distinct(dat$genus), big.mark = ",")
  )
}

# # Data downloading from PBDB --------------------------------------------
# if (isTRUE(params$download) || !file.exists("./data/raw/pbdb_data.RDS")) {
#   library(httr); library(readr); library(dplyr)
#   root_dir <- "./data/raw/pbdb_pages"
#   page_dir <- file.path(root_dir, format(Sys.Date(), "%Y-%m"))
#   dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
#   old_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
#   if (length(old_dirs))
#     unlink(old_dirs[basename(old_dirs) != basename(page_dir)], recursive = TRUE)
#   dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)
#   read_pbdb <- function(f) {
#     hdr <- which(grepl("^occurrence_", readLines(f, n = 200)))
#     read_csv(
#       f,
#       skip = ifelse(length(hdr), hdr[1] - 1, 0),
#       show_col_types = FALSE,
#       guess_max = 100000
#     )
#   }
#   q <- params$query
#   q$limit <- q$limit %||% 50000L
#   q$offset <- 0L
#   q$rowcount <- NULL
#   q$datainfo <- NULL
#   pages <- list()
#   total <- 0L
#   repeat {
#     url <- modify_url(params$base_url, query = q)
#     resp <- RETRY(
#       "GET", url,
#       user_agent("LDG_climate_state/1.0"),
#       config(http_version = 1.1, accept_encoding = "gzip, deflate"),
#       times = 5
#     )
#     f <- file.path(page_dir, sprintf("page_%05d.csv", q$offset))
#     writeBin(content(resp, "raw"), f)
#     dat <- read_pbdb(f)
#     pages[[length(pages) + 1]] <- dat
#     total <- total + nrow(dat)
#     message(sprintf("Fetched %d rows (total %d)", nrow(dat), total))
#     if (nrow(dat) < q$limit) break
#     q$offset <- q$offset + q$limit
#   }
#   occdf <- bind_rows(pages)
#   message("Rows: ", format(nrow(occdf), big.mark = ","))
#   message("Unique occurrence_no: ", format(n_distinct(occdf$occurrence_no), big.mark = ","))
#   message("Duplicates: ", format(sum(duplicated(occdf$occurrence_no)), big.mark = ","))
#   saveRDS(occdf, "./data/raw/pbdb_data.RDS")
# } else {
#   occdf <- readRDS("./data/raw/pbdb_data.RDS") # 910734
# }
# Full PBDB download using system curl as a test---------------------------------
if (isTRUE(params$download) || !file.exists("./data/raw/pbdb_data.RDS")) {
  library(httr); library(readr); library(dplyr)
  dir.create("./data/raw", recursive = TRUE, showWarnings = FALSE)
  url <- modify_url(params$base_url, query = params$query)
  date_tag <- format(Sys.Date(), "%Y%m%d")
  csv <- file.path("./data/raw", paste0("pbdb_data_", date_tag, ".csv"))
  tmp <- paste0(csv, ".tmp")
  curl_bin <- Sys.which("curl")
  if (!nzchar(curl_bin)) stop("System curl was not found.")
  status <- system2(
    curl_bin,
    c(
      "--location",
      "--fail",
      "--show-error",
      "--http1.1",
      "--retry", "20",
      "--retry-all-errors",
      "--retry-delay", "10",
      "--connect-timeout", "60",
      "--max-time", "7200",
      "--remove-on-error",
      "--output", shQuote(tmp),
      shQuote(url)
    )
  )
  if (status != 0 || !file.exists(tmp)) {
    stop("PBDB download failed; incomplete file was not retained.")
  }
  occdf <- read_csv(tmp, show_col_types = FALSE, progress = TRUE)
  if (!"occurrence_no" %in% names(occdf)) {
    unlink(tmp)
    stop("Invalid PBDB response: occurrence_no not found.")
  }
  file.rename(tmp, csv)
  message("Rows: ", format(nrow(occdf), big.mark = ","))
  message(
    "Unique occurrence_no: ",
    format(n_distinct(occdf$occurrence_no), big.mark = ",")
  )
  message(
    "Duplicates: ",
    format(sum(duplicated(occdf$occurrence_no)), big.mark = ",")
  )
  saveRDS(occdf, "./data/raw/pbdb_data_2.RDS")
} else {
  occdf <- readRDS("./data/raw/pbdb_data_2.RDS") 
  #occ 860,604 genus 30,899 # coll 144,777
}
check_pbdb(occdf, "Raw occurrence information")
#. marine invertebrate occ 774,659, genus 27,465, coll 125,655-------------
max_bin_age <- max(bins$max_ma, na.rm = TRUE)

occdf <- occdf %>%
  filter(
    !is.na(max_ma),
    max_ma <= 485.4
  )
occdf <- occdf %>%
  filter(str_detect(research_group, "marine invertebrate")) #774659 occurrence
check_pbdb(occdf, "After marine invertebrate filter")

#. Canonical phyla occ 740,660 genus 26,666 coll 120,989
canonical_phyla <- c("Mollusca", "Arthropoda", "Brachiopoda", "Bryozoa",
  "Echinodermata", "Cnidaria", "Porifera") # 780206

occdf <- occdf %>%
  mutate(
    focal_group = case_when(
      phylum %in% canonical_phyla ~ phylum,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(focal_group))
check_pbdb(occdf, "After selecting canonical phyla")

# Set up time bins -----------------------------------------------------
bins <- time_bins(interval = "Phanerozoic",
                  rank = params$rank,
                  scale = params$GTS)
# Collapse Holocene equivalent bins
vec <- which(bins$interval_name == "Greenlandian")
bins$interval_name[vec] <- "Holocene"
bins$abbr[vec] <- "H"
# Update min_ma
bins$min_ma[vec] <- 0.0000
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Drop rows
bins <- bins[-which(bins$interval_name %in% c("Meghalayan", "Northgrippian")), ]
# Collapse Pleistocene equivalent bins
# Drop bins
pleis <- c("Late Pleistocene", "Chibanian", "Calabrian")
bins <- bins[-which(bins$interval_name %in% pleis), ]
# update Gelasian to be all of the Pleistocene
vec <- which(bins$interval_name == "Gelasian")
bins$interval_name[vec] <- "Pleistocene"
bins$abbr[vec] <- "Ple"
# Update min_ma
bins$min_ma[vec] <- bins[which(bins$interval_name == "Holocene"), "max_ma"]
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Update bin numbers
bins$bin <- 1:nrow(bins)
row.names(bins) <- 1:nrow(bins)

## GTS 2023
GTS_2023 <- read.csv('./data/GTS_2023.csv')
bins <- bins %>%
  left_join(GTS_2023, by="bin") %>%
  mutate(
    max_ma = bottom,
    mid_ma = mid,
    min_ma = top,
    duration_myr = dur
  ) %>%
  select(bin, interval_name, rank, max_ma, mid_ma, min_ma,
         duration_myr, short, font, sys, system, series,
         systemCol, seriesCol, stageCol, stageRGB)
bins <- bins %>% filter(mid_ma < 485.4)
# Save time bins
saveRDS(object = bins, file = "./data/time_bins.RDS")

# # Data cleaning and processing -----------------------------------------

# Remove suffixes from genus names
occdf$genus <- sub(" .*", "", occdf$genus)

## Remove rows where the 'genus' column contains the value "NO_GENUS_SPECIFIED" 575880
occdf <- occdf %>%
filter(genus != "NO_GENUS_SPECIFIED")
## Remove rows where any of the columns 'genus', 'lat', or 'lng' have NA values
occdf <- occdf %>%
  filter(!is.na(genus) & !is.na(lng) & !is.na(lat)) # 780190

# Round off coordinates to stack collections 
occdf$lng <- round(occdf$lng, digits = params$n_decs)
occdf$lat <- round(occdf$lat, digits = params$n_decs)

# Temporal binning -----------------------------------------------------
# Use collections to speed up binning and palaeogeographic reconstruction
colldf <- unique(occdf[, c("collection_no", "lng", "lat", "max_ma", "min_ma")])
# Use the majority method
colldf <- bin_time(occdf = colldf, bins = bins, method = params$method)
# Remove data which do not hit the majority threshold (params$threshold)
colldf <- colldf[-which(colldf$overlap_percentage < params$threshold), ] # colldf 124061

# Palaeorotate collections ---------------------------------------------
colldf <- palaeorotate(occdf = colldf,
                       lng = params$lng,
                       lat = params$lat,
                       age = params$age,
                       model = params$models,
                       method = "grid",
                       uncertainty = FALSE,
                       round = NULL)

# Exclude collections which palaeocoordinates could not be estimated for
colldf <- subset(
  colldf,
  is.finite(p_lat) &
    is.finite(p_lng) &
    p_lat >= -90 & p_lat <= 90 &
    p_lng >= -180 & p_lng <= 180
) # colldf 110684

# Join datasets --------------------------------------------------------
# Retain collections present in colldf 
occdf <- occdf[which(occdf$collection_no %in% colldf$collection_no), ]
# Join datasets
m <- match(x = occdf$collection_no, table = colldf$collection_no)
# Add data
occdf[, colnames(colldf)] <- colldf[m, colnames(colldf)] #occdf 701846
# Filter for unique occurrences from stacked collections 之前461054 现在575,419 
occdf <- distinct(occdf, lat, lng, family, genus, bin_assignment, collection_no,
                  .keep_all = TRUE)
check_pbdb(occdf, "After removing duplicates")
# number of occ 575,419, genus 23,211, coll 110,684
# Save processed data
saveRDS(object = occdf, file = "./data/processed/pbdb_data.RDS")
# Notify
if (params$notify) {
  beepr::beep(4)
}
