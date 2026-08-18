# Function: incfreq
# Description:
#   This function calculates the frequency of occurrences for each unique genus
#   in the input dataset and returns a list containing the total number of rows 
#   in the dataset and the frequency of each genus in descending order.
# 
# Arguments:
#   dat - A dataframe that must include a column named 'genus', representing 
#         biological genus names or similar categories.
# 
# Returns:
#   A list containing:
#     1. The total number of rows in the input dataset (total occurrences).
#     2. A vector of frequencies for each genus, sorted in descending order.
# 
# Example:
#   Input:
#     dat <- data.frame(genus = c("Panthera", "Canis", "Panthera", 
#                                 "Felis", "Panthera", "Felis"))
#   Output:
#     list(c(6, 3, 2, 1)) 
#     # 6 = total rows, 3 = Panthera, 2 = Felis, 1 = Canis

incfreq <- function(dat) {
  incidence_dat <- dat %>%
    distinct(collection_no, genus)
  
  freq <- incidence_dat %>%
    count(genus, name = "frequency") %>%
    arrange(desc(frequency))
  
  nT <- n_distinct(dat$collection_no)
  
  as.numeric(c(nT, freq$frequency))
}


Fun.ince <- function(x) {
  nT <- x[1]  # The total sample size, which is the first element in the vector x
  x <- x[-1]  # The remaining elements in x are the counts of species or occurrences
  U <- sum(x)  # Total number of occurrences (sum of all remaining elements)
  # Calculate Qk: a vector where Qk[k] represents the number of species occurring exactly k times
  Qk <- sapply(1:10, function(k) sum(x == k))
  Q1 <- Qk[1]  # Number of species that occur exactly once
  Q2 <- Qk[2]  # Number of species that occur exactly twice
  Sobs <- sum(x > 0)  # Number of observed species (species with at least one occurrence)
  # Estimate the number of undetected species (Q0.hat)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1)/nT * Q1 * (Q1 - 1)/2, 
                   (nT - 1)/nT * Q1^2/2/Q2)
  # Calculate A, which is used to adjust the estimate of sample coverage
  A <- ifelse(Q1 > 0, nT * Q0.hat/(nT * Q0.hat + Q1), 1)
  # Calculate the estimated sample coverage (Chat), rounded to four decimal places
  Chat <- round(1 - Q1/U * A, 4)
  Good.u <- 1 - Q1/U
  multition.rate <- (Sobs - Q1) / Sobs # cited from John Alroy on 21 March 2017
  # Return a vector of results including nT, U, Sobs, Chat, and Qk (Q1 to Q10)
  out <- c(nT, U, Sobs, Chat, Good.u, multition.rate, Qk)
}


process_data <- function(dat,
                         q = c(0, 1, 2),
                         datatype = "incidence_freq",
                         base = "coverage",
                         level = NULL,
                         nboot = 50,
                         req_cols=c("occurrence_no", "collection_no", "reference_no")) {
  # Step 1: Calculate frequency per tsbin
  inc_dat <- lapply(split(dat, dat$tsbin), incfreq)
  
  # Step 2: Diversity estimation using iNEXT
  rich_est <- iNEXT::estimateD(x = inc_dat, q = q, datatype = datatype, 
                               base = base, level = level, nboot = nboot)
  colnames(rich_est)[colnames(rich_est) == "Assemblage"] <- "tsbin"
  # Step 3: Extract diversity info
  rich_info <- do.call(rbind, lapply(inc_dat, Fun.ince))
  rich_info <- as.data.frame(rich_info)
  colnames(rich_info) <- c("nT", "nOcc", "Sobs", "Chat", "GoodU", "Mrate", paste0("Q", 1:10))
  rich_info$tsbin <- rownames(rich_info)
  num_cols <- c("nT", "nOcc", "Sobs", "Chat", "GoodU", "Mrate", paste0("Q", 1:10))
  rich_info[num_cols] <- lapply(rich_info[num_cols], as.numeric)
  
  # Step 4: Unique counts for required columns
  ts_stats <- do.call(rbind, lapply(split(dat, dat$tsbin), function(g) {
    stats <- sapply(req_cols, function(col) length(unique(g[[col]])))
    c(tsbin = unique(g$tsbin), stats)
  }))
  ts_stats <- as.data.frame(ts_stats)
  ts_stats[req_cols] <- lapply(ts_stats[req_cols], as.numeric)
  
  # Step 5: Merge diversity info and unique stats
  # Step 5: Merge all results
  rich_est$tsbin <- as.character(rich_est$tsbin)
  rich_info$tsbin <- as.character(rich_info$tsbin)
  ts_stats$tsbin <- as.character(ts_stats$tsbin)
  rich_finalInfo <- Reduce(function(x, y) merge(x, y, by = "tsbin", all = TRUE), 
                           list(rich_est, rich_info, ts_stats))
  return(rich_finalInfo)
}
