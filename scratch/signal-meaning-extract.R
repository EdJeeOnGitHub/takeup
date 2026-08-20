#!/usr/bin/env Rscript

# Extract the endline item-meaning free-text responses into (i) a blind,
# deduplicated coding sheet and (ii) a respondent-level analysis file with the
# original randomized distance assignment.  The coding sheet deliberately
# excludes every treatment and distance field so that coders cannot condition
# on them.  Base R only, so the script runs without the full renv library.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

out_dir <- arg_value("--output-dir", "ref-reports/signal-meaning-coding")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
clean_id <- function(x) as.integer(trimws(as.character(x)))

message("Loading endline data")
endline <- readRDS("data/clean-data/clean-endline-data.rds")

items <- data.frame(
  item = c("bracelet", "ink", "calendar"),
  meaning_var = c("bracelet_meaning", "ink_meaning", "cal_meaning"),
  seen_var = c("seen_bracelet", "seen_ink", "seen_cal"),
  expected_n = c(950L, 620L, 669L),
  stringsAsFactors = FALSE
)

# One row per respondent in each item's own arm.  The meaning question was
# only administered to respondents who reported having seen the item, so
# non-seers carry an empty text and support the unconditional sensitivity.
respondent <- do.call(rbind, lapply(seq_len(nrow(items)), function(i) {
  meaning <- as.character(endline[[items$meaning_var[i]]])
  has_response <- !is.na(meaning) & nchar(trimws(meaning)) > 0
  keep <- has_response | as.character(endline$assigned.treatment) == items$item[i]
  data.frame(
    person_id = endline$KEY.individ[keep],
    cluster_id = clean_id(endline$cluster.id[keep]),
    county = as.character(endline$county[keep]),
    item = items$item[i],
    arm = as.character(endline$assigned.treatment[keep]),
    final_far = as.integer(endline$dist.pot.group[keep] == "far"),
    sms_treatment = as.character(endline$sms.treatment[keep]),
    seen_item = as.integer(endline[[items$seen_var[i]]][keep]),
    has_response = as.integer(has_response[keep]),
    age = endline$age[keep],
    text = ifelse(has_response[keep], trimws(meaning[keep]), NA_character_),
    stringsAsFactors = FALSE
  )
}))

for (i in seq_len(nrow(items))) {
  n_item <- sum(respondent$item == items$item[i] & respondent$has_response == 1L)
  assert(n_item == items$expected_n[i],
         sprintf("%s: expected %d responses, got %d",
                 items$item[i], items$expected_n[i], n_item))
  # Each meaning question was only administered in its own arm.
  assert(all(respondent$arm[respondent$item == items$item[i] &
                              respondent$has_response == 1L] == items$item[i]),
         sprintf("%s meaning responses found outside the %s arm",
                 items$item[i], items$item[i]))
}

message("Merging original distance assignment and covariates")
schools <- readRDS("data/rct_targetable_schools_2.0.rds")
original <- unique(data.frame(
  cluster_id = clean_id(schools$pot.cluster.id),
  original_far = as.integer(schools$assigned.dist.cat == "far")
))
respondent <- merge(respondent, original, by = "cluster_id", all.x = TRUE, sort = FALSE)
assert(!anyNA(respondent$original_far),
       "Some meaning respondents have no original distance assignment")

full <- readRDS("data/clean-data/full-takeup-data.rds")
person_cov <- unique(data.frame(
  person_id = full$KEY.individ,
  female = as.numeric(full$gender == "female"),
  stringsAsFactors = FALSE
))
respondent <- merge(respondent, person_cov, by = "person_id", all.x = TRUE, sort = FALSE)

takeup <- read.csv("temp-data/analysis-cluster-recentered-covariate-data.csv",
                   check.names = FALSE)
cluster_cov <- unique(data.frame(
  cluster_id = clean_id(takeup$cluster.id.x),
  mu_d = takeup$mu_d,
  distance_km = takeup$cluster.dist.to.pot / 1000
))
assert(!any(duplicated(cluster_cov$cluster_id)),
       "Cluster covariates are not constant within community")
respondent <- merge(respondent, cluster_cov, by = "cluster_id", all.x = TRUE, sort = FALSE)

# Blind, deduplicated coding sheet: one row per distinct (item, normalized text).
responder <- respondent$has_response == 1L
text_norm <- tolower(gsub("\\s+", " ", respondent$text))
dedup_key <- ifelse(responder, paste(respondent$item, text_norm, sep = "\r"), NA)
first_idx <- responder & !duplicated(dedup_key)
dedup <- respondent[first_idx, c("item", "text")]
dedup$n_respondents <- as.integer(table(dedup_key[responder])[dedup_key[first_idx]])
dedup <- dedup[order(dedup$item, tolower(dedup$text)), ]
dedup$response_id <- sprintf("%s-%03d", substr(dedup$item, 1, 3),
                             ave(seq_len(nrow(dedup)), dedup$item, FUN = seq_along))
dedup <- dedup[c("response_id", "item", "text", "n_respondents")]

dedup_lookup_key <- paste(dedup$item, tolower(gsub("\\s+", " ", dedup$text)), sep = "\r")
respondent$response_id <- dedup$response_id[match(dedup_key, dedup_lookup_key)]
assert(!anyNA(respondent$response_id[responder]),
       "Some responses failed to match a dedup row")
assert(sum(dedup$n_respondents) == sum(responder),
       "Dedup counts do not add up to the responder total")

write.csv(dedup, file.path(out_dir, "responses-dedup.csv"), row.names = FALSE)
saveRDS(respondent, file.path(out_dir, "respondent-level.rds"))
write.csv(respondent, file.path(out_dir, "respondent-level.csv"), row.names = FALSE)

message(sprintf("Coding sheet: %d distinct responses (%s)",
                nrow(dedup),
                paste(sprintf("%s %d", items$item,
                              tabulate(factor(dedup$item, items$item))), collapse = ", ")))
message(sprintf("Respondent file: %d rows, %d communities",
                nrow(respondent), length(unique(respondent$cluster_id))))
for (itm in items$item) {
  sub <- respondent[respondent$item == itm, ]
  message(sprintf("  %s: %d in arm (%d responded, %d seen==0, %d seen NA), %d communities, original Far share %.2f",
                  itm, nrow(sub), sum(sub$has_response),
                  sum(sub$seen_item == 0, na.rm = TRUE), sum(is.na(sub$seen_item)),
                  length(unique(sub$cluster_id)), mean(sub$original_far)))
}
