## Knowledge Table Attrition Investigation
## Tests the three mechanisms in docs/notes/knowledge-table-attrition.md
## Run after sourcing reduced-form-setup.R (as in balance.R)

library(tidyverse)
library(fixest)
library(magrittr)

# reduced-form-setup.R needs params even when run non-interactively
params <- list(
  table_output_path = "presentations/rf-tables/main-specs",
  show_probs        = FALSE,
  width             = 0.95,
  cache             = FALSE,
  fit               = FALSE,
  fit_version       = "",
  output_path       = "temp-data",
  stat              = "std.error"
)

source(file.path("scratch", "reduced-form-setup.R"))

# endline_data already has in_know_table from clean-endline-data.R
# survey.type = "individual" / "paired" (from endline_and_reconsent.csv pre-load)
# endline.type = "endline" / "endline.backup" / "reconsent" / "reconsent.backup"
# knowledge_id1a ... knowledge_id20a are the pre-loaded IDs

analysis_df <- endline_data %>%
  mutate(
    not_in_know_table = !in_know_table,
    is_backup         = endline.type %in% c("endline.backup", "reconsent", "reconsent.backup"),
    # Has at least one pre-loaded knowledge ID (proxy for being in the knowledge list)
    has_know_ids      = !is.na(knowledge_id1a),
    assigned.treatment = factor(assigned.treatment, levels = c("control", "bracelet", "ink", "calendar"))
  )


# ── 1. Are missing respondents disproportionately backups? ───────────────────

cat("\n=== 1. Backup status by knowledge table presence ===\n")

analysis_df %>%
  count(is_backup, not_in_know_table) %>%
  group_by(is_backup) %>%
  mutate(pct = n / sum(n)) %>%
  print()

# Is backup status itself treatment-correlated?
cat("\n=== 1b. Backup rate by treatment arm ===\n")

analysis_df %>%
  group_by(assigned.treatment, dist.pot.group) %>%
  summarise(
    n          = n(),
    n_backup   = sum(is_backup),
    pct_backup = mean(is_backup),
    .groups = "drop"
  ) %>%
  arrange(assigned.treatment, dist.pot.group) %>%
  print()

# Regression: does backup status predict attrition from know table?
backup_predicts_attrition <- feols(
  not_in_know_table ~ is_backup | county,
  data    = analysis_df,
  cluster = ~cluster.id
)

cat("\n=== 1c. Does backup status predict attrition? ===\n")
summary(backup_predicts_attrition)


# ── 2. Does conditioning on backup absorb the treatment effect? ──────────────

# Baseline (same as balance.R attrition_treat_dist)
treat_only <- feols(
  not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") | county,
  data    = analysis_df,
  cluster = ~cluster.id
)

# Add backup indicator
treat_plus_backup <- feols(
  not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") + is_backup | county,
  data    = analysis_df,
  cluster = ~cluster.id
)

# Add has_know_ids indicator (did the form have pre-loaded knowledge IDs?)
treat_plus_knowids <- feols(
  not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") + has_know_ids | county,
  data    = analysis_df,
  cluster = ~cluster.id
)

# Both together
treat_plus_both <- feols(
  not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") + is_backup + has_know_ids | county,
  data    = analysis_df,
  cluster = ~cluster.id
)
names(coef(treat_only))
cat("\n=== 2. Does conditioning on backup / know IDs absorb treatment effect? ===\n")
etable(treat_only, treat_plus_backup, treat_plus_knowids, treat_plus_both,
       headers = c("Treat only", "+ Backup", "+ Know IDs", "+ Both"),
      #  keep    = c("bracelet", "is_backup", "has_know_ids"),
       fitstat = c("n", "r2", "my"),
       dict = c(
        "assigned.treatmentbracelet"                         = "Bracelet",
        "assigned.treatmentink"                              = "Ink",
        "assigned.treatmentcalendar"                         = "Calendar",
        "dist.pot.groupfar"                                  = "Far distance",
        "assigned.treatment::bracelet:dist.pot.group::close" = "Bracelet $\\times$ Close",
        "assigned.treatment::ink:dist.pot.group::close"      = "Ink $\\times$ Close",
        "assigned.treatment::calendar:dist.pot.group::close" = "Calendar $\\times$ Close",
        "is_backupTRUE"                                      = "Backup respondent",
        "dewormed_selfTRUE"                                  = "Dewormed"
       )
       )


# ── 3. Is attrition predicted by own deworming within treatment arm? ─────────

# dewormed.reported from the endline (1 = yes, or similar)
# Check what values it takes
cat("\n=== 3. Deworming variable distribution ===\n")
analysis_df %>% count(dewormed.reported) %>% print()

# Recode to logical if needed
analysis_df <- analysis_df %>%
  mutate(dewormed_self = dewormed.reported == 1)

# Within each treatment arm: does own deworming predict not being in know table?
cat("\n=== 3. Within-arm: does own deworming predict attrition? ===\n")
analysis_df %>%
  group_by(assigned.treatment) %>%
  summarise(
    pct_missing_if_dewormed     = mean(not_in_know_table[dewormed_self],  na.rm = TRUE),
    pct_missing_if_not_dewormed = mean(not_in_know_table[!dewormed_self], na.rm = TRUE),
    n_dewormed                  = sum(dewormed_self, na.rm = TRUE),
    n_not_dewormed              = sum(!dewormed_self, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(diff = pct_missing_if_dewormed - pct_missing_if_not_dewormed) %>%
  print()

# Regression: own deworming (controlling for treatment) predicting attrition
dworm_attrition <- feols(
  not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") + dewormed_self | county,
  data    = analysis_df,
  cluster = ~cluster.id
)

# Additive: backup + own deworming
dworm_plus_backup_attrition <- feols(
  not_in_know_table ~ i(assigned.treatment, dist.pot.group, ref = "control") + is_backup + dewormed_self | county,
  data    = analysis_df,
  cluster = ~cluster.id
)

cat("\n=== 3b. Own deworming + treatment predicting attrition ===\n")
summary(dworm_attrition)


# ── 4. Summary: how much of the bracelet gap is each mechanism? ─────────────

cat("\n=== 4. Summary of bracelet close coefficient across specs ===\n")
etable(
  treat_only, treat_plus_backup, treat_plus_knowids, treat_plus_both,
  headers = c("Baseline", "+ Backup", "+ Know IDs", "+ Both"),
  keep    = c("%bracelet", "%is_backup", "%has_know"),
  fitstat = c("n", "r2")
)

# Also print backup rates by treatment arm explicitly
cat("\n=== 1b full. Backup rate by treatment arm ===\n")
analysis_df %>%
  group_by(assigned.treatment, dist.pot.group) %>%
  summarise(n = n(), n_backup = sum(is_backup), pct_backup = mean(is_backup), .groups = "drop") %>%
  print(n = Inf)


# ── 5. Save .tex table ───────────────────────────────────────────────────────
names(coef(treat_only))
etable(
  # treat_only,
  treat_plus_backup,
  dworm_plus_backup_attrition,
  headers      = c("Backup", "+ Own deworming"),
  keep         = c("%bracelet", "%is_backup", "%dewormed_self"),
  fitstat = c("n", "r2", "my"),
  dict = c(
    "assigned.treatment::bracelet:dist.pot.group::close" = "Bracelet $\\times$ Close",
    "assigned.treatment::bracelet:dist.pot.group::far"   = "Bracelet $\\times$ Far",
    "assigned.treatment::ink:dist.pot.group::close"      = "Ink $\\times$ Close",
    "assigned.treatment::ink:dist.pot.group::far"        = "Ink $\\times$ Far",
    "assigned.treatment::calendar:dist.pot.group::close" = "Calendar $\\times$ Close",
    "assigned.treatment::calendar:dist.pot.group::far"   = "Calendar $\\times$ Far",
    "is_backupTRUE"                                      = "Backup respondent",
    "dewormed_selfTRUE"                                  = "Dewormed"
  ),
  depvar       = FALSE,
  notes        = "",
  tex          = TRUE,
  file         = file.path("presentations", "tables", "know-attrition-backup-dworm-table.tex"),
  replace      = TRUE,
  digits = 3
)
