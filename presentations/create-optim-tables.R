#!/usr/bin/Rscript
script_options <- docopt::docopt(
    "Usage: create-optim-tables.R [options]

    Options:
    --optim-input-path=<path>   Path to optim data dir [default: optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots]
    --table-output-path=<path>  Path to write .tex files [default: presentations/tables/fit105]
    --model=<model>             Model name [default: STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP]
    ",
    args = if (interactive()) "" else commandArgs(trailingOnly = TRUE)
)

library(tidyverse)
library(knitr)
library(kableExtra)

optim_input_path  <- script_options$optim_input_path
table_output_path <- script_options$table_output_path
models_we_want    <- script_options$model

dir.create(table_output_path, showWarnings = FALSE, recursive = TRUE)

custom_save_latex_table <- function(table, table_name) {
    table_conn <- file(file.path(table_output_path, paste0(table_name, ".tex")))
    attr(table, "kable_meta")$contents <-
        str_replace_all(attr(table, "kable_meta")$contents, "removeme12345", " ")
    table[1] <- str_replace_all(table[1], "removeme12345", " ")
    clean_table <- table %>%
        str_remove(fixed("\\begin{table}")) %>%
        str_remove("\\\\caption\\{.*\\}") %>%
        str_remove("\\\\end\\{table\\}") %>%
        str_replace_all("\\\\vphantom\\{[^}]+\\}\\s*", "")
    writeLines(clean_table, table_conn)
    close(table_conn)
    invisible(table)
}

# ── Load data ─────────────────────────────────────────────────────────────────
all_optim_df <- read_csv(
    file.path(optim_input_path, "posterior-clean-summ-optim.csv")
) %>%
    mutate(fix_type = replace_na(fix_type, "no-fix"))

optim_group_vars <- c(
    "private_benefit_z", "visibility_z", "static_vstar",
    "allocation_type", "distance_constraint", "fix_type",
    "fix_distance", "rep_type", "model", "cutoff_type"
)

summ_optim_df <- all_optim_df %>%
    group_by(across(all_of(optim_group_vars))) %>%
    mutate(mean_dist = mean_dist / 1000) %>%
    filter(model %in% c(models_we_want, "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP")) %>%
    filter(cutoff_type == "cutoff") %>%
    summarise(
        across(
            .cols = c(mean_dist, mean_demand),
            list(
                estimate = ~ mean(.x, na.rm = TRUE) %>% round(3),
                CI = ~ paste0("(", round(quantile(.x, 0.025, na.rm = TRUE), 3),
                              ", ", round(quantile(.x, 0.975, na.rm = TRUE), 3), ")")
            )
        ),
        across(
            .cols = c(n_pot),
            list(
                estimate = ~ mean(.x, na.rm = TRUE) %>% round(0),
                CI = ~ paste0("(", round(quantile(.x, 0.025, na.rm = TRUE), 0),
                              ", ", round(quantile(.x, 0.975, na.rm = TRUE), 0), ")")
            )
        )
    ) %>%
    ungroup() %>%
    select(-model, -cutoff_type) %>%
    arrange(private_benefit_z, visibility_z, rep_type) %>%
    rename(B_z = private_benefit_z, mu_z = visibility_z)

wide_summ_df_rep <- summ_optim_df %>%
    filter(rep_type == "rep") %>%
    gather(variable, value, -any_of(optim_group_vars), -B_z, -mu_z) %>%
    mutate(CI = if_else(str_detect(variable, "CI"), "ci", "estim")) %>%
    mutate(variable = str_remove(variable, "_CI")) %>%
    mutate(variable = str_remove(variable, "_estimate")) %>%
    pivot_wider(
        id_cols = c(B_z, mu_z, variable, any_of(optim_group_vars)),
        names_from = CI, values_from = value
    ) %>%
    mutate(estim = if_else(
        variable %in% c("n_pot"), as.numeric(estim),
        as.numeric(estim) %>% round(2)
    )) %>%
    mutate(estim_value = linebreak(paste0(estim, "\n", ci), align = "c")) %>%
    select(-estim, -ci) %>%
    spread(variable, estim_value)

wide_summ_df_sup_rep <- summ_optim_df %>%
    filter(rep_type == "suppress_rep") %>%
    gather(variable, value, -any_of(optim_group_vars), -B_z, -mu_z) %>%
    mutate(CI = if_else(str_detect(variable, "CI"), "ci", "estim")) %>%
    mutate(variable = str_remove(variable, "_CI")) %>%
    mutate(variable = str_remove(variable, "_estimate")) %>%
    pivot_wider(
        id_cols = c(B_z, mu_z, variable, rep_type, static_vstar,
                    allocation_type, distance_constraint),
        names_from = CI, values_from = value
    ) %>%
    mutate(estim_value = linebreak(paste0(estim, "\n", ci), align = "c")) %>%
    select(-estim, -ci) %>%
    spread(variable, estim_value)

optim_summ_input_df <- bind_rows(
    wide_summ_df_rep %>% mutate(rep_type = "rep") %>% filter(B_z == "control"),
    wide_summ_df_sup_rep %>% mutate(rep_type = "suppress_rep") %>% filter(B_z == mu_z)
) %>%
    filter(fix_type == "no-fix" | is.na(fix_type)) %>%
    mutate(
        mu_z = if_else(rep_type == "suppress_rep", "No visibility", mu_z),
        mu_z = case_when(
            static_vstar == FALSE ~ str_to_title(mu_z),
            static_vstar == TRUE & mu_z == "control" ~ "Signal value fixed at control 0.5km",
            static_vstar == TRUE & mu_z == "bracelet" ~ "Signal value fixed at bracelet 0.5km",
            TRUE ~ mu_z
        )
    ) %>%
    mutate(mu_z = factor(mu_z, c(
        "Bracelet", "Calendar", "Ink", "Control", "No Visibility",
        "Signal value fixed at control 0.5km",
        "Signal value fixed at bracelet 0.5km"
    ))) %>%
    mutate(B_z = factor(B_z, c("bracelet", "calendar", "ink", "control"))) %>%
    mutate(B_z = fct_relabel(B_z, str_to_title)) %>%
    mutate(mu_z = fct_relevel(mu_z, c(
        "Control", "Ink", "Calendar",
        "Signal value fixed at control 0.5km",
        "Signal value fixed at bracelet 0.5km",
        "Bracelet", "No Visibility"
    ))) %>%
    mutate(mu_z = fct_recode(mu_z,
        "Signalling value at control 0.5km"  = "Signal value fixed at control 0.5km",
        "Signalling value at bracelet 0.5km" = "Signal value fixed at bracelet 0.5km"
    )) %>%
    mutate(
        mean_dist = if_else(
            allocation_type == "experimental",
            str_remove(mean_dist, "\\\\\\(([^)]+)\\)") %>% str_replace(., "\\\\\\}", "}"),
            mean_dist
        ),
        n_pot = if_else(
            allocation_type == "experimental",
            str_remove(n_pot, "\\\\\\(([^)]+)\\)") %>% str_replace(., "\\\\\\}", "}"),
            n_pot
        )
    )

old_layout_order <- c(
    "Control",
    "Bracelet",
    "Control social image returns at 0.5km",
    "Bracelet social image returns at 0.5km",
    "No social image returns"
)

old_layout_robust_order <- c(
    "Control",
    "Bracelet",
    "Bracelet social image returns at 0.5km",
    "No social image returns"
)

format_old_layout <- function(data) {
    data %>%
        mutate(
            Observability = fct_recode(
                mu_z,
                "Control social image returns at 0.5km" =
                    "Signalling value at control 0.5km",
                "Bracelet social image returns at 0.5km" =
                    "Signalling value at bracelet 0.5km",
                "No social image returns" = "No Visibility"
            )
        )
}

# ── Main table (optim-table-main) ─────────────────────────────────────────────
optim_summ_tbl <- optim_summ_input_df %>%
    filter(distance_constraint == 3500 | allocation_type == "experimental") %>%
    format_old_layout() %>%
    mutate(Observability = factor(Observability, old_layout_order)) %>%
    arrange(rep_type, Observability) %>%
    select(Observability, n_pot, mean_demand, mean_dist) %>%
    knitr::kable(
        col.names = c("Observability", "Assigned PoTs", "Mean take-up",
                      "Mean distance (km)"),
        align = "lccc", booktabs = TRUE, format = "latex",
        linesep = "", escape = FALSE
    ) %>%
    kableExtra::kable_styling(latex_options = c("scale_down")) %>%
    pack_rows(
        index = c("Panel A: Experimental allocation" = 1,
                  "Panel B: Policymaker allocation" = 5),
        italic = TRUE, escape = FALSE, hline_after = TRUE,
        hline_before = TRUE, bold = TRUE, indent = FALSE
    )

custom_save_latex_table(optim_summ_tbl, "optim-summ-table")
message("Written: ", file.path(table_output_path, "optim-summ-table.tex"))

# ── Robust table (optim-summ-table-robust) ────────────────────────────────────
optim_summ_robust_tbl <- optim_summ_input_df %>%
    filter(allocation_type == "experimental" | distance_constraint %in% c(3500, 4500, 5500, 10000)) %>%
    format_old_layout() %>%
    filter(allocation_type == "experimental" | Observability %in% old_layout_robust_order) %>%
    mutate(Observability = factor(Observability, old_layout_robust_order)) %>%
    arrange(allocation_type, distance_constraint, Observability) %>%
    select(Observability, n_pot, mean_demand, mean_dist) %>%
    knitr::kable(
        col.names = c("Observability", "\\makecell{Number of PoTs}",
                      "Mean take-up", "\\makecell{Mean distance (km)}"),
        align = "lccc", booktabs = TRUE, escape = FALSE, format = "latex"
    ) %>%
    pack_rows(
        index = c(
            "Panel A: Experimental allocation" = 1,
            "Panel B: Policymaker allocation, 3.5km" = 4,
            "Panel C: Policymaker allocation, 4.5km" = 4,
            "Panel D: Policymaker allocation, 5.5km" = 4,
            "Panel E: Policymaker allocation, 10km"  = 4
        ),
        italic = TRUE, escape = FALSE, hline_after = TRUE,
        hline_before = TRUE, bold = TRUE, indent = FALSE
    )

custom_save_latex_table(optim_summ_robust_tbl, "optim-summ-robust-table")
message("Written: ", file.path(table_output_path, "optim-summ-robust-table.tex"))
