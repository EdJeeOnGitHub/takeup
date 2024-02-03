library(tidyverse)
library(posterior)
library(kableExtra)




rfnr = "REDUCED_FORM_NO_RESTRICT"
rfnrngp = "REDUCED_FORM_NO_RESTRICT_NO_GP"
rfnrdb = "REDUCED_FORM_NO_RESTRICT_DIFFUSE_BETA"
struct_dbdc = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_DIFFUSE_BETA_DIFFUSE_CLUSTER"
struct_hier_s = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_SOB"
struct_hier_f = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FOB"
struct_hier_f_fixed = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FIXED_FOB"

models_we_want = tribble(
    ~model, ~fit_version,
    rfnr, 97,
    rfnr, 98,
    rfnrngp, 98,
    rfnrdb, 99,
    rfnr, 99,
    struct_dbdc, 98,
    struct_hier_s, 98,
    struct_hier_f, 95,
    struct_hier_f, 100,
    struct_hier_f, 101,
    struct_hier_f_fixed, 101
)

ate_rvar_df = models_we_want %>%
    mutate(
        fp = str_glue(
                "rvar_processed_dist_fit{fit_version}_ates_{model}_\\d-\\d.rds$"
            ),
        files = map(fp, ~fs::dir_ls("temp-data", regexp = .x))
        ) %>%
    unnest(files) %>%
    mutate(
        last_digit = str_extract(files, "\\d+(?=\\.rds$)") %>% as.numeric,
    ) %>%
    group_by(model, fit_version) %>%
    filter(last_digit == max(last_digit)) %>%
    select(-last_digit) %>%
    ungroup() %>%
    mutate(data = map(files, read_rds))  




ate_rvar_df = ate_rvar_df %>%
    group_by(fit_version, model) %>%
    group_split() %>%
    map(
        ~ .x %>%
            select(-model, -fit_version) %>%
            unnest(data) %>%
        mutate(
            model_type = if_else(str_detect(model, "STRUCT"), "structural", "reduced form") %>% factor(),
            fit_type = factor(fit_type),
            model = factor(model)
        )
    )




# ate_rvar_df$value = c(ate_rvar_rf$value)

create_cis = function(.data, .width = 0.95, latex = TRUE ) {
  med_fun = function(x) {
    mean_x = mean(x) %>% round(3)
    conf.low = quantile(x, (1 - .width)/2) %>% round(3)
    conf.high = quantile(x, 1 - (1 - .width)/2) %>% round(3)
    if (latex) {
        return_val = linebreak(
        paste0(
            mean_x, "\n", str_glue(
            "({conf.low}, {conf.high})"
            )
        ), align = "c"
        )

        } else {
        return_val = paste0(
            mean_x, "\n", str_glue(
            "({conf.low}, {conf.high})"
            )
        )
    }
    return(return_val)
  }
  .data %>%
    mutate(across(where(is_rvar), med_fun))
}


create_ate_table = function(.data, .estimand, group_var = treatment) {
  .data %>%
    filter(estimand == .estimand) %>%
    select(
      model,
      {{ group_var }},
      dist_group,
      value
    ) %>%
    pivot_wider(
      names_from = dist_group,
      values_from = value
    ) %>%
    select(model, {{ group_var }}, any_of(c("combined", "close", "far"))) %>%
    arrange({{ group_var }}) %>%
    bind_rows(
      # bracelet minus calendar row
      .data %>%
        filter({{ group_var }} %in% c("Bracelet", "Calendar")) %>%
        filter(estimand == .estimand) %>%
        pivot_wider(names_from = {{ group_var }}, values_from = value) %>%
        mutate(
          bracelet_minus_calendar = Bracelet - Calendar
        ) %>%
        select(model, dist_group, bracelet_minus_calendar) %>%
        pivot_wider(
          names_from = dist_group,
          values_from = bracelet_minus_calendar
        ) %>%
        mutate("{{ group_var }}" := "bracelet_minus_calendar")
    )
}

ate_rvar_df =  ate_rvar_df %>%
    map(
        ~ group_by(.x, model, fit_version) %>%
            group_nest() %>%
            mutate(
                data = map2(data, model, ~mutate(.x, model = .y)),
                ate_table = map(data, ~create_ate_table(.x, .estimand = "overall", group_var = treatment))
            )
    )


ate_tbl_levels = c(
    "Bracelet",
    "Calendar",
    "Ink",
    "bracelet_minus_calendar",
    "Control"
)

recode_control_mean = function(data) {
  data %>%
    mutate(
      treatment = fct_recode(treatment, "Control mean" = "Control", "Bracelet - Calendar" = "bracelet_minus_calendar")
    )
}

nice_kbl_table = function(tbl, cap, stat = "ci") {

  linesep_str = if_else(stat == "ci", "\\addlinespace", "")

  nice_kbl = tbl %>%
  kbl(
    col.names = c(
      # "Estimand", 
      # "Treatment", 
      "Dependent variable: Take-up",
      paste0("(", 1:4, ")")
    ), 
    # format = "latex", 
    linesep = linesep_str, 
    booktabs = TRUE, 
    escape = FALSE, 
    align = "lcccc", 
    caption = cap
  )  %>%
  kable_styling(
    latex_options = c("scale_down")
  ) %>%
  add_header_above(
    c(" ", 
      "Combined", 
      "Close", 
      "Far", 
      "Far - Close"
      ), 
    line = FALSE
  ) %>%
  add_header_above(
    c(
      " " = 1,
      "Reduced Form" = 4
      )
  ) %>%
  row_spec(c(3), hline_after = TRUE) 

}

spread_rf = function(data) {
  data %>%
    # mutate(model_type = "rf") %>%
    select(-model) %>%
    pivot_wider(
      names_from = model_type,
      id_cols = treatment,
      names_glue = "{model_type}_{.value}",
      values_from = any_of(c("combined", "close", "far", "far_minus_close"))
    ) %>%
    select(treatment, starts_with("rf"), starts_with("struct"))
}

incentive_ate_df = ate_rvar_df %>%
    map(
        ~
    mutate(.x,
        tbl = map(
            ate_table,
            ~mutate(
                .x,
                model_type = if_else(str_detect(model, "STRUCT"), "struct", "rf") %>% factor(),
                fit_type = "fit",
                model = factor(model)
            ) %>%
            mutate(treatment = factor(treatment, levels = 
                ate_tbl_levels
            )) %>%
            mutate(far_minus_close = far - close) %>%
            arrange(treatment)   %>%
            create_cis(.width = 0.95, latex = FALSE) %>%
            recode_control_mean()  %>%
            spread_rf() 
        )
    )
    )

incentive_ate_df = incentive_ate_df %>%
    map(
        ~mutate(.x,
        cap_string = str_glue(
            "Stan: {str_replace_all(model, '_', ' ')}, Version: {fit_version}"
        ),
        kbl_table = map2(
            tbl, cap_string,
            nice_kbl_table
        )
    )

    )

all_incentive_ate_df = incentive_ate_df %>%
    bind_rows()


all_incentive_ate_df %>%
    filter(str_detect(model, struct_hier_f)) %>%
    pull(kbl_table)



