#!/usr/bin/Rscript

script_options <- docopt::docopt(
  "Usage:
  create-gpt-reason-plot.R [options]

Options:
  --input=<path>            GPT-classified reason data [default: temp-data/second-order-know-notknow-gpt-reason-full-data.csv]
  --output-dir=<path>       Directory where figure files should be saved [default: figures]
  --output-basename=<name>  Output filename without extension [default: gpt-reason-plot]
  --formats=<formats>       Comma-separated formats for ggsave [default: pdf,png]
  --width=<width>           Figure width in inches [default: 8]
  --height=<height>         Figure height in inches [default: 6]
  --dpi=<dpi>               DPI for raster outputs [default: 300]
",
  args = if (interactive()) "
    --input=temp-data/second-order-know-notknow-gpt-reason-full-data.csv
    --output-dir=figures
    --output-basename=gpt-reason-plot
    --formats=pdf,png
  " else commandArgs(trailingOnly = TRUE)
)

library(tidyverse)
library(here)
library(RColorBrewer)

options(dplyr.show_progress = FALSE, digits = 4)

theme_set(
  theme_minimal() +
    theme(legend.position = "bottom")
)

path_from_root <- function(path) {
  if (grepl("^/|^[A-Za-z]:[/\\\\]", path)) path else here(path)
}

input_path <- path_from_root(script_options$input)
output_dir <- path_from_root(script_options$output_dir)
output_basename <- script_options$output_basename
output_formats <- str_split(script_options$formats, ",", simplify = TRUE) %>%
  as.character() %>%
  str_trim() %>%
  discard(~ .x == "")
fig_width <- as.numeric(script_options$width)
fig_height <- as.numeric(script_options$height)
fig_dpi <- as.numeric(script_options$dpi)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gpt_endline_know_notknow_table_data <- read_csv(
  input_path,
  show_col_types = FALSE
)

gpt_endline_know_notknow_table_data <- gpt_endline_know_notknow_table_data %>%
  mutate(
    category_sob_reason_short = case_when(
      str_detect(second.order.reason, "see") ~ "campaign",
      assigned.treatment == "control" & category_sob_reason_short == "signal" ~ "other",
      TRUE ~ category_sob_reason_short
    )
  )

gpt_translation_df <- tribble(
  ~gpt_name, ~clean_name,
  "campaign", "Direct Observation",
  "communication", "Communication",
  "relationship", "Social Proximity",
  "signal", "Item/signal",
  "type", "Personal Characteristics",
  "circumstances", "Personal Characteristics",
  "other", "Other",
  "Don't Know", "Don't Know"
)

reason_levels <- c(
  "Direct Observation",
  "Communication",
  "Social Proximity",
  "Item/signal",
  "Personal Characteristics",
  "Other",
  "Don't Know"
)

gpt_df <- gpt_endline_know_notknow_table_data %>%
  mutate(
    gpt_name = if_else(second.order == "don't know", "Don't Know", category_sob_reason_short)
  ) %>%
  left_join(
    gpt_translation_df,
    by = c("gpt_name" = "gpt_name")
  ) %>%
  filter(!is.na(clean_name)) %>%
  mutate(
    clean_name = factor(clean_name, levels = reason_levels)
  )

col_gpt_df <- gpt_translation_df %>%
  distinct(clean_name) %>%
  mutate(
    colors = RColorBrewer::brewer.pal(n = 8, name = "Accent")[seq_len(n())],
    colors = if_else(clean_name == "Don't Know", "#AEAEAE", colors)
  )

gpt_reason_plot <- gpt_df %>%
  select(
    knowledge,
    dist.pot.group,
    assigned.treatment,
    category_sob_reason_short,
    clean_name
  ) %>%
  mutate(
    assigned.treatment = str_to_title(assigned.treatment),
    dist.pot.group = str_to_title(dist.pot.group),
    assigned.treatment = factor(
      assigned.treatment,
      levels = c("Control", "Ink", "Calendar", "Bracelet")
    ),
    dist.pot.group = factor(dist.pot.group, levels = c("Close", "Far"))
  ) %>%
  ggplot(aes(
    x = assigned.treatment,
    fill = clean_name
  )) +
  geom_bar(
    position = "fill",
    colour = "black"
  ) +
  facet_grid(~dist.pot.group) +
  scale_fill_manual(
    "",
    values = col_gpt_df$colors,
    labels = col_gpt_df$clean_name
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Treatment",
    y = "Proportion"
  )

walk(
  output_formats,
  ~ggsave(
    filename = file.path(output_dir, str_glue("{output_basename}.{.x}")),
    plot = gpt_reason_plot,
    width = fig_width,
    height = fig_height,
    dpi = fig_dpi
  )
)

message(
  str_glue(
    "Saved {length(output_formats)} file(s) to {output_dir}: ",
    "{str_c(str_glue('{output_basename}.{output_formats}'), collapse = ', ')}"
  )
)
