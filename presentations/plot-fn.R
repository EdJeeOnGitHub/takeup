
plot_ord_beliefs_est = function(beliefs_results, 
                             top_title = NULL, 
                             width = 0.3, 
                             crossbar_width = 0.2, 
                             order = 1,
                             caption = default_plot_caption
                             ) {
  pos_dodge <- position_dodge(width = width)

  first_plot = beliefs_results %>%
    filter(param == "belief_probs") %>%
    select(tidy_draws) %>%
    unnest(tidy_draws) %>%
    filter(str_detect(variable, as.character(order))) %>%
    filter(dist_group %in% c("Close", "Far")) %>%
    ggplot(aes(
      y = treatment,
      group = dist_group,
      xmin = conf.low,
      xmax = conf.high,
      x = value
    )) +
      geom_linerange(aes(
        color = dist_group), 
        position = pos_dodge, 
        size = 0.4,
        data = . %>% filter(.width == 0.95)
        ) +
      geom_crossbar(aes(
        color = dist_group), 
        position = pos_dodge, 
        fatten = 2, 
        size = 0.4, 
        width = crossbar_width,
        data = . %>% filter(.width == 0.8)
        ) +
      geom_linerange(aes(
        color = dist_group), 
        position = pos_dodge, 
        alpha = 0.4, 
        size = 2.5,
        data = . %>% filter(.width == 0.5)
        ) +
      geom_point(aes(color = dist_group), position = pos_dodge, size = 1.8) +
      geom_point(position = pos_dodge, color = "white", size = 0.6) +
      scale_color_canva("", labels = str_to_title, palette = canva_palette_vibrant) + 
      scale_y_discrete("", labels = str_to_title) +
      labs(
        title = "",
        subtitle = "Panel A: Levels",
        x = "") +
      labs(x = "Proportion (%)") +
      scale_x_continuous(
        labels = scales::label_percent(suffix = "")
      )  +
      theme_minimal() +
      theme(legend.position = "top") +
      panel_border(remove = TRUE) +
      NULL


  if (!is.null(top_title)) {
    rel_heights = c(0.1, 1, 0.08)
  } else {
    rel_heights = c(0, 1, 0.08)
  }



  cowplot::plot_grid(
    if (!is_null(top_title)) { 
      cowplot::ggdraw() +
        cowplot::draw_label(top_title, size = 20, fontface = "italic")
    },
    
    cowplot::plot_grid(
      first_plot +
        theme(
          legend.position = "none"
        ) +
        NULL,



      second_plot = beliefs_results %>%
        filter(param == "belief_ates") %>%
        select(tidy_draws) %>%
        unnest(tidy_draws) %>%
        filter(str_detect(variable, as.character(order))) %>%
        filter(dist_group %in% c("Close", "Far")) %>%
        ggplot(aes(
          y = treatment,
          group = dist_group,
          xmin = conf.low,
          xmax = conf.high,
          x = value
          )) +
        geom_vline(xintercept = 0, linetype = "dotted") +
        geom_linerange(aes(
          color = dist_group), 
          position = pos_dodge, 
          size = 0.3,
          data = . %>% filter(.width == 0.95)
          ) +
        geom_crossbar(aes(
          color = dist_group), 
          position = pos_dodge, 
          fatten = 2, 
          size = 0.4, 
          width = crossbar_width,
          data = . %>% filter(.width == 0.8)
          ) +
        geom_linerange(aes(
          color = dist_group), 
          position = pos_dodge, 
          alpha = 0.4, 
          size = 2.25,
          data = . %>% filter(.width == 0.5)
          ) +
        geom_point(aes(
          color = dist_group), position = pos_dodge, size = 1.8) +
        geom_point(position = pos_dodge, color = "white", size = 0.6) +
        scale_y_discrete(drop = FALSE) +
        scale_color_canva("", labels = str_to_title, palette = canva_palette_vibrant) + 
        labs(
          title = "",
          subtitle = "Panel B: Treatment Effects",
          x = "", y = "") +
    labs(x = "Percentage Points") +
    caption +
    scale_x_continuous(
      labels = scales::label_percent(suffix = "")
    )  +
    theme_minimal() +
    panel_border(remove = TRUE) +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none"
        ) +
        NULL,
        ncol = 2, axis = "b", align = "h" 
    ),
    
    cowplot::get_legend(first_plot),
    ncol = 1, 
    rel_heights = rel_heights
  )
}
