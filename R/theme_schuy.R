theme_schuy <- function(
    graph_type     = "all",
    angle          = 0,
    base_size      = 10,
    base_family    = "",
    base_line_size = base_size / 22,
    base_rect_size = base_size / 22
) {
  schuy_theme <- theme_gray(
    base_size      = base_size,
    base_family    = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size) %+replace%
    theme(
      panel.background   = element_rect(fill = "white", colour = NA),
      panel.border       = element_blank(),
      panel.grid         = element_line(colour = "#F1F0F2"),
      panel.grid.minor.x = element_blank(),
      legend.title       = element_text(size = base_size * 1.2, face = "bold",
                            hjust = 0),
      legend.text        = element_text(size = base_size),
      legend.spacing.x   = unit(0.005, "npc"),
      legend.key         = element_rect(fill = NA, colour = NA),
      legend.key.size    = unit(4, "mm"),
      axis.ticks.y = element_blank(),
      axis.text.x        = element_text(size = base_size, angle = 0,
                            vjust = 0, hjust = 0.5),
      axis.text.y        = element_text(size = base_size),
      axis.title.x       = element_text(size = base_size * 1.2, face = "bold",
                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y       = element_text(size = base_size * 1.2, face = "bold",
                            margin = margin(t = 0, r = 10, b = 0, l = 0),
                            angle = 90),
      panel.spacing      = unit(.002, 'npc'),
      strip.text.x       = element_text(size = base_size * 0.9, face = "bold"),
      strip.background   = element_rect(colour = "black", fill = "#DCDDDF",
                            linewidth = 1.4),
      complete           = TRUE
    )
  if (angle > 0 & angle < 90) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        axis.text.x = element_text(size = base_size,
        angle = angle, vjust = 1, hjust = 1))
    }
  if (angle < 0 & angle > -90) {
    schuy_theme <- schuy_theme %+replace% 
      theme(
        axis.text.x = element_text(size = base_size,
        angle = angle, vjust = 1, hjust = 0.05))
  }
  if (abs(angle) == 90) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        axis.text.x = element_text(size = base_size,
        angle = angle, vjust = 1, hjust = 1))
  }
  if (grepl("bar", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
      )
  }
  if (grepl("box", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
      )
  }
  if (grepl("line", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
      )
  }
  if (grepl("heat", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(size = base_size, hjust = 1)
      )
  }
  if (grepl("ord|pc|rda|ca", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
      )
  }
  if (grepl("nmds|tsne", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
      )
  }
  if (grepl("graph|net", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
      )
  }
  if (grepl("dend", graph_type)) {
    schuy_theme <- schuy_theme %+replace%
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "left"
      )
  }
  return(schuy_theme)
}
