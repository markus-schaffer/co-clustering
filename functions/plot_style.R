theme_nice <- function() {
  theme_bw() +
    theme(
      text = element_text(size = 8, colour = "black"),
      axis.text = element_text(size = 6, colour = "black"),
      line = element_line(colour = "black", linewidth = 0.1),
      rect = element_rect(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(colour = "black", linewidth = 0.1),
      panel.border = element_rect(colour = "black", linewidth = 0.1),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white", linewidth = 0.1)
    )
}
