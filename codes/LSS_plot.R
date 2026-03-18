library(this.path)
setwd(this.dir())

library(showtext)
library(ggplot2)

font_add_google("Josefin Sans", "josefin")
showtext_auto()

my_ggplot_theme <- theme_bw() +
  theme(
    plot.title   = element_text(size = 20, family = "josefin"),
    axis.title.x = element_text(size = 20, family = "josefin"),
    axis.title.y = element_text(size = 20, family = "josefin"),
    axis.text.x  = element_text(size = 15, family = "josefin"),
    axis.text.y  = element_text(size = 15, family = "josefin"),
    legend.text  = element_text(size = 15, family = "josefin"),
    legend.title = element_text(size = 20, family = "josefin"),
    legend.box.background = element_rect()
  )

draw_plot <- function(file, col, mu, var, out) {
  x <- read.csv(file)[[col]]
  
  p <- ggplot(data.frame(x = (x - mu) / sqrt(var)), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)),
                   color = "blue", fill = "cyan") +
    geom_function(fun = dnorm, color = "red", linewidth = 1) +
    labs(x = "Eigenvalue", y = "Density") +
    my_ggplot_theme
  
  ggsave(out, plot = p, width = 6, height = 4)
}

draw_plot("../data/LSS/LSS-Chisq-2.csv", 2, -6, 24,       "../figure/LSS/LSS-Chisq-x1.pdf")
draw_plot("../data/LSS/LSS-Chisq-2.csv", 3, -24, 1600,    "../figure/LSS/LSS-Chisq-x2.pdf")
draw_plot("../data/LSS/LSS-Chisq-2.csv", 4, -80, 96000,   "../figure/LSS/LSS-Chisq-x3.pdf")

draw_plot("../data/LSS/LSS-Exp-2.csv", 2, -2, 4,          "../figure/LSS/LSS-Exp-x1.pdf")
draw_plot("../data/LSS/LSS-Exp-2.csv", 3, -4, 68,         "../figure/LSS/LSS-Exp-x2.pdf")
draw_plot("../data/LSS/LSS-Exp-2.csv", 4, -7, 1050,       "../figure/LSS/LSS-Exp-x3.pdf")