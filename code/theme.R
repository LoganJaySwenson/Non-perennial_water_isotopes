# theme for plots
plot_theme <- function(...){
  theme_bw(base_size=10)+ 
    theme(
      text=element_text(color='black'),
      axis.text=element_text(size=rel(1)),
      strip.text=element_text(size=rel(1)),
      legend.title=element_text(size=rel(1)),
      legend.text=element_text(size=rel(0.9)),
      legend.position = 'right',
      panel.grid=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin=unit(c(1,1,1,1), "mm"),
      strip.background=element_blank())
}
theme_set(plot_theme())