#install.packages("ggplot2")
#install.packages("grid")
library(ggplot2)
library(grid)

## Theme
my.theme <- theme(
    plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
    axis.text = element_text(colour="black", size=12),
    axis.title.x = element_text(colour="black", size=12),
    axis.title.y = element_text(colour="black", size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
)

## Comparison of the R^2,(1-mer + 1-shape) model v.s. (1-mer) model
model_1mer <- c(0.717, 0.702, 0.726)
model_1mer_1shape <- c(0.734, 0.720, 0.730)

## Ploting
bp <- ggplot() +
    geom_point(aes(x = model_1mer, y = model_1mer_1shape), color = "red", size=1) +
    geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    my.theme  
bp
bp + ggtitle("1-mer+1-shape v.sv 1-mer")
