library(devtools)
devtools::install_github('hyeonlee-stats/psacomp')
library(psacomp)

library(plyr)
library(dplyr)
library(ggplot2)
library(ggtern)
library(GGally)
library(viridis)
library(scales)
library(ggnewscale)
library(cowplot) # plot_grid
library(reshape) # parallel coordinate plot
library(latex2exp)
library(compositions)

invisible(lapply(list.files('utils', pattern = '.R', full.names = T), source))
