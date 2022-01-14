library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(tidybayes)
library(latex2exp)
library(tidyr)
library(cowplot)

enalapril <- read.csv('../../cleaned_data/Cleaned_mice_data_step_May23.csv', header = TRUE, sep = ",")

enalapril$age.bin <- as.factor(enalapril$age.bin)
enalapril$real.age <- enalapril$time + enalapril$baseline.age

enalapril$prepair <- (enalapril$repair/enalapril$n)/enalapril$delta.t
enalapril$pdamage <- (enalapril$damage/(enalapril$N - enalapril$n))/enalapril$delta.t

enalapril$age <- enalapril$real.age

img <- ggplot() +
    geom_step(data=enalapril[enalapril$mouse == 'E11.4',], mapping=aes(x=age, y = d0)) +
        theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0.02, 0.02, 0.02, 0.02), "cm")) + labs(x='Age (months)', y = 'Deficit value') +
    scale_y_continuous(breaks=c(0,1), label=c("0", "1")) +
    scale_x_continuous(breaks=c(17,19, 21), label=c("17", '19', '21'), expand = c(0, 0)) + coord_cartesian(xlim=c(16.5, 22.3)) +
    ggtitle('a) ')

ggsave("transitions_schematic.pdf", img, width=4, height=1.25)
