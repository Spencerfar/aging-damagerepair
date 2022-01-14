library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(tidybayes)
library(latex2exp)
library(tidyr)
library(cowplot)


x <- seq(0, 6, by = 6./1000)




img <- ggplot() +
    geom_line(aes(x=x, y=exp(x)-exp(1)+1)) +
    geom_line(aes(x=x, y=2.5*x-1.5), linetype='dotdash') +
    #geom_line(aes(x=x, y=1+log(x)*1.4), linetype='longdash') +
    #geom_line(aes(x=x, y=2-(2-x)**2), linetype='longdash') +
    geom_line(aes(x=x, y=2-(2-x)**2), linetype='longdash') +
    #geom_line(aes(x=x, y=1+((x-1)**0.25)*1.4), linetype='longdash') +
        theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title.y=element_text(size=7),
        axis.title.x=element_text(size=7,vjust=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0.02, 0.02, 0.02, 0.02), "cm")) + labs(x='Age', y = 'Frailty Index') + ggtitle('a) ') + coord_cartesian(ylim=c(0.95, 5),xlim=c(1,2)) + scale_x_continuous(breaks=c(1.0, 1.25, 1.5, 1.75, 2.0), labels=c("","","","",""))+
    scale_y_continuous(breaks=c(1, 2, 3, 4, 5), labels=c("","","","",""))+
    geom_text(aes(x=1.55, y=4.35, label='Positive curvature'),size=2) +
    geom_text(aes(x=1.87, y=2.5, label='Zero curvature'),size=2) +
    geom_text(aes(x=1.75, y=1.5, label='Negative curvature'),size=2)
    #scale_y_continuous(breaks=c(0,1), label=c("0", "1")) +
    #scale_x_continuous(breaks=c(17,19, 21), label=c("17", '19', '21'), expand = c(0, 0)) + coord_cartesian(xlim=c(16.5, 22.3)) +
    

ggsave("curvature_schematic.pdf", img, width=2, height=1.5)
