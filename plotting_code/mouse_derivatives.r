library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(bayesplot)
library(bayestestR)
library(tidybayes)
library(latex2exp)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(cowplot)
source("../utils/functions.r")
source("../utils/palettes.r")

##### mouse dataset 1
enalapril <- read.csv('../datasets/enalapril_data.csv', header = TRUE, sep = ",")

enalapril$age <- enalapril$time + enalapril$baseline.age
enalapril$prepair <- (enalapril$repair/enalapril$n)/enalapril$delta.t
enalapril$pdamage <- (enalapril$damage/(enalapril$N - enalapril$n))/enalapril$delta.t
enalapril$treatment <- as.factor(enalapril$treatment)
levels(enalapril$treatment) <- c('Control', 'Enalapril')
enalapril.test <- na.omit(enalapril)

# FI standard deviation
f.sd = sd(enalapril.test$f)

bins <- seq(14.5, 29, by=1)
binned <- cut(enalapril$age, bins, include.lowest = TRUE)
binned <- midpoints(binned)
enalapril$time.bin <- binned

# get stan fit
fit <- readRDS('../fits/mouse_1.rds')

bins <- seq(0, 20, by=1)
binned.test <- cut(enalapril.test$time, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)


full.repair.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_r_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(f = enalapril.test[n,'f']) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>% # need to unstandardize
    mutate(deriv_rr = deriv_r + dotf*deriv_r_f)  %>% # total derivative
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(lambda_r = mean(deriv_rr)) %>%
    ungroup() %>%
    group_by(sex, treatment, age) %>%
    median_hdci(lambda_r, .width=0.95)

full.damage.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_d[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(f = enalapril.test[n,'f']) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>% # need to unstandardize
    mutate(deriv_dd = deriv_d + dotf*deriv_d_f) %>% # total derivative
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(lambda_d = mean(deriv_dd)) %>%
    ungroup() %>%
    group_by(sex, treatment, age) %>%
    median_hdci(lambda_d, .width=0.95)


enalapril.repair <- ggplot() +
    geom_lineribbon(data=full.repair.plot.stats, mapping=aes(x=age, y=lambda_r, fill=treatment,
                                                             color=NA, ymin=.lower, ymax=.upper), alpha=0.5) +
    geom_line(data=full.repair.plot.stats, mapping=aes(x=age, y=lambda_r,
          fill=treatment, color=treatment,group=treatment), alpha=1, size=1.75, color='white') +
    geom_line(data=full.repair.plot.stats, mapping=aes(x=age, y=lambda_r,
         fill=treatment, color=treatment, group=treatment), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='Repair rate time slope', x='Time since intervention (months)',
         color='', fill='') +
    scale_fill_d3() +
    scale_color_d3() + ggtitle('c) Mouse dataset 1 (Keller et al. 2019)')



enalapril.damage <- ggplot() +
    geom_lineribbon(data=full.damage.plot.stats, mapping=aes(x=age, y=lambda_d, fill=treatment,
                                                             color=NA, ymin=.lower, ymax=.upper), alpha=0.5) +
    geom_line(data=full.damage.plot.stats, mapping=aes(x=age, y=lambda_d,
          fill=treatment, color=treatment,group=treatment), alpha=1, size=1.75, color='white') +
    geom_line(data=full.damage.plot.stats, mapping=aes(x=age, y=lambda_d,
         fill=treatment, color=treatment, group=treatment), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6), 
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='Damage rate time slope', x='Time since intervention (months)',
         color='', fill='') +
    scale_fill_d3() +
    scale_color_d3() + xlim(0, 10.25)



###### mouse dataset 2
exercise <- read.csv('../datasets/exercise_data.csv', header = TRUE, sep = ",")

exercise$sex <- as.factor(exercise$sex)
exercise$exercise <- as.factor(exercise$exercise)
exercise$age <- exercise$time + exercise$baseline.age
exercise$prepair <- (exercise$repair/exercise$n)/exercise$delta.t
exercise$pdamage <- (exercise$damage/(exercise$N - exercise$n))/exercise$delta.t
levels(exercise$exercise) <- c('Control', 'Exercise')

exercise.test <- na.omit(exercise)

# frailty index standard deviation
f.sd = sd(exercise.test$f)



bins <- seq(20.5, 27, by=0.5)
binned <- cut(exercise$age, bins, include.lowest = TRUE)
binned <- midpoints(binned)
exercise$time.bin <- binned

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))

# get stan fit
fit <- readRDS('../fits/mouse_2.rds')

bins <- seq(0, 15, by=0.5)
binned.test <- cut(exercise.test$time, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)


full.repair.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_r_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(f = exercise.test[n,'f']) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>% # unstandardize
    mutate(deriv_rr = deriv_r + dotf*deriv_r_f)  %>% # total derivative
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_r = mean(deriv_rr)) %>%
    ungroup() %>%
    group_by(sex, exercise, age) %>%
    median_hdci(lambda_r, .width=0.95)


full.damage.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_d[n],deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(f = exercise.test[n,'f']) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>% # unstandardize
    mutate(deriv_dd = deriv_d + dotf*deriv_d_f) %>% # total derivative
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_d = mean(deriv_dd)) %>%
    ungroup() %>%
    group_by(sex, exercise, age) %>%
    median_hdci(lambda_d, .width=0.95)





exercise.repair <- ggplot() +
    geom_lineribbon(data=full.repair.plot.stats, mapping=aes(x=age, y=lambda_r, fill=exercise,
                                                             color=NA, ymin=.lower, ymax=.upper), alpha=0.5) +
    geom_line(data=full.repair.plot.stats, mapping=aes(x=age, y=lambda_r,
          fill=exercise, color=exercise,group=exercise), alpha=1, size=1.75, color='white') +
    geom_line(data=full.repair.plot.stats, mapping=aes(x=age, y=lambda_r,
         fill=exercise, color=exercise, group=exercise), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='', x='Time since intervention (months)',
         color='', fill='') + 
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette)  +
    scale_x_continuous(breaks=c(0, 1, 2, 3),label = c("0", "1", "2", "3"))+ xlim(0, 3.25) + ggtitle('d) Mouse dataset 2 (Bisset et al. 2021)')

exercise.damage <- ggplot() +
    geom_lineribbon(data=full.damage.plot.stats, mapping=aes(x=age, y=lambda_d, fill=exercise,
                                                             color=NA, ymin=.lower, ymax=.upper), alpha=0.5) +
    geom_line(data=full.damage.plot.stats, mapping=aes(x=age, y=lambda_d,
          fill=exercise, color=exercise,group=exercise), alpha=1, size=1.75, color='white') +
    geom_line(data=full.damage.plot.stats, mapping=aes(x=age, y=lambda_d,
         fill=exercise, color=exercise, group=exercise), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='', x='Time since intervention (months)',
         color='', fill='') + 
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette) + scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3"))+ xlim(0, 3.25)


img.enalapril <- plot_grid(enalapril.repair, enalapril.damage, ncol=1,nrow=2, rel_heights = c(1.6/3, 1.4/3)) 
img.exercise <- plot_grid(exercise.repair, exercise.damage, ncol=1,nrow=2, rel_heights = c(1.6/3, 1.4/3)) 

img <- grid.arrange(img.enalapril, img.exercise, ncol=2,nrow=1)

ggsave("../plots/mice_rates_derivatives.pdf", img, width=8, height=3)
