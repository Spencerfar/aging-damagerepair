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

enalapril.test <- na.omit(enalapril)


f.sd = sd(enalapril.test$f)


bins <- seq(14.5, 29, by=1)
binned <- cut(enalapril$age, bins, include.lowest = TRUE)
binned <- midpoints(binned)
enalapril$time.bin <- binned



# get stan fit
fit <- readRDS('../fits/mouse_1.rds')

bins <- seq(14.5, 29, by=1)
binned.test.age <- cut(enalapril.test$age, bins, include.lowest = TRUE)
binned.test.age <- midpoints(binned.test.age)


fdot.plot <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(f = enalapril.test[n, 'f']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    mutate(age = binned.test.age[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r =  - f*deriv_r - dotf*lambda_r) %>%
    mutate(term.d =  (1-f)*deriv_d - dotf*lambda_d) %>%
    mutate(ratio = term.d - term.r) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(ratio = mean(ratio)) %>%
    ungroup() %>%
    group_by(sex, treatment, age) %>%
    median_hdci(ratio, .width = 0.95)


fdot.plot$sex <- as.factor(fdot.plot$sex)
fdot.plot$treatment <- as.factor(fdot.plot$treatment)
levels(fdot.plot$treatment) <- c('Control', 'Enalapril')

levels(fdot.plot$sex) <- c('Female', 'Male')

fdot.plot <- fdot.plot[fdot.plot$treatment=='Control',]

terms.enalapril <- ggplot() +
    geom_errorbar(data=fdot.plot, mapping=aes(x=age, y=ratio, ymin=.lower, ymax=.upper, fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=fdot.plot, mapping=aes(x=age, y=ratio), alpha=1, size=1.75, color='white') +
    geom_point(data=fdot.plot, mapping=aes(x=age, y=ratio, color=sex), alpha=1, size=1.25) +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex) +
    labs(y='FI curvature terms difference', x='Age (months)',
         color='', fill='') +
    scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
                       label = c("16", "18", "20", "22", "24")) +
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + ggtitle('a) Mouse dataset 1 (Keller et al. 2019)')


##### mouse dataset 2
exercise <- read.csv('../datasets/exercise_data.csv', header = TRUE, sep = ",")

exercise$sex <- as.factor(exercise$sex)
exercise$exercise <- as.factor(exercise$exercise)
exercise$age <- exercise$time + exercise$baseline.age

exercise$prepair <- (exercise$repair/exercise$n)/exercise$delta.t
exercise$pdamage <- (exercise$damage/(exercise$N - exercise$n))/exercise$delta.t

exercise.test <- na.omit(exercise)


f.sd = sd(exercise.test$f)


# get stan fit
fit <- readRDS('../fits/mouse_2.rds')

bins <- seq(20.5, 27, by=0.5)
binned.test <- cut(exercise.test$age, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

fdot.plot <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(f = exercise.test[n, 'f']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r =  - f*deriv_r - dotf*lambda_r) %>%
    mutate(term.d =  (1-f)*deriv_d - dotf*lambda_d) %>%
    mutate(ratio = term.d - term.r) %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(ratio = mean(ratio)) %>%
    ungroup() %>%
    group_by(sex, exercise, age) %>%
    median_hdci(ratio, .width = 0.95)

fdot.plot$exercise <- as.factor(fdot.plot$exercise)
fdot.plot$sex <- as.factor(fdot.plot$sex)
levels(fdot.plot$exercise) <- c('Control', 'Exercise')
levels(fdot.plot$sex) <- c('Female', 'Male')

fdot.plot <- fdot.plot[fdot.plot$exercise=='Control',]


terms.exercise <- ggplot() +
    geom_errorbar(data=fdot.plot, mapping=aes(x=age, y=ratio, ymin=.lower, ymax=.upper, fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=fdot.plot, mapping=aes(x=age, y=ratio), alpha=1, size=1.75, color='white') +
    geom_point(data=fdot.plot, mapping=aes(x=age, y=ratio, color=sex), alpha=1, size=1.25) +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex) +
    labs(y='', x='Age (months)',
         color='', fill='') +
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)+
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + ggtitle('b) Mouse dataset 2 (Bisset et al. 2021)')


##### mouse dataset 3
schultz <- read.csv('../datasets/schultz_data.csv', header = TRUE, sep = ",")

schultz$age <- schultz$time + schultz$baseline.age
schultz$prepair <- (schultz$repair/schultz$n)/schultz$delta.t
schultz$pdamage <- (schultz$damage/(schultz$N - schultz$n))/schultz$delta.t
schultz.test <- na.omit(schultz)

f.sd <- sd(schultz.test$f)


fit <- readRDS('../fits/mouse_3.rds')

bins <- seq(20, 38, by=1.5)
binned.test <- cut(schultz.test$age, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)


fdot.plot <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = schultz.test[n, 'sex']) %>%
    mutate(f = schultz.test[n, 'f']) %>%
    mutate(mouse = schultz.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r = -dotf*lambda_r- f*deriv_r) %>%
    mutate(term.d = -dotf*lambda_d + (1-f)*deriv_d) %>%
    mutate(ratio = term.d - term.r) %>%
    ungroup() %>%
    group_by(sex, age, .draw) %>%
    summarize(ratio = mean(ratio)) %>%
    ungroup() %>%
    group_by(sex, age) %>%
    median_hdci(ratio, .width = 0.95)


terms.schultz <- ggplot() +
    geom_errorbar(data=fdot.plot, mapping=aes(x=age, y=ratio, ymin=.lower, ymax=.upper, fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=fdot.plot, mapping=aes(x=age, y=ratio), alpha=1, size=1.75, color='white') +
    geom_point(data=fdot.plot, mapping=aes(x=age, y=ratio, color=sex), alpha=1, size=1.25) +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0, 0.06, 0, 0), "cm"))+
    facet_grid(~sex) +
    labs(y='', x='Age (months)',
         color='', fill='') +
    scale_color_nejm() +
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36),
                       label = c("22", "24", "26", "28", "30", "32", "34", "36")) + ggtitle('c) Mouse dataset 3 (Schultz et al. 2020)') +
    scale_color_manual(values=rev(sex.palette)) +
    scale_fill_manual(values=rev(sex.palette)) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1)))


img <- plot_grid(terms.enalapril, terms.exercise, terms.schultz, ncol=3,nrow=1, rel_widths = c(1.675/5, 1.675/5, 1.65/5))

#ggsave("../plots/mice_control_rates.pdf", img, width=8, height=4.5)
ggsave("../plots/mice_curvature_control_test.pdf", img, width=8, height=1.5)
