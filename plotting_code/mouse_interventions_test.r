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


f.sd = sd(enalapril.test$f)

# get stan fit
fit <- readRDS('../fits/mouse_1.rds')

bins <- seq(0, 20, by=1)
binned.test <- cut(enalapril.test$time, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

full.repair.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(f = enalapril.test[n,'f']) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_rr = deriv_r + dotf*deriv_r_f)  %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(lambda_r = mean(lambda_r), deriv_r = mean(deriv_rr)) %>%
    ungroup() %>% filter(age <= 8.5) %>% as.data.frame()


full.damage.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(f = enalapril.test[n,'f']) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_dd = deriv_d + dotf*deriv_d_f) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(lambda_d = mean(lambda_d), deriv_d = mean(deriv_dd)) %>%
    ungroup() %>% filter(age <= 8.5) %>% as.data.frame()


# rates
repair.diff.stats <- full.repair.plot.stats[full.repair.plot.stats$treatment == 'Control',]
repair.diff.stats$diff <- full.repair.plot.stats[full.repair.plot.stats$treatment == 'Enalapril',]$lambda_r - full.repair.plot.stats[full.repair.plot.stats$treatment == 'Control',]$lambda_r
repair.diff.stats <- repair.diff.stats %>% group_by(sex, treatment, age) %>% median_hdci(diff, .width=0.95)

damage.diff.stats <- full.damage.plot.stats[full.damage.plot.stats$treatment == 'Control',]
damage.diff.stats$diff <- full.damage.plot.stats[full.damage.plot.stats$treatment == 'Enalapril',]$lambda_d - full.damage.plot.stats[full.damage.plot.stats$treatment == 'Control',]$lambda_d
damage.diff.stats <- damage.diff.stats %>% group_by(sex, treatment, age) %>% median_hdci(diff, .width=0.95)


# derivs
repair.diff.deriv.stats <- full.repair.plot.stats[full.repair.plot.stats$treatment == 'Control',]
repair.diff.deriv.stats$diff <- full.repair.plot.stats[full.repair.plot.stats$treatment == 'Enalapril',]$deriv_r - full.repair.plot.stats[full.repair.plot.stats$treatment == 'Control',]$deriv_r
repair.diff.deriv.stats <- repair.diff.deriv.stats %>% group_by(sex, treatment, age) %>% median_hdci(diff, .width=0.95)

damage.diff.deriv.stats <- full.damage.plot.stats[full.damage.plot.stats$treatment == 'Control',]
damage.diff.deriv.stats$diff <- full.damage.plot.stats[full.damage.plot.stats$treatment == 'Enalapril',]$deriv_d - full.damage.plot.stats[full.damage.plot.stats$treatment == 'Control',]$deriv_d
damage.diff.deriv.stats <- damage.diff.deriv.stats %>% group_by(sex, treatment, age) %>% median_hdci(diff, .width=0.95)

enalapril.repair <- ggplot() +
    geom_errorbar(data=repair.diff.stats, mapping=aes(x=age, y=diff,
                                                             color=sex, ymin=.lower, ymax=.upper, fill=sex), alpha=0.5) +
    #geom_line(data=repair.diff.stats, mapping=aes(x=age, y=diff, fill=sex, color=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=repair.diff.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.25 )+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='Difference in repair rate', x='Time since intervention (months)',
         color='', fill='') +   scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) + ggtitle('c) Mouse dataset 1 (Keller et al. 2019)')


enalapril.damage <- ggplot() +
    #geom_line(data=full.damage.plot, mapping=aes(x=age, y=lambda_d, color=treatment, group=interaction(treatment,.draw)), alpha=0.05) +
    geom_errorbar(data=damage.diff.stats, mapping=aes(x=age, y=diff, ymin=.lower, ymax=.upper,fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=damage.diff.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=damage.diff.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0.01), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='Difference in damage rate', x='Time since intervention (months)',
         color='', fill='') +xlim(0, 10.25) +     scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette)



enalapril.deriv.repair <- ggplot() +
    geom_errorbar(data=repair.diff.deriv.stats, mapping=aes(x=age, y=diff, ymin=.lower, ymax=.upper, fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=repair.diff.deriv.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=repair.diff.deriv.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.position = "none",
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='Difference in repair rate slope', x='Time since intervention (months)',
         color='', fill='') +     scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette)


enalapril.deriv.damage <- ggplot() +
    geom_errorbar(data=damage.diff.deriv.stats, mapping=aes(x=age, y=diff, ymin=.lower, ymax=.upper, fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=damage.diff.deriv.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=damage.diff.deriv.stats, mapping=aes(x=age, y=diff, fill=sex,color=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.position = "none",
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='Difference in damage rate slope', x='Time since intervention (months)',
         color='', fill='') +xlim(0, 10.25) +     scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette)




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

bins <- seq(0, 15, by=0.5)
binned.test <- cut(exercise.test$time, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)


full.repair.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(f = exercise.test[n,'f']) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(age = binned.test[n]) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_rr = deriv_r + dotf*deriv_r_f)  %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_r = mean(lambda_r), deriv_r = mean(deriv_rr)) %>%
    ungroup()


full.damage.plot.stats <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(f = exercise.test[n,'f']) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(age = binned.test[n]) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(deriv_dd = deriv_d + dotf*deriv_d_f)  %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_d = mean(lambda_d), deriv_d = mean(deriv_dd)) %>%
    ungroup()

# rates
repair.diff.stats <- full.repair.plot.stats[full.repair.plot.stats$exercise == 'no',]
repair.diff.stats$diff <- full.repair.plot.stats[full.repair.plot.stats$exercise == 'yes',]$lambda_r - full.repair.plot.stats[full.repair.plot.stats$exercise == 'no',]$lambda_r
repair.diff.stats <- repair.diff.stats %>% group_by(sex, exercise, age) %>% median_hdci(diff, .width=0.95)

damage.diff.stats <- full.damage.plot.stats[full.damage.plot.stats$exercise == 'no',]
damage.diff.stats$diff <- full.damage.plot.stats[full.damage.plot.stats$exercise == 'yes',]$lambda_d - full.damage.plot.stats[full.damage.plot.stats$exercise == 'no',]$lambda_d
damage.diff.stats <- damage.diff.stats %>% group_by(sex, exercise, age) %>% median_hdci(diff, .width=0.95)


# derivs
repair.diff.deriv.stats <- full.repair.plot.stats[full.repair.plot.stats$exercise == 'no',]
repair.diff.deriv.stats$diff <- full.repair.plot.stats[full.repair.plot.stats$exercise == 'yes',]$deriv_r - full.repair.plot.stats[full.repair.plot.stats$exercise == 'no',]$deriv_r
repair.diff.deriv.stats <- repair.diff.deriv.stats %>% group_by(sex, exercise, age) %>% median_hdci(diff, .width=0.95)

damage.diff.deriv.stats <- full.damage.plot.stats[full.damage.plot.stats$exercise == 'no',]
damage.diff.deriv.stats$diff <- full.damage.plot.stats[full.damage.plot.stats$exercise == 'yes',]$deriv_d - full.damage.plot.stats[full.damage.plot.stats$exercise == 'no',]$deriv_d
damage.diff.deriv.stats <- damage.diff.deriv.stats %>% group_by(sex, exercise, age) %>% median_hdci(diff, .width=0.95)

exercise.repair <- ggplot() +
    geom_errorbar(data=repair.diff.stats, mapping=aes(x=age, y=diff, ymin=.lower, ymax=.upper, fill=sex,color=sex), alpha=0.5) +
    #geom_line(data=repair.diff.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=repair.diff.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0.01, 0, 0, 0), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='', x='Time since intervention (months)',
         color='', fill='') + 
    scale_x_continuous(breaks=c(0, 1, 2, 3),label = c("0", "1", "2", "3"))+ xlim(0, 3.25)+     scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) + ggtitle("d) Mouse dataset 2 (Bisset et al. 2021)")

exercise.damage <- ggplot() +
    geom_errorbar(data=damage.diff.stats, mapping=aes(x=age, y=diff, ymin=.lower, ymax=.upper,fill=sex, color=sex), alpha=0.5) +
    #geom_line(data=damage.diff.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=damage.diff.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='', x='Time since intervention (months)',
         color='', fill='') + scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3"))+ xlim(0, 3.25)+     scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette)



exercise.deriv.repair <- ggplot() +
    geom_errorbar(data=repair.diff.deriv.stats, mapping=aes(x=age, y=diff,ymin=.lower, ymax=.upper,fill=sex,color=sex), alpha=0.5) +
    #geom_line(data=repair.diff.deriv.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=repair.diff.deriv.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0.01, 0, 0, 0), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='', x='Time since intervention (months)',
         color='', fill='') + #Repair rate (/month)
    scale_x_continuous(breaks=c(0, 1, 2, 3),label = c("0", "1", "2", "3"))+ xlim(0, 3.25)+     scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette)

exercise.deriv.damage <- ggplot() +
    geom_errorbar(data=damage.diff.deriv.stats, mapping=aes(x=age, y=diff,ymin=.lower, ymax=.upper, fill=sex,color=sex), alpha=0.5) +
    #geom_line(data=damage.diff.deriv.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.75, color='white') +
    geom_point(data=damage.diff.deriv.stats, mapping=aes(x=age, y=diff, color=sex,fill=sex), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    facet_grid(~sex, labeller = as_labeller(c('F'='Female', 'M'='Male'))) +
    labs(y='', x='Time since intervention (months)',
         color='', fill='') + 
    scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3"))+ xlim(0, 3.25)+    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette)


img.enalapril <- plot_grid(enalapril.repair, enalapril.damage, enalapril.deriv.repair, enalapril.deriv.damage, ncol=1,nrow=4, rel_heights=c(1.1/4, 0.96666666666/4, 0.96666666666/4, 0.96666666666/4))

img.exercise <- plot_grid(exercise.repair, exercise.damage, exercise.deriv.repair, exercise.deriv.damage, ncol=1,nrow=4, rel_heights=c(1.1/4, 0.96666666666/4, 0.96666666666/4, 0.96666666666/4)) #exercise.f

img <- grid.arrange(img.enalapril, img.exercise, ncol=2,nrow=1)


ggsave("../plots/mice_combined_rates_tests.pdf", img, width=8, height=4.5)
