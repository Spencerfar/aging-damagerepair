library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(tidybayes)
library(latex2exp)
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


bins <- seq(0, 15, by=1)
binned <- cut(enalapril$time, bins, include.lowest = TRUE)
binned <- midpoints(binned)
enalapril$time.bin <- binned


binned.data.rates <- 
enalapril %>% 
  group_by(time.bin, treatment, sex) %>%
  summarise(mean.damage = mean(pdamage,na.rm=TRUE), 
            mean.repair = mean(prepair,na.rm=TRUE),
            mean.f = mean(f, na.rm=TRUE),
            se.repair = se(prepair),
            se.damage = se(pdamage),
            se.f = se(f),
	    repair.num = length(prepair),
	    damage.num = length(pdamage),
        f.num = length(f)) %>%
mutate(age = time.bin) %>% as.data.frame()

binned.data.rates$lower.repair <- binned.data.rates$mean.repair-binned.data.rates$se.repair
binned.data.rates$upper.repair <- binned.data.rates$se.repair+binned.data.rates$mean.repair

binned.data.rates$lower.damage <- binned.data.rates$mean.damage-binned.data.rates$se.damage
binned.data.rates$upper.damage <- binned.data.rates$se.damage+binned.data.rates$mean.damage

binned.data.rates$lower.f <- binned.data.rates$mean.f-binned.data.rates$se.f
binned.data.rates$upper.f <- binned.data.rates$se.f+binned.data.rates$mean.f

# get stan fit
fit <- readRDS('../fits/mouse_1.rds')

binned.test <- cut(enalapril.test$time, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

full.repair.plot <- fit %>%
    spread_draws(lambda_r[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(lambda_r = mean(lambda_r)) %>%
    ungroup()


full.damage.plot <- fit %>%
    spread_draws(lambda_d[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(lambda_d = mean(lambda_d)) %>%
    ungroup()


binned <- cut(enalapril$time, bins, include.lowest = TRUE)
binned <- midpoints(binned)

full.f.plot <- fit %>%
    spread_draws(sampled_n[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = enalapril[n, 'sex']) %>%
    mutate(treatment = enalapril[n, 'treatment']) %>%
    mutate(age = binned[n]) %>%
    mutate(mean_f = sampled_n[n]/124.0) %>%
    mutate(mouse = enalapril[n,'mouse']) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(mean_f = mean(mean_f, na.rm=TRUE)) %>%
    ungroup()

levels(full.repair.plot$treatment) <- c('Control', 'Enalapril')
levels(full.damage.plot$treatment) <- c('Control', 'Enalapril')
levels(full.f.plot$treatment) <- c('Control', 'Enalapril')
levels(binned.data.rates$treatment) <- c('Control', 'Enalapril')

# plot repair
enalapril.repair <- ggplot() +
    geom_point(binned.data.rates %>% filter(repair.num >= 5), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                              ymax=upper.repair, group=treatment, color=treatment, shape=treatment), size=2) +
    geom_errorbar(binned.data.rates %>% filter(repair.num >= 5), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                             ymax=upper.repair, group=treatment, color=treatment), width=0.1) +
    geom_line(data=full.repair.plot, mapping=aes(x=age, y=lambda_r, color=treatment, group=interaction(.draw, treatment)), alpha=0.05) +
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
    labs(y='Repair rate (/month)', x='Time since intervention (months)',
         color='', fill='', shape='') + 
    scale_fill_d3() +
    scale_color_d3() + ggtitle("a) Mouse dataset 1 (Keller et al. 2019)")  + xlim(0, 10.25)

# plot damage
enalapril.damage <- ggplot() +
    geom_point(binned.data.rates %>% filter(damage.num >= 5), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                              ymax=upper.damage, group=treatment, color=treatment, shape=treatment), size=2) +
    geom_errorbar(binned.data.rates %>% filter(damage.num >= 5), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                             ymax=upper.damage, group=treatment, color=treatment), width=0.1) +
    geom_line(data=full.damage.plot, mapping=aes(x=age, y=lambda_d, color=treatment, group=interaction(.draw, treatment)), alpha=0.05) +
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
    labs(y='Damage rate (/month)', x='Time since intervention (months)',
         color='', fill='', shape='') + 
    scale_fill_d3() +
    scale_color_d3() + xlim(0, 10.25)

# plot Frailty Index
enalapril.f <- ggplot() +
    geom_point(binned.data.rates %>% filter(f.num >= 5), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                              ymax=upper.f, group=treatment, color=treatment, shape=treatment), size=2) +
    geom_errorbar(binned.data.rates %>% filter(f.num >= 5), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                        ymax=upper.f, group=treatment, color=treatment), width=0.1) +
    geom_line(data=full.f.plot, mapping=aes(x=age, y=mean_f, color=treatment, group=interaction(.draw, treatment)), alpha=0.05) +
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
    labs(y='Frailty Index', x='Time since intervention (months)',
         color='', fill='', shape='') +
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

bins <- seq(0, 15, by=0.5)
binned <- cut(exercise$time, bins, include.lowest = TRUE)
binned <- midpoints(binned)
exercise$time.bin <- binned

binned.data.rates <- 
exercise %>% 
  group_by(time.bin, exercise, sex) %>%
  summarise(mean.damage = mean(pdamage,na.rm=TRUE), 
            mean.repair = mean(prepair,na.rm=TRUE),
            mean.f = mean(f, na.rm=TRUE),
            se.repair = se(prepair),
            se.damage = se(pdamage),
            se.f = se(f),
	    repair.num = length(prepair),
	    damage.num = length(pdamage),
        f.num = length(f)) %>%
  mutate(age = time.bin)


binned.data.rates$lower.repair <- binned.data.rates$mean.repair-binned.data.rates$se.repair
binned.data.rates$upper.repair <- binned.data.rates$se.repair+binned.data.rates$mean.repair

binned.data.rates$lower.damage <- binned.data.rates$mean.damage-binned.data.rates$se.damage
binned.data.rates$upper.damage <- binned.data.rates$se.damage+binned.data.rates$mean.damage

binned.data.rates$lower.f <- binned.data.rates$mean.f-binned.data.rates$se.f
binned.data.rates$upper.f <- binned.data.rates$se.f+binned.data.rates$mean.f


binned.data.rates <- na.omit(binned.data.rates)


# get stan fit
fit <- readRDS('../fits/mouse_2.rds')

bins <- seq(0, 15, by=0.5)
binned.test <- cut(exercise.test$time, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)


full.repair.plot <- fit %>%
    spread_draws(lambda_r[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_r = mean(lambda_r)) %>%
    ungroup()  %>% as.data.frame()

full.damage.plot <- fit %>%
    spread_draws(lambda_d[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    #filter(survival >= 0.5) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    #mutate(age = exercise.test[n, 'real.age']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_d = mean(lambda_d)) %>%
    ungroup()

bins <- seq(0, 16, by=0.5)
binned <- cut(exercise$time, bins, include.lowest = TRUE)
binned <- midpoints(binned)

full.f.plot <- fit %>%
    spread_draws(sampled_n[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise[n, 'sex']) %>%
    mutate(exercise = exercise[n, 'exercise']) %>%
    mutate(mean_f = sampled_n[n]/124.0) %>%
    mutate(age = binned[n]) %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(mean_f = mean(mean_f)) %>%
    ungroup()


levels(full.repair.plot$exercise) <- c('Control', 'Exercise')
levels(full.damage.plot$exercise) <- c('Control', 'Exercise')
levels(full.f.plot$exercise) <- c('Control', 'Exercise')


# plot repair
exercise.repair <- ggplot() +
    geom_point(data=binned.data.rates %>% filter(repair.num >= 5), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                              ymax=upper.repair, group=exercise, color=exercise, shape=exercise), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(repair.num >= 5), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                                  ymax=upper.repair, group=exercise, color=exercise), width=0.1) +
    geom_line(data=full.repair.plot, mapping=aes(x=age, y=lambda_r, color=exercise, group=interaction(.draw, exercise)), alpha=0.05)+
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
         color='', fill='', shape='') + 
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual( values = exercise.palette) + ggtitle('b) Mouse dataset 2 (Bisset et al. 2021)')  +
    scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3")) + xlim(0, 3.25)

# plot damage
exercise.damage <- ggplot() +
    geom_point(binned.data.rates %>% filter(damage.num >= 5), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                              ymax=upper.damage, group=exercise, color=exercise, shape=exercise), size=2) +
    geom_errorbar(binned.data.rates %>% filter(damage.num >= 5), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                             ymax=upper.damage, group=exercise, color=exercise), width=0.1) +
    geom_line(data=full.damage.plot, mapping=aes(x=age, y=lambda_d, color=exercise,group=interaction(.draw, exercise)), alpha=0.05)+
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
         color='', fill='',shape='') + 
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette)  +
    scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3"))+ xlim(0, 3.25)

# plot frailty Index
exercise.f <- ggplot() +
    geom_point(binned.data.rates %>% filter(f.num >= 5), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                              ymax=upper.f, group=exercise, color=exercise, shape=exercise), size=2) +
    geom_errorbar(binned.data.rates %>% filter(f.num >= 5), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                             ymax=upper.f, group=exercise, color=exercise), width=0.1) +
    geom_line(data=full.f.plot, mapping=aes(x=age, y=mean_f, color=exercise, group=interaction(.draw, exercise)), alpha=0.05)+
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
         color='', fill='', shape='') + 
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette)  +
    scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3"))+ xlim(0, 3.25)

img.enalapril <- grid.arrange(enalapril.repair, enalapril.damage, enalapril.f, ncol=1,nrow=3)
img.exercise <- grid.arrange(exercise.repair, exercise.damage, exercise.f, ncol=1,nrow=3)

img <- grid.arrange(img.enalapril, img.exercise, ncol=2,nrow=1)

ggsave("../plots/mice_interventions.pdf", img, width=8, height=4)
