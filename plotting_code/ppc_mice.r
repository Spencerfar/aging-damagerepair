library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(tidybayes)
library(tidyr)
library(cowplot)
library(survival)
library(survminer)
library(ggfortify)
library(pammtools)
source("../utils/functions.r")
source("../utils/palettes.r")


mice <- read.csv('../datasets/enalapril_data.csv', header = TRUE, sep = ",")
mice.surv <- read.csv('../datasets/enalapril_surv_data.csv', header = TRUE, sep = ",")
mice.surv$event.time <- mice.surv$death.age-mice.surv$baseline.age

mice$age <- mice$time + mice$baseline.age
mice$treatment <- as.factor(mice$treatment)
levels(mice$treatment) <- c('Control', 'Enalapril')
mice$sex <- as.factor(mice$sex)
levels(mice$sex) <- c('Female', 'Male')

mice$prepair <- (mice$repair/mice$n)/mice$delta.t
mice$pdamage <- (mice$damage/(mice$N - mice$n))/mice$delta.t

mice.test <- na.omit(mice)


fit <- readRDS('../fits/mouse_1.rds')

sampled.data <- fit %>%
    spread_draws(sampled_repair[n], sampled_damage[n], sampled_n[n])

bins <- seq(0,124)

enalapril.repair <- ggplot() +
    geom_histogram(data=mice.test, aes(x=repair,y=..density..), breaks=seq(0,124, by = 2), fill = rates.palette[2]) + xlim(0, 25) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_repair, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_repair, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Repair count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=7,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle('a) Mouse dataset 1 (Keller et al. 2019)')


enalapril.damage <- ggplot() +
    geom_histogram(data=mice.test, aes(x=damage,y=..density..), breaks=seq(0,124, by = 2), fill = rates.palette[1]) + xlim(0, 25) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_damage, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_damage, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Damage count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle(' ')


enalapril.count <- ggplot() +
    geom_histogram(data=mice, aes(x=n,y=..density..), breaks=seq(0,124, by = 2), alpha=0.7) + xlim(0, 50) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_n, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_n, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Deficit count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle(' ')


km <- fortify(survfit(Surv(event.time, status) ~ 1, data=mice.surv))[,c('time', 'surv')]

bins <- seq(0, 15, by=0.5)
data <- fit %>%
    spread_draws(output_age[n], output_status[n]) %>%
    group_by(.draw) %>%
    mutate(event.time = output_age * sd(mice$time) + mean(mice$time))%>%
    mutate(status = output_status * mice.surv[n, 'status']) %>% 
    summarize(survival = fortify(survfit(Surv(event.time, status) ~ 1))$surv, time = fortify(survfit(Surv(event.time, status) ~ 1))$time) %>%
    filter(time < 12) %>%
    mutate(binned.time =  midpoints(cut(time, bins, include.lowest = TRUE))) %>%
    ungroup() %>%
    group_by(binned.time) %>%
    median_hdci(survival, .width=0.95) %>%
    as.data.frame()


enalapril.surv <- ggplot() + geom_stepribbon(data=data, mapping=aes(x=binned.time, y=survival, ymin=.lower, ymax=.upper), alpha=0.5) + geom_step(data=data, mapping=aes(x=binned.time, y=survival),size=0.75) + geom_step(data=km, mapping=aes(x=time, y=surv), color = rates.palette[3],size=0.75) + labs(y = 'Survival', x = 'Time') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=7,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')

enalapril.img <- grid.arrange(enalapril.repair, enalapril.damage, enalapril.count, enalapril.surv, nrow=1, ncol=4)


mice <- read.csv('../datasets/exercise_data.csv', header = TRUE, sep = ",")
mice.surv <- read.csv('../datasets/exercise_surv_data.csv', header = TRUE, sep = ",")
mice.surv$event.time <- mice.surv$death.age#-mice.surv$baseline.age

mice$age <- mice$time + mice$baseline.age
mice$sex <- as.factor(mice$sex)
mice$exercise <- as.factor(mice$exercise)
levels(mice$exercise) <- c('Control', 'Exercise')
levels(mice$sex) <- c('Female', 'Male')
mice$prepair <- (mice$repair/mice$n)/mice$delta.t
mice$pdamage <- (mice$damage/(mice$N - mice$n))/mice$delta.t

mice.test <- na.omit(mice)


fit <- readRDS('../fits/mouse_2.rds')

sampled.data <- fit %>%
    spread_draws(sampled_repair[n], sampled_damage[n], sampled_n[n])

bins <- seq(0,124)

exercise.repair <- ggplot() +
    geom_histogram(data=mice.test, aes(x=repair,y=..density..), breaks=seq(0,124, by = 2), fill = rates.palette[2]) + xlim(0, 25) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_repair, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_repair, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Repair count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=7,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle('b) Mouse dataset 2 (Bisset et al. 2021)')


exercise.damage <- ggplot() +
    geom_histogram(data=mice.test, aes(x=damage,y=..density..), breaks=seq(0,124, by = 2), fill = rates.palette[1]) + xlim(0, 25) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_damage, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_damage, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Damage count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')


exercise.count <- ggplot() +
    geom_histogram(data=mice, aes(x=n,y=..density..), breaks=seq(0,124, by = 2), alpha=0.7) + xlim(0, 50) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_n, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_n, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Deficit count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')



#mice <- na.omit(read.csv('../../cleaned_data/Cleaned_mice_exercise_data_step_May28.csv', header = TRUE, sep = ","))
#mice.surv <- read.csv('../../cleaned_data/Cleaned_mice_exercise_surv_data_step_May28.csv', header = TRUE, sep = ",")



km <- fortify(survfit(Surv(event.time, status) ~ 1, data=mice.surv))[,c('time', 'surv')]

bins <- seq(0, 4, by=0.25)
data <- fit %>%
    spread_draws(output_age[n], output_status[n]) %>%
    group_by(.draw) %>%
    mutate(event.time = output_age * sd(mice$time) + mean(mice$time))%>%
    mutate(status = output_status * mice.surv[n, 'status']) %>% 
    summarize(survival = fortify(survfit(Surv(event.time, status) ~ 1))$surv, time = fortify(survfit(Surv(event.time, status) ~ 1))$time) %>%
    filter(time < 4) %>%
    mutate(binned.time =  midpoints(cut(time, bins, include.lowest = TRUE))) %>%
    ungroup() %>%
    group_by(binned.time) %>%
    median_hdci(survival, .width=0.95) %>%
    as.data.frame()

exercise.surv <- ggplot() + geom_stepribbon(data=data, mapping=aes(x=binned.time, y=survival, ymin=.lower, ymax=.upper), alpha=0.5) + geom_step(data=data, mapping=aes(x=binned.time, y=survival),size=0.75)  + geom_step(data=km, mapping=aes(x=time, y=surv), color = rates.palette[3],size=0.75) + labs(y = 'Survival', x = 'Time') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=7,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')

exercise.img <- grid.arrange(exercise.repair, exercise.damage, exercise.count, exercise.surv, nrow=1, ncol=4)


mice <- read.csv('../datasets/schultz_data.csv', header = TRUE, sep = ",")
mice.surv <- read.csv('../datasets/schultz_surv_data.csv', header = TRUE, sep = ",")
mice.surv$event.time <- mice.surv$death.age-mice.surv$baseline.age

mice$age <- mice$time + mice$baseline.age
mice$prepair <- (mice$repair/mice$n)/mice$delta.t
mice$pdamage <- (mice$damage/(mice$N - mice$n))/mice$delta.t
mice.test <- na.omit(mice)

fit <- readRDS('../fits/mouse_3.rds')


sampled.data <- fit %>%
    spread_draws(sampled_repair[n], sampled_damage[n], sampled_n[n])

bins <- seq(0,124)

schultz.repair <- ggplot() +
    geom_histogram(data=mice.test, aes(x=repair,y=..density..), breaks=seq(0,124, by = 2), fill = rates.palette[2]) + xlim(0, 25) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_repair, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_repair, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Repair count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=7,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle('c) Mouse dataset 3 (Schultz et al. 2020)')


schultz.damage <- ggplot() +
    geom_histogram(data=mice.test, aes(x=damage,y=..density..), breaks=seq(0,124, by = 2), fill = rates.palette[1]) + xlim(0, 25) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_damage, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_damage, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Damage count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')

bins <- seq(0,124, by=2)

schultz.count <- ggplot() +
    geom_histogram(data=mice, aes(x=n,y=..density..), breaks=seq(0,124, by = 4), alpha=0.7) + xlim(10, 70) +
    geom_pointrange(data=sampled.data %>% group_by(.draw) %>%
                        summarize(prop = hist(sampled_n, breaks=bins, plot=FALSE)$density,
                                  bins = hist(sampled_n, breaks=bins, plot=FALSE)$mids) %>%
                        ungroup() %>%
                        group_by(bins) %>%
                        median_hdci(prop, .width=0.95),
                    aes(x=bins, y=prop, ymin=.lower, ymax=.upper), size=0.2) +
    labs(y = '', x = 'Deficit count') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')



#mice <- na.omit(read.csv('../../cleaned_data/Cleaned_mice_Schultz_data_step_May28.csv', header = TRUE, sep = ","))
#mice.surv <- read.csv('../../cleaned_data/Cleaned_mice_surv_Schultz_data_step_May28.csv', header = TRUE, sep = ",")


km <- fortify(survfit(Surv(event.time, status) ~ 1, data=mice.surv))[,c('time', 'surv')]

bins <- seq(0, 15, by=0.5)
data <- fit %>%
    spread_draws(output_age[n], output_status[n]) %>%
    group_by(.draw) %>%
    mutate(event.time = output_age * sd(mice$time) + mean(mice$time))%>%
    mutate(status = output_status * mice.surv[n, 'status']) %>% 
    summarize(survival = fortify(survfit(Surv(event.time, status) ~ 1))$surv, time = fortify(survfit(Surv(event.time, status) ~ 1))$time) %>%
    filter(time < 15) %>%
    mutate(binned.time =  midpoints(cut(time, bins, include.lowest = TRUE))) %>%
    ungroup() %>%
    group_by(binned.time) %>%
    median_hdci(survival, .width=0.95) %>%
    as.data.frame()


schultz.surv <- ggplot() + geom_stepribbon(data=data, mapping=aes(x=binned.time, y=survival, ymin=.lower, ymax=.upper), alpha=0.5) + geom_step(data=data, mapping=aes(x=binned.time, y=survival),size=0.75)  + geom_step(data=km, mapping=aes(x=time, y=surv), color = rates.palette[3], size=0.75) + labs(y = 'Survival', x = 'Time') + theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=7,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7))  + ggtitle(' ')

schultz.img <- grid.arrange(schultz.repair, schultz.damage, schultz.count, schultz.surv, nrow=1, ncol=4)

img <- grid.arrange(enalapril.img, exercise.img, schultz.img, nrow=3, ncol=1)


ggsave("../plots/ppc_mice.pdf", img, width=6, height=3)
