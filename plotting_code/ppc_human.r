library(rstan)
library(cmdstanr)
library(posterior)
set_cmdstan_path("~/Downloads/cmdstan-2.26.0/")
library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(tidybayes)
library(tidyr)
library(cowplot)
source("../utils/functions.r")
source("../utils/palettes.r")

elsa <- read.csv('../datasets/human_data.csv', header = TRUE, sep = ",")
elsa$id <- as.character(elsa$id)
elsa$time <- elsa$age - elsa$baseline.age
elsa <- elsa %>% group_by(id) %>% filter( (baseline.age >= 50)  & (baseline.age < 90))
elsa.test <- na.omit(elsa)


elsa$f <- elsa$n/elsa$N
elsa$wealth <- log(elsa$wealth + mean(elsa$wealth))

fit <- as_draws_df(readRDS("../fits/long_human_fit.RDS")$draws(c("sampled_repair", "sampled_damage", "sampled_n")))

sampled.data <- fit %>%
    spread_draws(sampled_repair[n], sampled_damage[n], sampled_n[n])

bins <- seq(0,50,by=1)

repair <- ggplot() +
    geom_histogram(data=elsa.test, aes(x=repair,y=..density..), breaks=seq(0,23, by = 1), fill = rates.palette[2])  +
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
        plot.title=element_text(size=8,hjust=0.5),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.1, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle('d) ELSA humans (Phelps et al. 2020)') +
    scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(0, 10, by = 2))+ xlim(0, 10)


damage <- ggplot() +
    geom_histogram(data=elsa.test, aes(x=damage,y=..density..), breaks=seq(0,23, by = 1), fill = rates.palette[1])  +
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
        plot.margin = unit(c(0.1, 0, 0, 0.), "cm"),
        legend.title = element_text(size = 7)) + ggtitle(' ') +
    scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(0, 10, by = 2)) + xlim(0, 10)


count <- ggplot() +
    geom_histogram(data=elsa, aes(x=n,y=..density..), breaks=seq(0,23, by = 1), alpha=0.7) + xlim(0, 15) +
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
        plot.margin = unit(c(0.1, 0, 0, 0.), "cm"),
        legend.title = element_text(size = 7)) + ggtitle(' ')

img <- grid.arrange(repair, damage, count, nrow=1, ncol=3)

ggsave("../plots/ppc_human.pdf", img, width=6, height=1)
