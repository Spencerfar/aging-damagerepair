library(posterior)
library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(tidybayes)
library(latex2exp)
library(cowplot)
source("../utils/functions.r")
source("../utils/palettes.r")

elsa <- read.csv('../datasets/human_data.csv', header = TRUE, sep = ",")

elsa$id <- as.character(elsa$id)
elsa$time <- elsa$age - elsa$baseline.age

elsa <- elsa %>% group_by(id) %>% filter( (baseline.age >= 50)  & (baseline.age < 90))
elsa$f <- elsa$n/elsa$N
elsa$wealth <- log(elsa$wealth + mean(elsa$wealth))

bins <- seq(50, 90, by=10)
elsa <- 
elsa %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE))
elsa$baseline.bin <- as.factor(elsa$baseline.bin)
elsa$sex <- as.factor(elsa$sex)
levels(elsa$sex) <- c('Male', 'Female')

elsa <- na.omit(elsa)

elsa$prepair <- (elsa$repair/elsa$n)/elsa$delta.t
elsa$pdamage <- (elsa$damage/(elsa$N - elsa$n))/elsa$delta.t

bins <- seq(20, 105, by=3)

fit <- as_draws_df(readRDS("../fits/long_human_fit.RDS")$draws(c("lambda_r", "lambda_d")))

baseline.bins <- seq(50, 90, by=4)
binned.age <- cut(elsa$age, bins, include.lowest=TRUE)
binned.age <- midpoints(binned.age)+1.5

repair.wealth.corr.med <- fit %>%
    spread_draws(lambda_r[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa[n, 'sex']$sex) %>%
    mutate(id = elsa[n, 'id']$id) %>%
    mutate(age = binned.age[n]) %>%
    mutate(wealth = elsa[n, 'wealth']$wealth) %>%
    mutate(baseline.bin = elsa[n,'baseline.bin']$baseline.bin) %>%
    filter(elsa[n,'n'] > 0) %>%
    ungroup() %>%
    group_by(id, sex, age, baseline.bin, .draw) %>%
    summarize(lambda_r = mean(lambda_r, na.rm=TRUE),
              wealth = mean(wealth)) %>%
    ungroup() %>%
    group_by(sex, age, baseline.bin, .draw) %>%
    filter(length(id) >= 3) %>%
    summarize(corr_r = cor(wealth, lambda_r, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, age, baseline.bin) %>%
    median_hdci(corr_r, .width=0.95)

repair.wealth.corr.med %>%as.data.frame()

repair.corr.plot <- ggplot() +
    geom_lineribbon(data = repair.wealth.corr.med, mapping = aes(x = age, y = corr_r, ymin = .lower, ymax=.upper),
                    alpha=0.5, fill = rates.palette[2]) +
    geom_line(data=repair.wealth.corr.med, mapping=aes(x=age, y=corr_r),
              alpha=1, size=1.25, color='white') + 
    geom_line(data=repair.wealth.corr.med, mapping=aes(x=age, y=corr_r),
              alpha=1, size=0.75, color = rates.palette[2]) +
    geom_hline(yintercept = 0, linetype="dotted") +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    labs(x="Age (years)", y="Repair rate wealth-correlation") +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + ggtitle("e) ELSA humans")

rm(repair.wealth.corr.med)
gc()


damage.wealth.corr.med <- fit %>%
    spread_draws(lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(sex = elsa[n, 'sex']$sex) %>%
    mutate(id = elsa[n, 'id']$id) %>%
    mutate(wealth = elsa[n, 'wealth']$wealth) %>%
    mutate(age = binned.age[n]) %>%    
    mutate(baseline.bin = elsa[n,'baseline.bin']$baseline.bin) %>%
    filter(elsa[n,'n'] > 0) %>%
    ungroup() %>%
    group_by(id, sex, age, baseline.bin, .draw) %>%
    summarize(lambda_d = mean(lambda_d, na.rm=TRUE),
              wealth = mean(wealth)) %>%
    ungroup() %>%
    group_by(sex, age, baseline.bin, .draw) %>%
    filter(length(id) >= 3) %>%
    summarize(corr_d = cor(wealth, lambda_d, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, age, baseline.bin) %>%
    median_hdci(corr_d, .width=0.95)

damage.corr.plot <- ggplot() +
    geom_lineribbon(data = damage.wealth.corr.med, mapping = aes(x = age, y = corr_d, ymin = .lower, ymax=.upper),
                    alpha=0.5, fill = rates.palette[1]) +
    geom_line(data=damage.wealth.corr.med, mapping=aes(x=age, y=corr_d),
              alpha=1, size=1.25, color='white') + 
    geom_line(data=damage.wealth.corr.med, mapping=aes(x=age, y=corr_d),
              alpha=1, size=0.75, color = rates.palette[1]) +
    geom_hline(yintercept = 0, linetype="dotted") +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    labs(x="Age (years)", y="Damage rate wealth-correlation") +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

img <- grid.arrange(repair.corr.plot, damage.corr.plot, nrow=1)

ggsave("../plots/human_wealth_rates_corr.pdf", img, width=8, height=2.)
