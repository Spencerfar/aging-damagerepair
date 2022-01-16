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
wealth.bins <- seq(min(elsa$wealth), 1.01*max(elsa$wealth), length.out = 7)
elsa <- 
elsa %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE)) %>%
    mutate(wealth.plotbin = midpoints(cut(wealth, wealth.bins, right=FALSE)))
elsa$wealth.bin <- as.factor(ntile(elsa$wealth, 3))
elsa$baseline.bin <- as.factor(elsa$baseline.bin)
levels(elsa$wealth.bin) <- c('poor', 'middle', 'rich')
elsa$sex <- as.factor(elsa$sex)
levels(elsa$sex) <- c('Male', 'Female')

elsa.test <- na.omit(elsa)
elsa.na <- na.omit(elsa.test)


elsa.na$prepair <- (elsa.na$repair/elsa.na$n)/elsa.na$delta.t
elsa.na$pdamage <- (elsa.na$damage/(elsa.na$N - elsa.na$n))/elsa.na$delta.t


bins <- seq(20, 105, by=2)
bins.time <- seq(0, 20, by=2)

binned.data <- elsa.na %>% 
    mutate(age = cut(age, bins, include.lowest=TRUE)) %>%
    group_by(baseline.bin, wealth.bin, age, sex) %>%
    summarise(mean.damage = mean(pdamage,na.rm=TRUE), 
              mean.repair = mean(prepair,na.rm=TRUE),
              se.repair = se(prepair),
              se.damage = se(pdamage),
              num.repair = length(prepair),
              num.damage = length(pdamage),
              count = length(unique(id))) %>%
    mutate(age = midpoints(age)+1) %>%
    as.data.frame()


binned.data$lower.repair <- binned.data$mean.repair-binned.data$se.repair
binned.data$upper.repair <- binned.data$se.repair+binned.data$mean.repair

binned.data$lower.damage <- binned.data$mean.damage-binned.data$se.damage
binned.data$upper.damage <- binned.data$se.damage+binned.data$mean.damage

binned.data <- na.omit(binned.data)

# read stan fits and combine chains
fit1 <- as_draws_df(readRDS("../fits/long_human_fit_chain1.RDS")$draws(c("lambda_r", "lambda_d", "sampled_n")))
fit2 <- as_draws_df(readRDS("../fits/long_human_fit_chain2.RDS")$draws(c("lambda_r", "lambda_d", "sampled_n")))
fit2$.chain <- 2
fit2$.draw <- max(fit1$.draw) + fit2$.draw + 1
fit <- rbind(fit1, fit2)
rm(fit1)
rm(fit2)
gc()

baseline.bins <- seq(50, 90, by=4)
binned.age <- cut(elsa.test$age, bins, include.lowest=TRUE)
binned.age <- midpoints(binned.age)+1

repair <- fit %>%
    spread_draws(lambda_r[n], n=100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']$sex) %>%
    mutate(wealth.bin = elsa.test[n, 'wealth.bin']$wealth.bin) %>%
    mutate(time = elsa.test[n, 'time']$time) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']$baseline.bin) %>%
    mutate(id = elsa.test[n, 'id']$id) %>%
    mutate(num = elsa.test[n,'n']) %>%
    filter(num > 0) %>%
    ungroup() %>%
    group_by(baseline.bin, wealth.bin, sex, age, .draw) %>%
    summarize(lambda_r=mean(lambda_r)) %>%
    ungroup() %>% as.data.frame()


damage <- fit %>%
    spread_draws(lambda_d[n], n=100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']$sex) %>%
    mutate(wealth.bin = elsa.test[n, 'wealth.bin']$wealth.bin) %>%
    mutate(time = elsa.test[n, 'time']$time) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']$baseline.bin) %>%
    mutate(id = elsa.test[n, 'id']$id) %>%
    ungroup() %>%
    group_by(baseline.bin, wealth.bin, sex, age, .draw) %>%
    summarize(lambda_d=mean(lambda_d)) %>%
    ungroup() %>% as.data.frame()

levels(binned.data$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(repair$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(damage$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')


repair.rates <- ggplot() +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    geom_point(data = binned.data, mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
               ymax=upper.repair, color=wealth.bin), size = 1) +
    geom_errorbar(data = binned.data, mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                  ymax=upper.repair, color=wealth.bin), width=0.25) +
    geom_line(data=repair, mapping=aes(x=age, y=lambda_r, color=wealth.bin, group=interaction(.draw, wealth.bin)), alpha=0.05) +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Repair rate (/year)", fill='Household wealth', color='Household wealth') +
    scale_color_manual(values=wealth.palette, 
                       labels = c('Lower tecile', 'Middle tercile', 'Upper tercile')) +
    scale_fill_manual(values=wealth.palette, 
                       labels = c('Lower tecile', 'Middle tercile', 'Upper tercile')) +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.title = element_text(size = 7)) + ggtitle("a)")


damage.rates <- ggplot() +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    geom_point(data = binned.data, mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
               ymax=upper.damage, color=wealth.bin), size = 1) +
    geom_errorbar(data = binned.data, mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                  ymax=upper.damage, color=wealth.bin), width=0.25) +
    geom_line(data=damage, mapping=aes(x=age, y=lambda_d, color=wealth.bin, group=interaction(.draw, wealth.bin)), alpha=0.05) +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Damage rate (/year)", fill='Household wealth', color='Household wealth') +
    scale_color_manual(values=wealth.palette, 
                       labels = c('Lower tercile', 'Middle tercile', 'Upper tercile')) +
    scale_fill_manual(values=wealth.palette, 
                       labels = c('Lower tercile', 'Middle tercile', 'Upper tercile')) + 
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7)) + coord_cartesian(ylim = c(0, 0.335))


img <- grid.arrange(repair.rates,damage.rates,nrow=2)
ggsave("../plots/human_wealth_rates.pdf", img, width=8, height=4)

rm(fit)
gc()

fit <- as_draws_df(readRDS("../fits/long_human_fit.RDS")$draws(c("sampled_n")))


binned.data <- elsa %>% 
    mutate(age = cut(age, bins, include.lowest=TRUE)) %>%
    group_by(baseline.bin, wealth.bin, age, sex) %>%
    summarise(mean.f = mean(f, na.rm=TRUE),
              se.f = se(f),
              count = length(unique(id))) %>%
    mutate(age = midpoints(age)+1) %>%
    as.data.frame()

binned.data$lower.f <- binned.data$mean.f-binned.data$se.f
binned.data$upper.f <- binned.data$se.f+binned.data$mean.f


binned.data <- na.omit(binned.data)

binned.age <- cut(elsa$age, bins, include.lowest=TRUE)
binned.age <- midpoints(binned.age)+1

f <- fit %>%
    spread_draws(sampled_n[n], n=100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa[n, 'sex']$sex) %>%
    mutate(wealth.bin = elsa[n, 'wealth.bin']$wealth.bin) %>%
    mutate(time = elsa[n, 'time']$time) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa[n,'baseline.bin']$baseline.bin) %>%
    mutate(id = elsa[n, 'id']$id) %>%
    ungroup() %>%
    group_by(baseline.bin, wealth.bin, sex, age, .draw) %>%
    summarize(mean_f=mean(sampled_n)/23.0) %>%
    ungroup() %>% as.data.frame()


levels(f$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(binned.data$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')


f.plot <- ggplot() +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    geom_point(data = binned.data, mapping=aes(x=age, y=mean.f, ymin=lower.f, 
               ymax=upper.f, color=wealth.bin), size = 1) +
    geom_errorbar(data = binned.data, mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                  ymax=upper.f, color=wealth.bin), width=0.25) +
    geom_line(data=f, mapping=aes(x=age, y=mean_f, color=wealth.bin, group=interaction(.draw, wealth.bin)), alpha=0.05) +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Frailty Index", fill='Household wealth', color='Household wealth') +
    scale_color_manual(values=wealth.palette, 
                       labels = c('Lower tercile', 'Middle tercile', 'Upper tercile')) +
    scale_fill_manual(values=wealth.palette, 
                       labels = c('Lower tercile', 'Middle tercile', 'Upper tercile')) + 
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        plot.margin = unit(c(0, 0, 0, 0.01), "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))


img <- grid.arrange(repair.rates, damage.rates, f.plot, nrow=3)
ggsave("../plots/human_wealth_rates_f.pdf", img, width=8, height=6)
