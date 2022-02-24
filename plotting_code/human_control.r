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
elsa$sex <- as.factor(elsa$sex)
elsa <- elsa %>% group_by(id) %>% filter( (baseline.age >= 50)  & (baseline.age < 90))

elsa$f <- elsa$n/elsa$N
elsa$wealth <- log(elsa$wealth + mean(elsa$wealth))

bins <- seq(50, 90, by=10)

elsa <- 
elsa %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE)) %>% as.data.frame()
elsa$baseline.bin <- as.factor(elsa$baseline.bin)
elsa$sex <- as.factor(elsa$sex)
levels(elsa$sex) <- c('Male', 'Female')

elsa.test <- na.omit(elsa)
elsa.na <- na.omit(elsa.test)


elsa.na$prepair <- (elsa.na$repair/elsa.na$n)/elsa.na$delta.t
elsa.na$pdamage <- (elsa.na$damage/(elsa.na$N - elsa.na$n))/elsa.na$delta.t


bins <- seq(20, 105, by=2)

binned.data <- 
    elsa.na %>% 
    mutate(age = cut(age, bins, include.lowest=TRUE)) %>%
    group_by(baseline.bin, age, sex) %>%
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
    spread_draws(lambda_r[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']) %>%
    mutate(time = elsa.test[n, 'time']) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']) %>%
    mutate(id = elsa.test[n, 'id']) %>%
    mutate(num = elsa.test[n,'n']) %>%
    filter(age < 97) %>%
    filter(num > 0) %>%
    ungroup() %>%
    group_by(baseline.bin, sex, age, .draw) %>%
    summarize(lambda_r=mean(lambda_r)) %>%
    ungroup()


damage <- fit %>%
    spread_draws(lambda_d[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']) %>%
    mutate(time = elsa.test[n, 'time']) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']) %>%
    mutate(id = elsa.test[n, 'id']) %>%
    filter(age < 97) %>%
    ungroup() %>%
    group_by(baseline.bin, sex, age, .draw) %>%
    summarize(lambda_d=mean(lambda_d)) %>%
    ungroup()


corr.repair <- fit %>%
    spread_draws(lambda_r[n]) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']) %>%
    mutate(time = elsa.test[n, 'time']) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']) %>%
    mutate(id = elsa.test[n, 'id']) %>%
    ungroup() %>%
    group_by(id, baseline.bin, .draw) %>%
    summarize(corr_r = cor(time, lambda_r, method="spearman")) %>%
    ungroup() %>%
    group_by(baseline.bin, .draw) %>%
    summarize(corr_r = mean(corr_r, na.rm=TRUE))

corr.damage <- fit %>%
    spread_draws(lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']) %>%
    mutate(time = elsa.test[n, 'time']) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']) %>%
    mutate(id = elsa.test[n, 'id']) %>%
    ungroup() %>%
    group_by(id, baseline.bin, .draw) %>%
    summarize(corr_d = cor(time, lambda_d, method="spearman")) %>%
    ungroup() %>%
    group_by(baseline.bin, .draw) %>%
    summarize(corr_d = mean(corr_d, na.rm=TRUE))

corr.repair.text <- corr.repair  %>%
    ungroup() %>%
    group_by(baseline.bin) %>%
    mean_hdci(corr_r, .width=0.95) %>%
    ungroup() %>%
    group_by(baseline.bin) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_r*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$')) %>%
    mutate(x = midpoints(baseline.bin))


corr.damage.text <- corr.damage  %>%
    ungroup() %>%
    group_by(baseline.bin) %>%
    mean_hdci(corr_d, .width=0.95) %>%
    ungroup() %>%
    group_by(baseline.bin) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_d*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$')) %>%
    mutate(x = midpoints(baseline.bin))


repair$color <- '0'
repair$fill <- '0'
damage$color <- '0'
damage$fill <- '0'
binned.data$color <- '0'
binned.data$fill <- '0'

levels(binned.data$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(repair$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(damage$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(corr.repair.text$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(corr.damage.text$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')


# save figure source data
full.plot.data <- repair
full.plot.data$lambda_d <- damage$lambda_d
full.plot.data <- subset(full.plot.data, select = -c(color, fill))
write.csv(subset(binned.data, select = -c(color, fill)), '../figure_data/figure2/human_binned_rates.csv')
write.csv(full.plot.data, '../figure_data/figure2/human_model_rates.csv')
write.csv(corr.repair.text, '../figure_data/figure2/human_rates_correlation_repair.csv')
write.csv(corr.damage.text, '../figure_data/figure2/human_rates_correlation_damage.csv')



repair.rates <- ggplot() +
    facet_grid(.~baseline.bin,scales='free') + theme_cowplot() +
    geom_point(data = binned.data, mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
               ymax=upper.repair, color=sex, shape=sex), size = 2) +
    geom_errorbar(data = binned.data, mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                  ymax=upper.repair, color=sex), width=0.25) +
    geom_line(data=repair, mapping=aes(x=age, y=lambda_r, color=sex, group=interaction(sex,.draw)), alpha=0.075) +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Repair rate (/year)", shape="Sex") +
    scale_color_manual(name="Sex", values=sex.palette) +
    scale_fill_manual(name="Sex", values=sex.palette) +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.title = element_text(size=8)) +
    ggtitle('d) ELSA humans (Phelps et al. 2020)') +
    geom_text(
        data    = corr.repair.text,
        mapping = aes(x = x, y = 0.4, label = TeX(label,output='character')),size=2.5,
        hjust   = 0,
        vjust   = 0, parse=TRUE)


damage.rates <- ggplot() +
    facet_grid(.~baseline.bin,scales='free') + theme_cowplot() +
    geom_point(data = binned.data, mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
               ymax=upper.damage, color=sex,shape=sex), size = 2) +
    geom_errorbar(data = binned.data, mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                  ymax=upper.damage, color=sex), width=0.25) +
    geom_line(data=damage, mapping=aes(x=age, y=lambda_d, color=sex, group=interaction(sex,.draw)), alpha=0.075) +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Damage rate (/year)", shape="Sex") +
    scale_color_manual(name="Sex", values=sex.palette) +
    scale_fill_manual(name="Sex", values=sex.palette) +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7)) +
    geom_text(
        data    = corr.damage.text,
        mapping = aes(x = x-5, y = 0.25, label = TeX(label,output='character')),size=2.5,
        hjust   = 0,
        vjust   = 0, parse=TRUE)


# compute binned frailty index
binned.data <- 
    elsa %>% 
    mutate(age = cut(age, bins, include.lowest=TRUE)) %>%
    group_by(baseline.bin, age, sex) %>%
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
    spread_draws(sampled_n[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa[n, 'sex']) %>%
    mutate(time = elsa[n, 'time']) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa[n,'baseline.bin']) %>%
    mutate(id = elsa[n, 'id']) %>%
    filter(age <= 98) %>%
    ungroup() %>%
    group_by(baseline.bin, sex, age, .draw) %>%
    summarize(mean_f=mean(sampled_n)/23.0) %>%
    ungroup() %>% as.data.frame()

levels(binned.data$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')
levels(f$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')

# save figure source data
write.csv(binned.data, '../figure_data/supplemental_figure2/human_binned_f.csv')
write.csv(f, '../figure_data/supplemental_figure2/human_model_f.csv')


f.plot <- ggplot() +
    facet_grid(.~baseline.bin,scales='free') + theme_cowplot() +
    geom_point(data = binned.data, mapping=aes(x=age, y=mean.f, ymin=lower.f, 
               ymax=upper.f, color=sex), size = 2) +
    geom_errorbar(data = binned.data, mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                  ymax=upper.f, color=sex), width=0.25) +
    geom_line(data=f, mapping=aes(x=age, y=mean_f, color=sex, group=interaction(sex,.draw)), alpha=0.05) +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Frailty Index", shape="Sex") +
    scale_color_manual(name="Sex", values=sex.palette) +
    scale_fill_manual(name="Sex", values=sex.palette) +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))# + ggtitle('d) ELSA humans (Phelps et al. 2020)')

img <- plot_grid(repair.rates, NULL, damage.rates,nrow=3, rel_heights=c(1.6/3, -0.2/3, 1.4/3))
ggsave("../plots/human_control_rates_sex.pdf", img, width=8, height=3.5)

ggsave("../plots/human_control_f_sex.pdf", f.plot, width=8, height=1.5)
