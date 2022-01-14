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

enalapril.test <- na.omit(enalapril)

bins <- seq(14.5, 29, by=1)
binned <- cut(enalapril$age, bins, include.lowest = TRUE)
binned <- midpoints(binned)
enalapril$time.bin <- binned

binned.data.rates <- 
enalapril %>% 
  group_by(time.bin, treatment, sex) %>%
  summarise(mean.damage = mean(pdamage,na.rm=TRUE), 
            mean.repair = mean(prepair,na.rm=TRUE),
            mean.f = mean(f, na.rm=TRUE),
            se.f = se(f),
            se.repair = se(prepair),
            se.damage = se(pdamage),
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

binned.data.rates <- binned.data.rates[binned.data.rates$treatment=='control',]

# get stan fit
fit <- readRDS('../fits/mouse_1.rds')

binned.test <- cut(enalapril.test$age, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

full.repair.plot <- fit %>%
    spread_draws(lambda_r[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
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


binned <- cut(enalapril$age, bins, include.lowest = TRUE)
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
    summarize(mean_f = mean(mean_f)) %>%
    ungroup() %>% as.data.frame()

# compute rank correlations
corr.repair <- fit %>%
    spread_draws(lambda_r[n]) %>%
    group_by(.draw) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(age = enalapril.test[n, 'age']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    ungroup() %>%
    group_by(sex, mouse, treatment, .draw) %>%
    summarize(corr_r = cor(age, lambda_r, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, treatment, .draw) %>%
    summarize(corr_r = mean(corr_r, na.rm=TRUE))

corr.damage <- fit %>%
    spread_draws(lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(age = enalapril.test[n, 'age']) %>%
    mutate(mouse = enalapril.test[n,'mouse']) %>%
    ungroup() %>%
    group_by(sex, mouse, treatment, .draw) %>%
    summarize(corr_d = cor(age, lambda_d, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, treatment, .draw) %>%
    summarize(corr_d = mean(corr_d, na.rm=TRUE))


full.repair.plot <- full.repair.plot[full.repair.plot$treatment=='control',]
full.damage.plot <- full.damage.plot[full.damage.plot$treatment=='control',]
full.f.plot <- full.f.plot[full.f.plot$treatment=='control',]

corr.repair <- corr.repair[corr.repair$treatment=='control',]
corr.damage <- corr.damage[corr.damage$treatment=='control',]


corr.repair.text <- corr.repair  %>%
    ungroup() %>%
    group_by(sex) %>%
    mean_hdci(corr_r, .width=0.95) %>%
    ungroup() %>%
    group_by(sex) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_r*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$')) %>%ungroup() %>%
    as.data.frame()


corr.damage.text <- corr.damage  %>%
    ungroup() %>%
    group_by(sex) %>%
    mean_hdci(corr_d, .width=0.95) %>%
    ungroup() %>%
    group_by(sex) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_d*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$')) %>%ungroup()


# fix labels
binned.data.rates$sex <- as.factor(binned.data.rates$sex)
full.repair.plot$sex <- as.factor(full.repair.plot$sex)
full.damage.plot$sex <- as.factor(full.damage.plot$sex)
full.f.plot$sex <- as.factor(full.f.plot$sex)
corr.repair.text$sex <- as.factor(corr.repair.text$sex)
corr.damage.text$sex <- as.factor(corr.damage.text$sex)
levels(binned.data.rates$sex) <- c('Female', 'Male')
levels(full.repair.plot$sex) <- c('Female', 'Male')
levels(full.damage.plot$sex) <- c('Female', 'Male')
levels(full.f.plot$sex) <- c('Female', 'Male')
levels(corr.repair.text$sex) <- c('Female', 'Male')
levels(corr.damage.text$sex) <- c('Female', 'Male')

corr.repair.text$y <- 0
corr.repair.text[corr.repair.text$sex=='Male',]$y <- 0.6
corr.repair.text[corr.repair.text$sex=='Female',]$y <- 0.48

corr.damage.text$y <- 0
corr.damage.text[corr.damage.text$sex=='Male',]$y <- 0.18
corr.damage.text[corr.damage.text$sex=='Female',]$y <- 0.145

# plot repair rates
enalapril.repair <- ggplot() +
    geom_point(binned.data.rates %>% filter(repair.num >= 0), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                              ymax=upper.repair, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(repair.num >= 0), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                             ymax=upper.repair, group=sex, color=sex), width=0.1) +
    geom_line(data=full.repair.plot, mapping=aes(x=age, y=lambda_r, color=sex, group=interaction(.draw, sex)), alpha=0.075) +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "none",
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) + 
    labs(y='Repair rate (/month)', x='',
         color='', fill='',shape='') +
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
                       label = c("16", "18", "20", "22", "24"))  + ggtitle("a) Mouse dataset 1 (Keller et al. 2019)")   +
    geom_text(
        data    = corr.repair.text,
        mapping = aes(x = 18, y = y, label = TeX(label,output='character'), color=sex, fontface = "bold"),size=2.75,
        hjust   = 0,
        vjust   = 0, parse=TRUE)

# plot damage rates
enalapril.damage <- ggplot() +
    geom_point(binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                              ymax=upper.damage, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                             ymax=upper.damage, group=sex, color=sex), width=0.1) +
    geom_line(data=full.damage.plot, mapping=aes(x=age, y=lambda_d, color=sex, group=interaction(.draw, sex)), alpha=0.075) +
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
        plot.margin = unit(c(0, 0, 0, 0.01), "cm")) + 
    labs(y='Damage rate (/month)', x='Age (months)',
         color='', fill='',shape='') +
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
                       label = c("16", "18", "20", "22", "24"))+
    scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
                       label = c("16", "18", "20", "22", "24")) +
    geom_text(
        data    = corr.damage.text,
        mapping = aes(x = 15.7, y = y, label = TeX(label,output='character'), color=sex, fontface = "bold"),size=2.75,
        hjust   = 0,
        vjust   = 0, parse=TRUE)

# plot frailty index
enalapril.f <- ggplot() +
    geom_point(binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                              ymax=upper.f, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                        ymax=upper.f, group=sex, color=sex), width=0.1) +
    geom_line(data=full.f.plot, mapping=aes(x=age, y=mean_f, color=sex, group=interaction(.draw, sex)), alpha=0.075) +
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
        plot.margin = unit(c(0, 0, 0, 0.01), "cm")) +
    labs(y='Frailty Index', x='Age (months)',
         color='', fill='',shape='') +
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_y_continuous(breaks=c(0.15, 0.2, 0.25, 0.3),
                       label = c("0.15", "0.20", "0.25", "0.30")) + 
    scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
                       label = c("16", "18", "20", "22", "24"))+ ggtitle('a) Mouse dataset 1 (Keller et al. 2019)')



##### mouse dataset 2
exercise <- read.csv('../datasets/exercise_data.csv', header = TRUE, sep = ",")

exercise$sex <- as.factor(exercise$sex)
exercise$exercise <- as.factor(exercise$exercise)
exercise$age <- exercise$time + exercise$baseline.age

exercise$prepair <- (exercise$repair/exercise$n)/exercise$delta.t
exercise$pdamage <- (exercise$damage/(exercise$N - exercise$n))/exercise$delta.t

exercise.test <- na.omit(exercise)

bins <- seq(20.5, 27, by=0.5)
binned <- cut(exercise$age, bins, include.lowest = TRUE)
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


levels(binned.data.rates$exercise) <- c('No', 'Yes')
binned.data.rates <- binned.data.rates[binned.data.rates$exercise=='No',]

# get stan fit
fit <- readRDS('../fits/mouse_2.rds')

bins <- seq(20.5, 27, by=0.5)
binned.test <- cut(exercise.test$age, bins, include.lowest = TRUE)
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
    ungroup() %>%
    filter(age <= 25.75) %>%
    filter(age >= 21.25) %>% as.data.frame()

full.damage.plot <- fit %>%
    spread_draws(lambda_d[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(lambda_d = mean(lambda_d)) %>%
    ungroup() %>%
    filter(age <= 25.75) %>%
    filter(age >= 21.25)


bins <- seq(20.5, 27, by=0.5)
binned <- cut(exercise$age, bins, include.lowest = TRUE)
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
    ungroup() %>%
    filter(age >= 21.25)

# compute rank correlations
corr.repair <- fit %>%
    spread_draws(lambda_r[n]) %>%
    group_by(.draw) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(age = exercise.test[n, 'age']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    ungroup() %>%
    group_by(sex, mouse, exercise, .draw) %>%
    summarize(corr_r = cor(age, lambda_r, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, exercise, .draw) %>%
    summarize(corr_r = mean(corr_r))

corr.damage <- fit %>%
    spread_draws(lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(age = exercise.test[n, 'age']) %>%
    mutate(mouse = exercise.test[n,'mouse']) %>%
    ungroup() %>%
    group_by(sex, mouse, exercise, .draw) %>%
    summarize(corr_d = cor(age, lambda_d, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, exercise, .draw) %>%
    summarize(corr_d = mean(corr_d))

# fix labels
levels(full.repair.plot$exercise) <- c('No', 'Yes')
levels(full.damage.plot$exercise) <- c('No', 'Yes')
levels(full.f.plot$exercise) <- c('No', 'Yes')
levels(corr.repair$exercise) <- c('No', 'Yes')
levels(corr.damage$exercise) <- c('No', 'Yes')

full.repair.plot <- full.repair.plot[full.repair.plot$exercise=='No',]
full.damage.plot <- full.damage.plot[full.damage.plot$exercise=='No',]
full.f.plot <- full.f.plot[full.f.plot$exercise=='No',]

corr.repair <- corr.repair[corr.repair$exercise=='No',]
corr.damage <- corr.damage[corr.damage$exercise=='No',]



corr.repair.text <- corr.repair  %>%
    ungroup() %>%
    group_by(sex) %>%
    mean_hdci(corr_r, .width=0.95) %>%
    ungroup() %>%
    group_by(sex) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_r*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$'))


corr.damage.text <- corr.damage  %>%
    ungroup() %>%
    group_by(sex) %>%
    mean_hdci(corr_d, .width=0.95) %>%
    ungroup() %>%
    group_by(sex) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_d*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$'))


# fix labels
binned.data.rates$sex <- as.factor(binned.data.rates$sex)
full.repair.plot$sex <- as.factor(full.repair.plot$sex)
full.damage.plot$sex <- as.factor(full.damage.plot$sex)
full.f.plot$sex <- as.factor(full.f.plot$sex)
corr.repair.text$sex <- as.factor(corr.repair.text$sex)
corr.damage.text$sex <- as.factor(corr.damage.text$sex)
levels(binned.data.rates$sex) <- c('Female', 'Male')
levels(full.repair.plot$sex) <- c('Female', 'Male')
levels(full.damage.plot$sex) <- c('Female', 'Male')
levels(full.f.plot$sex) <- c('Female', 'Male')
levels(corr.repair.text$sex) <- c('Female', 'Male')
levels(corr.damage.text$sex) <- c('Female', 'Male')

corr.repair.text$y <- 0
corr.repair.text[corr.repair.text$sex=='Male',]$y <- 1.12
corr.repair.text[corr.repair.text$sex=='Female',]$y <- 0.91

corr.damage.text$y <- 0
corr.damage.text[corr.damage.text$sex=='Male',]$y <- 0.455
corr.damage.text[corr.damage.text$sex=='Female',]$y <- 0.405

# plot repair
exercise.repair <- ggplot() +
    geom_point(data=binned.data.rates %>% filter(repair.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                              ymax=upper.repair, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(repair.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                                  ymax=upper.repair, group=sex, color=sex), width=0.1) +
    geom_line(data=full.repair.plot, mapping=aes(x=age, y=lambda_r, color=sex, group=interaction(.draw, sex)), alpha=0.075)+
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0), "cm")) +
    labs(y='', x='',
         color='', fill='',shape='') + 
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) + ggtitle('b) Mouse dataset 2 (Bisset et al. 2021)') +
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)+
    geom_text(
        data    = corr.repair.text,
        mapping = aes(x = 22, y = y, label = TeX(label,output='character'), color=sex, fontface = "bold"),size=2.75,
        hjust   = 0,
        vjust   = 0,parse=TRUE)

# plot damage
exercise.damage <- ggplot() +
    geom_point(binned.data.rates %>% filter(damage.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                              ymax=upper.damage, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(damage.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                             ymax=upper.damage, group=sex, color=sex), width=0.1) +
    geom_line(data=full.damage.plot, mapping=aes(x=age, y=lambda_d, color=sex, group=interaction(.draw, sex)), alpha=0.075)+
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
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
    labs(y='', x='Age (months)',
         color='', fill='',shape='') + 
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)+
    geom_text(
       data    = corr.damage.text,
        mapping = aes(x = 21.5, y = y, label = TeX(label,output='character'), color=sex, fontface = "bold"),size=2.75,
        hjust   = 0,
       vjust   = 0, parse=TRUE)

# plot frailty index
exercise.f <- ggplot() +
    geom_point(binned.data.rates %>% filter(f.num >= 0) %>% filter((age >= 21)), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                              ymax=upper.f, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(f.num >= 0) %>% filter((age >= 21)), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                             ymax=upper.f, group=sex, color=sex), width=0.1) +
    geom_line(data=full.f.plot, mapping=aes(x=age, y=mean_f, color=sex, group=interaction(.draw, sex)), alpha=0.075)+
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
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
    labs(y='', x='Age (months)',
         color='', fill='',shape='') +
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)+ ggtitle('b) Mouse dataset 2 (Bisset et al. 2021)')



##### mouse dataset 3
schultz <- read.csv('../datasets/schultz_data.csv', header = TRUE, sep = ",")

schultz$age <- schultz$time + schultz$baseline.age
schultz$prepair <- (schultz$repair/schultz$n)/schultz$delta.t
schultz$pdamage <- (schultz$damage/(schultz$N - schultz$n))/schultz$delta.t
schultz.test <- na.omit(schultz)


bins <- seq(20, 38, by=1.5)
binned <- cut(schultz$age, bins, include.lowest = TRUE)
binned <- midpoints(binned)
schultz$time.bin <- binned

binned.data.rates <- 
schultz %>% 
  group_by(sex, time.bin) %>%
  summarise(mean.damage = mean(pdamage,na.rm=TRUE), 
            mean.repair = mean(prepair,na.rm=TRUE),
            mean.f = mean(f, na.rm=TRUE),
            se.repair = se(prepair),
            se.damage = se(pdamage),
            se.f = se(f),
            var.repair = var(prepair,na.rm=TRUE),
            var.damage = var(pdamage,na.rm=TRUE),
            var.f = var(f,na.rm=TRUE),
	    repair.num = length(prepair),
	    damage.num = length(pdamage),
        f.num = length(f)) %>%
mutate(age = time.bin) %>%
filter(age <=36) %>%
as.data.frame()

binned.data.rates$lower.repair <- binned.data.rates$mean.repair-binned.data.rates$se.repair
binned.data.rates$upper.repair <- binned.data.rates$se.repair+binned.data.rates$mean.repair

binned.data.rates$lower.damage <- binned.data.rates$mean.damage-binned.data.rates$se.damage
binned.data.rates$upper.damage <- binned.data.rates$se.damage+binned.data.rates$mean.damage

binned.data.rates$lower.f <- binned.data.rates$mean.f-binned.data.rates$se.f
binned.data.rates$upper.f <- binned.data.rates$se.f+binned.data.rates$mean.f

# get stan fit
fit <- readRDS('../fits/mouse_3.rds')


bins <- seq(20, 38, by=1.5)
binned.test <- cut(schultz.test$age, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

full.repair.plot <- fit %>%
    spread_draws(lambda_r[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = schultz.test[n, 'sex']) %>%
    mutate(mouse = schultz.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, age, .draw) %>%
    summarize(lambda_r = mean(lambda_r)) %>%
    ungroup() %>%
    filter(age <=36)

full.damage.plot <- fit %>%
    spread_draws(lambda_d[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = schultz.test[n, 'sex']) %>%
    mutate(mouse = schultz.test[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    ungroup() %>%
    group_by(sex, age, .draw) %>%
    summarize(lambda_d = mean(lambda_d)) %>%
    ungroup() %>%
    filter(age <=36)


bins <- seq(20, 38, by=1.5)
binned <- cut(schultz$age, bins, include.lowest = TRUE)
binned <- midpoints(binned)

full.f.plot <- fit %>%
    spread_draws(sampled_n[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = schultz[n, 'sex']) %>%
    mutate(mouse = schultz[n,'mouse']) %>%
    mutate(age = binned[n]) %>%
    mutate(mean_f = sampled_n[n]/116.0) %>%
    ungroup() %>%
    group_by(sex, age, .draw) %>%
    summarize(mean_f = mean(mean_f, na.rm=TRUE)) %>%
    ungroup() %>%
    filter(age <=36)

# compute rank correlations
corr.repair <- fit %>%
    spread_draws(lambda_r[n]) %>%
    group_by(.draw) %>%
    mutate(sex = schultz.test[n, 'sex']) %>%
    mutate(age = schultz.test[n, 'age']) %>%
    mutate(mouse = schultz.test[n,'mouse']) %>%
    ungroup() %>%
    group_by(mouse, sex,  .draw) %>%
    summarize(corr_r = cor(age, lambda_r, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, .draw) %>%
    summarize(corr_r = mean(corr_r, na.rm=TRUE))

corr.damage <- fit %>%
    spread_draws(lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(sex = schultz.test[n, 'sex']) %>%
    mutate(age = schultz.test[n, 'age']) %>%
    mutate(mouse = schultz.test[n,'mouse']) %>%
    ungroup() %>%
    group_by(mouse, sex, .draw) %>%
    summarize(corr_d = cor(age, lambda_d, method="spearman")) %>%
    ungroup() %>%
    group_by(sex, .draw) %>%
    summarize(corr_d = mean(corr_d, na.rm=TRUE))

corr.repair.text <- corr.repair  %>%
    ungroup() %>%
    group_by(sex) %>%
    mean_hdci(corr_r, .width=0.95) %>%
    ungroup() %>%
    group_by(sex) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_r*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$'))

corr.damage.text <- corr.damage  %>%
    ungroup() %>%
    group_by(sex) %>%
    mean_hdci(corr_d, .width=0.95) %>%
    ungroup() %>%
    group_by(sex) %>%
    mutate(label = paste0('$\\rho = ', trunc(corr_d*100)/100, '  (',trunc(100*.lower)/100,', ',trunc(100*.upper)/100, ')$'))


df.legend <- binned.data.rates
df.legend$mean.repair <- NaN
df.legend$mean.damage <- NaN
df.legend$mean.f <- NaN
df.legend$sex <- 'Male'
df.legend[1,'sex'] <- 'Female'
levels(df.legend) <- c('Male', 'Female')

# plot repair
schultz.rep <- ggplot() + geom_point(data=df.legend, mapping=aes(x=age, y=mean.repair, color=sex, shape=sex), size=2) +
    geom_point(data=binned.data.rates %>% filter(repair.num >= 0) , mapping=aes(x=age, y=mean.repair, color=sex,shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(repair.num >= 0) , mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                                   ymax=upper.repair, color=sex,shape=sex), width=0.1) +
        geom_line(data=full.repair.plot, mapping=aes(x=age, y=lambda_r, group=.draw, color=sex), alpha=0.075) + 
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0.06, 0, 0), "cm")) +
    labs(y='', x='',
         color='', fill='',shape='') + 
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36), label = c("22", "24", "26", "28", "30", "32", "34", "36")) + ggtitle('c) Mouse dataset 3 (Schultz et al. 2020)')+
    geom_text(
       data    = corr.repair.text,
        mapping = aes(x = 26, y = 0.02, label = TeX(label,output='character'), color=sex, fontface = "bold"),size=2.75,
        hjust   = 0,
        vjust   = 0, parse=TRUE)

# plot damage
schultz.dam <- ggplot()  + geom_point(data=df.legend, mapping=aes(x=age, y=mean.repair, color=sex, shape=sex), size=2) +
    geom_point(data=binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, color=sex,shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                                  ymax=upper.damage, color=sex,shape=sex), width=0.1) +
        geom_line(data=full.damage.plot, mapping=aes(x=age, y=lambda_d, group=.draw, color=sex), alpha=0.075) + 
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0.06, 0, 0), "cm")) +
    labs(y='', x='Age (months)',
         color='', fill='',shape='') + 
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36),
                       label = c("22", "24", "26", "28", "30", "32", "34", "36"))+
    geom_text(
       data    = corr.damage.text,
        mapping = aes(x = 26, y = 0.15, label = TeX(label,output='character'), color=sex, fontface = "bold"),size=2.75,
        hjust   = 0,
        vjust   = 0, parse=TRUE)

# plot frailty index
schultz.f <- ggplot() + geom_point(data=df.legend, mapping=aes(x=age, y=mean.repair, color=sex, shape=sex), size=2) +
    geom_point(data=binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, color=sex,shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                             ymax=upper.f, color=sex,shape=sex), width=0.1) +
    geom_line(data=full.f.plot, mapping=aes(x=age, y=mean_f, group=.draw, color=sex), alpha=0.075) + 
    theme_cowplot()  + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0, 0.06, 0, 0), "cm")) +
    labs(y='', x='Age (months)',
         color='', fill='',shape='') + 
    scale_color_manual(values=sex.palette) +
    scale_fill_manual(values=sex.palette) +
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36),
                       label = c("22", "24", "26", "28", "30", "32", "34", "36"))+ ggtitle('c) Mouse dataset 3 (Schultz et al. 2020)')


img.enalapril <- grid.arrange(enalapril.repair, enalapril.damage,  ncol=1,nrow=2)
img.exercise <- grid.arrange(exercise.repair, exercise.damage,  ncol=1,nrow=2)
img.schultz <- grid.arrange(schultz.rep, schultz.dam,  ncol=1,nrow=2)

img <- plot_grid(img.enalapril, img.exercise, img.schultz, nrow = 1, ncol=3, rel_widths = c(1.52/5, 1.52/5, 1.96/5))

ggsave("../plots/mice_control_rates_sex.pdf", img, width=8, height=2.75)


img.enalapril <- grid.arrange(enalapril.repair, enalapril.damage,  ncol=1,nrow=2)
img.exercise <- grid.arrange(exercise.repair, exercise.damage,  ncol=1,nrow=2)
img.schultz <- grid.arrange(schultz.rep, schultz.dam,  ncol=1,nrow=2)

img <- plot_grid(enalapril.f, exercise.f, schultz.f, nrow = 1, ncol=3, rel_widths = c(1.52/5, 1.52/5, 1.96/5))

ggsave("../plots/mice_control_f_sex.pdf", img, width=8, height=1.25)
