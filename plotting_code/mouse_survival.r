library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(survminer)
library(tidybayes)
library(latex2exp)
library(tidyr)
library(ggfortify)
library(cowplot)
library(posterior)
source("../utils/functions.r")
source("../utils/palettes.r")

enalapril <- read.csv('../datasets/enalapril_data.csv', header = TRUE, sep = ",")
enalapril.surv <- read.csv('../datasets/enalapril_surv_data.csv', header = TRUE, sep = ",")

enalapril$sex <- as.factor(enalapril$sex)
enalapril$treatment <- as.factor(enalapril$treatment)
enalapril$age <- enalapril$time + enalapril$baseline.age
enalapril$prepair <- (enalapril$repair/enalapril$n)/enalapril$delta.t
enalapril$pdamage <- (enalapril$damage/(enalapril$N - enalapril$n))/enalapril$delta.t
enalapril$sex.numeric <- as.numeric(enalapril$sex)-1
enalapril$treatment.numeric <- as.numeric(enalapril$treatment)-1

enalapril.test <- na.omit(enalapril)

sex.mean <- mean(enalapril.test$sex.numeric)
treatment.mean <- mean(enalapril.test$treatment.numeric)
sex.sd <- sd(enalapril.test$sex.numeric)
treatment.sd <- sd(enalapril.test$treatment.numeric)

enalapril.test$event.time <- enalapril.test$death.age - enalapril.test$baseline.age
enalapril.surv$event.time <- enalapril.surv$death.age - enalapril.surv$baseline.age

enalapril.surv2 <- tmerge(enalapril.surv, enalapril.surv, id=mouse, endpt = event(event.time, status))

enalapril.surv.long <- tmerge(enalapril.surv2, enalapril.test, id=mouse, f = tdc(time, f), n = tdc(time, n), age = tdc(time, age))
enalapril.surv.long$tstart <- enalapril.surv.long$tstart
enalapril.surv.long$tstop.orig <- enalapril.surv.long$tstop
enalapril.surv.long$tstop <- enalapril.surv.long$tstop + enalapril.surv.long$baseline.age
enalapril.surv.long$event.time <- enalapril.surv.long$event.time
enalapril.surv.long$sex <- as.factor(enalapril.surv.long$sex)


fit <- readRDS('../fits/mouse_1.rds')



bins <- seq(14.5, 35, by=1)
binned.test <- cut(enalapril.surv.long$tstop, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

sd.enalapril <- fit %>% spread_draws(lambda_r[n], lambda_d[n]) %>%
    group_by(.draw) %>%
    summarize(sd.r = sd(invsoftplus(lambda_r)),  sd.d = sd(invsoftplus(lambda_d))) %>%
    ungroup() %>%
    summarize(sd.r = median(sd.r), sd.d = median(sd.d)) %>%
    as.data.frame()
sd.enalapril <- as.numeric(sd.enalapril[1,])


sd.enalapril.nologit <- fit %>% spread_draws(lambda_r[n], lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(mouse = enalapril.surv.long[n,'mouse']) %>%
    mutate(f = enalapril.surv.long[n,'f']) %>%
    group_by(mouse) %>%
    summarize(sd.r = sd((lambda_r)),  sd.d = sd((lambda_d)), sd.f = sd(f)) %>%
    ungroup() %>%
    summarize(sd.r = median(sd.r), sd.d = median(sd.d), sd.f = median(sd.f)) %>%
    as.data.frame()
sd.enalapril.nologit <- as.numeric(sd.enalapril.nologit[1,])


mean.enalapril <- fit %>% spread_draws(lambda_r[n], lambda_d[n]) %>%
    group_by(.draw) %>%
    mutate(mouse = enalapril.surv.long[n,'mouse']) %>%
    mutate(f = enalapril.surv.long[n,'f']) %>%
    group_by(mouse) %>%
    summarize(mean.r = mean((lambda_r)),  mean.d = mean((lambda_d)), mean.f = mean(f)) %>%
    ungroup() %>%
    summarize(mean.r = median(mean.r), mean.d = median(mean.d), mean.f = median(mean.f)) %>%
    as.data.frame()
mean.enalapril <- as.numeric(mean.enalapril[1,])



enalapril.group <- fit %>% spread_draws(lambda_r[n], lambda_d[n]) %>%
    mutate(sex = enalapril.surv.long[n,'sex']) %>%
    mutate(mouse = enalapril.surv.long[n,'mouse']) %>%
    mutate(f = enalapril.surv.long[n,'f']) %>%
    group_by(sex, mouse, .draw) %>%
    summarize(mean.repair = mean(lambda_r),
              mean.damage = mean(lambda_d),
              mean.f = mean(f))%>%
    ungroup() %>%
    group_by(sex, mouse) %>%
    summarize(mean.repair = median(mean.repair),
              mean.damage = median(mean.damage),
              mean.f = median(mean.f)) %>%
    ungroup() %>%
    group_by(sex, mouse) %>%
    mutate(group.repair = cut(mean.repair, c(-Inf, mean.enalapril[1] - sd.enalapril.nologit[1], mean.enalapril[1] + sd.enalapril.nologit[1], Inf))) %>%
    mutate(group.damage = cut(mean.damage, c(-Inf, mean.enalapril[2] - sd.enalapril.nologit[2], mean.enalapril[2] + sd.enalapril.nologit[2], Inf))) %>%
    mutate(group.f = cut(mean.damage, c(-Inf, mean.enalapril[3] - sd.enalapril.nologit[3], mean.enalapril[3] + sd.enalapril.nologit[3], Inf)))

enalapril.surv$group.repair <- enalapril.group$group.repair
enalapril.surv$group.damage <- enalapril.group$group.damage
enalapril.surv$group.f <- enalapril.group$group.damage


enalapril.female.control <- fit %>%
    gather_draws(female_control_repair, female_control_damage) %>%
    mutate(sex = 'Female') %>%
    mutate(treatment = 'Control') %>%
    as.data.frame()

enalapril.female.drug <- fit %>%
    gather_draws(female_drug_repair, female_drug_damage) %>%
    mutate(sex = 'Female') %>%
    mutate(treatment = 'Enalapril') %>%
    as.data.frame()

enalapril.male.control <- fit %>%
    gather_draws(male_control_repair,  male_control_damage) %>%
    mutate(sex = 'Male') %>%
    mutate(treatment = 'Control') %>%
    as.data.frame()

enalapril.male.drug <- fit %>%
    gather_draws(male_drug_repair, male_drug_damage) %>%
    mutate(sex = 'Male') %>%
    mutate(treatment = 'Enalapril') %>%
    as.data.frame()


enalapril.hazard <- rbind(enalapril.male.control, enalapril.male.drug, enalapril.female.control, enalapril.female.drug)

enalapril.hazard[enalapril.hazard$.variable=='female_control_repair',]$.variable = 'Repair rate'
enalapril.hazard[enalapril.hazard$.variable=='female_drug_repair',]$.variable = 'Repair rate'
enalapril.hazard[enalapril.hazard$.variable=='male_control_repair',]$.variable = 'Repair rate'
enalapril.hazard[enalapril.hazard$.variable=='male_drug_repair',]$.variable = 'Repair rate'

enalapril.hazard[enalapril.hazard$.variable=='female_control_damage',]$.variable = 'Damage rate'
enalapril.hazard[enalapril.hazard$.variable=='female_drug_damage',]$.variable = 'Damage rate'
enalapril.hazard[enalapril.hazard$.variable=='male_control_damage',]$.variable = 'Damage rate'
enalapril.hazard[enalapril.hazard$.variable=='male_drug_damage',]$.variable = 'Damage rate'

enalapril.hazard[enalapril.hazard$.variable=='Repair rate',]$.value = enalapril.hazard[enalapril.hazard$.variable=='Repair rate',]$.value * sd.enalapril[1]
enalapril.hazard[enalapril.hazard$.variable=='Damage rate',]$.value = enalapril.hazard[enalapril.hazard$.variable=='Damage rate',]$.value * sd.enalapril[2]


enalapril.hazard <- enalapril.hazard[enalapril.hazard$treatment=='Control', ]


# save figure source data
write.csv(enalapril.hazard %>% group_by(.variable, sex) %>% median_hdci(.value), '../figure_data/figure2/mouse1_hazards.csv')


enalapril.hazard <- enalapril.hazard %>%
    mutate(color = '0') %>%
    ggplot() +
    stat_eye(aes(x=.variable,y=.value,fill=sex),alpha=0.75, .width=c(0.95), position = position_dodge(width = 0.75), point_interval=median_hdci) +
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
        plot.margin = unit(c(0.01, 0.0, 0, 0.02), "cm")) + ggtitle('e) Mouse dataset 1\n    (Keller et al. 2019)') + labs(y='Log hazard ratio', x='')+ scale_fill_manual(values=sex.palette) + scale_y_continuous(breaks=c(-2,-1,0,1,2),label = c("-2", "-1", "0", "1", "2"), limits=c(-2.75,2.75))


###### mouse dataset 2
exercise <- read.csv('../datasets/exercise_data.csv', header = TRUE, sep = ",")
exercise.surv <- read.csv('../datasets/exercise_surv_data.csv', header = TRUE, sep = ",")

           
exercise$sex <- as.factor(exercise$sex)
exercise$exercise <- as.factor(exercise$exercise)
exercise$age <- exercise$time + exercise$baseline.age
exercise$prepair <- (exercise$repair/exercise$n)/exercise$delta.t
exercise$pdamage <- (exercise$damage/(exercise$N - exercise$n))/exercise$delta.t
exercise$sex.numeric <- as.numeric(exercise$sex)-1
exercise$exercise.numeric <- as.numeric(exercise$exercise)-1

exercise.test <- na.omit(exercise)

sex.mean <- mean(exercise.test$sex.numeric)
exercise.mean <- mean(exercise.test$exercise.numeric)
sex.sd <- sd(exercise.test$sex.numeric)
exercise.sd <- sd(exercise.test$exercise.numeric)


exercise.test$event.time <- exercise.test$death.age
exercise.surv$event.time <- exercise.surv$death.age
exercise.surv <- exercise.surv[exercise.surv$mouse %in% unique(exercise.test$mouse),]

exercise.surv2 <- tmerge(exercise.surv, exercise.surv, id=mouse, endpt = event(event.time, status))

exercise.surv.long <- tmerge(exercise.surv2, exercise.test, id=mouse, f = tdc(time, f), n = tdc(time, n), age = tdc(time, age))
exercise.surv.long$tstart <- exercise.surv.long$tstart
exercise.surv.long$tstop.orig <- exercise.surv.long$tstop
exercise.surv.long$tstop <- exercise.surv.long$tstop + exercise.surv.long$baseline.age
exercise.surv.long$event.time <- exercise.surv.long$event.time
exercise.surv.long$sex <- as.factor(exercise.surv.long$sex)

fit <- readRDS('../fits/mouse_2.rds')

bins <- seq(20.5, 35, by=0.75)
binned.test <- cut(exercise.surv.long$tstop, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

sd.exercise <- fit %>% spread_draws(lambda_r[n], lambda_d[n]) %>%
    group_by(.draw) %>%
    summarize(sd.r = sd(invsoftplus(lambda_r)), sd.d = sd(invsoftplus(lambda_d))) %>%
    ungroup() %>%
    summarize(sd.r = median(sd.r), sd.d = median(sd.d)) %>%
    as.data.frame()
sd.exercise <- as.numeric(sd.exercise[1,])



exercise.female.control <- fit %>%
    gather_draws(female_control_repair, female_control_damage) %>%
    mutate(sex = 'Female') %>%
    mutate(exercise = 'Control') %>%
    as.data.frame()

exercise.female.drug <- fit %>%
    gather_draws(female_drug_repair, female_drug_damage) %>%
    mutate(sex = 'Female') %>%
    mutate(exercise = 'Exercise') %>%
    as.data.frame()

exercise.male.control <- fit %>%
    gather_draws(male_control_repair,  male_control_damage) %>%
    mutate(sex = 'Male') %>%
    mutate(exercise = 'Control') %>%
    as.data.frame()

exercise.male.drug <- fit %>%
    gather_draws(male_drug_repair, male_drug_damage) %>%
    mutate(sex = 'Male') %>%
    mutate(exercise = 'Exercise') %>%
    as.data.frame()

exercise.hazard <- rbind(exercise.male.control, exercise.male.drug, exercise.female.control, exercise.female.drug)

exercise.hazard[exercise.hazard$.variable=='female_control_repair',]$.variable = 'Repair rate'
exercise.hazard[exercise.hazard$.variable=='female_drug_repair',]$.variable = 'Repair rate'
exercise.hazard[exercise.hazard$.variable=='male_control_repair',]$.variable = 'Repair rate'
exercise.hazard[exercise.hazard$.variable=='male_drug_repair',]$.variable = 'Repair rate'
exercise.hazard[exercise.hazard$.variable=='female_control_damage',]$.variable = 'Damage rate'
exercise.hazard[exercise.hazard$.variable=='female_drug_damage',]$.variable = 'Damage rate'
exercise.hazard[exercise.hazard$.variable=='male_control_damage',]$.variable = 'Damage rate'
exercise.hazard[exercise.hazard$.variable=='male_drug_damage',]$.variable = 'Damage rate'


exercise.hazard[exercise.hazard$.variable=='Repair rate',]$.value = exercise.hazard[exercise.hazard$.variable=='Repair rate',]$.value * sd.exercise[1]
exercise.hazard[exercise.hazard$.variable=='Damage rate',]$.value = exercise.hazard[exercise.hazard$.variable=='Damage rate',]$.value * sd.exercise[2]


exercise.hazard <- exercise.hazard[exercise.hazard$exercise=='Control', ]


# save figure source data
write.csv(exercise.hazard %>% group_by(.variable, sex) %>% median_hdci(.value), '../figure_data/figure2/mouse2_hazards.csv')


exercise.hazard <- exercise.hazard %>%
    ggplot() +
    stat_eye(aes(x=.variable,y=.value,fill=sex),alpha=0.75, .width=c(0.95), position = position_dodge(width = 0.75), point_interval=median_hdci)+
    geom_hline(yintercept = 0, linetype="dotted")  + 
    scale_fill_manual(values = sex.palette) +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0.0, 0.0, 0, 0), "cm")) + ggtitle('f) Mouse dataset 2\n    (Bisset et al. 2021)')  + labs(y='', x='') +  scale_y_continuous(breaks=c(-2,-1,0,1,2),label = c("-2", "-1", "0", "1", "2"),limits=c(-2.75,2.75))

##### mouse dataset 3
schultz <- read.csv('../datasets/schultz_data.csv', header = TRUE, sep = ",")
schultz.surv <- read.csv('../datasets/schultz_surv_data.csv', header = TRUE, sep = ",")

schultz$prepair <- (schultz$repair/schultz$n)/schultz$delta.t
schultz$pdamage <- (schultz$damage/(schultz$N - schultz$n))/schultz$delta.t
schultz.test <- na.omit(schultz)

f.mean <- mean(schultz.test$f)
f.sd <- sd(schultz.test$f)
prepair.sd <- sd(invsoftplus(schultz.test[schultz.test$prepair>0,]$prepair))
pdamage.sd <- sd(invsoftplus(schultz.test[schultz.test$pdamage>0,]$pdamage))

schultz$age <- schultz$time + schultz$baseline.age

schultz$prepair <- (schultz$repair/schultz$n)/schultz$delta.t
schultz$pdamage <- (schultz$damage/(schultz$N - schultz$n))/schultz$delta.t

schultz.test <- na.omit(schultz)

schultz.test$event.time <- schultz.test$death.age - schultz.test$baseline.age
schultz.surv$event.time <- schultz.surv$death.age - schultz.surv$baseline.age
schultz.surv <- schultz.surv[schultz.surv$mouse %in% unique(schultz.test$mouse),]

schultz.surv2 <- tmerge(schultz.surv, schultz.surv, id=mouse, endpt = event(event.time, status))

schultz.surv.long <- tmerge(schultz.surv2, schultz.test, id=mouse, f = tdc(time, f), n = tdc(time, n), age = tdc(time, age))
schultz.surv.long$tstart <- schultz.surv.long$tstart
schultz.surv.long$tstop.orig <- schultz.surv.long$tstop
schultz.surv.long$tstop <- schultz.surv.long$tstop + schultz.surv.long$baseline.age
schultz.surv.long$event.time <- schultz.surv.long$event.time
schultz.surv.long$sex <- as.factor(schultz.surv.long$sex)


fit <- readRDS('../fits/mouse_3.rds')

bins <- seq(20, 45, by=1.5)
binned.test <- cut(schultz.surv.long$tstop, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)


sd.schultz <- fit %>% spread_draws(lambda_r[n], lambda_d[n]) %>%
    group_by(.draw) %>%
    summarize(sd.r = sd(invsoftplus(lambda_r)), sd.d = sd(invsoftplus(lambda_d))) %>%
    ungroup() %>%
    summarize(sd.r = median(sd.r), sd.d = median(sd.d)) %>%
    as.data.frame()
sd.schultz <- as.numeric(sd.schultz[1,])


schultz.hazard <- fit %>%
    spread_draws(gamma_rand[m]) %>%
    mutate(sex = 'Male') %>%
    mutate(.value = gamma_rand * sd.schultz[m]) %>%
    mutate(.variable = c('Repair rate', 'Damage rate')[m])


# save figure source data
write.csv(schultz.hazard %>% group_by(.variable, sex) %>% median_hdci(.value), '../figure_data/figure2/mouse3_hazards.csv')

schultz.hazard <- schultz.hazard %>% ggplot() +
    stat_eye(aes(.variable, .value),alpha=0.75, .width=c(0.95), fill = sex.palette[2], point_interval=median_hdci) + geom_hline(yintercept = 0, linetype="dotted")   +
    theme_cowplot() + theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = "none",
        plot.margin = unit(c(0.0 , 0.05, 0, 0), "cm"))  +
    ggtitle('g) Mouse dataset 3\n    (Schultz et al. 2020)') +
    labs(y='', x='') +
    scale_fill_manual(values=sex.palette) +
    scale_y_continuous(breaks=c(-2,-1,0,1,2),label = c("-2", "-1", "0", "1", "2"), limits=c(-2.75,2.75))


hazard.img <- plot_grid(enalapril.hazard, exercise.hazard, schultz.hazard, nrow=1, ncol=3, rel_widths = c(1.95/5, 1.95/5, 1.1/5))

ggsave("../plots/mice_survival_hazards.pdf", hazard.img, width=8, height=1.5)
