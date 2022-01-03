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


bins <- seq(0, 15, by=1)
binned <- cut(enalapril.test$time, bins, include.lowest = TRUE)
binned <- midpoints(binned)


# get stan fit
fit <- readRDS('../fits/mouse_1.rds') 


fdotdot.plot <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = enalapril.test[n, 'sex']) %>%
    mutate(treatment = enalapril.test[n, 'treatment']) %>%
    mutate(f = enalapril.test[n, 'f']) %>%
    mutate(age = binned[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(dotdotf=mean(dotdotf)) %>%
    ungroup() %>%
    group_by(sex, treatment, age) %>%
    median_hdci(dotdotf, .width=0.95)


fdotdot.plot$sex <- as.factor(fdotdot.plot$sex)
fdotdot.plot$treatment <- as.factor(fdotdot.plot$treatment)

levels(fdotdot.plot$sex) <- c('Female', 'Male')

enalapril.plot <- ggplot() +
    geom_lineribbon(data=fdotdot.plot, mapping=aes(x=age, y=dotdotf, fill=treatment,
                                                             color=NA, ymin=.lower, ymax=.upper), alpha=0.5) +
    geom_line(data=fdotdot.plot, mapping=aes(x=age, y=dotdotf,
          fill=treatment, color=treatment,group=treatment), alpha=1, size=1.75, color='white') +
    geom_line(data=fdotdot.plot, mapping=aes(x=age, y=dotdotf,
         fill=treatment, color=treatment, group=treatment), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
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
    facet_grid(~sex) +
    labs(y='Frailty Index curvature terms', x='Age (months)',
         color='', fill='') +
    scale_fill_d3() +
    scale_color_d3() + xlim(0, 10.25)+
    ggtitle('a) Mouse dataset 1 (Keller et al. 2019)')
    #scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
    #                   label = c("16", "18", "20", "22", "24")) +
    #scale_color_manual(values=c('Damage rate\n terms' = '#CC79A7', 'Repair rate\n terms'='#009E73')) +
    #guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + 



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


fdotdot.plot <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = exercise.test[n, 'sex']) %>%
    mutate(exercise = exercise.test[n, 'exercise']) %>%
    mutate(f = exercise.test[n, 'f']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    ungroup() %>%
    group_by(sex, exercise, age, .draw) %>%
    summarize(dotdotf=mean(dotdotf)) %>%
    ungroup() %>%
    group_by(sex, exercise, age) %>%
    median_hdci(dotdotf, .width=0.95) %>% as.data.frame()

fdotdot.plot$exercise <- as.factor(fdotdot.plot$exercise)
fdotdot.plot$sex <- as.factor(fdotdot.plot$sex)
levels(fdotdot.plot$exercise) <- c('Control', 'Exercise')
levels(fdotdot.plot$sex) <- c('Female', 'Male')



exercise.plot <- ggplot() +
    geom_lineribbon(data=fdotdot.plot, mapping=aes(x=age, y=dotdotf, fill=exercise,
                                                             color=NA, ymin=.lower, ymax=.upper), alpha=0.5) +
    geom_line(data=fdotdot.plot, mapping=aes(x=age, y=dotdotf,
          fill=exercise, color=exercise,group=exercise), alpha=1, size=1.75, color='white') +
    geom_line(data=fdotdot.plot, mapping=aes(x=age, y=dotdotf,
         fill=exercise, color=exercise, group=exercise), alpha=1, size=1.25)+
    geom_hline(yintercept = 0, linetype="dotted") +
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
    facet_grid(~sex) +
    labs(y='', x='Age (months)',
         color='', fill='') +  
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette) + scale_x_continuous(breaks=c(0, 1, 2, 3), label = c("0", "1", "2", "3"))+ xlim(0, 3.25) +
    ggtitle('b) Mouse dataset 2 (Bisset et al. 2021)')
    #scale_color_manual(values=c('Damage rate\n terms' = '#e04d52', 'Repair rate\n terms'='#54b14e'))+
    #scale_color_manual(values=c('Damage rate\n terms' = '#CC79A7', 'Repair rate\n terms'='#009E73')) +
    #guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + 

img <- grid.arrange(enalapril.plot, exercise.plot, ncol=2,nrow=1)

#ggsave("../plots/mice_control_rates.pdf", img, width=8, height=4.5)
ggsave("../plots/mice_curvature_interventions.pdf", img, width=8, height=1.5)
