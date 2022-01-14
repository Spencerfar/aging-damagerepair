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

enalapril <- na.omit(enalapril)

f.sd = sd(enalapril$f)

# get stan fit
fit <- readRDS('../fits/mouse_1.rds')

bins <- seq(14.5, 29, by=1)
binned.age <- cut(enalapril$age, bins, include.lowest = TRUE)
binned.age <- midpoints(binned.age)

# compute mouse dataset 1 terms
terms.enalapril <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = enalapril[n, 'sex']) %>%
    mutate(treatment = enalapril[n, 'treatment']) %>%
    mutate(f = enalapril[n, 'f']) %>%
    mutate(mouse = enalapril[n,'mouse']) %>%
    mutate(age = binned.age[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r = -f*deriv_r - dotf*lambda_r) %>%
    mutate(term.d = (1-f)*deriv_d - dotf*lambda_d) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(term.r=mean(term.r), term.d=mean(term.d)) %>%
    ungroup()

# compute mouse dataset 1 significance
enalapril.diff.test <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.draw) %>%
    mutate(sex = enalapril[n, 'sex']) %>%
    mutate(treatment = enalapril[n, 'treatment']) %>%
    mutate(f = enalapril[n, 'f']) %>%
    mutate(mouse = enalapril[n,'mouse']) %>%
    mutate(age = binned.age[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r =  - f*deriv_r - dotf*lambda_r) %>%
    mutate(term.d =  (1-f)*deriv_d - dotf*lambda_d) %>%
    mutate(diff = term.d - term.r) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(diff = mean(diff)) %>%
    ungroup() %>%
    group_by(sex, treatment, age) %>%
    summarize(p = mean(diff > 0)) %>% as.data.frame()
enalapril.diff.test
enalapril.text <- data.frame(sex=c('Female', 'Male'), label = c('', '$p<0.05$ for age$<22$'))

terms.enalapril$sex <- as.factor(terms.enalapril$sex)
terms.enalapril$treatment <- as.factor(terms.enalapril$treatment)
levels(terms.enalapril$sex) <- c('Female', 'Male')
levels(terms.enalapril$treatment) <- c('Control', 'Enalapril')

# select only control
terms.enalapril <- terms.enalapril[terms.enalapril$treatment=='Control',]

plot.enalapril <- ggplot() +
    geom_line(data=terms.enalapril, mapping=aes(x=age, y=term.d, group=interaction(treatment,.draw), color='Damage rate\n terms'), alpha=0.1) +
    geom_line(data=terms.enalapril, mapping=aes(x=age, y=term.r, group=interaction(treatment,.draw), color='Repair rate\n terms'), alpha=0.1) +
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
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex) +
    labs(y='Frailty Index curvature terms', x='Age (months)',
         color='', fill='') +
    scale_x_continuous(breaks=c(16, 18, 20, 22, 24),
                       label = c("16", "18", "20", "22", "24")) +
    scale_color_manual(values=c('Damage rate\n terms' = rates.palette[1], 'Repair rate\n terms'=rates.palette[2])) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + ggtitle('b) Mouse dataset 1 (Keller et al. 2019)') +
    geom_text(
        data    = enalapril.text,
        mapping = aes(x = 16, y = 0.4, label = TeX(label,output='character'), fontface = "bold"),size=2.25,color='black', 
        hjust   = 0,
        vjust   = 0, parse=TRUE)

##### mouse dataset 2
exercise <- read.csv('../datasets/exercise_data.csv', header = TRUE, sep = ",")

exercise$sex <- as.factor(exercise$sex)
exercise$exercise <- as.factor(exercise$exercise)
exercise$age <- exercise$time + exercise$baseline.age

exercise$prepair <- (exercise$repair/exercise$n)/exercise$delta.t
exercise$pdamage <- (exercise$damage/(exercise$N - exercise$n))/exercise$delta.t

exercise <- na.omit(exercise)


f.sd = sd(exercise$f)


# get stan fit
fit <- readRDS('../fits/mouse_2.rds')

bins <- seq(20.5, 27, by=0.5)
binned.test <- cut(exercise$age, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

# compute mouse dataset 2 terms
terms.exercise <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n], n = 100) %>%
    group_by(.draw) %>%
    mutate(sex = exercise[n, 'sex']) %>%
    mutate(treatment = exercise[n, 'exercise']) %>%
    mutate(f = exercise[n, 'f']) %>%
    mutate(mouse = exercise[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r =  - f*deriv_r - dotf*lambda_r) %>%
    mutate(term.d =  (1-f)*deriv_d - dotf*lambda_d) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(term.r=mean(term.r), term.d=mean(term.d)) %>%
    ungroup() %>%
    filter(age >= 21.25)

# compute mouse dataset 2 significance
exercise.diff.test <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.draw) %>%
    mutate(sex = exercise[n, 'sex']) %>%
    mutate(treatment = exercise[n, 'exercise']) %>%
    mutate(f = exercise[n, 'f']) %>%
    mutate(mouse = exercise[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r =  - f*deriv_r - dotf*lambda_r) %>%
    mutate(term.d =  (1-f)*deriv_d - dotf*lambda_d) %>%
    mutate(diff = term.d - term.r) %>%
    ungroup() %>%
    group_by(sex, treatment, age, .draw) %>%
    summarize(diff = mean(diff)) %>%
    ungroup() %>%
    group_by(sex, treatment, age) %>%
    summarize(p = mean(diff > 0)) %>% as.data.frame()

exercise.diff.test
exercise.text <- data.frame(sex=c('Female', 'Male'), label = c('$p<0.05$ for all ages', ''))

terms.exercise$treatment <- as.factor(terms.exercise$treatment)
terms.exercise$sex <- as.factor(terms.exercise$sex)
levels(terms.exercise$sex) <- c('Female', 'Male')
levels(terms.exercise$treatment) <- c('Control', 'Exercise')

# select only control
terms.exercise <- terms.exercise[terms.exercise$treatment=='Control',]

plot.exercise <- ggplot() +
    geom_line(data=terms.exercise, mapping=aes(x=age, y=term.d, group=.draw, color='Damage rate\n terms'), alpha=0.1) +
    geom_line(data=terms.exercise, mapping=aes(x=age, y=term.r, group=.draw, color='Repair rate\n terms'), alpha=0.1) +
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
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    facet_grid(~sex) +
    labs(y='', x='Age (months)',
         color='', fill='') +
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5) +
    scale_color_manual(values=c('Damage rate\n terms' = rates.palette[1], 'Repair rate\n terms'=rates.palette[2])) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + ggtitle('c) Mouse dataset 2 (Bisset et al. 2021)') +
    geom_text(
        data    = exercise.text,
        mapping = aes(x = 22, y = 0.59, label = TeX(label,output='character'), fontface = "bold"),size=2.25,color='black', 
        hjust   = 0,
        vjust   = 0, parse=TRUE)


##### mouse dataset 3
schultz <- read.csv('../datasets/schultz_data.csv', header = TRUE, sep = ",")

schultz$age <- schultz$time + schultz$baseline.age
schultz$prepair <- (schultz$repair/schultz$n)/schultz$delta.t
schultz$pdamage <- (schultz$damage/(schultz$N - schultz$n))/schultz$delta.t
schultz <- na.omit(schultz)

f.sd <- sd(schultz$f)


fit <- readRDS('../fits/mouse_3.rds')

bins <- seq(20, 38, by=1.5)
binned.test <- cut(schultz$age, bins, include.lowest = TRUE)
binned.test <- midpoints(binned.test)

# compute mouse dataset 3 terms
terms.schultz <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n], n = 100) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = schultz[n, 'sex']) %>%
    mutate(f = schultz[n, 'f']) %>%
    mutate(mouse = schultz[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r = -dotf*lambda_r- f*deriv_r) %>%
    mutate(term.d = -dotf*lambda_d + (1-f)*deriv_d) %>%
    ungroup() %>%
    group_by(sex, age, .draw) %>%
    summarize(term.r = mean(term.r), term.d = mean(term.d)) %>%
    ungroup() %>% as.data.frame()

# compute mouse dataset 3 significance
schultz.diff.test <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_d[n], deriv_r_f[n], deriv_d_f[n]) %>%
    group_by(.iteration, .chain) %>%
    mutate(sex = schultz[n, 'sex']) %>%
    mutate(f = schultz[n, 'f']) %>%
    mutate(mouse = schultz[n,'mouse']) %>%
    mutate(age = binned.test[n]) %>%
    mutate(deriv_r_f = deriv_r_f * f.sd) %>%
    mutate(deriv_d_f = deriv_d_f * f.sd) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(dotdotf = -dotf*(lambda_d+lambda_r) - f*(deriv_r + dotf*deriv_r_f) + (1-f)*(deriv_d + dotf*deriv_d_f))  %>%
    mutate(term.r = -dotf*lambda_r- f*deriv_r) %>%
    mutate(term.d = -dotf*lambda_d + (1-f)*deriv_d) %>%
    mutate(diff = term.d - term.r) %>%
    ungroup() %>%
    group_by(sex, age, .draw) %>%
    summarize(diff = mean(diff)) %>%
    ungroup() %>%
    group_by(sex, age) %>%
    median_hdci(diff, .width = 0.95) %>%
    ungroup() %>%
    group_by(sex, age) %>%
    summarize(p = mean(diff > 0)) %>% as.data.frame()

schultz.diff.test
schultz.text <- data.frame(sex=c('Male'), label = c('$p<0.05$ for ages$<35$'))


plot.schultz <- ggplot() +
    geom_line(data=terms.schultz, mapping=aes(x=age, y=term.d, group=.draw, color='Damage rate\n terms'), alpha=0.1) +
    geom_line(data=terms.schultz, mapping=aes(x=age, y=term.r, group=.draw, color='Repair rate\n terms'), alpha=0.1) +
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
        plot.margin = unit(c(0, 0.06, 0, 0), "cm"))+
    facet_grid(~sex) +
    labs(y='', x='Age (months)',
         color='', fill='') +
    scale_color_nejm() +
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36),
                       label = c("22", "24", "26", "28", "30", "32", "34", "36")) + ggtitle('d) Mouse dataset 3 (Schultz et al. 2020)') +
    scale_color_manual(values=c('Damage rate\n terms' = rates.palette[1], 'Repair rate\n terms'=rates.palette[2])) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) +
    geom_text(
        data    = schultz.text,
        mapping = aes(x = 27, y = 0.47, label = TeX(label,output='character'), fontface = "bold"),size=2.25,color='black', 
        hjust   = 0,
        vjust   = 0, parse=TRUE)

img <- plot_grid(plot.enalapril, plot.exercise, plot.schultz, ncol=3,nrow=1, rel_widths = c(1.675/5, 1.675/5, 1.65/5))

ggsave("../plots/mice_curvature_control_terms.pdf", img, width=8, height=1.5)
