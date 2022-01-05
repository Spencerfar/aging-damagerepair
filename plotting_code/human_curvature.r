library(posterior)
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
library(cowplot)
source("../utils/functions.r")
source("../utils/palettes.r")
    
# read elsa data
elsa <- read.csv('../datasets/human_data.csv', header = TRUE, sep = ",")
elsa$id <- as.character(elsa$id)
elsa$time <- elsa$age - elsa$baseline.age
elsa$f <- elsa$n/elsa$N
elsa$wealth <- log(elsa$wealth + mean(elsa$wealth))

# create baseline age bins
bins <- seq(50, 90, by=10)
elsa <- elsa %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE))

elsa$baseline.bin <- as.factor(elsa$baseline.bin)
elsa$sex <- as.factor(elsa$sex)
levels(elsa$sex) <- c('Male', 'Female')

elsa <- na.omit(elsa)

f.mean <- mean(elsa$f)
f.sd <- sd(elsa$f)

# create text for significance
labels <- c('', '','', '$p<0.05$ for age$<97$', '','','$p<0.05$ for age$>70$','$p<0.05$ for age$<95$')
xs <- c(55, 55, 75, 85, 55,55, 75,85)
human.text <- data.frame(baseline.bin=c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)', 'Baseline age: [50,60)','[60,70)','[70,80)','[80,90)'),sex=c('Female', 'Female', 'Female', 'Female', 'Male', 'Male','Male','Male'), label = labels, x= xs)
human.text$baseline.bin <- factor(human.text$baseline.bin, levels=c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)'))

# read stan fit
fit <- as_draws_df(readRDS("../fits/long_human_fit.RDS")$draws(c("lambda_r", "lambda_d", "deriv_r", "deriv_d", "deriv_r_f", "deriv_d_f")))

bins <- seq(20, 105, by=2)
binned.age <- cut(elsa$age, bins, include.lowest=TRUE)
binned.age <- midpoints(binned.age)+1

# compute terms
terms <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_r_f[n], deriv_d[n], deriv_d_f[n], n=100) %>%
    group_by(.draw) %>%
    mutate(sex = elsa[n, 'sex']$sex) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa[n,'baseline.bin']$baseline.bin) %>%
    mutate(f = elsa[n,'f']$f) %>%
    mutate(num = elsa[n,'n']$n) %>%
    mutate(deriv_r_f = deriv_r_f*f.sd) %>%
    mutate(deriv_d_f = deriv_d_f*f.sd) %>%
    filter(num > 0) %>%
    mutate(dotf = (1 - f)*lambda_d - f*lambda_r )  %>%
    mutate(term.r = -f*deriv_r-dotf*lambda_r) %>%
    mutate(term.d = (1-f)*deriv_d-dotf*lambda_d) %>% 
    ungroup() %>%
    group_by(baseline.bin, sex, age, .draw) %>%
    summarize(term.r=mean(term.r), term.d=mean(term.d)) %>%
    ungroup() %>%
    as.data.frame()

levels(terms$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')

plot.terms <- ggplot() +
    geom_hline(yintercept = 0, linetype="dotted") +
    geom_line(data=terms, mapping=aes(x=age, y=term.d, group=.draw, color='Damage rate\n terms'), alpha=0.1) +
    geom_line(data=terms, mapping=aes(x=age, y=term.r, group=.draw, color='Repair rate\n terms'), alpha=0.1) +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="Frailty Index curvature terms", color='') +
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
        legend.title = element_text(size = 7)) + ggtitle('e) ELSA humans (Phelps et al. 2020)') +scale_color_manual(values=c('Damage rate\n terms'=rates.palette[1], 'Repair rate\n terms'=rates.palette[2]))+ guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + ylim(-2.5,4.5) +
    geom_text(
        data    = human.text,
        mapping = aes(x = x, y = 3.75, label = TeX(label,output='character'), fontface = "bold"),size=2.25,color='black', 
        hjust   = 0,
        vjust   = 0, parse=TRUE)


ggsave("../plots/human_curvature_terms.pdf", plot.terms, width=8, height=2)
