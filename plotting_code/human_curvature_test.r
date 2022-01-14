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

f.mean <- mean(elsa.test$f)
f.sd <- sd(elsa.test$f)

elsa.na <- na.omit(elsa.test)

elsa.na$prepair <- (elsa.na$repair/elsa.na$n)/elsa.na$delta.t
elsa.na$pdamage <- (elsa.na$damage/(elsa.na$N - elsa.na$n))/elsa.na$delta.t


bins <- seq(20, 105, by=2)
bins.time <- seq(0, 20, by=2)
elsa.na$baseline.num <- as.numeric(elsa.na$baseline.bin)


fit <- as_draws_df(readRDS("../fits/long_human_fit.RDS")$draws(c("lambda_r", "lambda_d", "deriv_r", "deriv_d", "deriv_r_f", "deriv_d_f")))


baseline.bins <- seq(50, 90, by=4)
binned.age <- cut(elsa.test$age, bins, include.lowest=TRUE)
binned.age <- midpoints(binned.age)+1

curv.plot.terms.control <- fit %>%
    spread_draws(lambda_r[n], lambda_d[n], deriv_r[n], deriv_r_f[n], deriv_d[n], deriv_d_f[n]) %>%
    group_by(.draw) %>%
    mutate(sex = elsa.test[n, 'sex']$sex) %>%
    mutate(age = binned.age[n]) %>%
    mutate(baseline.bin = elsa.test[n,'baseline.bin']$baseline.bin) %>%
    mutate(diff = ((1-(elsa.test[n,'f']$f))*(deriv_d+deriv_d_f*f.sd) - ((1 - (elsa.test[n,'f']$f))*lambda_d - (elsa.test[n,'f']$f)*lambda_r)*lambda_d) - (-(elsa.test[n,'f']$f)*(deriv_r +deriv_r_f*f.sd) -  ((1 - (elsa.test[n,'f']$f))*lambda_d - (elsa.test[n,'f']$f)*lambda_r)*lambda_r) ) %>%
    filter(elsa.test[n,'n']$n > 0) %>%
    ungroup() %>%
    group_by(baseline.bin, sex, age, .draw) %>%
    summarize(diff=mean(diff, na.rm=TRUE)) %>%
    ungroup() %>%
    group_by(baseline.bin, sex, age) %>%
    median_hdci(diff, .width = 0.95)

levels(curv.plot.terms.control$baseline.bin) <- c('Baseline age: [50,60)','[60,70)','[70,80)','[80,90)')

fdotdot.terms.control <- ggplot() +
    geom_errorbar(data=curv.plot.terms.control, mapping=aes(x=age, y=diff, color=sex, ymin=.lower, ymax=.upper, fill=sex)) +
    geom_point(data=curv.plot.terms.control, mapping=aes(x=age, y=diff, color=sex), alpha=1) +
    geom_hline(yintercept = 0, linetype="dotted") +
    facet_grid(sex~baseline.bin,scales='free') + theme_cowplot() +
    theme(strip.background = element_blank()) +
    labs(x="Age (years)", y="FI curvature terms difference", color='') +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0.01), "cm")) + ggtitle('f) ELSA humans (Phelps et al. 2020)') +
    scale_color_manual(values=c('Male' = sex.palette[2], 'Female'=sex.palette[1])) +
    scale_fill_manual(values=c('Male' = sex.palette[2], 'Female'=sex.palette[1])) +
    guides(colour = guide_legend(override.aes = list(size=1., alpha=1))) + ylim(-2.5,4.5)



ggsave("../plots/human_curvature_test.pdf", fdotdot.terms.control, width=8, height=2)
