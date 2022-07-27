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
enalapril <- read.csv('../datasets/enalapril_data_pruned.csv', header = TRUE, sep = ",")

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
        f.num = length(f))  %>%
  mutate(age = time.bin) %>% as.data.frame()


binned.data.rates$lower.repair <- binned.data.rates$mean.repair-binned.data.rates$se.repair
binned.data.rates$upper.repair <- binned.data.rates$se.repair+binned.data.rates$mean.repair

binned.data.rates$lower.damage <- binned.data.rates$mean.damage-binned.data.rates$se.damage
binned.data.rates$upper.damage <- binned.data.rates$se.damage+binned.data.rates$mean.damage

binned.data.rates$lower.f <- binned.data.rates$mean.f-binned.data.rates$se.f
binned.data.rates$upper.f <- binned.data.rates$se.f+binned.data.rates$mean.f

binned.data.rates <- binned.data.rates[binned.data.rates$treatment=='control',]




# fix labels
binned.data.rates$sex <- as.factor(binned.data.rates$sex)
levels(binned.data.rates$sex) <- c('Female', 'Male')
enalapril$sex <- as.factor(enalapril$sex)
levels(enalapril$sex) <- c('Female', 'Male')

write.csv(binned.data.rates, '../figure_data/figure5_supplement3/mouse1_binned_rates.csv')


# plot repair rates
enalapril.repair <- ggplot() +
    geom_point(binned.data.rates %>% filter(repair.num >= 0), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                              ymax=upper.repair, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(repair.num >= 0), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                             ymax=upper.repair, group=sex, color=sex), width=0.1) +
    geom_smooth(data=enalapril, mapping=aes(x=age, y=prepair, group=sex, color=sex), width=0.1, alpha=0, method='lm') +
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
                       label = c("16", "18", "20", "22", "24"))  + ggtitle("a) Mouse dataset 1 (Keller et al. 2019)") 


# plot damage rates
enalapril.damage <- ggplot() +
    geom_point(binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                              ymax=upper.damage, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                             ymax=upper.damage, group=sex, color=sex), width=0.1) +
    geom_smooth(data=enalapril, mapping=aes(x=age, y=pdamage, group=sex, color=sex), width=0.1, alpha=0, method='lm') +
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
                       label = c("16", "18", "20", "22", "24"))


# plot frailty index
enalapril.f <- ggplot() +
    geom_point(binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                              ymax=upper.f, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                        ymax=upper.f, group=sex, color=sex), width=0.1) +
    geom_smooth(data=enalapril, mapping=aes(x=age, y=f, group=sex, color=sex), width=0.1, alpha=0, method='lm') +
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
exercise <- read.csv('../datasets/exercise_data_pruned.csv', header = TRUE, sep = ",")

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


# fix labels
binned.data.rates$sex <- as.factor(binned.data.rates$sex)
levels(binned.data.rates$sex) <- c('Female', 'Male')
levels(exercise$sex) <- c('Female', 'Male')

write.csv(binned.data.rates, '../figure_data/figure5_supplement3/mouse2_binned_rates.csv')

# plot repair
exercise.repair <- ggplot() +
    geom_point(data=binned.data.rates %>% filter(repair.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                              ymax=upper.repair, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(repair.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                                  ymax=upper.repair, group=sex, color=sex), width=0.1) +
    geom_smooth(data=exercise %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=prepair, group=sex, color=sex), width=0.1, alpha=0, method='lm') +
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
    scale_fill_manual(values=sex.palette) + ggtitle('b) Mouse dataset 2 (Bisset et al. 2022)') +
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)
exercise.repair


# plot damage
exercise.damage <- ggplot() +
    geom_point(binned.data.rates %>% filter(damage.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                              ymax=upper.damage, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(damage.num >= 0) %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                             ymax=upper.damage, group=sex, color=sex), width=0.1) +
    geom_smooth(data=exercise %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=pdamage, group=sex, color=sex), width=0.1, alpha=0., method='lm') +
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
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)
exercise.damage

# plot frailty index
exercise.f <- ggplot() +
    geom_point(binned.data.rates %>% filter(f.num >= 0) %>% filter((age >= 21)), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                              ymax=upper.f, group=sex, color=sex, shape=sex), size=2) +
    geom_errorbar(binned.data.rates %>% filter(f.num >= 0) %>% filter((age >= 21)), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                             ymax=upper.f, group=sex, color=sex), width=0.1) +
    geom_smooth(data=exercise %>% filter((age >= 21) & (age <= 26)), mapping=aes(x=age, y=f, group=sex, color=sex), width=0.1, alpha=0., method='lm') +
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
    scale_x_continuous(breaks=c(21,22,23,24,25,26), label=c("21","22", "23", "24", "25", "26"))+ xlim(21, 26.5)+ ggtitle('b) Mouse dataset 2 (Bisset et al. 2022)')


##### mouse dataset 3
schultz <- read.csv('../datasets/schultz_data_pruned.csv', header = TRUE, sep = ",")

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


write.csv(binned.data.rates, '../figure_data/figure5_supplement3/mouse3_binned_rates.csv')


df.legend <- binned.data.rates
df.legend$mean.repair <- NaN
df.legend$mean.damage <- NaN
df.legend$mean.f <- NaN
df.legend$sex <- 'Male'
df.legend[1,'sex'] <- 'Female'
levels(df.legend) <- c('Male', 'Female')
levels(schultz$sex) <- c('Male', 'Female')

# plot repair
schultz.rep <- ggplot() + geom_point(data=df.legend, mapping=aes(x=age, y=mean.repair, color=sex, shape=sex), size=2) +
    geom_point(data=binned.data.rates %>% filter(repair.num >= 0) , mapping=aes(x=age, y=mean.repair, color=sex,shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(repair.num >= 0) , mapping=aes(x=age, y=mean.repair, ymin=lower.repair, 
                                                                                   ymax=upper.repair, color=sex,shape=sex), width=0.1) +
    geom_smooth(data=schultz, mapping=aes(x=age, y=prepair, group=sex, color=sex), width=0.1, alpha=0., method='lm') +
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
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36), label = c("22", "24", "26", "28", "30", "32", "34", "36")) + ggtitle('c) Mouse dataset 3 (Schultz et al. 2020)')

# plot damage
schultz.dam <- ggplot()  + geom_point(data=df.legend, mapping=aes(x=age, y=mean.repair, color=sex, shape=sex), size=2) +
    geom_point(data=binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, color=sex,shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(damage.num >= 0), mapping=aes(x=age, y=mean.damage, ymin=lower.damage, 
                                                                                  ymax=upper.damage, color=sex,shape=sex), width=0.1) +
    geom_smooth(data=schultz, mapping=aes(x=age, y=pdamage, group=sex, color=sex), width=0.1, alpha=0., method='lm') +
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
                       label = c("22", "24", "26", "28", "30", "32", "34", "36"))

# plot frailty index
schultz.f <- ggplot() + geom_point(data=df.legend, mapping=aes(x=age, y=mean.repair, color=sex, shape=sex), size=2) +
    geom_point(data=binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, color=sex,shape=sex), size=2) +
    geom_errorbar(data=binned.data.rates %>% filter(f.num >= 0), mapping=aes(x=age, y=mean.f, ymin=lower.f, 
                                                                             ymax=upper.f, color=sex,shape=sex), width=0.1) +
    geom_smooth(data=schultz, mapping=aes(x=age, y=f, group=sex, color=sex), width=0.1, alpha=0., method='lm') +
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

ggsave("mice_control_rates_sex_pruned.pdf", img, width=8, height=2.25)


img.enalapril <- grid.arrange(enalapril.repair, enalapril.damage,  ncol=1,nrow=2)
img.exercise <- grid.arrange(exercise.repair, exercise.damage,  ncol=1,nrow=2)
img.schultz <- grid.arrange(schultz.rep, schultz.dam,  ncol=1,nrow=2)

img <- plot_grid(enalapril.f, exercise.f, schultz.f, nrow = 1, ncol=3, rel_widths = c(1.52/5, 1.52/5, 1.96/5))

#ggsave("../plots/mice_control_f_sex.pdf", img, width=8, height=1.25)
