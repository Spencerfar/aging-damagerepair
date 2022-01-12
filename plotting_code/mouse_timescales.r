library(rstan)
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
library(splines2)
library(pammtools)
library(glrt)
source("../utils/functions.r")
source("../utils/palettes.r")

###### mouse dataset 1
mice <- read.csv('../datasets/enalapril_repair_deficits.csv', header = TRUE, sep = ",")
mice$diff <- mice$repair.age - mice$damage.age
deficits <- unique(mice$deficit)
mice$sex <- as.factor(mice$sex)
mice$treatment <- as.factor(mice$treatment)
mice$deficit <- as.factor(mice$deficit)

mice$damage.time <- mice$damage.age - mice$baseline.age

rstan_options(auto_write = TRUE)


levels(mice$sex) <- c('F'='Female', 'M'='Male')
levels(mice$treatment) <- c('control'='Control', 'drug'='Enalapril')


mice$group <- NA
counter <- 1
for(s in levels(mice$sex)) {
    for(t in levels(mice$treatment)) {
	mice[(mice$sex == s) & (mice$treatment == t), 'group'] <- counter
	counter <- counter + 1
    }	 
}


mice.censored <- mice[mice$status == 0,]
mice.dead <- mice[mice$status == 1,]


input <- mice[,c('diff.left', 'diff.right', 'treatment', 'sex')]
input = filter(input, !((mice$diff.left < 0.0000001) & (is.infinite(mice$diff.right))))
input$treatment <- as.numeric(input$treatment)-1

male.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Male',c('diff.left', 'diff.right', 'treatment')])$p*1000)/1000)
female.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Female',c('diff.left', 'diff.right', 'treatment')])$p*1000)/1000)

pvals <- data.frame(pval = c(male.pval, female.pval), sex=c('Male', 'Female'))


counter <- 1
combined.list <- list()
for (g in unique(mice$group)) {

    selected = mice[mice$group == g,]
    
    selected.censored <- mice.censored[mice.censored$group == g,]
    selected.dead <- mice.dead[mice.dead$group == g,]
    
    test.times <- sort(unique(append(unique(selected$diff.left), unique(selected$diff.right[is.finite(selected$diff.right)]))))
    test.times <- test.times[test.times>0]
    
    num_knots <- 30 # number of knots for fitting
    spline_degree <- 3
    
    knots <- unname(quantile( selected$diff.left, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
    knots <- knots[knots > 0]
    boundary.knots <- c(0, max(selected$diff.right[is.finite(selected$diff.right)])+0.01)
    
    splines.lower <- iSpline(selected.dead$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.upper <- iSpline(selected.dead$diff.right,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.right <- iSpline(selected.censored$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.test <- iSpline(test.times,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)

    num_basis <- dim(splines.lower)[2]
    Ntest = length(test.times)
    
    stan.data <- list(n_censored=dim(selected.censored)[1],
                      n_dead=dim(selected.dead)[1],
                      spline_alpha=1,
                      isplines_lower = splines.lower,
                      isplines_upper = splines.upper,
                      isplines_right = splines.right,
                      isplines_test = splines.test,
                      num_basis=num_basis,
                      Ntest=Ntest)

    fit <- stan(file = '../models/bayes_spline_surv.stan', 
                data = stan.data, chains=2, 
                iter = 4000, cores = 2, warmup=1000,
                control=list(adapt_delta=0.9, max_treedepth=15),
                refresh=-1)
    
    test.times <- append(0, test.times)
    
    combined.fit <- fit %>%
        spread_draws(survival[n]) %>%
        group_by(.chain, .iteration) %>%
        mutate(time = test.times[n]) %>%
        mutate(survival = survival) %>%
        ungroup() %>%
        group_by(time) %>%
        median_hdci(survival, .width=0.95) %>%
        as.data.frame()
    
    combined.fit$sex <- unique(selected$sex)
    combined.fit$treatment <- unique(selected$treatment)
    
    combined.list[[counter]] <- combined.fit
    
    counter <- counter + 1
    
}


combined.fit <- do.call(rbind, combined.list)

enalapril.repair <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y= survival, 
                                       color = treatment),
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper, fill=treatment),
                           size=0.75, alpha = 0.5) +
	  facet_grid(.~sex) +
	  labs(x = 'Months since damage occured', y = 'Probability of remaining damaged', color='', fill='') +
    	  theme_cowplot() + scale_fill_d3() + scale_color_d3() +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position=c(.68,0.9),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) + ggtitle("b) Mouse dataset 1\n    (Keller et al. 2019)") +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       label = c("0.0", "0.25", "0.5", "0.75", "1.0"), limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8),
                       label = c("0", "2", "4", "6", "8")) +
    coord_cartesian(xlim = c(0, 9)) +
    geom_text(
        data    = pvals,
        mapping = aes(x = 0, y = 0.05, label = pval), size=2.25,
        hjust   = 0,
        vjust   = 0)






mice <- read.csv('../datasets/enalapril_damage_deficits.csv', header = TRUE, sep = ",")
mice$diff <- mice$damage.age - mice$repair.age
deficits <- unique(mice$deficit)
mice$sex <- as.factor(mice$sex)
mice$treatment <- as.factor(mice$treatment)
mice$deficit <- as.factor(mice$deficit)

levels(mice$sex) <- c('F'='Female', 'M'='Male')
levels(mice$treatment) <- c('control'='Control', 'drug'='Enalapril')

mice$group <- NA
counter <- 1
for(s in levels(mice$sex)) {
    for(t in levels(mice$treatment)) {
	mice[(mice$sex == s) & (mice$treatment == t), 'group'] <- counter
	counter <- counter + 1
    }	 
}

mice.censored <- mice[mice$status == 0,]
mice.dead <- mice[mice$status == 1,]


input <- mice[,c('diff.left', 'diff.right', 'treatment', 'sex')]
input = filter(input, !((mice$diff.left < 0.0000001) & (is.infinite(mice$diff.right))))
input$treatment <- as.numeric(input$treatment)-1

male.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Male',c('diff.left', 'diff.right', 'treatment')])$p*1000)/1000)
female.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Female',c('diff.left', 'diff.right', 'treatment')])$p*1000)/1000)

pvals <- data.frame(pval = c(male.pval, female.pval), sex=c('Male', 'Female'))



counter <- 1
combined.list <- list()
for (g in unique(mice$group)) {
    
    selected = mice[mice$group == g,]
    selected.censored <- mice.censored[mice.censored$group == g,]
    selected.dead <- mice.dead[mice.dead$group == g,]
    
    test.times <- sort(unique(append(unique(selected$diff.left), unique(selected$diff.right[is.finite(selected$diff.right)]))))
    test.times <- test.times[test.times > 0]
    
    num_knots <- 30 # number of knots for fitting
    spline_degree <- 3
    
    knots <- unname(quantile( selected$diff.left, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
    knots <- knots[knots > 0]
    boundary.knots <- c(0, max(selected$diff.right[is.finite(selected$diff.right)])+0.01)
    
    splines.lower <- iSpline(selected.dead$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.upper <- iSpline(selected.dead$diff.right,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.right <- iSpline(selected.censored$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.test <- iSpline(test.times,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)

    num_basis <- dim(splines.lower)[2]
    Ntest = length(test.times)
    
    stan.data <- list(n_censored=dim(selected.censored)[1],
                      n_dead=dim(selected.dead)[1],
                      spline_alpha=1,
                      isplines_lower = splines.lower,
                      isplines_upper = splines.upper,
                      isplines_right = splines.right,
                      isplines_test = splines.test,
                      num_basis=num_basis,
                      Ntest=Ntest)
    
    fit <- stan(file = '../models/bayes_spline_surv.stan', 
                data = stan.data, chains=2, 
                iter = 4000, cores = 2, warmup=1000,
                control=list(adapt_delta=0.9, max_treedepth=15),
                refresh=-1)
    
    
    test.times <- append(0, test.times)
    
    combined.fit <- fit %>%
        spread_draws(survival[n]) %>%
        group_by(.chain, .iteration) %>%
        mutate(time = test.times[n]) %>%
        mutate(survival = survival) %>%
        ungroup() %>%
        group_by(time) %>%
        median_hdci(survival, .width=0.95) %>% as.data.frame()

    
    combined.fit$sex <- unique(selected$sex)
    combined.fit$treatment <- unique(selected$treatment)
    
    combined.list[[counter]] <- combined.fit
    
    counter <- counter + 1
    
}


combined.fit <- do.call(rbind, combined.list)


enalapril.damage <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y= survival, 
                                       color = treatment),
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper, fill=treatment),
                           size=0.75, alpha = 0.5) +
	  facet_grid(.~sex) +
	  labs(x = 'Months since repair occured', y = 'Probability of remaining undamaged', color='', fill='') +
    	  theme_cowplot() + scale_fill_d3() + scale_color_d3() +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position=c(.74,0.9),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       label = c("0.0", "0.25", "0.5", "0.75", "1.0"), limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8),
                       label = c("0", "2", "4", "6", "8")) +# +
    coord_cartesian(xlim = c(0, 9)) +
    geom_text(
        data    = pvals,
        mapping = aes(x = 0, y = 0.05, label = pval), size=2.25,
        hjust   = 0, vjust   = 0)




##### mouse dataset 2
mice <- read.csv('../datasets/exercise_repair_deficits.csv', header = TRUE, sep = ",")
mice$diff <- mice$repair.age - mice$damage.age

deficits <- unique(mice$deficit)
mice$sex <- as.factor(mice$sex)
mice$exercise <- as.factor(mice$exercise)

levels(mice$sex) <- c('F'='Female', 'M'='Male')
levels(mice$exercise) <- c('no'='Control', 'yes'='Exercise')



mice$group <- NA
counter <- 1
for(s in levels(mice$sex)) {
    for(t in levels(mice$exercise)) {
	mice[(mice$sex == s) & (mice$exercise == t), 'group'] <- counter
	counter <- counter + 1
    }	 
}

mice.censored <- mice[mice$status == 0,]
mice.dead <- mice[mice$status == 1,]


input <- mice[,c('diff.left', 'diff.right', 'exercise', 'sex')]
input = filter(input, !((mice$diff.left < 0.0000001) & (is.infinite(mice$diff.right))))
input$exercise <- as.numeric(input$exercise)-1

male.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Male',c('diff.left', 'diff.right', 'exercise')])$p*1000)/1000)
female.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Female',c('diff.left', 'diff.right', 'exercise')])$p*1000)/1000)

pvals <- data.frame(pval = c(male.pval, female.pval), sex=c('Male', 'Female'))



counter <- 1
combined.list <- list()
for (g in unique(mice$group)) {
    
    selected = mice[mice$group == g,]
    selected.censored <- mice.censored[mice.censored$group == g,]
    selected.dead <- mice.dead[mice.dead$group == g,]
    
    test.times <- sort(unique(append(unique(selected$diff.left), unique(selected$diff.right[is.finite(selected$diff.right)]))))
    test.times <- test.times[test.times > 0]
    
    num_knots <- 30 # number of knots for fitting
    spline_degree <- 3

    knots <- unname(quantile( selected$diff.left, probs=seq(from=0.1, to=0.90, length.out = num_knots)))
    knots <- knots[knots > 0]
boundary.knots <- c(0, max(selected$diff.right[is.finite(selected$diff.right)])+0.01)

    splines.lower <- iSpline(selected.dead$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.upper <- iSpline(selected.dead$diff.right,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    
    splines.right <- iSpline(selected.censored$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.test <- iSpline(test.times,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)

    num_basis <- dim(splines.lower)[2]
    Ntest = length(test.times)
    
    stan.data <- list(n_censored=dim(selected.censored)[1],
                      n_dead=dim(selected.dead)[1],
                      spline_alpha=1,
                      isplines_lower = splines.lower,
                      isplines_upper = splines.upper,
                      isplines_right = splines.right,
                      isplines_test = splines.test,
                      num_basis=num_basis,
                      Ntest=Ntest)
    
    fit <- stan(file = '../models/bayes_spline_surv.stan', 
                data = stan.data, chains=2, 
                iter = 4000, cores = 2, warmup=1000,
                control=list(adapt_delta=0.9, max_treedepth=15),
                refresh=-1)
    
    test.times <- append(0, test.times)

    combined.fit <- fit %>%
        spread_draws(survival[n]) %>%
        group_by(.chain, .iteration) %>%
        mutate(time = test.times[n]) %>%
        mutate(survival = survival) %>%
        ungroup() %>%
        group_by(time) %>%
        median_hdci(survival, .width=0.95) %>%
        as.data.frame()
    
    
    combined.fit$sex <- unique(selected$sex)
    combined.fit$exercise <- unique(selected$exercise)
    
    combined.list[[counter]] <- combined.fit
    
    counter <- counter + 1
    
}

combined.fit <- do.call(rbind, combined.list)

exercise.repair <- ggplot() +
	  geom_step(data = combined.fit,
                mapping=aes(x = time, y=survival, 
                            color = exercise),
                size=0.75, alpha = 1) +
    geom_stepribbon(data = combined.fit,
                    mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper, fill=exercise),
                    size=0.75, alpha = 0.5) +
    facet_grid(.~sex) +
    labs(x = 'Months since damage occured', y = '', color='', fill='') +
    theme_cowplot() +
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette) +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position=c(.74,0.9),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) + ggtitle('c) Mouse dataset 2\n    (Bisset et al. 2021)') +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       label = c("0.0", "0.25", "0.5", "0.75", "1.0"), limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3")) +
    coord_cartesian(xlim = c(0, 3.5)) +
    geom_text(
        data    = pvals,
        mapping = aes(x = 0, y = 0.05, label = pval), size=2.25,
        hjust   = 0,
        vjust   = 0)




mice <- read.csv('../datasets/exercise_damage_deficits.csv', header = TRUE, sep = ",")
mice$diff <- mice$damage.age - mice$repair.age

deficits <- unique(mice$deficit)
mice$sex <- as.factor(mice$sex)
mice$exercise <- as.factor(mice$exercise)

survdiff(Surv(diff, status) ~ exercise, data=mice, subset=sex=='M')
survdiff(Surv(diff, status) ~ exercise, data=mice, subset=sex=='F')

levels(mice$sex) <- c('F'='Female', 'M'='Male')
levels(mice$exercise) <- c('no'='Control', 'yes'='Exercise')

mice$group <- NA
counter <- 1
for(s in levels(mice$sex)) {
    for(t in levels(mice$exercise)) {
	mice[(mice$sex == s) & (mice$exercise == t), 'group'] <- counter
	counter <- counter + 1
    }	 
}


mice.censored <- mice[mice$status == 0,]
mice.dead <- mice[mice$status == 1,]


input <- mice[,c('diff.left', 'diff.right', 'exercise', 'sex')]
input = filter(input, !((mice$diff.left < 0.0000001) & (is.infinite(mice$diff.right))))
input$exercise <- as.numeric(input$exercise)-1

male.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Male',c('diff.left', 'diff.right', 'exercise')])$p*1000)/1000)
female.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Female',c('diff.left', 'diff.right', 'exercise')])$p*1000)/1000)

pvals <- data.frame(pval = c(male.pval, female.pval), sex=c('Male', 'Female'))


counter <- 1
combined.list <- list()
for (g in unique(mice$group)) {
    
    selected = mice[mice$group == g,]
    selected.censored <- mice.censored[mice.censored$group == g,]
    selected.dead <- mice.dead[mice.dead$group == g,]
    
    test.times <- sort(unique(append(unique(selected$diff.left), unique(selected$diff.right[is.finite(selected$diff.right)]))))
    test.times <- test.times[test.times > 0]
    
    num_knots <- 30 # number of knots for fitting
    spline_degree <- 3
    
    knots <- unname(quantile( selected$diff.left, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
    knots <- knots[knots > 0]
    boundary.knots <- c(0, max(selected$diff.right[is.finite(selected$diff.right)])+0.01)
    
    splines.lower <- iSpline(selected.dead$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.upper <- iSpline(selected.dead$diff.right,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.right <- iSpline(selected.censored$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.test <- iSpline(test.times,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)

    num_basis <- dim(splines.lower)[2]
    Ntest = length(test.times)
    
    stan.data <- list(n_censored=dim(selected.censored)[1],
                      n_dead=dim(selected.dead)[1],
                      spline_alpha=1,
                      isplines_lower = splines.lower,
                      isplines_upper = splines.upper,
                      isplines_right = splines.right,
                      isplines_test = splines.test,
                      num_basis=num_basis,
                      Ntest=Ntest)
    
    fit <- stan(file = '../models/bayes_spline_surv.stan', 
                data = stan.data, chains=2, 
                iter = 4000, cores = 2, warmup=1000,
                control=list(adapt_delta=0.9, max_treedepth=15),
                refresh=-1)
    
    test.times <- append(0, test.times)
    
    combined.fit <- fit %>%
        spread_draws(survival[n]) %>%
        group_by(.chain, .iteration) %>%
        mutate(time = test.times[n]) %>%
        mutate(survival = survival) %>%
        ungroup() %>%
        group_by(time) %>%
        median_hdci(survival, .width=0.95) %>% as.data.frame()
    
    
    combined.fit$sex <- unique(selected$sex)
    combined.fit$exercise <- unique(selected$exercise)
    
    combined.list[[counter]] <- combined.fit
    
    counter <- counter + 1
    
}

combined.fit <- do.call(rbind, combined.list)

exercise.damage <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y=survival, 
                                       color = exercise),
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper, fill=exercise),
                           size=0.75, alpha = 0.5) +
	  facet_grid(.~sex) +
	  labs(x = 'Months since repair occured', y = '', color='', fill='') +
    theme_cowplot() +
    scale_fill_manual(values = exercise.palette) +
    scale_color_manual(values = exercise.palette) +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position=c(.72,0.9),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm"))  +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       label = c("0.0", "0.25", "0.5", "0.75", "1.0"), limits = c(0.0, 1.0)) +# +
    scale_x_continuous(breaks=c(0, 1, 2, 3),
                       label = c("0", "1", "2", "3")) +
    coord_cartesian(xlim = c(0, 3.5)) +
    geom_text(
        data    = pvals,
        mapping = aes(x = 0, y = 0.05, label = pval), size=2.25,
        hjust   = 0,
        vjust   = 0)





###### mouse dataset 3
mice <- read.csv('../datasets/schultz_repair_deficits.csv', header = TRUE, sep = ",")
mice$diff <- mice$repair.age - mice$damage.age

deficits <- unique(mice$deficit)
mice$sex <- as.factor(mice$sex)

mice$group <- NA
counter <- 1
for(s in levels(mice$sex)) {
	mice[(mice$sex == s), 'group'] <- counter
	counter <- counter + 1 
}

mice.censored <- mice[mice$status == 0,]
mice.dead <- mice[mice$status == 1,]

counter <- 1
combined.list <- list()
for (g in unique(mice$group)) {
    
    selected = mice[mice$group == g,]
    selected.censored <- mice.censored[mice.censored$group == g,]
    selected.dead <- mice.dead[mice.dead$group == g,]
    
    test.times <- sort(unique(append(unique(selected$diff.left), unique(selected$diff.right[is.finite(selected$diff.right)]))))
    test.times <- test.times[test.times > 0]
    
    num_knots <- 30 # number of knots for fitting
    spline_degree <- 3
    
    knots <- unname(quantile( selected$diff.left, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
    knots <- knots[knots > 0]
boundary.knots <- c(0, max(selected$diff.right[is.finite(selected$diff.right)])+0.01)

    splines.lower <- iSpline(selected.dead$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.upper <- iSpline(selected.dead$diff.right,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.right <- iSpline(selected.censored$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.test <- iSpline(test.times,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)

    num_basis <- dim(splines.lower)[2]
    Ntest = length(test.times)
    
    
    stan.data <- list(n_censored=dim(selected.censored)[1],
                      n_dead=dim(selected.dead)[1],
                      spline_alpha=1,
                      isplines_lower = splines.lower,
                      isplines_upper = splines.upper,
                      isplines_right = splines.right,
                      isplines_test = splines.test,
                      num_basis=num_basis,
                      Ntest=Ntest)

    fit <- stan(file = '../models/bayes_spline_surv.stan', 
                data = stan.data, chains=2, 
                iter = 4000, cores = 2, warmup=1000,
                control=list(adapt_delta=0.9, max_treedepth=15),
                refresh=-1)

    test.times <- append(0, test.times)
    
    combined.fit <- fit %>%
        spread_draws(survival[n]) %>%
        group_by(.chain, .iteration) %>%
        mutate(time = test.times[n]) %>%
        mutate(survival = survival) %>%
        ungroup() %>%
        group_by(time) %>%
        median_hdci(survival, .width=0.95) %>%
        as.data.frame()
    

    combined.fit$sex <- unique(selected$sex)
    combined.list[[counter]] <- combined.fit
    counter <- counter + 1
    
}

combined.fit <- do.call(rbind, combined.list)

schultz.repair <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y=survival),color='#1c9099',
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper),
                           size=0.75, alpha = 0.5, fill='#1c9099') +
	  facet_grid(.~sex) +
	  labs(x = 'Months since damage occured', y = '', color='', fill='') +
    theme_cowplot() +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position=c(.68,0.9),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) + ggtitle('d) Mouse dataset 3\n    (Schultz et al. 2020)') +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       label = c("0.0", "0.25", "0.5", "0.75", "1.0"), limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8),
                       label = c("0", "2", "4", "6", "8")) +
    coord_cartesian(xlim = c(0, 8.5))









mice <- read.csv('../datasets/schultz_damage_deficits.csv', header = TRUE, sep = ",")
mice$diff <- mice$damage.age - mice$repair.age

deficits <- unique(mice$deficit)
mice$sex <- as.factor(mice$sex)

mice$group <- NA
counter <- 1
for(s in levels(mice$sex)) {
	mice[(mice$sex == s), 'group'] <- counter
	counter <- counter + 1 
}

mice.censored <- mice[mice$status == 0,]
mice.dead <- mice[mice$status == 1,]

counter <- 1
combined.list <- list()
for (g in unique(mice$group)) {
    
    selected = mice[mice$group == g,]
    selected.censored <- mice.censored[mice.censored$group == g,]
    selected.dead <- mice.dead[mice.dead$group == g,]
    
    test.times <- sort(unique(append(unique(selected$diff.left), unique(selected$diff.right[is.finite(selected$diff.right)]))))
    test.times <- test.times[test.times > 0]
    
    num_knots <- 30 # number of knots for fitting
    spline_degree <- 3
    
    knots <- unname(quantile( selected$diff.left, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
    knots <- knots[knots > 0]
boundary.knots <- c(0, max(selected$diff.right[is.finite(selected$diff.right)])+0.01)

    splines.lower <- iSpline(selected.dead$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.upper <- iSpline(selected.dead$diff.right,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.right <- iSpline(selected.censored$diff.left,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)
    splines.test <- iSpline(test.times,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = TRUE)

    num_basis <- dim(splines.lower)[2]
    Ntest = length(test.times)
    
    stan.data <- list(n_censored=dim(selected.censored)[1],
                      n_dead=dim(selected.dead)[1],
                      spline_alpha=1,
                      isplines_lower = splines.lower,
                      isplines_upper = splines.upper,
                      isplines_right = splines.right,
                      isplines_test = splines.test,
                      num_basis=num_basis,
                      Ntest=Ntest)
    
    fit <- stan(file = '../models/bayes_spline_surv.stan', 
            data = stan.data, chains=2, 
            iter = 4000, cores = 2, warmup=1000,
            control=list(adapt_delta=0.9, max_treedepth=15),
            refresh=-1)
    
    test.times <- append(0, test.times)
    
    combined.fit <- fit %>%
        spread_draws(survival[n]) %>%
        group_by(.chain, .iteration) %>%
        mutate(time = test.times[n]) %>%
        mutate(survival = survival) %>%
        ungroup() %>%
        group_by(time) %>%
        median_hdci(survival, .width=0.95) %>%
        as.data.frame()
    
    
    combined.fit$sex <- unique(selected$sex)
    combined.list[[counter]] <- combined.fit
    counter <- counter + 1
    
}

combined.fit <- do.call(rbind, combined.list)


schultz.damage <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y=survival),color='#1c9099',
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper),
                           size=0.75, alpha = 0.5, fill='#1c9099') +
	  facet_grid(.~sex) +
	  labs(x = 'Months since repair occured', y = '', color='', fill='') +
    theme_cowplot() +
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=6.),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position=c(.68,0.9),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(0.01, 0, 0, 0.01), "cm")) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       label = c("0.0", "0.25", "0.5", "0.75", "1.0"), limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8),
                       label = c("0", "2", "4", "6", "8")) +
    coord_cartesian(xlim = c(0, 8.5))



img <- plot_grid(enalapril.repair, exercise.repair, schultz.repair,
                  nrow = 1, ncol=3, rel_widths = c(1.95/5, 1.95/5, 1.1/5))

ggsave("../plots/mouse_timescales_resilience.pdf", img, , width=8, height=2)


enalapril <- plot_grid(enalapril.repair, enalapril.damage, nrow=2, ncol=1, rel_heights=c(1.65/3, 1.35/3))
exercise <- plot_grid(exercise.repair, exercise.damage, nrow=2, ncol=1, rel_heights=c(1.65/3, 1.35/3))
schultz <- plot_grid(schultz.repair, schultz.damage, nrow=2, ncol=1, rel_heights=c(1.65/3, 1.35/3))

img <- plot_grid(enalapril, exercise, schultz,
                 nrow = 1, ncol=3, rel_widths = c(1.95/5, 1.95/5, 1.1/5))

ggsave("../plots/mouse_timescales_both.pdf", img, , width=8, height=3.)
