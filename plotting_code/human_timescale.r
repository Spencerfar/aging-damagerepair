library(rstan)
library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(survival)
library(survminer)
library(tidybayes)
library(latex2exp)
library(tidyr)
library(cowplot)
library(ggfortify)
library(splines2)
library(pammtools)
#library(forcats)
library(glrt)
source("../utils/functions.r")
source("../utils/palettes.r")

elsa <- read.csv('../datasets/human_data.csv', header = TRUE, sep = ",")

elsa$id <- as.character(elsa$id)
elsa$time <- elsa$age - elsa$baseline.age

elsa <- elsa %>% group_by(id) %>% filter( (baseline.age >= 45)  & (age < 95))

elsa <- na.omit(elsa)

elsa$f <- elsa$n/elsa$N
mean.wealth <- mean(elsa$wealth)
elsa$wealth <- log(elsa$wealth + mean.wealth)

elsa <- elsa %>% group_by(id) %>% filter(all(delta.t<=4))
elsa <- elsa %>% group_by(id) %>% filter(n()>5) %>% as.data.frame()


bins <- seq(50, 90, by=10)
elsa <- 
elsa %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE))
elsa$wealth.bin <- as.factor(ntile(elsa$wealth, 3))
elsa$baseline.bin <- as.factor(elsa$baseline.bin)
levels(elsa$wealth.bin) <- c('poor', 'middle', 'rich')
elsa$sex <- as.factor(elsa$sex)
levels(elsa$sex) <- c('M', 'F')
elsa <- na.omit(elsa)



data <- read.csv('../datasets/human_repair_deficits.csv', header = TRUE, sep = ",")
data$wealth <- log(data$wealth + mean.wealth)


data <- data[data$baseline.age >= 50,]
data <- data[data$baseline.age < 90,]

bins <- seq(50, 90, by=10)
data <- 
data %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE)) %>%
    as.data.frame()

deficits <- c('d0'='Difficulty walking 100 yads','d1'='Difficulty sitting for 2 hours','d2'='Difficulty getting up from a chair',
              'd3'='Difficulty climbing several flights of stairs', 'd4'='Difficulty climbing one flight of stairs',
	      'd5'='Difficulty stooping, kneeling, or crouching',
	      'd6'='Difficulty reaching/extending arms', 'd7'='Difficulty pulling/pushing objects',
	      'd8'='Difficulty lifting/carrying 10 lbs weights', 'd9'='Difficulty picking up a coin',
	      'd10'='Difficulty dressing','d11'='Difficulty walking across a room','d12'='Difficulty bathing','d13'='Difficulty Eating',
	      'd14'='Difficulty getting in/out of bed','d15'='Difficulty using the toilet','d16'='Difficulty using a map',
	      'd17'='Difficulty preparing a hot meal','d18'='Difficulty shopping','d19'='Difficulty using telephone',
	      'd20'='Difficulty taking medications','d21'='Difficulty doing work around the house','d22'='Difficulty managing money')



names.deficits <- c()
for(i in 1:length(deficits)) {
      names.deficits <- append(names.deficits, paste('d', i-1, sep=""))
}
data$deficit <- as.factor(data$deficit)
levels(data$deficit) <- deficits

data[data$sex==0,]$sex <- 'M'
data[data$sex==1,]$sex <- 'F'
data$sex <- as.factor(data$sex)
data$diff <- data$repair.age - data$damage.age

data <- na.omit(data)

data$wealth.bin <- 
  as.factor(cut(data$wealth, breaks=c(-Inf, quantile(elsa$wealth, c(1/3, 2/3)), Inf)))


levels(data$wealth.bin) <- c("Bottom tercile","Middle tercile","Upper tercile")
levels(data$sex) <- c("Female", "Male")
    
rstan_options(auto_write = TRUE)

data$group <- NA
counter <- 1
for(w in levels(data$wealth.bin)) {
    for(s in levels(data$sex)) {
	    data[(data$sex == s) & (data$wealth.bin == w), 'group'] <- counter
	    counter <- counter + 1
    }	 
}


data.censored <- data[data$status == 0,]
data.dead <- data[data$status == 1,]


input <- data[,c('diff.left', 'diff.right', 'wealth.bin', 'sex')]
input = filter(input, !((data$diff.left < 0.0000001) & (is.infinite(data$diff.right))))
input$wealth.bin <- as.numeric(input$wealth.bin)-1


male.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Male',c('diff.left', 'diff.right', 'wealth.bin')],k=3)$p*1000)/1000)
female.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Female',c('diff.left', 'diff.right', 'wealth.bin')],k=3)$p*1000)/1000)

pvals <- data.frame(pval = c(male.pval, female.pval), sex=c('Male', 'Female'))


counter <- 1
combined.list <- list()
for (g in unique(data$group)) {
    
    selected = data[data$group == g,]
    selected.censored <- data.censored[data.censored$group == g,]
    selected.dead <- data.dead[data.dead$group == g,]
    
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
                data = stan.data, chains=4, 
                iter = 5000, cores = 4, warmup=2000, control=list(adapt_delta=0.9, max_treedepth=15), refresh=-1)

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
    combined.fit$wealth.bin <- unique(selected$wealth.bin)
    
    combined.list[[counter]] <- combined.fit
    
    counter <- counter + 1
    
}

combined.fit <- do.call(rbind, combined.list)

wealth.repair <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y=survival, 
                                       color = wealth.bin),
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper, fill=wealth.bin),
                           size=0.75, alpha = 0.5) +
	  facet_grid(.~sex) +
    theme_cowplot()+
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position=c(.77, 0.65),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = unit(c(0.05, 0, 0, 0.03), "cm")) + 
	  scale_fill_manual(values = c("#2F8DFA","#1FD0BF","#EB648B")) +
	  scale_color_manual(values = c("#2F8DFA","#1FD0BF","#EB648B")) +
    labs(x="Years since damage occured", y="Probability of remaining damaged", color = "Household wealth", fill ="Household wealth") +
    ggtitle('e) ELSA humans (Phelps et al. 2020)') +
    geom_text(
        data    = pvals,
        mapping = aes(x = 0, y = 0.05, label = pval), size=2.25,
        hjust   = 0,
        vjust   = 0)


data <- read.csv('../datasets/human_damage_deficits.csv', header = TRUE, sep = ",")
data$wealth <- log(data$wealth + mean.wealth)


data <- data[data$baseline.age >= 50,]
data <- data[data$baseline.age < 90,]

bins <- seq(50, 90, by=10)
data <- 
data %>% 
    group_by(id) %>%
    mutate(baseline.bin = cut(baseline.age, bins, right=FALSE)) %>%
    as.data.frame()


deficits <- c('d0'='Difficulty walking 100 yads','d1'='Difficulty sitting for 2 hours','d2'='Difficulty getting up from a chair',
              'd3'='Difficulty climbing several flights of stairs', 'd4'='Difficulty climbing one flight of stairs',
	      'd5'='Difficulty stooping, kneeling, or crouching',
	      'd6'='Difficulty reaching/extending arms', 'd7'='Difficulty pulling/pushing objects',
	      'd8'='Difficulty lifting/carrying 10 lbs weights', 'd9'='Difficulty picking up a coin',
	      'd10'='Difficulty dressing','d11'='Difficulty walking across a room','d12'='Difficulty bathing','d13'='Difficulty Eating',
	      'd14'='Difficulty getting in/out of bed','d15'='Difficulty using the toilet','d16'='Difficulty using a map',
	      'd17'='Difficulty preparing a hot meal','d18'='Difficulty shopping','d19'='Difficulty using telephone',
	      'd20'='Difficulty taking medications','d21'='Difficulty doing work around the house','d22'='Difficulty managing money')



names.deficits <- c()
for(i in 1:length(deficits)) {
      names.deficits <- append(names.deficits, paste('d', i-1, sep=""))
}
data$deficit <- as.factor(data$deficit)
levels(data$deficit) <- deficits

data[data$sex==0,]$sex <- 'M'
data[data$sex==1,]$sex <- 'F'
data$sex <- as.factor(data$sex)
data$diff <- data$damage.age - data$repair.age

data <- na.omit(data)

data$wealth.bin <- 
  as.factor(cut(data$wealth, breaks=c(-Inf, quantile(elsa$wealth, c(1/3, 2/3)), Inf)))


levels(data$wealth.bin) <- c("Bottom tercile","Middle tercile","Upper tercile")
levels(data$sex) <- c("Female", "Male")

rstan_options(auto_write = TRUE)

data$group <- NA
counter <- 1
for(w in levels(data$wealth.bin)) {
    for(s in levels(data$sex)) {
	    data[(data$sex == s) & (data$wealth.bin == w), 'group'] <- counter
	    counter <- counter + 1
    }	 
}

data.censored <- data[data$status == 0,]
data.dead <- data[data$status == 1,]


input <- data[,c('diff.left', 'diff.right', 'wealth.bin', 'sex')]
input = filter(input, !((data$diff.left < 0.0000001) & (is.infinite(data$diff.right))))
input$wealth.bin <- as.numeric(input$wealth.bin)-1


male.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Male',c('diff.left', 'diff.right', 'wealth.bin')],k=3)$p*1000)/1000)
female.pval <- paste0("p=",trunc(gLRT3(input[input$sex=='Female',c('diff.left', 'diff.right', 'wealth.bin')],k=3)$p*1000)/1000)

pvals <- data.frame(pval = c(male.pval, female.pval), sex=c('Male', 'Female'))

pvals

counter <- 1
combined.list <- list()
for (g in unique(data$group)) {

    selected = data[data$group == g,]
    selected.censored <- data.censored[data.censored$group == g,]
    selected.dead <- data.dead[data.dead$group == g,]
    
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
                data = stan.data, chains=4, 
                iter = 5000, cores = 4, warmup=2000, control=list(adapt_delta=0.9, max_treedepth=15), refresh=-1)

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
    combined.fit$wealth.bin <- unique(selected$wealth.bin)
    
    combined.list[[counter]] <- combined.fit
    
    counter <- counter + 1
    
}

combined.fit <- do.call(rbind, combined.list)

wealth.damage <- ggplot() +
	  geom_step(data = combined.fit,
                           mapping=aes(x = time, y=survival, 
                                       color = wealth.bin),
                           size=0.75, alpha = 1) +
	  geom_stepribbon(data = combined.fit,
                           mapping=aes(x = time, y=survival, ymin=.lower, ymax=.upper, fill=wealth.bin),
                           size=0.75, alpha = 0.5) +
	  facet_grid(.~sex) +
    theme_cowplot()+
    theme(
        strip.background = element_blank(),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x=element_text(size=7),
        strip.text.y=element_text(size=7),
        plot.title=element_text(size=8),
        legend.position="none", #c(.77, 0.75)
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = unit(c(0.05, 0, 0, 0.03), "cm")) + 
	  scale_fill_manual(values = c("#2F8DFA","#1FD0BF","#EB648B")) +
	  scale_color_manual(values = c("#2F8DFA","#1FD0BF","#EB648B")) +
    labs(x="Years since repair occured", y="Probability of remaining undamaged", color = "Household wealth", fill ="Household wealth") + ggtitle('') +
    geom_text(
        data    = pvals,
        mapping = aes(x = 0, y = 0.05, label = pval), size=2.25,
        hjust   = 0,
        vjust   = 0)


img <- grid.arrange(wealth.repair, wealth.damage, ncol=2, nrow=1)

ggsave("../plots/Human_wealth_timescale.pdf", img, , width=8, height=1.65)
