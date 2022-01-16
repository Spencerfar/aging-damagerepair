library(dplyr)
library(caret)
library(splines2)
library(cmdstanr)
set_cmdstan_path("~/Downloads/cmdstan-2.26.0/")

# create folder for fits
if (!dir.exists('fits')) {dir.create('fits')}

elsa <- read.csv('datasets/human_data.csv', header = TRUE, sep = ",")

elsa$id <- as.character(elsa$id)
elsa$time <- elsa$age - elsa$baseline.age
elsa <- elsa %>% group_by(id) %>% filter( (baseline.age >= 50)  & (baseline.age < 90))
elsa <- na.omit(elsa)
elsa$f <- elsa$n/elsa$N
elsa$wealth <- log(elsa$wealth + mean(elsa$wealth))


elsa.X <- elsa[,c('sex', 'baseline.age', 'time', 'f', 'wealth')]
elsa.X.slope <- elsa[,c('sex', 'baseline.age', 'time', 'f', 'wealth')]
elsa.Z <- elsa[,c('time', 'f')]
elsa.Z$intercept <- 1
elsa.Z <- elsa.Z[,c('intercept', 'time')]


# create 2nd and 3rd order interactions
elsa.X$sexXtime <- elsa.X$sex * elsa.X$time
elsa.X$baselineXtime <- elsa.X$baseline.age * elsa.X$time
elsa.X$sexXbaseline <- elsa.X$sex * elsa.X$baseline.age
elsa.X$sexXbaselineXtime <- elsa.X$sex * elsa.X$baseline.age * elsa.X$time
elsa.X$timeXwealth<- elsa.X$time * elsa.X$wealth
elsa.X$sexXf <- elsa.X$sex * elsa.X$f
elsa.X$wealthXf <- elsa.X$wealth * elsa.X$f
elsa.X$baselineXf <- elsa.X$baseline.age * elsa.X$f
elsa.X$sexXwealth <- elsa.X$sex * elsa.X$wealth
elsa.X$sexXtimeXwealth <- elsa.X$sex * elsa.X$time * elsa.X$wealth

normELSA <- preProcess(elsa.X)
elsa.X <- predict(normELSA, elsa.X)

num_knots <- 3 # number of knots for fitting
spline_degree <- 3

baseline.knots <- unname(quantile(elsa.X$baseline.age, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
baseline.boundary.knots <- c(min(elsa.X$baseline.age), max(elsa.X$baseline.age))
baseline.splines <- bSpline(elsa.X$baseline.age,  knots = baseline.knots, degree = spline_degree, Boundary.knots = baseline.boundary.knots, intercept = FALSE)

wealth.knots <- unname(quantile(elsa.X$wealth, probs=seq(from=0.1, to=0.9, length.out = num_knots)))
wealth.boundary.knots <- c(min(elsa.X$wealth), max(elsa.X$wealth))
wealth.splines <- bSpline(elsa.X$wealth,  knots = wealth.knots, degree = spline_degree, Boundary.knots = wealth.boundary.knots, intercept = FALSE)
wealth.splines.deriv <- bSpline(elsa.X$wealth,  knots = wealth.knots, degree = spline_degree, Boundary.knots = wealth.boundary.knots, intercept = FALSE, derivs=1)


tensor_spline <- array(rep(NaN, dim(wealth.splines)[1]*dim(wealth.splines)[2]*dim(baseline.splines)[2]), c(dim(wealth.splines)[1], dim(wealth.splines)[2], dim(baseline.splines)[2]))
tensor_spline.deriv <- array(rep(NaN, dim(wealth.splines)[1]*dim(wealth.splines)[2]*dim(baseline.splines)[2]), c(dim(wealth.splines)[1], dim(wealth.splines)[2], dim(baseline.splines)[2]))

for(i in 1:dim(wealth.splines)[1]) {
    tensor_spline[i,,] <- outer(wealth.splines[i,], baseline.splines[i,])
    tensor_spline.deriv[i,,] <- outer(wealth.splines.deriv[i,], baseline.splines[i,])
}


elsa.repair <- elsa$repair
elsa.damage <- elsa$damage
elsa.numdef <- elsa[,'n']
elsa.offsets <- elsa[,c('n', 'N', 'delta.t')]
elsa.index <- as.integer(as.numeric(factor(elsa$id)))


elsa.X <- elsa.X[,c('sex', 'time', 'sexXtime', 'wealth', 'baseline.age', 'sexXwealth', 'sexXbaseline', 'baselineXtime', 'timeXwealth', 'sexXbaselineXtime', 'sexXtimeXwealth', 'f', 'sexXf', 'wealthXf', 'baselineXf')]


data <- list(m = length(unique(elsa$id)),
             n = length(elsa$repair),
             p = length(names(elsa.X)),
             X = elsa.X,
             index = elsa.index,
             offsets = elsa.offsets,
             repair = elsa.repair,
             damage = elsa.damage,
             num_basis = dim(baseline.splines)[2],
             spline_2d = tensor_spline, spline_2d_deriv = tensor_spline.deriv)




file <- file.path('models/long_human.stan')
model <- cmdstan_model(file, cpp_options = list(stan_threads = TRUE))
model$compile(force_recompile = TRUE, cpp_options=list(stan_threads=TRUE))

model$check_syntax(pedantic = TRUE) 

fit <- model$sample(
  data = data,
  seed = 2,
  chains = 1,
  parallel_chains = 1,
  refresh = 5,
  iter_warmup = 1000,
  iter_sampling = 3000,
  threads_per_chain=40
)

fit$save_object(file = "fits/long_human_fit_chain2.RDS")
fit$save_output_files(dir = "fits/")

#fit$save_object(file = "~/resilience/long_human_fit_100.RDS")
#fit$save_output_files(dir = "~/resilience/")
