library(rstan)
library(dplyr)
library(survival)
library(splines2)
library(statmod)
library(abind)
library(caret)

### setup data 
mice <- read.csv('datasets/exercise_data.csv', header = TRUE, sep = ",")
mice.surv <- read.csv('datasets/exercise_surv_data.csv', header = TRUE, sep = ",")


mice$sex <- as.factor(mice$sex)
mice$exercise <- as.factor(mice$exercise)
mice$real.age <- mice$time + mice$baseline.age

mice$prepair <- (mice$repair/mice$n)/mice$delta.t
mice$pdamage <- (mice$damage/(mice$N - mice$n))/mice$delta.t

mice <- na.omit(mice)
mice.surv <- na.omit(mice.surv)

mice$event.time <- mice$death.age
mice.surv$event.time <- mice.surv$death.age
mice.surv <- mice.surv[mice.surv$mouse %in% unique(mice$mouse),]

mice.surv2 <- tmerge(mice.surv, mice.surv, id=mouse, endpt = event(event.time, status))

mice.surv.long <- tmerge(mice.surv2, mice, id=mouse, f = tdc(time, f), n = tdc(time, n), 
                        real.age = tdc(time, real.age), time = tdc(time, time))

mice.surv.long$sex <- as.factor(mice.surv.long$sex)
mice.surv.long$exercise <- as.factor(mice.surv.long$exercise)

### create stan input

mice.surv.X <- mice.surv.long[,c('sex', 'exercise', 'time', 'f', 
                                 'baseline.age', 'tstart', 'tstop')]

mice.X <- mice[,c('sex', 'exercise', 'time', 'f', 'baseline.age')]
mice.Z <- mice[,c('time', 'f')]
mice.Z$intercept <- 1
mice.Z <- mice.Z[,c('intercept', 'time')]

# scale
time.mean <- mean(mice.X$time)
time.sd <- sd(mice.X$time)
age.mean <- mean(mice.surv.X$baseline.age)
age.sd <- sd(mice.surv.X$baseline.age)
f.mean <- mean(mice.X$f)
f.sd <- sd(mice.X$f)
std_scale <- function(x, m, s) return ((x-m)/s)
inv_scale <- function(x, m, s) return (s*x + m)


# fix factors -> integers
mice.X$sex <- as.numeric(mice.X$sex)-1
mice.X$exercise <- as.numeric(mice.X$exercise)-1

mice.X_T <- mice.X[,c('sex', 'exercise', 'time', 'f', 'baseline.age')]
mice.X_T$time <- mice.surv.long$tstop

# create interactions
mice.X$sexXtime <- mice.X$sex * mice.X$time
mice.X$exerciseXtime <- mice.X$exercise * mice.X$time
mice.X$sexXexercise <- mice.X$sex * mice.X$exercise
mice.X$sexXexerciseXtime <- mice.X$sex * mice.X$exercise * mice.X$time
mice.X$baselineXsex <- mice.X$baseline.age * mice.X$sex

mice.X$sexXf <- mice.X$sex * mice.X$f
mice.X$exerciseXf <- mice.X$exercise * mice.X$f
mice.X$sexXexerciseXf <- mice.X$sex * mice.X$exercise * mice.X$f

mice.X_T$sexXtime <- mice.X_T$sex * mice.X_T$time
mice.X_T$exerciseXtime <- mice.X_T$exercise * mice.X_T$time
mice.X_T$sexXexercise <- mice.X_T$sex * mice.X_T$exercise
mice.X_T$sexXexerciseXtime <- mice.X_T$sex * mice.X_T$exercise * mice.X_T$time
mice.X_T$baselineXsex <- mice.X_T$baseline.age * mice.X_T$sex

mice.X_T$sexXf <- mice.X_T$sex * mice.X_T$f
mice.X_T$exerciseXf <- mice.X_T$exercise * mice.X_T$f
mice.X_T$sexXexerciseXf <- mice.X_T$sex * mice.X_T$exercise * mice.X_T$f

# standardize
normMice <- preProcess(mice.X)
mice.X <- predict(normMice, mice.X)
mice.Z$time <- (mice.Z$time - (normMice$mean)['time'])/(normMice$std)['time']

X.means <- as.vector(normMice$mean)
X.stds <- as.vector(normMice$std)

mice.X_T <- predict(normMice, mice.X_T)

mice.repair <- mice[,'repair']
mice.damage <- mice[,'damage']
mice.numdef <- mice[,'n']
mice.offsets <- mice[,c('n', 'N', 'delta.t')]
mice.index <- as.integer(as.numeric(factor(mice$mouse)))



mice.surv.X$sex <- as.numeric(mice.surv.X$sex)-1
mice.surv.X$exercise <- as.numeric(mice.surv.X$exercise)-1

mice.start.age <- mice.surv.long$tstart
mice.death.age <- mice.surv.long$tstop
mice.status <- as.numeric(mice.surv.long$endpt)
dt <- mice.surv.X$tstop - mice.surv.X$tstart
mice.surv.X <- mice.surv.X[,c('sex', 'exercise', 'f', 'baseline.age')]

mice.surv.X$sexXexercise <- mice.surv.X$sex * mice.surv.X$exercise

normSurvMice <- preProcess(mice.surv.X)
mice.surv.X <- predict(normSurvMice, mice.surv.X)

# create splines
mice.start.age <- std_scale(mice.start.age, time.mean, time.sd)
mice.death.age <- std_scale(mice.death.age, time.mean, time.sd)

num_knots <- 15 # number of knots for fitting
spline_degree <- 3

knots <- unname(quantile(mice.death.age, probs=seq(from=0.05, to=0.95, length.out = num_knots)))
boundary.knots <- c(min(mice.start.age), max(mice.death.age))

msMat <- mSpline(mice.death.age,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = FALSE)

num_basis <- dim(msMat)[2]

quadMat <- c()

quad <- gauss.quad(5, kind="legendre",alpha=0,beta=0)
quad.weights <- quad$weights
quad.nodes <- quad$nodes

for (i in 1:dim(msMat)[1]) {
    delta = (mice.death.age[i] - mice.start.age[i])/2
    average = (mice.death.age[i] + mice.start.age[i])/2.0
    x = delta * quad.nodes + average
    quadMat[[i]] <- mSpline(x,  knots = knots, degree = spline_degree, Boundary.knots = boundary.knots, intercept = FALSE)
}

quadMat <- abind( quadMat, along=0 )


# compute size of each mouse time-series
index.size <- rep(0, length(unique(mice$mouse)))
for(i in 1:length(unique(mice$mouse))) {
    index.size[i] <- length(mice.index[mice.index == i])
}






data <- list(m = length(unique(mice.index)),
             n = length(mice$repair),
             p = length(names(mice.X)),
             ps = length(names(mice.surv.X)),
             q = length(names(mice.Z)),
             N = 124,
             X = mice.X,
             surv_X = mice.surv.X,
             Z = mice.Z,
             index = mice.index,
             index_size = index.size,
             offsets = mice.offsets,
             repair = mice.repair,
             damage = mice.damage,
             Tstart = mice.start.age,
             T = mice.death.age,
             status = mice.status,
             num_basis = num_basis,
             num_quad = 5,
             quad_splines = quadMat,
             quad_nodes = quad.nodes,
             quad_weights = quad.weights,
             msplines = msMat,
             means = X.means,
             stds = X.stds,
             X_T = mice.X_T)

rstan_options(auto_write = TRUE)
fit <- stan(file = 'models/joint_long_interventions.stan', data = data, chains=4, iter=3000, cores = 4, verbose=TRUE, warmup=1500, control=list(adapt_delta=0.95, max_treedepth=20))
saveRDS(fit, 'fits/mouse_2.rds')
