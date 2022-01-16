library(rstan)
library(dplyr)
library(survival)
library(splines2)
library(statmod)
library(abind)
library(caret)

# create folder for fits
if (!dir.exists('fits')) {dir.create('fits')}

### setup data 
mice <- read.csv('datasets/schultz_data.csv', header = TRUE, sep = ",")
mice.surv <- read.csv('datasets/schultz_surv_data.csv', header = TRUE, sep = ",")


mice$sex <- as.factor(mice$sex)
mice$real.age <- mice$time + mice$baseline.age

mice$prepair <- (mice$repair/mice$n)/mice$delta.t
mice$pdamage <- (mice$damage/(mice$N - mice$n))/mice$delta.t

mice <- na.omit(mice)
mice.surv <- na.omit(mice.surv)

mice$event.time <- mice$death.age - mice$baseline.age
mice.surv$event.time <- mice.surv$death.age - mice.surv$baseline.age
mice.surv <- mice.surv[mice.surv$mouse %in% unique(mice$mouse),]

mice.surv2 <- tmerge(mice.surv, mice.surv, id=mouse, endpt = event(event.time, status))


mice.surv.long <- tmerge(mice.surv2, mice, id=mouse, f = tdc(time, f), n = tdc(time, n), 
                        real.age = tdc(time, real.age))
mice.surv.long$tstart <- mice.surv.long$tstart
mice.surv.long$tstop <- mice.surv.long$tstop
mice.surv.long$event.time <- mice.surv.long$event.time


mice.surv.long$sex <- as.factor(mice.surv.long$sex)


### create stan input

mice.X <- mice[,c('time', 'f', 'baseline.age')]
mice.Z <- mice[,c('time', 'f')]
mice.Z$intercept <- 1
mice.Z <- mice.Z[,c('intercept', 'time')]

# scale
age.mean <- mean(mice.X$baseline.age)
age.sd <- sd(mice.X$baseline.age)
time.mean <- mean(mice.X$time)
time.sd <- sd(mice.X$time)
f.mean <- mean(mice.X$f)
f.sd <- sd(mice.X$f)
std_scale <- function(x, m, s) return ((x-m)/s)
inv_scale <- function(x, m, s) return (s*x + m)

# standardize
normMice <- preProcess(mice.X)
mice.X <- predict(normMice, mice.X)
mice.Z$time <- (mice.Z$time - (normMice$mean)['time'])/(normMice$std)['time']

mice.repair <- mice[,'repair']
mice.damage <- mice[,'damage']
mice.numdef <- mice[,'n']
mice.offsets <- mice[,c('n', 'N', 'delta.t')]
mice.index <- as.integer(as.numeric(factor(mice$mouse)))


mice.surv.X <- mice.surv.long[,c('f', 
                                 'baseline.age', 'tstart', 'tstop')]


mice.start.age <- mice.surv.long$tstart
mice.death.age <- mice.surv.long$tstop
mice.status <- as.numeric(mice.surv.long$endpt)
dt <- mice.surv.X$tstop - mice.surv.X$tstart
mice.surv.X <- mice.surv.X[,c('f', 'baseline.age')]


normSurvMice <- preProcess(mice.surv.X)
mice.surv.X.normed <- predict(normSurvMice, mice.surv.X)
mice.surv.X <- mice.surv.X.normed


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
             N = 116,
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
             msplines = msMat)

rstan_options(auto_write = TRUE)
fit <- stan(file = 'models/joint_long.stan', data = data, chains=4, iter = 10000, cores = 4, verbose=TRUE, warmup=4000)
saveRDS(fit, 'fits/mouse_3.rds')
