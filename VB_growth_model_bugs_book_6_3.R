library(R2OpenBUGS)
library(R2jags)

# Von Bertalanffy growth models in section 6.3.1 of Lunn et al. (2013) "The BUGS Book"
# note this contains 4 slightly different versions of the VB model
# 1 and 2 are the same parameterization, with different priors
# 3 and 4 have use the same alternate parameterization, but different priors 
sink("dugongVB.txt")
cat("
model{
  for(j in 1:N) {
    for(i in 1:2) {
  y[1,j] ~ dnorm(mu[1,j], tau[i])
  }
  mu[1,j] <- Linf[1] - (Linf[1] - L0[1]) * exp(-K[1]*x[1,j])
  mu[2,j] <- Linf[2] - (Linf[2] - L0[2]) * exp(-K[2]*x[2,j])
  mu[3,j] <- alpha[3] - beta[3] * pow(gamma[3], x[3,j])
  mu[4,j] <- alpha[4] - beta[4] * pow(gamma[4], x[4,j])
  }
  L0[1] ~ dunif(0, 100)
  L0[2] ~ dnorm(0, 0.0001)I(0, Linf[2]) # truncated normal distribution
  Linf[1] <- L0[1] + beta[1]
  Linf[2] ~ dnorm(0, 0.0001)I(L0[2], ) # right truncated
  K[1] ~ dunif(0, 100)
  K[2] ~ dunif(0, 100)
  for(i in 1:2) {alpha[i] <- Linf[i]}
  for(i in 3:4) {alpha[i] ~ dunif(0, 100)}
  beta[1] ~ dunif(0, 100)
  beta[2] <- Linf[2] - L0[2]
  for(i in 3:4) {beta[i] ~ dunif(0, 100)}
  for(i in 1:2) {gamma[i] <- exp(-K[i])}
  gamma[3] ~ dunif(0, 1)
  gamma[4] ~ dgamma(0.001, 0.001)I(0, 1)
  for (i in 1:4){
    tau[i] <- 1/sigma2[i]
    log(sigma2[i]) <- 2*log.sigma[i]
    log.sigma[i] ~ dunif(-10, 10)
  }
}
", fill=TRUE)
sink()

# Data
dugong.data <- list(x = structure(
  .Data = c(1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0,
            8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
            13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5,
            1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0,
            8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
            13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5,
            1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0,
            8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
            13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5,
            1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0,
            8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
            13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5),
  .Dim = c(4, 27)),
  y = structure(
    .Data = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
              2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
              2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57,
              1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
              2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
              2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57,
              1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
              2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
              2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57,
              1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
              2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
              2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57),
    .Dim = c(4, 27)), N = 27)

params <- c("alpha", "beta", "gamma", "sigma2")

inits <- function() list(Linf[2] = 3)

ni <- 10000
nc <- 3
nb <- 2000
nt <- 4

out.dugong <- bugs(data=dugong.data, inits=NULL, parameters.to.save=params, 
                   model.file="dugongVB.txt", n.iter=ni, n.burnin = nb, n.chains = nc,
                   n.thin = nt, debug=TRUE) #, OpenBUGS.pgm = bugs.dir, working.directory=getwd())

####################
### Model 1 only ###
####################

sink("dugongVB1.txt")
cat("
    model{
    for(j in 1:N) {
    
    y[j] ~ dnorm(mu[j], tau)
    
    mu[j] <- Linf - (Linf - L0) * exp(-K*x[j])

    }
    L0 ~ dunif(0, 100)
    Linf <- L0 + beta
    K ~ dunif(0, 100)
    alpha <- Linf
    beta ~ dunif(0, 100)
    gamma <- exp(-K)
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-10, 10)
    
    }
    ", fill=TRUE)
sink()

# Data
dugong.data.m1 <- list(x = dugong.data$x[1,], y = dugong.data$y[1,], N = 27)


params <- c("alpha", "beta", "gamma", "sigma2")

# model runs poorly if you do not submit initial values
inits <- function() list(L0=2, beta=0.8, K=0.1, log.sigma=-2)

ni <- 60000
nc <- 3
nb <- 10000
nt <- 5

out.dugong.m1 <- bugs(data=dugong.data.m1, inits=inits, 
                      parameters.to.save=params, model.file="dugongVB1.txt", 
                      n.iter=ni, n.burnin = nb, n.chains = nc, n.thin = nt, 
                      debug=TRUE)

out.dugong.m1.jags <- jags(data=dugong.data.m1, inits=inits, 
                           parameters.to.save=params, model.file="dugongVB1.txt",
                           n.iter=ni, n.burnin = nb, n.chains = nc, n.thin = nt,
                           working.directory = getwd())

print(out.dugong.m1, digits=4)
print(out.dugong.m1.jags, digits=4)

#######################################################
#### VB model in terms of size and time increments ####
#######################################################

# Based on Eq. 18 in Hart and Chute (2009) ICES J. Mar. Sci.

sink("ggsVB2.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau) # y[j] is size increment grown in row j
    mu[j] <- (Linf - Lt[j]) * (1 - exp(-K*x[j])) # x is time increment of row j
  } # j
    
    #Lt ~ dunif(0, 1000) # Lt is data for each individual, not a parameter to be estimated
    #Linf <- L0 + beta
    #Linf ~ dunif(0, 2000)
    
    Linf ~ dnorm(1000, 0.0001)
    K ~ dunif(0, 100)
    
    #alpha <- Linf
    #beta ~ dunif(0, 2000)
    #gamma <- exp(-K)
    
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-2000, 2000)
    
}
", fill=TRUE)
sink()

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, 
                 y = ggs.selfjoin.nodupes.30.pos$SVL.diff,   
                 Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , 
                 N = nrow(ggs.selfjoin.nodupes.30.pos))

ggs.data <- list(x = filter2.ggs.pos$Date.diff, y = filter2.ggs.pos$SVL.diff, 
                 Lt=filter2.ggs.pos$SVL.x , N = nrow(filter2.ggs.pos))


# model runs poorly if you do not submit initial values
inits <- function() list(Linf=1000, K=0.0009, log.sigma=3)

params <- c("Linf", "K", "sigma2")

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

out.ggs.m2 <- bugs(data=ggs.data, inits=inits, parameters.to.save=params, 
                   model.file="ggsVB2.txt", n.iter=ni, n.burnin = nb, 
                   n.chains = nc, n.thin = nt, debug=TRUE)

# run in jags
ptm <- proc.time()
out.ggs.m2.jags <- jags(data = ggs.data, inits=inits, parameters.to.save = params,
                        model.file="ggsVB2.txt", n.chains = nc, n.iter = ni, 
                        n.burnin = nb, n.thin = nt, working.directory = getwd())
proc.time() - ptm

# run jags in parallel
ptm <- proc.time()
out.ggs.m2.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, parameters.to.save = params,
                        model.file="ggsVB2.txt", n.chains = nc, n.iter = 20000, 
                        n.burnin = 4000, n.thin = 1, n.cluster = nc,
                        working.directory = getwd())
proc.time() - ptm


print(out.ggs.m2.jags.parallel, digits=6)
hist(out.ggs.m2.jags$BUGSoutput$sims.list$K)

# plot von Bertalanffy curve using jags output
age.vec <- seq(0, 3650, by = 365)
traj.jags <- vonBert(K = 0.000618, w = 1 - (200/1000), Linf = 1000, A = age.vec)
plot(age.vec/365, traj.jags, ylab = "SVL (mm)", xlab = "Age (years)", ylim=c(0,1000))
# size trajectory does not make sense, value of K too large

# compare to estimation of a and b using linear regression
traj.lm <- vonBert(K = -b, Linf=a/-b, w=1-200/(a/-b), A = seq(0,3650,by=365))
points(seq(0, 3650, by = 365)/365, traj.lm, pch=2) #ylab = "SVL (mm)", xlab = "Age (years)")

## now with individual level random effect for Asymptotic size (Linf)

sink("ggsVB2RE.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau)                        # y[j] is size increment grown in row j
    mu[j] <- (Linf[j] - Lt[j]) * (1 - exp(-K*x[j])) # x is time increment of row j
    Linf[j] <- alpha + lambda[j]                    # Linf varies among individuals
    
  } # j
    
    # Priors
    # loop over all rows - some individuals duplicated though
    for(j in 1:N) {                                 
     lambda[j] ~ dnorm(0, inv.omega.lambda.squared)
    } # j
    alpha ~ dnorm(1000, 0.0001)
    K ~ dunif(0, 100)
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-2000, 2000)
    
    # Hyperprior
    sigma.lambda ~ dunif(0, 100)
    sigma2.lambda <- sigma.lambda * sigma.lambda
    inv.omega.lambda.squared <- 1 / pow(sigma.lambda, 2)
    
}
", fill=TRUE)
sink()

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.diff, 
                 Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , N = nrow(ggs.selfjoin.nodupes.30.pos))

ggs.data <- list(x = filter2.ggs.pos$Date.diff, y = filter2.ggs.pos$SVL.diff, 
                 Lt=filter2.ggs.pos$SVL.x , N = nrow(filter2.ggs.pos))


# model runs poorly if you do not submit initial values
inits <- function() list(alpha=1000, K=0.0009, log.sigma=3)

params <- c("alpha", "K", "sigma2", "sigma2.lambda")

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

ptm <- proc.time()
out.ggs.m2RE.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, parameters.to.save = params,
                                          model.file="ggsVB2RE.txt", n.chains = nc, n.iter = 20000, 
                                          n.burnin = 4000, n.thin = 1, n.cluster = nc,
                                          working.directory = getwd())
proc.time() - ptm

print(out.ggs.m2RE.jags.parallel, digits=6)

age.vec <- seq(0, 3650, by = 365)
traj.jags2 <- vonBert(K = 0.000762, w = 1 - (200/972.7), Linf = 972.7, A = age.vec)
plot(age.vec/365, traj.jags, ylab = "SVL (mm)", xlab = "Age (years)", ylim=c(0,1000))


## add sex effect on asymptotic size (Linf) to model

# create dummy variable for male status
ggs.selfjoin.nodupes.30.pos$male <- ifelse(ggs.selfjoin.nodupes.30.pos$Sex.y == "M", 1, 0)

sink("ggsVB2REsex.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau)                        # y[j] is size increment grown in row j
    mu[j] <- (Linf[j] - Lt[j]) * (1 - exp(-K*x[j])) # x is time increment of row j
    Linf[j] <- alpha + beta0 * sex[j] + lambda[j]                    # Linf varies among individuals
    
  } # j
    
    # Priors
    # loop over all rows - some individuals duplicated though
    for(j in 1:N) {                                 
     lambda[j] ~ dnorm(0, inv.omega.lambda.squared)
    } # j
    alpha ~ dnorm(1000, 0.0001)
    beta0 ~ dnorm(0, 0.0001)
    K ~ dunif(0, 100)
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-2000, 2000)
    
    # Hyperprior
    sigma.lambda ~ dunif(0, 100)
    sigma2.L <- sigma.lambda * sigma.lambda
    inv.omega.lambda.squared <- 1 / pow(sigma.lambda, 2)
    
}
", fill=TRUE)
sink()

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.diff, 
                 Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , N = nrow(ggs.selfjoin.nodupes.30.pos), 
                 sex = ggs.selfjoin.nodupes.30.pos$male)

ggs.data <- list(x = filter2.ggs.pos$Date.diff, y = filter2.ggs.pos$SVL.diff, 
                 Lt=filter2.ggs.pos$SVL.x , N = nrow(filter2.ggs.pos), 
                 sex = filter2.ggs.pos$male)

# model runs poorly if you do not submit initial values
inits <- function() list(alpha=1000, K=0.0009, log.sigma=3, beta0=rnorm(1))

params <- c("alpha", "K", "sigma2", "sigma2.L", "beta0")

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

ptm <- proc.time()
out.ggs.m2REsex.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, parameters.to.save = params,
                                            model.file="ggsVB2REsex.txt", n.chains = nc, n.iter = 20000, 
                                            n.burnin = 4000, n.thin = 1, n.cluster = nc,
                                            working.directory = getwd())
proc.time() - ptm

print(out.ggs.m2REsex.jags.parallel, digits=6)

hist(out.ggs.m2RE.jags.parallel$BUGSoutput$sims.list$beta0)

## now add sex effect on K, the growth coefficient

sink("ggsVB2REsex2.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau)                              # y[j] is size increment grown in row j
    mu[j] <- (Linf[j] - Lt[j]) * (1 - exp(-K[j]*x[j]))    # x is time increment of row j
    Linf[j] <- alpha + beta0 * sex[j] + lambda[j]         # Linf varies among individuals
    K[j] <- alphaK + beta1 * sex[j]                       # K now depends on sex
    } # j
    
    # Priors
    # loop over all rows - some individuals duplicated though
    for(j in 1:N) {                                 
    lambda[j] ~ dnorm(0, inv.omega.lambda.squared)
    } # j
    alpha ~ dnorm(1000, 0.0001)
    beta0 ~ dnorm(0, 0.0001)
    alphaK ~ dnorm(0, 0.0001)
    beta1 ~ dnorm(0, 0.0001)
    #K ~ dunif(0, 100)
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-2000, 2000)
    
    # Hyperprior
    sigma.lambda ~ dunif(0, 100)
    sigma2.L <- sigma.lambda * sigma.lambda
    inv.omega.lambda.squared <- 1 / pow(sigma.lambda, 2)
    
    }
    ", fill=TRUE)
sink()

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.diff, 
                 Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , N = nrow(ggs.selfjoin.nodupes.30.pos), 
                 sex = ggs.selfjoin.nodupes.30.pos$male)

# model runs poorly if you do not submit initial values
inits <- function() list(alpha = 1000, alphaK = 0.0009, log.sigma = 3,
                         beta0 = rnorm(1), beta1 = rnorm(1))

params <- c("alpha", "alphaK", "sigma2", "sigma2.L", "beta0", "beta1")

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

ptm <- proc.time()
out.ggs.m2REsex2.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, parameters.to.save = params,
                                               model.file="ggsVB2REsex2.txt", n.chains = nc, n.iter = 20000, 
                                               n.burnin = 4000, n.thin = 1, n.cluster = nc,
                                               working.directory = getwd())
proc.time() - ptm

print(out.ggs.m2REsex2.jags.parallel, digits=6)
traceplot(out.ggs.m2REsex2.jags.parallel)

# convert model output to an MCMC object
out.ggs.mcmc <- as.mcmc(out.ggs.m2REsex2.jags.parallel)
summary(out.ggs.mcmc)

hist(out.ggs.m2REsex2.jags.parallel$BUGSoutput$sims.list$beta1)

###############################
### add random effect for K ###
###############################
sink("ggsVB2RE2.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau)                            # y[j] is size increment grown in row j
    mu[j] <- (Linf[j] - Lt[j]) * (1 - exp(-K[j]*x[j]))  # x is time increment of row j
    Linf[j] <- alphaL + lambda[j]                       # Linf varies among individuals
    K[j] <- alphaK + gamma[j]
    
  } # j
    
    # Priors
    # loop over all rows - some individuals duplicated though
    for(j in 1:N) {                                 
     lambda[j] ~ dnorm(0, inv.lambda.squared)
     gamma[j] ~ dnorm(0, inv.gamma.squared)
    } # j
    alphaL ~ dnorm(1000, 0.0001)
    alphaK ~ dnorm(0, 0.0001)
    
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    sigma <- pow(sigma2, 0.5)
    log.sigma ~ dunif(-2000, 2000)
    
    # Hyperpriors
    sigma.L ~ dunif(0, 100)
    #sigma2.L <- sigma.L* sigma.L
    inv.lambda.squared <- 1 / pow(sigma.L, 2)

    sigma.K ~ dunif(0, 100)
    #sigma2.K <- sigma.K * sigma.K
    inv.gamma.squared <- 1 / pow(sigma.K, 2)
    
}
", fill=TRUE)
sink()

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.diff, 
                 Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , N = nrow(ggs.selfjoin.nodupes.30.pos))

ggs.data <- list(x = filter2.ggs.pos$Date.diff, y = filter2.ggs.pos$SVL.diff, 
                 Lt=filter2.ggs.pos$SVL.x , N = nrow(filter2.ggs.pos))

# model runs poorly if you do not submit initial values
inits <- function() list(alphaL=1000, alphaK=0.0009, log.sigma=3)

params <- c("alphaL", "alphaK", "sigma", "sigma.L", "sigma.K")

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

ptm <- proc.time()
out.ggs.m2RE2.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, parameters.to.save = params,
                                            model.file="ggsVB2RE2.txt", n.chains = nc, n.iter = 20000, 
                                            n.burnin = 4000, n.thin = 1, n.cluster = nc,
                                            working.directory = getwd())
proc.time() - ptm

print(out.ggs.m2RE2.jags.parallel, digits = 6)
traceplot(out.ggs.m2RE2.jags.parallel)


## add fixed effect of sex on Linf and K

sink("ggsVB2RE2sex.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau)                            # y[j] is size increment grown in row j
    mu[j] <- (Linf[j] - Lt[j]) * (1 - exp(-K[j]*x[j]))  # x is time increment of row j
    Linf[j] <- alphaL + beta0 * sex[j] + lambda[j]      # Linf varies among individuals and by sex
    K[j] <- alphaK + beta1 * sex[j] + gamma[j]          # K varies among individuals and by sex
    
    } # j
    
    # Priors
    # loop over all rows - some individuals duplicated though
    for(j in 1:N) {                                 
    lambda[j] ~ dnorm(0, inv.lambda.squared)
    gamma[j] ~ dnorm(0, inv.gamma.squared)
    } # j

    alphaL ~ dnorm(1000, 0.0001)
    alphaK ~ dnorm(0, 0.0001)
    beta0 ~ dnorm(0, 0.0001)
    beta1 ~ dnorm(0, 0.0001)
    
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    sigma <- pow(sigma2, 0.5)
    log.sigma ~ dunif(-2000, 2000)
    
    # Hyperpriors
    sigma.L ~ dunif(0, 100)
    #sigma2.L <- sigma.L* sigma.L
    inv.lambda.squared <- 1 / pow(sigma.L, 2)
    
    sigma.K ~ dunif(0, 100)
    #sigma2.K <- sigma.K * sigma.K
    inv.gamma.squared <- 1 / pow(sigma.K, 2)
    
    }
    ", fill=TRUE)
sink()

#ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.diff, 
                # Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , N = nrow(ggs.selfjoin.nodupes.30.pos), 
                # sex = ggs.selfjoin.nodupes.30.pos$male)

# model runs poorly if you do not submit initial values
inits <- function() list(alphaL = 1000, alphaK = 0.0009, log.sigma = 3,
                         beta0 = rnorm(1), beta1 = rnorm(1))

params <- c("alphaL", "alphaK", "sigma", "sigma.L", "sigma.K", "beta0", "beta1")
            #, "lambda", "gamma")

ni <- 10000
nc <- 3
nb <- 2000
nt <- 1

ptm <- proc.time()
out.ggs.m2RE2sex.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, 
                                                parameters.to.save = params,
                                                model.file="ggsVB2RE2sex.txt", 
                                                n.chains = nc, n.iter = 20000, 
                                                n.burnin = 4000, n.thin = 1, 
                                                n.cluster = nc,
                                                working.directory = getwd())
proc.time() - ptm

out.ggs.RE.sex <- as.mcmc(out.ggs.m2RE2sex.jags.parallel)

print(out.ggs.m2RE2sex.jags.parallel, digits = 6)
traceplot(out.ggs.m2RE2sex.jags.parallel)

K.M <- out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$alphaK + out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$beta1
      #+ out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$gamma

K.F <- out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$alphaK 
       #+ out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$gamma

Linf.F <- out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$alphaL
          #+ out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$lambda
 
Linf.M <- out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$alphaL + out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$beta0
          #+ out.ggs.m2RE2sex.jags.parallel$BUGSoutput$sims.list$lambda

L0 <- 188.5 # data from Julia's neonate study (sd = 9.36)

age.vec <- seq(0, 3650, by = 365)
traj.jags.F <- vonBert(K = mean(K.F), w = 1 - (L0/mean(Linf.F)),
                       Linf = mean(Linf.F), A = age.vec)

traj.jags.M <- vonBert(K = mean(K.M), w = 1 - (L0/mean(Linf.M)), 
                       Linf = mean(Linf.M), A = age.vec)

plot(age.vec/365, traj.jags.F, ylab = "SVL (mm)", xlab = "Age (years)", 
     ylim=c(0,1000), type = "l", bty="n")
lines(age.vec/365, traj.jags.M, lty = 2)

# add credible intervals - only represent uncertainty around means - no ind RE
HPDinterval(out.ggs.RE.sex, prob=0.95)

library(rethinking)
KF.HPDI <- HPDI(K.F, prob = 0.95)
KM.HPDI <- HPDI(K.M, prob = 0.95)

LF.HPDI <- HPDI(Linf.F, prob = 0.95)
LM.HPDI <- HPDI(Linf.M, prob = 0.95)

traj.jags.F.lo <- vonBert(K = HPDI(K.F, prob=0.95)[1], 
                    w = 1 - (L0/HPDI(Linf.F, prob=0.95)[1]),
                    Linf = HPDI(Linf.F, prob=0.95)[1], A = age.vec)

traj.jags.F.hi <- vonBert(K = HPDI(K.F, prob=0.95)[2], 
                          w = 1 - (L0/HPDI(Linf.F, prob=0.95)[2]),
                          Linf = HPDI(Linf.F, prob=0.95)[2], A = age.vec)

traj.jags.M.lo <- vonBert(K = HPDI(K.M, prob=0.95)[1], 
                          w = 1 - (L0/HPDI(Linf.M, prob=0.95)[1]),
                          Linf = HPDI(Linf.M, prob=0.95)[1], A = age.vec)

traj.jags.M.hi <- vonBert(K = HPDI(K.M, prob=0.95)[2], 
                          w = 1 - (L0/HPDI(Linf.M, prob=0.95)[2]),
                          Linf = HPDI(Linf.M, prob=0.95)[2], A = age.vec)


lines(age.vec/365, traj.jags.F.lo, lty=1, col="gray")
lines(age.vec/365, traj.jags.F.hi, lty=1, col="gray")

lines(age.vec/365, traj.jags.M.lo, lty=2, col="gray")
lines(age.vec/365, traj.jags.M.hi, lty=2, col="gray")

# add legend
legend(x = 6, y = 550, lty=c(1, 2), legend = c("Female", "Male"), bty = "n")


##################################
### with filtered/cleaned data ###
##################################
sink("ggsVB2RE2sex.txt")
cat("
    model{
    for(j in 1:N) {
      for(i in 1:length(PIT.vec)) {
    y[j] ~ dnorm(mu[j], tau)                            # y[j] is size increment grown in row j
    mu[j] <- (Linf[j] - Lt[j]) * (1 - exp(-K[j]*x[j]))  # x is time increment of row j
    Linf[j] <- alphaL + beta0 * sex[j] + lambda[i]      # Linf varies among individuals and by sex
    K[j] <- alphaK + beta1 * sex[j] + gamma[i]          # K varies among individuals and by sex
      } #i
    } # j
    
    # Priors
    # loop over all PIT values, not over all rows 
    for(i in 1:length(PITvec)) {                                 
    lambda[i] ~ dnorm(0, inv.lambda.squared)
    gamma[i] ~ dnorm(0, inv.gamma.squared)
    } # i
    
    alphaL ~ dnorm(1000, 0.0001)
    alphaK ~ dnorm(0, 0.0001)
    beta0 ~ dnorm(0, 0.0001)
    beta1 ~ dnorm(0, 0.0001)
    
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    sigma <- pow(sigma2, 0.5)
    log.sigma ~ dunif(-2000, 2000)
    
    # Hyperpriors
    sigma.L ~ dunif(0, 100)
    #sigma2.L <- sigma.L* sigma.L
    inv.lambda.squared <- 1 / pow(sigma.L, 2)
    
    sigma.K ~ dunif(0, 100)
    #sigma2.K <- sigma.K * sigma.K
    inv.gamma.squared <- 1 / pow(sigma.K, 2)
    
    }
    ", fill=TRUE)
sink()

PIT.vec <- unique(filter2.ggs$PIT)

ggs.data <- list(x = filter2.ggs.pos$Date.diff, y = filter2.ggs.pos$SVL.diff, 
                 Lt=filter2.ggs.pos$SVL.x , N = nrow(filter2.ggs.pos), 
                 sex = filter2.ggs.pos$male)#, PIT = filter2.ggs.pos$PIT,
                 #PITvec <- unique(filter2.ggs$PIT))

# model runs poorly if you do not submit initial values
inits <- function() list(alphaL = 1000, alphaK = 0.0009, log.sigma = 3,
                         beta0 = rnorm(1), beta1 = rnorm(1))

params <- c("alphaL", "alphaK", "sigma", "sigma.L", "sigma.K", "beta0", "beta1")
#, "lambda", "gamma")

ni <- 30000
nc <- 3
nb <- 5000
nt <- 1

ptm <- proc.time()
out.ggs.m2RE2sex.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, 
                                                parameters.to.save = params,
                                                model.file="ggsVB2RE2sex.txt", 
                                                n.chains = nc, n.iter = 20000, 
                                                n.burnin = 4000, n.thin = 1, 
                                                n.cluster = nc,
                                                working.directory = getwd())
proc.time() - ptm

print(out.ggs.m2RE2sex.jags.parallel, digits=6)
#############################################################
##### Alternative parameterization based on Lt and Lt+1 #####
#############################################################

# Based on Eq. 1 in Hart and Chute (2009) ICES J. Mar. Sci.

sink("ggsVB3.txt")
cat("
    model{
    for(j in 1:N) {
    y[j] ~ dnorm(mu[j], tau) # y[j] is size at time t + 1
    mu[j] <- Lt[j] + (Linf - Lt[j]) * (1 - exp(-K*x[j] / 365)) 
    # x[j] is time increment of row j
    # Lt[j] is size at time t for individual j
  } # j
    
    #L0 ~ dunif(0, 1000) # L0 is data for each individual, not a parameter to be estimated
    #Linf <- L0 + beta
    #Linf ~ dunif(0, 2000)
    
    Linf ~ dnorm(1000, 0.0001)
    K ~ dunif(0, 100)
    
    #alpha <- Linf
    #beta ~ dunif(0, 2000)
    #gamma <- exp(-K)
    
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-2000, 2000)
    
}
", fill=TRUE)
sink()

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.y, 
                 Lt=ggs.selfjoin.nodupes.30.pos$SVL.x , N = nrow(ggs.selfjoin.nodupes.30.pos))

# model runs poorly if you do not submit initial values
inits <- function() list(Linf=1000, K=0.0009, log.sigma=3)

params <- c("Linf", "K", "sigma2")

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

out.ggs.m3.jags <- jags(data = ggs.data, inits=inits, parameters.to.save = params,
                                          model.file="ggsVB3.txt", n.chains = nc, n.iter = ni, 
                                          n.burnin = nb, n.thin = nt,
                                          working.directory = getwd())

out.ggs.m3.jags.parallel <- jags.parallel(data = ggs.data, inits=inits, parameters.to.save = params,
                                          model.file="ggsVB3.txt", n.chains = nc, n.iter = ni, 
                                          n.burnin = nb, n.thin = nt, n.cluster = nc,
                                          working.directory = getwd())

print(out.ggs.m3.jags, digits=6)
hist(out.ggs.m3.jags$BUGSoutput$sims.list$K)

traj.jags3 <- vonBert(K = 0.225, w = 1 - (200/1000), Linf = 1000, A = 0:10)
plot(0:10, traj.jags3, ylab = "SVL (mm)", xlab = "Age (years)", ylim=c(0,1000))

############################################
#### BUGS book VB model fit to GGS data ####
############################################

# problem - we do not have age data, only data on length of time between measurements
# rather than x = age, we will have t-t0, the growth increment?

sink("ggsVB1.txt")
cat("
    model{
    for(j in 1:N) {
    
    y[j] ~ dnorm(mu[j], tau) # y[j] is size measurement of individual j
    
    mu[j] <- Linf - (Linf - L0) * exp(-K*x[j]) # x is age of individual j
    
    }
    L0 ~ dunif(0, 1000)
    Linf <- L0 + beta
    K ~ dunif(0, 100)
    alpha <- Linf
    beta ~ dunif(0, 2000)
    gamma <- exp(-K)
    tau <- 1/sigma2
    log(sigma2) <- 2*log.sigma
    log.sigma ~ dunif(-10, 10)
    
    }
    ", fill=TRUE)
sink()


# Data

ggs.data <- list(x = ggs.selfjoin.nodupes.30.pos$Date.diff, y = ggs.selfjoin.nodupes.30.pos$SVL.mean, N = nrow(ggs.selfjoin.nodupes.30.pos))

# model runs poorly if you do not submit initial values
inits <- function() list(L0=200, beta=800, K=0.0009732578, log.sigma=-2)

ni <- 20000
nc <- 3
nb <- 4000
nt <- 1

params <- c("K", "Linf", "L0", "sigma2", "beta")

out.ggs.m1 <- bugs(data=ggs.data, inits=inits, parameters.to.save=params, 
                      model.file="ggsVB1.txt", n.iter=ni, n.burnin = nb, n.chains = nc,
                      n.thin = nt, debug=TRUE)

# try running in jags instead - K, Linf, and deviance seem to be unidentifiable
out.ggs.m1.jags <- jags(data=ggs.data, inits=inits, parameters.to.save=params, 
                   model.file="ggsVB1.txt", n.iter=ni, n.burnin = nb, 
                   n.chains = nc, n.thin = nt, working.directory = getwd())

print(out.ggs.m1.jags, digits=4)


