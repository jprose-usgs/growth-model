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



