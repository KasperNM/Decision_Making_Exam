install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm, VGAM, reshape2)

set.seed(1982)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
groupSize <- 4
ntrials <- 10
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)

niterations <- 20 # 100

# THIS WILL MATCH NUMBER OF GROUPS IN THE NO PUNISHMENT CONDITION IN THE ORIG DATA
ngroups <- 269

true_mu_alpha_log <- array(NA,c(niterations))
true_mu_rho_probit <- array(NA,c(niterations))
true_mu_omega_probit <- array(NA,c(niterations))

infer_mu_alpha_log <- array(NA,c(niterations))
infer_mu_rho_probit <- array(NA,c(niterations))
infer_mu_omega_probit <- array(NA,c(niterations))

source('CC_simple.R')

start_time = Sys.time()
for (i in 1:niterations) {
  
  # let's see how robust the model is. Does it recover all sorts of values?
  #-------------------------------------------------------------------
  #-------------------  Group-level priors -------------
  #-------------------------------------------------------------------
  mu_alpha_log <- runif(1,1.5,3) # close estimate to getting most values between 0 and 20
  
  mu_rho_probit <- runif(1,-3,3)
  
  mu_omega_probit <- runif(1,-3,3)
  
  #------ Individual priors ---------------------------------------
  #-----------------------------------------------------------------
  #------ variance priors ------------------------------------------
  # standard deviation for initial belief (inverse gamma priors)
  tau_alpha <- rgamma(1,100,1)
  sigma_alpha <- 1/sqrt(tau_alpha) 
  
  # concentration (precision) of rate parameters (for beta priors)
  sigma_rho <- runif(1,1,100)
  sigma_omega <- runif(1,1,100)
  
  # reparameterising gamma prior for beliefs about others on first trial 
  mu_alpha <- loglink(mu_alpha_log, inverse=TRUE) # standardisation of initial belief mode - log link
  rate_alpha <- ( mu_alpha + sqrt( mu_alpha^2 + 4*sigma_alpha^2 ))/
    (2*sigma_alpha^2) 
  shape_alpha <- 1 + mu_alpha * rate_alpha
  
  #reparamaterising beta prior for slope of preference sin CC model
  mu_rho <- probitlink(mu_rho_probit, inverse=TRUE) # standardisation of rate estimate mean - probit link
  shape1_rho <- (mu_rho) * sigma_rho
  shape2_rho <- (1 - mu_rho) * sigma_rho 
  
  #reparamaterising beta prior for belief updating in CC model
  mu_omega <- probitlink(mu_omega_probit, inverse=TRUE) # standardisation of rate estimate mean - probit link
  shape1_omega <- (mu_omega) * sigma_omega
  shape2_omega <- (1 - mu_omega) * sigma_omega
  
  source('CC_simple.R')
  CC_sims <- CC_single(groupSize, ntrials, ngroups, rate_alpha, shape_alpha,
                       shape1_rho, shape2_rho, shape1_omega, shape2_omega)
  
  c <- CC_sims$c
  Ga <- CC_sims$Ga
  
  data <- list("groupSize", "ngroups", "ntrials","c","Ga") 
  params <- c("mu_alpha_log", "mu_rho_probit", "mu_omega_probit") 
  
  # - run jags code
  CC.samples <- jags.parallel(data, inits=NULL, params,
                              model.file ="CC_corr_simple.txt",
                              n.chains=4, n.iter=100000, n.burnin=20000, n.thin=1, n.cluster=64)
  
  true_mu_alpha_log[i] <- mu_alpha_log
  true_mu_rho_probit[i] <- mu_rho_probit
  true_mu_omega_probit[i] <- mu_omega_probit
  
  # find maximum a posteriori
  Y <- CC.samples$BUGSoutput$sims.list
  infer_mu_alpha_log[i] <- MPD(Y$mu_alpha_log)
  infer_mu_rho_probit[i] <- MPD(Y$mu_rho_probit)
  infer_mu_omega_probit[i] <- MPD(Y$mu_omega_probit)
  
  print(i)
  
}
end_time = Sys.time()
end_time - start_time


# plotting code courtesy of Lasse
source('recov_plot.R')
pl1 <- recov_plot(true_mu_alpha_log, infer_mu_alpha_log, c("true mu alpha", "infer mu alpha"), 'smoothed linear fit')
pl2 <- recov_plot(true_mu_rho_probit, infer_mu_rho_probit, c("true mu rho", "infer mu rho"), 'smoothed linear fit')
pl3 <- recov_plot(true_mu_omega_probit, infer_mu_omega_probit, c("true mu omega", "infer mu omega"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3)

# traceplot(CC.samples, mfrow=c(4,2))