#Setting a seed
seed_id = 1000
set.seed(seed_id)

#Installing packages
install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


##### Preprocessing #####

#Setting groupsize to four as it was four in the experiment
groupSize <- 4
#Played for 10 rounds (trials)
ntrials <- 10
#What the the common pool of money will be multiplied with before redistribution
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)
#vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

rawDat <- read.csv("HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

#- create covariates in raw data matrix
# Indoex for the different nations - note two cities from Austria (=9)
rawDat$nation <- c()
rawDat$nation[rawDat$city=="Melbourne"]=1
rawDat$nation[rawDat$city=="Minsk"]=2
rawDat$nation[rawDat$city=="Chengdu"]=3
rawDat$nation[rawDat$city=="Copenhagen"]=4
rawDat$nation[rawDat$city=="Bonn"]=5
rawDat$nation[rawDat$city=="Athens"]=6
rawDat$nation[rawDat$city=="Seoul"]=7
rawDat$nation[rawDat$city=="Samara"]=8
rawDat$nation[rawDat$city=="Zurich"]=9
rawDat$nation[rawDat$city=="St. Gallen"]=9
rawDat$nation[rawDat$city=="Istanbul"]=10
rawDat$nation[rawDat$city=="Nottingham"]=11
rawDat$nation[rawDat$city=="Dnipropetrovs'k"]=12
rawDat$nation[rawDat$city=="Boston"]=13
#Saudi
rawDat$nation[rawDat$city=="Riyadh"]=14
#Oman
rawDat$nation[rawDat$city=="Muscat"]=15

# create variable for Vaccine compliance. 
#Data gathered from the global covid vaccine database 
#from http://www.nature.com/articles/s41562-021-01122-8
#Above all cities it is stated which date the datapoint is from
#The index is number of people vaccinated per 100 

rawDat$vax <- c()
#Australia 2022-02-18
rawDat$vax[rawDat$city=="Melbourne"]=83.92
#Belarus 2022-02-13
rawDat$vax[rawDat$city=="Minsk"]=57.27
#China #Denmark 2022-02-18
rawDat$vax[rawDat$city=="Chengdu"]=88.94
#Denmark 2022-02-18
rawDat$vax[rawDat$city=="Copenhagen"]=81.33
#Germany 2022-02-18
rawDat$vax[rawDat$city=="Bonn"]=77.15
#Greece 2022-02-18
rawDat$vax[rawDat$city=="Athens"]=75.47
#SouthKorea 2022-02-18
rawDat$vax[rawDat$city=="Seoul"]=86.04
#Russia 2022-02-14
rawDat$vax[rawDat$city=="Samara"]=53.89
#Austria 2022-02-18
rawDat$vax[rawDat$city=="Zurich"]=76.64
rawDat$vax[rawDat$city=="St. Gallen"]=76.64
#Turkey 2022-02-18
rawDat$vax[rawDat$city=="Istanbul"]=67.51
#England 2022-02-18
rawDat$vax[rawDat$city=="Nottingham"]=78.14
#Ukraine 2022-02-18
rawDat$vax[rawDat$city=="Dnipropetrovs'k"]=39.46
#USA 2022-02-18
rawDat$vax[rawDat$city=="Boston"]=76.30
#Saudi 2022-02-18
rawDat$vax[rawDat$city=="Riyadh"]=70.97
#Oman 2022-02-14
rawDat$vax[rawDat$city=="Muscat"]=70.03

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

missing <- array(0,ngroups)

#
for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    #singlechild[s,g] <- redDat$singlechild[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][s+(s-1)*10]
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish

Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish


Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

c_choice_index <- c

vax <- array(0, ngroups)
Nation <- array(0, ngroups)
for (g in 1:ngroups) {
  vax[g] <- mean(redDat$vax[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
  Nation[g] <- mean(redDat$nation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
}

c_win <- c_no_punish[,,!is.na(vax)]
c_keep <- rep()-c_win

Gga_punish <- Gga_punish[,!is.na(vax)]
Gga_no_punish <- Gga_no_punish[,!is.na(vax)]

c <- c[,,!is.na(vax),]
Gga <- Gga[,!is.na(vax),]
Ggas <- Ggas[,,!is.na(vax),]
Gc <- Gc[,,!is.na(vax),]
Ga <- Ga[,,!is.na(vax),]
vax <- vax[!is.na(vax)]
Nation <- Nation[!is.na(Nation)]

#redefine number of groups after removing those without civic scores
ngroups <- length(vax)

# aggregate vax to just 1 number per Nation-index (using the mean here should be unproblematic since all groups within a given nation should have been given the same vax-coefficient)
vax <- aggregate(vax~Nation, FUN=mean)[,2]

nnations <- length(vax)

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

################################################################################
########################### Conditional cooperation model ######################
################################################################################

# JZS priors for partial correlation. Method described here
# https://link.springer.com/article/10.3758/s13423-012-0295-x
# Code available here
# https://github.com/MicheleNuijten/BayesMed/blob/master/R/jzs_corSD.R
# Paper where code is used here (mediation paper)
# https://link.springer.com/article/10.3758/s13428-014-0470-2

#################################################################
#------------------ Winnings analysis ---------------------------
#################################################################

#-------------------  Regress vax on winnings ---------------

# standardise variables

X <- vax
X <- (X-mean(X))/sd(X)

invSigma <- solve(t(X)%*%X) # required for JZS priors

Y <- (winnings-mean(winnings))/sd(winnings)

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX") 

# - run jags code
win.samples <- jags.parallel(data, inits=NULL, params,
                             model.file ="win_corr.txt",
                             n.chains=4, n.iter=100000, n.burnin=20000, n.thin=1, n.cluster=64)

# #----------------------- Control for all four national variables using partial correlation _------------
# fullControl.win <- function (X1,X2,X3,X4,X5,ngroups,nnations,winnings,Nation) {
#   
#   # standardise covariates
#   X1 <- (X1-mean(X1))/sd(X1)
#   X2 <- (X2-mean(X2))/sd(X2)
#   X3 <- (X3-mean(X3))/sd(X3)
#   X4 <- (X4-mean(X4))/sd(X4)
#   X5 <- (X5-mean(X5))/sd(X5)
#   
#   X <- cbind(X1,X2,X3,X4,X5)
#   V <- solve(t(X)%*%X)
#   
#   Y <- (winnings-mean(winnings))/sd(winnings)
#   
#   data <- list("ngroups", "Y", "nnations","X","V","Nation") #data inputted into jags
#   params <- c("beta0","betaX","prior_T") #parameters we'll track in jags
#   
#   # - run jags code
#   win.samples <- jags(data, inits=NULL, params,
#                       model.file ="win_partcor_4control.txt",
#                       n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1)
#   
#   return(list(win.samples))
#   
# }
# 
# Full_control.win <- fullControl.win(GDP,Indiv,Trust,Civic,vax,ngroups,nnations,winnings,Nation)

#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  Regress vax on belief weights and slope of prefs in CC model ---------------

# standardise covariate
X <- vax
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

# Ga_old <- Ga
# Ga <- Ggas

data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                            model.file ="CC_corr.txt",
                            n.chains=4, n.iter=100000, n.burnin=20000, n.thin=1, n.cluster=64)
end_time = Sys.time()
end_time - start_time

# #----------------------- Control for all four variables using partial correlation _------------
# fullControl.CC <- function (X1,X2,X3,X4,X5,groupSize,ngroups,ntrials,nnations,c,Ga,Nation) {
#   
#   # standardise covariates
#   X1 <- (X1-mean(X1))/sd(X1)
#   X2 <- (X2-mean(X2))/sd(X2)
#   X3 <- (X3-mean(X3))/sd(X3)
#   X4 <- (X4-mean(X4))/sd(X4)
#   X5 <- (X5-mean(X5))/sd(X5)
#     
#   X <- cbind(X1,X2,X3,X4,X5)
#   V <- solve(t(X)%*%X)
#   
#   data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","V","Nation") #data inputted into jags
#   params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") #parameters we'll track in jags
#   
#   # - run jags code
#   CC.samples <- jags(data, inits=NULL, params,
#                      model.file ="CC_partcor_4control.txt",
#                      n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1)
#   
#   return(list(CC.samples))
#   
# }
# 
# Full_control.CC <- fullControl.CC(GDP,Indiv,Trust,Civic,vax,groupSize,ngroups,ntrials,nnations,c,Ga,Nation)

#################################################################
#------------------ Plotting ---------------------------
#################################################################

# ------ Create empirical and parameter arrays for plots -------
# empirical group winnings data - means and standard deviations
empirical.win <- array(0,c(3,length(vax)))
for (i in 1:length(vax)) {
  empirical.win[1,i] <- mean(winnings[Nation==i]) - sd(winnings[Nation==i]) 
  empirical.win[2,i] <- mean(winnings[Nation==i]) 
  empirical.win[3,i] <- mean(winnings[Nation==i]) + sd(winnings[Nation==i])
}

empirical.initial <- array(0,c(3,length(vax)))
for (i in 1:length(vax)) {
  empirical.initial[1,i] <- mean(colMeans(c[,1,,1])[Nation==i]) - sd(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[2,i] <- mean(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[3,i] <- mean(colMeans(c[,1,,1])[Nation==i]) + sd(colMeans(c[,1,,1])[Nation==i])
}

# # summarise posterior estimates for controlled winnings model
# win.model.summary <- cbind(Full_control.win[[1]]$BUGSoutput$summary[2:6,3],
#                             Full_control.win[[1]]$BUGSoutput$summary[2:6,1],
#                             Full_control.win[[1]]$BUGSoutput$summary[2:6,7])
# 
# # summarise posterior estimates for controlled CC model - initial beliefs
# CC.model.summary.alpha <- cbind(Full_control.CC[[1]]$BUGSoutput$summary[4:8,3],
#                                 Full_control.CC[[1]]$BUGSoutput$summary[4:8,1],
#                                 Full_control.CC[[1]]$BUGSoutput$summary[4:8,7])
# 
# # summarise posterior estimates for controlled CC model - belief weighting
# CC.model.summary.omega <- cbind(Full_control.CC[[1]]$BUGSoutput$summary[9:13,3],
#                               Full_control.CC[[1]]$BUGSoutput$summary[9:13,1],
#                               Full_control.CC[[1]]$BUGSoutput$summary[9:13,7])
# 
# # summarise posterior estimates for controlled CC model - preference slopes
# CC.model.summary.rho <- cbind(Full_control.CC[[1]]$BUGSoutput$summary[14:18,3],
#                               Full_control.CC[[1]]$BUGSoutput$summary[14:18,1],
#                               Full_control.CC[[1]]$BUGSoutput$summary[14:18,7])
# 

#-------------- create high res png for saving ---------------------------------

png(filename="PGG_figure_new.png",
    height = 1500,
    width = 1800)

# set layout for graphic
layout(rbind(c(1,1,1,2,2,2,4,4,4), 
             c(1,1,1,2,2,2,4,4,4),
             c(1,1,1,3,3,3,4,4,4), 
             c(1,1,1,3,3,3,4,4,4),
             
             c(5,5,5,7,7,7,9,9,9), 
             c(5,5,5,7,7,7,9,9,9),
             c(6,6,6,8,8,8,10,10,10), 
             c(6,6,6,8,8,8,10,10,10)))

# text scale plot parameter
par(cex=1.5)

# set margins - will change when there is no title required
par(mar=c(5,3,5,2))

#----- empirical correlation - vax and group winnings --------------------
plot(c(40,100), c(0,900), type = "n", main = "A: Mean group contribution in relation to vaccine compliance", 
     xlab = "People vaccinated per 100 for each nation", ylab = "Group Contribution",axes=FALSE)
for (i in 1:length(vax)) {
  lines(c(vax[i],vax[i]),c(empirical.win[1,i],empirical.win[3,i]))
  points(vax[i],empirical.win[2,i])
}
axis(1)
axis(2)

#----- posterior of effect - vax and group winnings --------------------
plot(density(win.samples$BUGSoutput$sims.list$betaX),frame=FALSE,lwd=2,ylim=c(0,7),
     cex=2,xlab = "Standardized Vaccine Estimate",ylab = " ",main="B: Contribution")
COL <- adjustcolor(c("red"))
lines(c(win.samples$BUGSoutput$summary[2,3],win.samples$BUGSoutput$summary[2,7]),c(-.1,-.1),col=COL,lwd=2)
points(win.samples$BUGSoutput$summary[2,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title

# #----- Controls - group winnings --------------------
# plot(c(-1,1), c(-.5,-5.5), type = "n", main = " ", xlab = "Standardised Estimates", 
#      ylab = " ",axes=FALSE)
# for (i in 1:5) {
#   ifelse(sign(win.model.summary[i,1]) == sign(win.model.summary[i,3]),
#          COL <- adjustcolor(c("red")),
#          COL <- adjustcolor(c("black")))
#   lines(c(win.model.summary[i,1],win.model.summary[i,3]),c(-i,-i),col=COL,lwd=2)
#   points(win.model.summary[i,2],c(-i),pch=19,col=COL)
#   abline(h=0,v=0,col="gray60")
# }
# axis(1)
# axis(2, at = c(-1, -2, -3, -4, -5),
#      labels = c("GDP","Indiv.","Trust","Civic",
#                 "vax"),las=1,pos=-1)  
# 
# par(mar=c(5,3,5,2))

#----- empirical correlation - vax and initial contribution --------------------
plot(c(40,100), c(0,20), type = "n", main = "C: Data - Initial Contribution", 
     xlab = "People vaccinated per 100 for each nation", ylab = "Initial Contribution",axes=FALSE)
for (i in 1:length(vax)) {
  lines(c(vax[i],vax[i]),c(empirical.initial[1,i],empirical.initial[3,i]))
  points(vax[i],empirical.initial[2,i])
}
axis(1)
axis(2)

#---------- Initial Belief -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_alpha),frame=FALSE,lwd=2,ylim=c(0,35),
     cex=2,xlab = "Standardized Vaccine Estimate",ylab = " ",main=" Initial Belief (Alpha)")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[4,3],CC.samples$BUGSoutput$summary[4,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[4,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title

# plot(c(-.6,.6), c(-.5,-5.5), type = "n", main = " ", xlab = "Standardised Estimates", 
#      ylab = " ",axes=FALSE)
# for (i in 1:5) {
#   ifelse(sign(CC.model.summary.alpha[i,1]) == sign(CC.model.summary.alpha[i,3]),
#          COL <- adjustcolor(c("red")),
#          COL <- adjustcolor(c("black")))
#   lines(c(CC.model.summary.alpha[i,1],CC.model.summary.alpha[i,3]),c(-i,-i),col=COL,lwd=2)
#   points(CC.model.summary.alpha[i,2],c(-i),pch=19,col=COL)
#   abline(h=0,v=0,col="gray60")
# }
# axis(1)
# axis(2, at = c(-1, -2, -3, -4, -5),
#      labels = c("GDP","Indiv.","Trust","Civic",
#                 "vax"),las=1,pos=-.6)  
# 
# par(mar=c(5,3,5,2))

#---------- Belief Learning -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_omega),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardized vaccine Estimate",ylab = " ",main="Belief Learning Weight (Omega)")
COL <- adjustcolor(c("black"))
lines(c(CC.samples$BUGSoutput$summary[5,3],CC.samples$BUGSoutput$summary[5,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[5,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))

# plot(c(-.6,.6), c(-.5,-5.5), type = "n", main = NULL, xlab = "Standardised Estimates", 
#      ylab = " ",axes=FALSE)
# for (i in 1:5) {
#   ifelse(sign(CC.model.summary.omega[i,1]) == sign(CC.model.summary.omega[i,3]),
#          COL <- adjustcolor(c("red")),
#          COL <- adjustcolor(c("black")))
#   lines(c(CC.model.summary.omega[i,1],CC.model.summary.omega[i,3]),c(-i,-i),col=COL,lwd=2)
#   points(CC.model.summary.omega[i,2],c(-i),pch=19,col=COL)
#   abline(h=0,v=0,col="gray60")
# }
# axis(1)
# axis(2, at = c(-1, -2, -3, -4, -5),
#      labels = c("GDP","Indiv.","Trust","Civic",
#                 "vax"),las=1,pos=-.6)  
# 
# par(mar=c(5,3,5,2))

#---------- Preference Slope -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_rho),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardized vaccine Estimate",ylab = " ",main="Conditional Preferences (Rho)")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[6,3],CC.samples$BUGSoutput$summary[6,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[6,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))

# plot(c(-.6,.6), c(-.5,-5.5), type = "n", main = " ", xlab = "Standardised Estimates", 
#      ylab = " ",axes=FALSE)
# for (i in 1:5) {
#   ifelse(sign(CC.model.summary.rho[i,1]) == sign(CC.model.summary.rho[i,3]),
#          COL <- adjustcolor(c("red")),
#          COL <- adjustcolor(c("black")))
#   lines(c(CC.model.summary.rho[i,1],CC.model.summary.rho[i,3]),c(-i,-i),col=COL,lwd=2)
#   points(CC.model.summary.rho[i,2],c(-i),pch=19,col=COL)
#   abline(h=0,v=0,col="gray60")
# }
# axis(1)
# axis(2, at = c(-1, -2, -3, -4, -5),
#      labels = c("GDP","Indiv.","Trust","Civic",
#                 "vax"),las=1,pos=-.6)  
# 
traceplot(CC.samples, mfrow=c())
CC.samples$BUGSoutput
dev.off()
