rm(list=ls())
library(readrba)
library(readabs)
library(xts)
library(tseries)
library(urca)
library(FinTS)
library(rmarkdown)
library(parallel)
library(MASS)
library(coda)
library(sads)
library(truncnorm)

# Data downloads
########################################################################### 

## Monetary policy variables
# Cash interest rate
cash_rate.dl   <- read_rba(series_id = "FIRMMCRTD")
cashrate <- to.quarterly(xts(cash_rate.dl$value, cash_rate.dl$date), 
                         OHLC = FALSE)
# M3 Money supply
M3.dl   <- read_rba(series_id = "DMAM3N")
M3 <- to.quarterly(xts(M3.dl$value, M3.dl$date), 
                   OHLC = FALSE)


# Unemployment rate: A84423092X
unemp_rate.dl <- read_abs(series_id = "A84423092X")
unemp <- to.quarterly(xts(unemp_rate.dl$value, unemp_rate.dl$date), 
                      OHLC = FALSE)

# real GDP: A2302459A
realGDP.dl   <- read_abs(series_id = "A2302459A")
realGDP <- to.quarterly(xts(realGDP.dl$value, realGDP.dl$date), 
                        OHLC = FALSE)

## Fiscal variables
# Total tax: A2301963V
totaltax.dl <- read_abs(series_id = "A2301963V")
totaltax <- to.quarterly(xts(totaltax.dl$value, totaltax.dl$date), 
                         OHLC = FALSE)
# Non-tax revenue: Gross income (A2302106V) - total tax 
govgrossinc.dl <- read_abs(series_id = "A2302106V")
govgrossinc <- to.quarterly(xts(govgrossinc.dl$value, govgrossinc.dl$date), 
                            OHLC = FALSE)
nontax <- govgrossinc - totaltax

# Public gross fixed capital formation: A2302555A
pubinv.dl <- read_abs(series_id = "A2302555A")
pubinv <- to.quarterly(xts(pubinv.dl$value, pubinv.dl$date), 
                       OHLC = FALSE)

# Social assistance benefits payments: A2301976F
pubtrans.dl <- read_abs(series_id = "A2301976F")
pubtrans <- to.quarterly(xts(pubtrans.dl$value, pubtrans.dl$date), 
                         OHLC = FALSE)
# Government final consumption: A2302527T
pubcons.dl <- read_abs(series_id = "A2302527T")
pubcons <- to.quarterly(xts(pubcons.dl$value, pubcons.dl$date),
                        OHLC = FALSE)

# Merge data into one matrix
Y.df  <- na.omit(merge(unemp, realGDP , totaltax, nontax, pubinv, 
                       pubtrans, pubcons, cashrate, M3))

varname_vec <- c("Unemployment rate", "Real GDP", "Tax revenue", "Non-tax revenue", "Gov't GFCF",
                 "Social benefits payments", "Gov't consumption", "Cash rate target", "M3 supply")
colnames(Y.df) <- varname_vec

Y.df <- Y.df[1:132,]
# Transform into natural logs
lnY.df <- log(Y.df)

date <- as.vector(index(cashrate))[1:132]
T <- length(date)


# data for input to estimation function 
lnY.df_num <- c()

for (i in 1:ncol(lnY.df)){
  lnY.df_num <- cbind(lnY.df_num, as.numeric(lnY.df[,i]))
}

# Create Y and X
N = ncol(lnY.df_num)

# p = no. of lags
p = 4
K = 1 + p*N

Y       = ts(lnY.df_num[(p+1):nrow(lnY.df_num),], start=c(1991,1), frequency=4)
# Y       = lnY.df_num[5:nrow(lnY.df_num),]
# nrow(Y)
X       = matrix(1,nrow(Y),1)
# nrow(X)
for (i in 1:p){
  X     = cbind(X,lnY.df_num[5:nrow(lnY.df_num)-i,])
}


# Metropolis-Hastings
########################################################################### 
S <- 1000
# initialize theta values
# v0, v1, v2, rho
theta_old <- c(17, 65, 20, 0.8)
c <- 0.000525
W <- diag(length(theta_old))
# draw S thetas 

Theta <- matrix(NA,S,4)
set.seed(1)
for (s in 1:S){
  
  # covid volatility likelihood
  v.lik <- function(V, N) det(diag(V^2))^(-N/2)
  
  # covid volatility prior: pareto(1,1) for v's and beta(3,1.5) for rho
  v.prior <- function(theta, pareto.a=1, pareto.b=1, beta.a=3, beta.b=1.5){
    
    beta.cons <- 1/beta(beta.a,beta.b)
    
    (pareto.a*pareto.b^pareto.a)/(theta[1]^(pareto.a+1))*
      (pareto.a*pareto.b^pareto.a)/(theta[2]^(pareto.a+1))*
      (pareto.a*pareto.b^pareto.a)/(theta[3]^(pareto.a+1))*
      beta.cons*theta[4]^(beta.a-1)*(1-theta[4])^(beta.b-1)
    #dnorm(theta[4],0,1)
  }
  
  # creates a vector of v's for the covid period (Q2 2020 til latest)
  covid.vec <- function(theta){
    vec <- theta[1:3]
    for (i in 4:12){
      vec <- c(vec, 1 + (theta[3]-1)*theta[4]^(i-1))
    }
    
    return(vec)
  }
  
  # creates vector of 1's
  v_ones <- ts(rep(1, nrow(Y)-12), c(1991,1), frequency = 4) 
  # appends covid period vector at the end
  V.old <- c(v_ones, covid.vec(theta_old))
  
  # new candidate parameters values
  theta_new <- mvrnorm(1, theta_old, c*W)
  # rho <- rtruncnorm(1, a=0, b=1, mean = theta_old[4], sd = sqrt(c))
  # theta_new[4] <- rho  
  
  V.new <- c(v_ones, covid.vec(theta_new))
  
  # calculate posteriors 
  v.posterior_old <- v.lik(V.old, N)*v.prior(theta_old)
  v.posterior_new <- v.lik(V.new, N)*v.prior(theta_new)
  
  # alpha for accpetance/rejection
  alpha <- min(1, v.posterior_new/v.posterior_old)
  
  u_star <- runif(1)
  
  if (alpha > u_star){
    Theta[s,] <- theta_new
  } else {Theta[s,] <- theta_old}
  
  # set theta_old to latest accepted theta
  theta_old <- Theta[s,]  
}

plot.ts(Theta)
# Acceptance rate
1 - rejectionRate(as.mcmc(Theta[,1]))

# Last estimated parameter values
Theta[nrow(Theta),]
# Covariance matrix
cov(Theta)


```{r}
mh.mcmc <- function(data, p=4, start_date = c(1991,1),
                    S.mh = 1000, c = 0.00015, W = diag(4), theta.init = c(13, 65, 20, 0.8)){
  N = ncol(data)
  # p = no. of lags
  K = 1 + p*N
  # forecast horizon
  # h       = 8
  
  Y       = ts(data[(p+1):nrow(data),], start=start_date, frequency=4)
  X       = matrix(1,nrow(Y),1)
  # nrow(X)
  for (i in 1:p){
    X     = cbind(X,data[(p+1):nrow(data)-i,])
  }
  
  # v0, v1, v2, rho
  Theta <- matrix(NA,S.mh,4)
  theta_old <- theta.init
  #theta_old <- Theta[nrow(Theta),]
  
  set.seed(1)
  pb = txtProgressBar(min = 0, max = S.mh, initial = 0) 
  for (s in 1:S.mh){
    setTxtProgressBar(pb,s)
    # Covid volatility likelihood
    v.lik <- function(V, N) det(diag(V^2))^(-N/2)
    
    # Covid volatility prior
    v.prior <- function(theta, pareto.a=1, pareto.b=1, beta.a=3, beta.b=1.5){
      beta.cons <- 1/beta(beta.a,beta.b)
      
      (pareto.a*pareto.b^pareto.a)/(theta[1]^(pareto.a+1))*
        (pareto.a*pareto.b^pareto.a)/(theta[2]^(pareto.a+1))*
        (pareto.a*pareto.b^pareto.a)/(theta[3]^(pareto.a+1))*
        beta.cons*theta[4]^(beta.a-1)*(1-theta[4])^(beta.b-1)
      #dnorm(theta[4],0,1)
    }
    
    covid.vec <- function(theta){
      vec <- theta[1:3]
      for (i in 4:12){
        vec <- c(vec, 1 + (theta[3]-1)*theta[4]^(i-3))
      }
      
      return(vec)
    }
    
    
    v_ones <- ts(rep(1, nrow(Y)-12), c(1991,1), frequency = 4) 
    V.old <- c(v_ones, covid.vec(theta_old))
    
    # New candidate parameters values
    theta_new <- mvrnorm(1, theta_old, c*W)
    
    V.new <- c(v_ones, covid.vec(theta_new))
    
    # calculate posteriors 
    v.posterior_old <- v.lik(V.old, N)*v.prior(theta_old)
    v.posterior_new <- v.lik(V.new, N)*v.prior(theta_new)
    
    alpha <- min(1, v.posterior_new/v.posterior_old)
    
    u_star <- runif(1)
    
    if (alpha > u_star){
      Theta[s,] <- theta_new
    } else {Theta[s,] <- theta_old}
    
    theta_old <- Theta[s,]  
  }
  
  colnames(Theta) <- c("v0", "v1" , "v2", "rho")
  #plot.ts(Theta)
  1 - rejectionRate(as.mcmc(Theta[,1]))
  #Theta[nrow(Theta),]
  re <- list(Theta=Theta, AcceptRate = 1 - rejectionRate(as.mcmc(Theta[,1])))
}


mh1 <- mh.mcmc(lnY.df_num, c = 0.00015, W = diag(4), theta.init = c(13, 65, 20, 0.8), S.mh = 1000)
plot.ts(mh1$Theta)
mh1$AcceptRate
mh1$Theta[nrow(mh1$Theta),]
```
