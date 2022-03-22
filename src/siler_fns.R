"
A bunch of functions to help plot and explore Siler mortality
"
require(numDeriv)

# A function to generate a mortality curve from siler function parameters
siler <- function(pars, ages){
  mu = exp(- (pars$b*ages+pars$B)) + exp(pars$c*ages-pars$C) + pars$d
  lmort = log(mu)
  mort = exp(lmort)
  
  mort[which(mort > 1)] <- 1
  return(mort)
}

# Generate a survival curve from siler function parameters
siler_survival <- function(pars, ages){
  lS = -pars$d*ages + (1/pars$b)*(exp(-(pars$b*ages + pars$B)) - exp(-pars$B)) -
    (1/pars$c)*(exp(pars$c*ages - pars$C) - exp(-pars$C))
  S = exp(lS)
  return(S) 
}

## Functions for derivative of S wrt B
siler_SB <- function(pars, ages){
  S_B <- (exp(-pars$b*pars$B) - exp(-pars$b*(ages + pars$B)))*siler_survival(pars, ages)
  return(S_B) 
}

siler_SB_num <- function(pars, ages){
  S_grad <- rep(NA, length(ages))
  for (ii in 1:length(ages)){
    func0 <- function(x){pars$B <- x; siler(pars, ages[ii])}
    S_grad[ii] <- grad(func0, 1.5)
  }
  return(S_grad)
}


## Functions for derivative of S wrt b
siler_Sb <- function(pars, ages){
  S_b <- (1/pars$b)*(pars$B*exp(-pars$b*pars$B) - 
                       (ages + pars$B)*exp(-pars$b*(ages + pars$B)))*siler_survival(pars, ages)
  return(S_b) 
}

## Functions for derivative of S wrt C
siler_SC <- function(pars, ages){
  S_C <- (exp(pars$c*(ages - pars$C)) - exp(-pars$c*pars$C))*siler_survival(pars, ages)
  return(S_C) 
}

## Functions for derivative of S wrt c
siler_Sc <- function(pars, ages){
  S_c <- (1/pars$c)*(pars$C*exp(-pars$c*pars$C) - 
                       (ages - pars$C)*exp(pars$c*(ages - pars$C)))*siler_survival(pars, ages)
  return(S_c) 
}

## Functions for derivative of S wrt d
siler_Sd <- function(pars, ages){
  S_d <- - ages*siler_survival(pars, ages)
  return(S_d) 
}


# Function for life expectancy
siler_LE <- function(pars, ages){
  #LE <-  siler_survival(pars, ages)/ 
  #  (pars$d + exp(-pars$b*(ages+pars$B)) + exp(pars$c*(ages-pars$C)))
  LE <- as.vector(rep(NA, length(ages)))
  for (ii in 1:length(ages)){
    #integrand <- function(x) {siler_survival(pars, x)}
    #LE[ii] <- integrate(integrand, lower = ages[ii], upper = Inf)$value
    LE[ii] <- sum(siler_survival(pars, ages[ii]:200))
    LE[ii] <- LE[ii]/siler_survival(pars, ages[ii])
  }
  return(LE)
}

# Function for lifespan inequality
siler_ineq <- function(pars, ages){
  #LE <-  siler_survival(pars, ages)/ 
  #  (pars$d + exp(-pars$b*(ages+pars$B)) + exp(pars$c*(ages-pars$C)))
  H <- as.vector(rep(NA, length(ages)))
  for (ii in 1:length(ages)){
    LE <- siler_LE(pars, ages[ii]) 
    S <- siler_survival(pars, ages[ii]:200)
    S <- S[S>0]
    H[ii] <- -sum(S*log(S))/LE
  }
  return(H)
}

siler_ineq_theta <- function(LE, LE_theta, H, S_theta, S){
  temp <- S_theta*log(S) 
  temp <- temp[!is.na(temp)]
  H_theta <- -(1/LE)*(LE_theta*(1+H) + sum(temp))
  return(H_theta)
} 

