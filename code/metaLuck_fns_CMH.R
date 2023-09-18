#####################################################
# Functions for luck analyses
#####################################################    

## some toy matrices for testing things out:
Umat<- matrix(c(0.5, 0.1, 0.1, 0, 0.5, 0.3, 0, 0, 0.8), ncol=3)
Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
Amat<- Umat+Fmat

## Fmat2<- matrix(0,3,3); Fmat2[1,]=1; 
## Amat2<- Umat + Fmat2; 

fundamental_matrix = function(P){
    solve( diag(ncol(P)) - P)
}

# Calculate the expected value of lifetime reproductive success
mean_LRO<- function(Umat, Fmat, mixdist=NULL, offspring_weight=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }
  
  ## Calculate Ex(R | birth size z)
  N<- exactLTRE::fundamental_matrix(Umat)
  expRCond_z<- rep(1,Nclasses)%*%Fmat%*%N
  
  if(!is.null(mixdist)){
    expR<- expRCond_z%*%mixdist
    return(expR)
  } else{
    return(expRCond_z)
  }
}

# Calculate the variance in lifetime reproductive success for a given population
# matrix, given a number of options: based on the Markov Chain method and moments
# repro_var can be: poisson, bernoulli, or fixed
# mixdist, if provided, should be a vector for combining across offspring types
# offspring_weight, if provided, should be a vector with a weight for each offspring type
var_LRO_mcr<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL, offspring_weight=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  
  # Add offspring weights:
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }
  
  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it
  
  # Z matrix operator:
  Zmat<- cbind(diag(1, Nclasses, Nclasses, names=FALSE), 0)
  # column matrix of 1's:
  OneVec<- rep(1, Nclasses+1)
  # The fundamental matrix: 
  Nmat<- exactLTRE::fundamental_matrix(Umat)
  # Take the sum of offspring:
  eff<- colSums(Fmat) #eff is the stage-specific offspring production
  
  # Calculate the raw moments of the reward matrix:
  R1<- OneVec %*% t(c(eff, 0)) # first moment is the stage-specific reproductive output
  if (repro_var %in% c("poisson", "Poisson")){
    R2<- R1 + (R1*R1)
    R3<- R1 + 3*(R1*R1) + (R1*R1*R1)
  } else if (repro_var %in% c("fixed", "Fixed")){
    R2<- R1*R1
    R3<- R1*R1*R1
  } else if (repro_var %in% c("bernoulli", "Bernoulli")){
    R2<- R1
    R3<- R1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }
  # Rtilde1 is the first raw moment of the reward matrix for only transient states
  Rtilde1<- Zmat%*%R1%*%t(Zmat)
  # Rtilde2 is the second raw moment of the reward matrix for only transient states
  Rtilde2<- Zmat%*%R2%*%t(Zmat)
  
  # Calculate the raw moments of LRO conditional on starting state:
  mu_prime1<- t(Nmat)%*%Zmat %*% t(Pmat*R1)%*%OneVec
  mu_prime2<- t(Nmat)%*% (Zmat %*% t(Pmat*R2)%*%OneVec + 2*t(Umat*Rtilde1) %*% mu_prime1)
  # mu_prime3<- t(Nmat)%*% (Zmat %*% t(Pmat*R3)%*%OneVec + 3*t(Umat*Rtilde2) %*% mu_prime1 + 3*t(Umat*Rtilde1) %*% mu_prime2)
  
  ## The first raw moment of LRO conditional on starting state is the expected value:
  expRCond_z<- mu_prime1
  
  ## Var(R | birth size z)
  varRCond_z<- mu_prime2 - (mu_prime1*mu_prime1)
  
  if(is.null(mixdist)){
    return(varRCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varR_within<- t(mixdist)%*%varRCond_z
    # variance in LRO due to differences among starting states:
    varR_between<- t(mixdist)%*%(expRCond_z^2) - (t(mixdist)%*%expRCond_z)^2
    # total variance in LRO, given the mixing distribution:
    varR<- varR_within + varR_between
    return(varR)
  }
}

var_LRO_ipmbook<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL, offspring_weight=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
  
  # Add offspring weights:
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }
  
  # Take the sum of offspring:
  betabar<- colSums(Fmat) #betabar is the average offspring production
  
  # Variance in per-capita number of offspring per time step:
  if (repro_var %in% c("poisson", "Poisson")){
    sigsqb<- betabar # poisson
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    sigsqb<- betabar*(1 - betabar) # bernoulli
  } else if (repro_var %in% c("fixed", "Fixed")){
    sigsqb<- 0
  }
  
  ## Calculate Ex(R | birth size z)
  expRCond_z<- mean_LRO(Umat, Fmat)
  
  ## Var(R | birth size z) (see ch. 3 of the IPM book, eq. 3.2.15 (p. 82 of pdf))
  N<- exactLTRE::fundamental_matrix(Umat)
  rbarPib = expRCond_z %*% Umat
  r2 = (sigsqb + (betabar)^2 + 2*betabar*rbarPib) %*% N
  varRCond_z = (r2 - expRCond_z^2)
  
  if(is.null(mixdist)){
    return(varRCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varR_within<- varRCond_z%*%mixdist 
    # variance in LRO due to differences among starting states:
    varR_between<- t(mixdist)%*%t(expRCond_z^2) - (t(mixdist)%*%t(expRCond_z))^2
    # total variance in LRO, given the mixing distribution:
    varR<- varR_within + varR_between
    return(varR)
  }
}

# Calculate skew in lifetime reproductive success:
skew_LRO<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL, offspring_weight=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  Nclasses<- dim(Umat)[1]
    
  # Add offspring weights:
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }
    
  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it
    
  # Z matrix operator:
  Zmat<- cbind(diag(1, Nclasses, Nclasses, names=FALSE), 0)
  # column matrix of 1's:
  OneVec<- rep(1, Nclasses+1)
  # The fundamental matrix: 
  Nmat<- exactLTRE::fundamental_matrix(Umat)
  # Take the sum of offspring:
  eff<- colSums(Fmat) #eff is the stage-specific offspring production
    
  # Calculate the raw moments of the reward matrix:
  R1<- OneVec %*% t(c(eff, 0)) # first moment is the stage-specific reproductive output
  if (repro_var %in% c("poisson", "Poisson")){
    R2<- R1 + (R1*R1)
    R3<- R1 + 3*(R1*R1) + (R1*R1*R1)
  } else if (repro_var %in% c("fixed", "Fixed")){
    R2<- R1*R1
    R3<- R1*R1*R1
  } else if (repro_var %in% c("bernoulli", "Bernoulli")){
    R2<- R1
    R3<- R1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }
  
  # Rtilde1 is the first raw moment of the reward matrix for only transient states
  Rtilde1<- Zmat%*%R1%*%t(Zmat)
  # Rtilde2 is the second raw moment of the reward matrix for only transient states
  Rtilde2<- Zmat%*%R2%*%t(Zmat)
    
  # Calculate the raw moments of LRO conditional on starting state:
  mu_prime1<- t(Nmat)%*%Zmat %*% t(Pmat*R1)%*%OneVec
  mu_prime2<- t(Nmat)%*% (Zmat %*% t(Pmat*R2)%*%OneVec + 2*t(Umat*Rtilde1) %*% mu_prime1)
  mu_prime3<- t(Nmat)%*% (Zmat %*% t(Pmat*R3)%*%OneVec + 3*t(Umat*Rtilde2) %*% mu_prime1 + 3*t(Umat*Rtilde1) %*% mu_prime2)
    
  ## The first raw moment of LRO conditional on starting state is the expected value:
  expRCond_z<- mu_prime1
    
  ## Second central moment: mu2 = Var(R | birth size z)
  varRCond_z<- mu_prime2 - (mu_prime1*mu_prime1)
  
  # third central moment of LRO, conditional on Z:
  mu3Cond_z<- mu_prime3 - 3*(mu_prime2*mu_prime1) + 2*(mu_prime1*mu_prime1*mu_prime1)
  # Skewness, conditional on Z:
  skewCond_z<- varRCond_z^(-3/2)*mu3Cond_z
  
  if(is.null(mixdist)){
    return(skewCond_z)
  } else{ # law of total cumulance for third central moment
    # expected value of third central moment:
    mu3_within<- t(mu3Cond_z) %*% mixdist
    # third central moment of the expected value:
    mu3_between<- t((mu_prime1-(t(mu_prime1)%*%mixdist)[1])^3) %*% mixdist
    # covariance of expected value and variance:
    covExpVar<- t(mu_prime1*varRCond_z)%*%mixdist - (t(mu_prime1)%*%mixdist) * (t(varRCond_z)%*%mixdist)
    # total third central moment:
    mu3_total<- mu3_within + mu3_between + 3*covExpVar
    
    # total variance:
    var_total<- t(varRCond_z)%*%mixdist + t(mu_prime1^2)%*%mixdist - (t(mu_prime1)%*%mixdist)^2
    
    # total skewness:
    skew_total<- (var_total)^(-3/2)*mu3_total
    return(skew_total)
  }
}

# Function for mean lifespan. I already wrote essentially the same function in 
# the exactLTRE package, but I want to add the mixing distribution option.
mean_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  ## Calculate Ex(R | current state)
  N<- exactLTRE::fundamental_matrix(Umat)
  expLCond_z<- rep(1,Nclasses)%*%N
  
  if(!is.null(mixdist)){
    expL<- expLCond_z%*%mixdist
    return(expL)
  } else{
    return(expLCond_z)
  }
}

# Calculate the variance in lifespan:
# note: this calculates the variance in the number of time steps!
var_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  # calculate the fundamental matrix
  N<- exactLTRE::fundamental_matrix(Umat)
  
  ## Calculate Ex(R | current state)
  expLCond_z<- mean_lifespan(Umat, mixdist = NULL)
  
  ## Var(L | current state) using eqn. 5.12 from Hal's book:
  eT<- matrix(data=1, ncol=Nclasses, nrow=1) # column vector of 1's
  varLCond_z<- eT %*% (2*N%*%N - N) - (expLCond_z)^2
  
  if(is.null(mixdist)){
    return(varLCond_z)
  } else{
    # variance in LRO due to differences along trajectories:
    varL_within<- varLCond_z %*% mixdist 
    # variance in LRO due to differences among starting states:
    varL_between<- t(mixdist)%*%t(expLCond_z^2) - (t(mixdist)%*%t(expLCond_z))^2
    # total variance in lifespan, given the mixing distribution:
    varL<- varL_within + varL_between
    return(varL)
  }
}

# Function for skew in lifespan:
skew_lifespan<- function(Umat, mixdist=NULL){
  # quick check that Umat is square:
  if (dim(Umat)[1]!=dim(Umat)[2]){
    warning('Umat and Fmat are not square matrices.')
  }
  Nclasses<- dim(Umat)[1]
  
  # calculate the fundamental matrix
  N<- exactLTRE::fundamental_matrix(Umat)
  # row vector of 1's
  eT<- matrix(data=1, ncol=Nclasses, nrow=1)
  # identity matrix
  eye<- diag(x=1, nrow=Nclasses, ncol=Nclasses, names=FALSE)
  
  # calculate the moments of lifespan
  eta1<- mean_lifespan(Umat)
  eta2<- eta1%*%(2*N-eye)
  eta3<- eta1%*%(6*N%*%N-6*N+eye)
  
  # calculate the central moments:
  eta_hat2<- eta2-eta1*eta1
  eta_hat3<- eta3 - 3*eta1*eta2 + 2*eta1*eta1*eta1
  
  # calculate skew:
  skewCond_z<- eta_hat2^(-3/2)*eta_hat3
  
  if(is.null(mixdist)){
    return(skewCond_z)
  } else{ # law of total cumulance for third central moment
    # expected value of third central moment:
    eta_hat3_within<- eta_hat3 %*% mixdist
    # third central moment of the expected value:
    eta_hat3_between<- ((eta1-(eta1%*%mixdist)[1])^3) %*% mixdist
    # covariance of expected value and variance:
    covExpVar<- (eta1*eta_hat2)%*%mixdist - (eta1%*%mixdist) * (eta_hat2%*%mixdist)
    # total third central moment:
    eta_hat3_total<- eta_hat3_within + eta_hat3_between + 3*covExpVar
    
    # total variance:
    var_total<- eta_hat2%*%mixdist + (eta1)^2%*%mixdist - (eta1%*%mixdist)^2
    
    # total skewness:
    skew_total<- (var_total)^(-3/2)*eta_hat3_total
    return(skew_total)
  }
}

# Function for offspring weighting by their probability of reproducing at least once:
prob_repro_wt<- function(Umat, Fmat, repro_var="poisson"){
  
  # calculate probability of reproducing:
  offspring_weight<- prob_repro(Umat, Fmat, repro_var)
  
  return(offspring_weight)
}

# Function for calculating the stable distribution:
stable_dist<- function(Amat){
  # Calculate the eigen values and vectors:
  eigz<- eigen(Amat)
  # Find the index of lambda:
  I<- which(Re(eigz$values)==max(Re(eigz$values)))
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,I]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
}

# Function for calculating the mixing distribution:
mixing_distro<- function(Amat, Fmat){
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,1]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
  # Multiply the stable distribution by Fmat to get a cohort of offspring:
  offspring<- Fmat%*%wmean
  # Rescale to sum to 1 (proportions):
  offspring<- offspring/sum(offspring)
  
  return(offspring)
}

# Function for coefficient of variation:
coeff_var<- function(mean, var){
  sqrt(var)/mean
}

# Function for the probability of reproducing at all:
prob_repro<- function(Umat, Fmat, repro_var = 'poisson'){
  # take the column sum of reproduction:
  betabar<- colSums(Fmat) # the average offspring production
  
  # Breeding probability values, p_b:
  if (repro_var %in% c("poisson", "Poisson")){
    p_b<- 1-ppois(0, lambda=betabar) # poisson: prob(x>0)
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    p_b<- betabar # bernoulli: prob(occurrence)
  } else if (repro_var %in% c("fixed", "Fixed")){
    p_b<- ifelse(betabar>0, 1, 0) # all individuals in that size reproduce the exact number of offspring shown
  }
  
  P_0<- sweep(Umat, MARGIN = 2, (1-p_b), '*') # survival matrix with reproduction as an absorbing state
  N_0<- fundamental_matrix(P_0)
  
  # probability of reproducing at least once:
  B<- p_b%*%N_0
  return(B)
}

# Function for expected age at first reproduction:
mean_age_repro<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL){
  # take the column sum of reproduction:
  betabar<- colSums(Fmat) #betabar is the average offspring production
  
  # Breeding probability values, p_b:
  if (repro_var %in% c("poisson", "Poisson")){
    p_b<- 1-ppois(0, lambda=betabar) # poisson: prob(x>0)
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    p_b<- betabar # bernoulli: prob(occurrence)
  } else if (repro_var %in% c("fixed", "Fixed")){
    p_b<- betabar # all individuals in that size reproduce the exact number of offspring shown
  }
  
  P_0<- sweep(Umat, MARGIN = 2, (1-p_b), '*') # survival matrix with reproduction as an absorbing state
  N_0<- fundamental_matrix(P_0) # matrix_inverse(I-P_0)
  # probability of reproducing at least once:
  B<- p_b%*%N_0
  
  # survival kernel P_0 conditioned on breeding at least once:
  operatormat<- matrix(B, nrow=nrow(P_0), ncol=ncol(P_0), byrow = FALSE)/matrix(B, nrow=nrow(P_0), ncol=ncol(P_0), byrow=TRUE)
  operatormat[is.nan(operatormat)]<- 0
  P_0<- P_0*operatormat
  
  # calculate the expected age at first reproduction:
  e<- rep(1,nrow(P_0))
  N_0<- fundamental_matrix(P_0) # matrix_inverse(I-P_0)
  expArCond_z<- e%*%N_0
  
  if(is.null(mixdist)){
    return(expArCond_z)
  } else{
    ## From Page 71 of IPM book, note that the mean age of repro for a cohort 
    ## has to be weighted by the different probabilities of reproducing.
    mixdist<- t(mixdist)*B/sum(t(mixdist)*B)
    expAr<- sum(expArCond_z*mixdist)
    return(expAr)
  }
}

lambda<- function(Amat){
  lambda<- Re(eigen(Amat, only.values = T)$values[1])
  return(lambda)
}

# Function to calculate the damping ratio (for any matrix):
damping_ratio<- function(Amat){
  eigenz<- Re(eigen(Amat, only.values = T)$values)
  eigenz<- eigenz[order(eigenz, decreasing = TRUE)]
  
  DR<- abs(eigenz[2])/eigenz[1]
  return(DR)
}

# Function for adding transparency to colors:
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# A function for calculating S, as an iteroparity score, from Demetrius' entropy:
# This function will calculate the age-based lx and mx for the provided mixing 
# distribution (i.e. cohort structure), and then calculates iteroparity following
# the calculation provided in Rage package: https://github.com/jonesor/Rage/blob/main/R/entropy_d.R
iteroparity_d<- function(Umat, Fmat, mixdist, nage=300){
  Amat<- Umat+Fmat
  e<- matrix(1, nrow=1, ncol=ncol(Umat))

  # age-specific survival and fertility:
  surv_age<- matrix(ncol=ncol(Umat), nrow=nage)
  surv_age[1,]<- 1
  repro_age<- matrix(ncol=ncol(Umat), nrow=nage)
  Umati<- diag(ncol(Umat)); # not Umat, because surv_age[2,] results from Umat, not Umat%*%Umat. 
  
  # repro_age[1] is calculated outside of the loop. The loop below starts at i=2. 
  repro_age[1,] = e%*%Fmat; 
  
  for (i in 2:nage){
    Umati<- Umat%*%Umati 
    surv_age[i,]<- e%*%Umati
    repro_age[i,]<- e%*%Fmat%*%Umati/surv_age[i,]
   }
  lx<- surv_age%*%mixdist
  repro_age[is.na(repro_age)]<- 0
  mx<- repro_age%*%mixdist
  
  lxmx<- lx*mx; px = lxmx/sum(lxmx); ## SPE: the scaling happens HERE.
  # if lxmx == 0, log(lxmx) == -Inf; for entropy calc below, these -Inf can be
  #  converted to 0, because lim(x * log(x)) as x->0 is 0
  log_px <- log(px)
  log_px[px==0]<- 0
  
  this_iteroparity<- -sum(px * log_px)
  return(list(iteroparity=this_iteroparity,lx=lx,mx=mx)) 
}

# Function to calculate typical "clutch size" (i.e., reproductive output per 
# time step) at the stable stage distribution:
clutch_size<- function(Amat, Fmat, offspring_weight=NULL){
  # If we are using offspring weights, rescale F:
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }
  
  # Calculate the stable population distribution:
  wmean<- Re(eigen(Amat)$vectors[,1]) # right eigenvector of the input matrix
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
  # Multiply the stable distribution by Fmat to get a cohort of offspring:
  offspring<- Fmat%*%wmean
  # Calculate the proportion of the population that is reproductively active:
  repro_active<- which(colSums(Fmat)>0)
  # The clutch size is the sum of the offspring cohort, divided by the 
  # reproductively active individuals:
  clutchsize<- sum(offspring)/sum(wmean[repro_active])
  
  return(clutchsize)
}

## Function to simulate an individual life history for a given Markov Chain
# Assumes Poisson distributed reproduction
SimulateMarkovIndiv<- function(Umat, Fmat, maxage=100, mixdist=NULL){
  # quick check that Umat and Fmat are the same size, and that they are both square:
  if (dim(Umat)[1]==dim(Fmat)[1] & dim(Umat)[2]==dim(Fmat)[2]){
    if (dim(Umat)[1]!=dim(Umat)[2]){
      warning('Umat and Fmat are not square matrices.')
    }
  } else {warning('Umat and Fmat are not the same size.')}
  
  Nclasses<- dim(Umat)[1]
  
  # Build the Markov Chain model:
  mortrow<- 1-colSums(Umat) # probability of mortality
  Pmat<- rbind(Umat, mortrow, deparse.level = 0) # add the mortality row
  Pmat<- cbind(Pmat, 0) # add a column of 0's for the dead individuals
  Pmat[Nclasses+1, Nclasses+1]<- 1 # Death is an absorbing state, individuals cannot leave it
  
  B = apply(Pmat,2,cumsum); ## Cumulative sums of each column, used in coin-toss for state transitions
  # For simulating reproduction, average repro output by state:
  effVec<- c(colSums(Fmat), 0) ## column sums of Fmat, plus no reproduction for dead individuals
  
  ## Why cumulative sums? To decide your next state, generate a Unif(0,1) random number
  ## and compare it to the i^th column of B, where i is your current state. If the random
  ## number exceeds none of the entries in B, you go to state 1. If the random number
  ## exceeds exactly one of the entries in B, you go to state 2. And so on -- this gives 
  ## correct transition probabilities conditional on your current state. 
  
  nt = maxage             ## Length of the simulation 
  states=numeric(nt+1);   ## Vector to hold states over time 
  repro<- numeric(nt+1)   ## Vector to hold reproductive output over time
  rd=runif(nt);           ## Do all the Unif(0,1) coin tosses at once 
  if(is.null(mixdist)){
    states[1] = 1;          ## Start in state 1 
  } else{
    birth_prob<- cumsum(mixdist)
    states[1]<- sum(birth_prob<runif(1))+1  ## Coin toss for starting state
  }
  
  for(i in 1:nt) {
    repro[i]<- rpois(1, effVec[states[i]])         ## Reproduce based on current state, Poisson distro 
    b=B[,states[i]];             ## Relevant cumulative probabilities, given current state  
    states[i+1]=sum(rd[i]>b)+1   ## ``Coin toss'' based on current state 
  }
  output<- list(states=states, repro=repro)
  return(output)
}

# Function for the age-decomposition of variance and skewness in LRO into
# contributions from birth state, survival trajectory, growth trajectory, and
# fecundity. This function assumes Poisson-distributed reproduction Inputs are:
# the survival matrix, P; the fertility matrix, Fmat; the distribution of
# starting/birth states, c0; the variance structure of reproductive output in
# each event, repro_var; the optional offspring weighting function,
# offspring_wt; and the maximum age to use for the decomposition, maxAge. The
# default for repro_var is "poisson" but it can also be "fixed." The default
# maxAge is 100.
getVarSkewnessPartitionsNoEnvVar = function (P, Fmat, c0, repro_var="poisson",
                                             offspring_weight=NULL,
                                             maxAge=100,debugging=FALSE) {
  # Input parameters:
  # P is the survival/growth transition matrix (called U in COMPADRE)
  # Fmat is the matrix of fertility transitions
  # c0 is the birth state distribution (also called mixing distribution)
  
  if (debugging) {
    ## some toy matrices for testing things out:
    P<- matrix(c(0.5, 0.07, 0.05, 0, 0.5, 0.21, 0, 0, 0.8), ncol=3)
    Fmat<- matrix(c(0.2, 0, 0, 1, 0.5, 0, 5, 3, 0), ncol=3)
    
    c0 = c(0.9, 0.1, 0)
    maxAge = 100
  }
  
  if (!is.null(offspring_weight)){
    Fmat<- sweep(Fmat, MARGIN = 1, offspring_weight, FUN="*")
  }
  
  mz = dim(P)[1]
  
  ## Make survival, growth, and fecundity bullet matrices
  Sbullet = diag (colSums (P))
  
  ## Does not work if any stage has zero probability of survival. 
  ## Replaced by code below that seems to deal with that issue. 
  #Sinv = diag (1 / diag(Sbullet))
  #Gbullet = P %*% Sinv
  
  Gbullet = P;  
  for(ic in 1:ncol(Gbullet)){
    if(Sbullet[ic,ic]>0) Gbullet[,ic]=Gbullet[,ic]/Sbullet[ic,ic]
  }       
  cat(range(Gbullet%*%Sbullet-P), "should equal 0", "\n"); 
  
  ## #kids isn't part of the state definition, so Fbullet = I
  Fbullet = diag (mz)
  
  ## create matrices for extended state space, which I will denote with
  ## the prefix "es".
  
  esmz = 3*mz
  esF = esP = matrix (0, esmz, esmz)
  
  esP[1:mz, 2*mz + 1:mz] = Gbullet
  esP[mz + 1:mz, 1:mz] = Fbullet
  esP[2*mz + 1:mz, mz + 1:mz] = Sbullet
  
  esF[mz + 1:mz, 1:mz] = Fmat
  
  ################################################################
  ## Extended state model that includes ur-stage alpha_z.  Used for
  ## calculating contributions from birth stage.  The "bs" before
  ## variable names stands for "birth state."
  ################################################################
  
  bsmz = 1 + mz
  ## initial state --- everybody starts off in the ur-state.
  bsC0 = rep(0, bsmz)
  bsC0[1] = 1
  
  bsF = bsP = matrix (0, bsmz, bsmz)
  ## From ur-state to birth state
  bsP[1 + (1:mz),1] = c0
  ## From normal state to normal state
  bsP[1 + 1:mz, 1 + 1:mz] = P
  
  ## You only start reproducing once you're past the ur-state.  Note
  ## that all we need are the column sums of bsF.
  bsF[1, 1 + 1:mz] = colSums (Fmat)
  
  ## set up reward matrices ###############################
  
  esR1 = esR2 = esR3 = matrix (0, esmz+1, esmz+1)
  bsR1 = bsR2 = bsR3 = matrix (0, bsmz+1, bsmz+1)
  
  ## First moment of clutch size (Poisson)
  for (j in 1:esmz) 
    esR1[,j] = sum(esF[,j])
  for (j in 1:bsmz) 
    bsR1[,j] = sum(bsF[,j])
  
  if (repro_var %in% c("Poisson", "poisson")){
    ## Second moment of clutch size  
    esR2 = esR1 + esR1^2
    bsR2 = bsR1 + bsR1^2
    ## Third moment of clutch size
    esR3 = esR1+ 3*esR1*esR1 + esR1*esR1*esR1
    bsR3 = bsR1+ 3*bsR1*bsR1 + bsR1*bsR1*bsR1
  } else if (repro_var %in% c("fixed", "Fixed")){
    ## Second moment of clutch size  
    esR2 = esR1^2
    bsR2 = bsR1^2
    ## Third moment of clutch size
    esR3 = esR1*esR1*esR1
    bsR3 = bsR1*bsR1*bsR1
  } else {
    stop("Reproductive random variable type not recognized. The available options are Poisson, Fixed, and Bernoulli")
  }
  
  ## Get moments conditional on init. state #########################
  
  out = calcMoments (esP, esR1, esR2, esR3)
  esRho1Vec = out$rho1Vec
  esMu2Vec = out$mu2Vec
  esSkewnessVec = out$skewnessVec
  
  out = calcMoments (bsP, bsR1, bsR2, bsR3)
  bsRho1Vec = out$rho1Vec
  bsMu2Vec = out$mu2Vec
  bsMu3Vec = out$mu3Vec
  bsSkewnessVec = out$skewnessVec
  
  ## Start calculating luck terms ###################################
  
  survTrajecSkewness = growthTrajecSkewness = 
    fecSkewness = numeric (maxAge)
  fecUpdateSkewness = survUpdateSkewness =
    growthUpdateSkewness =  numeric (maxAge)
  survUpdateVar = growthUpdateVar = fecUpdateVar =
    numeric (maxAge) 
  survTrajecVar = growthTrajecVar = fecVar =
    numeric (maxAge)
  surv = numeric (maxAge)
  
  ## for ur-state var.
  ## expectation wrt z of Var(success | z)
  expZVarCondZ = bsMu2Vec %*% bsC0
  ## variance wrt z of Exp(success | z)
  varZExpCondZ = sum(bsC0*bsRho1Vec^2) - (sum(bsC0*bsRho1Vec)^2)
  ## Law of total variance
  urVar = expZVarCondZ + varZExpCondZ
  
  justBirthStateVar = bsMu2Vec %*% bsP %*% bsC0
  birthStateVar = urVar - justBirthStateVar
  
  ## for ur-state skewness
  expZMu3CondZ = bsMu3Vec %*% bsC0
  ## 3rd moment wrt z of Exp(success | z)
  rho3ZExpCondZ = sum(bsC0*(bsRho1Vec)^3)
  ## 2nd moment wrt z of Exp(success | z)
  rho2ZExpCondZ = sum(bsC0*(bsRho1Vec)^2)
  ## 1st moment wrt z of Exp(success | z)
  expZExpCondZ = sum(bsC0*bsRho1Vec)
  ## 3rd central moment wrt z of Exp(success | z)
  mu3ZExpCondZ = rho3ZExpCondZ -
    3*expZExpCondZ*rho2ZExpCondZ +
    2*expZExpCondZ^3
  ## covariance wrt z of Exp(success | z) and Var(success | z)
  covZExpVarCondZ = sum(bsC0*bsRho1Vec*bsMu2Vec) -
    sum(bsC0*bsRho1Vec)*sum(bsC0*bsMu2Vec)
  ## Law of total cumulance
  urMu3 = expZMu3CondZ + mu3ZExpCondZ +
    3*covZExpVarCondZ
  urSkewness = urMu3 / urVar^1.5
  
  justBirthStateSkewness = bsSkewnessVec %*% bsP %*% bsC0
  birthStateSkewness = urSkewness - justBirthStateSkewness
  
  ## Starting the age loop
  PaC0 = c0  ## P^a c_0
  mzZero = rep(0, mz)
  
  ## Sanity check: passes
  if (debugging) {
    N = solve(diag(mz) - P)
    rho1Vec = colSums (Fmat %*% N)
    cat (esRho1Vec %*% c(c0, mzZero, mzZero),
         "should = ", rho1Vec %*% c0, "\n")
  }
  
  fecUpdateSkewness[1] = esSkewnessVec %*% c(mzZero, PaC0, mzZero)
  survUpdateSkewness[1] = esSkewnessVec %*% 
    c(mzZero, mzZero, Sbullet %*% PaC0)
  growthUpdateSkewness[1] = esSkewnessVec %*% 
    c(P %*% PaC0, mzZero, mzZero)
  
  fecUpdateVar[1] = esMu2Vec %*% c(mzZero, PaC0, mzZero)
  survUpdateVar[1] = esMu2Vec %*% c(mzZero, mzZero, Sbullet %*% PaC0)
  growthUpdateVar[1] = esMu2Vec %*% c(P %*% PaC0, mzZero, mzZero)
  
  surv[1] = sum (PaC0)
  ## a is actually age + 1, since the first year of life is age 0
  for (a in 2:maxAge) {
    PaC0 = P %*% PaC0
    
    fecUpdateSkewness[a] = esSkewnessVec %*% c(mzZero, PaC0, mzZero)
    survUpdateSkewness[a] = esSkewnessVec %*% 
      c(mzZero, mzZero, Sbullet %*% PaC0)
    growthUpdateSkewness[a] = esSkewnessVec %*% 
      c(P %*% PaC0, mzZero, mzZero)
    
    fecUpdateVar[a] = esMu2Vec %*% c(mzZero, PaC0, mzZero)
    survUpdateVar[a] = esMu2Vec %*% c(mzZero, mzZero, Sbullet %*% PaC0)
    growthUpdateVar[a] = esMu2Vec %*% c(P %*% PaC0, mzZero, mzZero)
    
    ## fecundity skewness and var.
    fecSkewness[a] = growthUpdateSkewness[a-1] - fecUpdateSkewness[a]
    fecVar[a] = growthUpdateVar[a-1] - fecUpdateVar[a]
    
    ## survival trajectory skewness and var.
    survTrajecSkewness[a-1] = fecUpdateSkewness[a-1] -
      survUpdateSkewness[a-1]
    survTrajecVar[a-1] = fecUpdateVar[a-1] - survUpdateVar[a-1]
    
    ## growth trajectory skewness and var.
    growthTrajecSkewness[a-1] = survUpdateSkewness[a-1] -
      growthUpdateSkewness[a-1] 
    growthTrajecVar[a-1] = survUpdateVar[a-1] - growthUpdateVar[a-1]
    
    surv[a] = sum(PaC0)
  }
  
  ## And get the first value of fecSkewness and fecVar
  fecSkewness[1] = justBirthStateSkewness - fecUpdateSkewness[1]
  fecVar[1] = justBirthStateVar - fecUpdateVar[1]
  
  ## By what age are 99% of individuals dead?
  lifespan = which(surv < 0.01)[1]
  
  totSurvTrajecSkewness = sum(survTrajecSkewness)
  totGrowthTrajecSkewness = sum(growthTrajecSkewness)
  totFecSkewness = sum(fecSkewness)
  totSkewness = bsSkewnessVec %*% bsC0
  ## Sanity check: passes
  if (debugging) 
    cat (totSurvTrajecSkewness + totGrowthTrajecSkewness +
           totFecSkewness + birthStateSkewness, "should = ",
         totSkewness, "\n")
  
  totSurvTrajecVar = sum(survTrajecVar)
  totGrowthTrajecVar = sum(growthTrajecVar)
  totFecVar = sum(fecVar)
  totVar = bsMu2Vec %*% bsC0
  ## Sanity check: passes
  if (debugging)
    cat (totSurvTrajecVar + totGrowthTrajecVar +
           totFecVar + birthStateVar, "should = ",
         totVar, "\n")
  
  if (debugging) {
    par (cex.lab=1.4, lwd=2)
    
    ymax = max(survTrajecSkewness, growthTrajecSkewness)
    ymin = min(survTrajecSkewness, growthTrajecSkewness)
    plot (0:25, survTrajecSkewness[1:26], col="black",
          xlab="Age", ylab="Skewness contributions", type="l",
          ylim=c(ymin, ymax))
    lines (0:25, growthTrajecSkewness[1:26], col="orange")
    lines (0:25, fecSkewness[1:26], col="blue")
    lines (rep(lifespan, 2), c(0, 0.5*ymax), col="black",
           lty=2, lwd=2)
    
    dev.new ()
    par (cex.lab=1.4, lwd=2)
    ymax = max(survTrajecVar, growthTrajecVar)
    ymin = min(survTrajecVar, growthTrajecVar)
    plot (0:25, survTrajecVar[1:26], col="black",
          xlab="Age", ylab="Var contributions", type="l",
          ylim=c(ymin, ymax))
    lines (0:25, growthTrajecVar[1:26], col="orange")
    lines (0:25, fecVar[1:26], col="blue")
    lines (rep(lifespan, 2), c(0, 0.45*ymax), col="black",
           lty=2, lwd=2)
  }
  
  return (out=list(birthStateVar=birthStateVar,
                   survTrajecVar=survTrajecVar,
                   growthTrajecVar=growthTrajecVar,
                   fecVar=fecVar,
                   birthStateSkewness=birthStateSkewness,
                   survTrajecSkewness=survTrajecSkewness,
                   growthTrajecSkewness=growthTrajecSkewness,
                   fecSkewness=fecSkewness))
}

# Function for the lifespan of individuals who reproduce at least once:
lifespan_reproducers<- function(Umat, Fmat, repro_var = 'poisson', mixdist=NULL){
  # take the column sum of reproduction:
  betabar<- colSums(Fmat) #betabar is the average offspring production
  
  # Breeding probability values, p_b:
  if (repro_var %in% c("poisson", "Poisson")){
    p_b<- 1-ppois(0, lambda=betabar) # poisson: prob(x>0)
  } else if (repro_var %in% c("Bernoulli", "bernoulli")){
    p_b<- betabar # bernoulli: prob(occurrence)
  } else if (repro_var %in% c("fixed", "Fixed")){
    p_b<- betabar # all individuals in that size reproduce the exact number of offspring shown
  }
  
  # sub-matrices for the markov chain w/ 2 absorbing states:
  P_T<- sweep(Umat, MARGIN = 2, (1-p_b), '*') # survival matrix with reproduction as an absorbing state
  P_R<- sweep(Umat, MARGIN = 2, (p_b), '*') # probability of transitioning to the reproduction state
  mort<- 1-colSums(Umat)
  
  # build the big markov chain:
  zeroz<- matrix(0, nrow=nrow(Umat), ncol=ncol(Umat))
  Q_0<- cbind(rbind(P_T, P_R), rbind(zeroz, Umat)) # transient states
  
  N_0<- fundamental_matrix(Q_0) # matrix_inverse(I-Q_0)
  
  # probability of being absorbed into omega2 on the next step:
  a_L<- c(rep(0, ncol(Umat)), mort)
  # probability of being absorved into omega2, conditional on starting state:
  q_L<- a_L%*%N_0
  # check that things are correct:
  q_s<- c(mort, rep(0, ncol(Umat)))%*%N_0
  if(sum((1-(q_L+q_s))>1e-12)>0){
    warning("There is an issue with calculating lifespan conditional on reproducing at least once.")
  }
  
  # survival kernel Q_0 conditioned on absorbing into omega2:
  operatormat<- matrix(q_L, nrow=nrow(Q_0), ncol=ncol(Q_0), byrow = FALSE)/matrix(q_L, nrow=nrow(Q_0), ncol=ncol(Q_0), byrow=TRUE)
  operatormat[is.nan(operatormat)]<- 0
  Q_0<- Q_0*operatormat
  
  # calculate the expected age at absorption into omega2:
  e<- rep(1,nrow(Q_0))
  N_0<- fundamental_matrix(Q_0) # matrix_inverse(I-P_b)
  expLCond_z<- e%*%N_0
  
  if(is.null(mixdist)){
    return(expLCond_z)
  } else{
    mixdist<- c(mixdist, rep(0, ncol(Umat))) # individuals cannot start in the post-reproductive transient stages
    if (sum(mixdist>0)==1){
      expL<- sum(expLCond_z*mixdist)
    } else {
      mixdist<- t(mixdist)*q_L/sum(t(mixdist)*q_L)
      expL<- sum(expLCond_z*mixdist)
    }
    return(expL)
  }
}

################################################################################
## Alternative measures for generation time:
################################################################################
gen_time_R0<- function(Umat, Fmat){
  Amat <- Umat + Fmat
  Nmat <- fundamental_matrix(Umat)
  R0 <- max(Re(eigen(Fmat %*% Nmat, only.values = T)$values))
  
  this_lambda<- max(Re(eigen(Amat, only.values=T)$values))
  
  gen_time<- log(R0)/log(this_lambda)
  return(gen_time)
}

gen_time_mu1_v<- function(Umat, Fmat){
  Nclasses<- ncol(Umat)
  e<- matrix(1, nrow=1, ncol=Nclasses)

  # Calculate the fundamental matrix for survival:
  N_0<- fundamental_matrix(Umat)
  # Calculate the stable stage distribution:
  eigz <- eigen(Umat+Fmat)
  ilambda <- which(Re(eigz$values) == max(Re(eigz$values)))
  w <- Re(eigz$vectors[, ilambda]); w = w/sum(w); 
  
  eigtz <- eigen(t(Umat+Fmat))
  ilambda <- which(Re(eigtz$values) == max(Re(eigtz$values)))
  v <- Re(eigtz$vectors[, ilambda]); v = v/sum(v); 
  v = matrix(v,nrow=1); 
  
  
  # calculation from Ellner 2018 (Am Nat)
  numerator<- (v%*%Fmat)%*%(N_0%*%N_0)%*%(Fmat%*%w)
  denominator<- (v%*%Fmat)%*%N_0%*%(Fmat%*%w)
  gen_time<- numerator/denominator
  
  return(gen_time)
}

gen_time_Ta<- function(Umat, Fmat){
    Amat<- Umat + Fmat
    eigz <- eigen(Amat)
    ilambda <- which(Re(eigz$values) == max(Re(eigz$values)))
    lambda <- max(Re(eigz$values))
    w <- Re(eigz$vectors[, ilambda])
    v <- Re(eigen(t(Amat))$vectors[, ilambda])
    gentime <- lambda*sum(v*w)/sum(v*(Fmat %*% w)) 
    return(gentime)
}

################################################################################
## Plotting functions
################################################################################

plot_var_components_quartiles<- function(exes, meanz, lower_bar, upper_bar, ax.labz, legend_loc){
  colz<- viridis::turbo(8) # black, (2) blue, (4) green, (6) orange, (8) brown
  
  plot(exes, meanz[,2], type='l', col=colz[6], ylim=c(0, 1), lwd=2,
       xlab=ax.labz[1], ylab=ax.labz[2])
  arrows(x0=exes, y0=lower_bar[,2], x1=exes, y1=upper_bar[,2], code=3, angle=90, length=0.05, col=colz[6])
  lines(exes, meanz[,5], type='l', col=colz[8], lwd=2)
  arrows(x0=exes, y0=lower_bar[,5], x1=exes, y1=upper_bar[,5], code=3, angle=90, length=0.05, col=colz[8])
  lines(exes, meanz[,3], type='l', col=colz[4], lwd=2)
  arrows(x0=exes, y0=lower_bar[,3], x1=exes, y1=upper_bar[,3], code=3, angle=90, length=0.05, col=colz[4])
  lines(exes, meanz[,4], type='l', col=colz[2], lwd=2)
  arrows(x0=exes, y0=lower_bar[,4], x1=exes, y1=upper_bar[,4], code=3, angle=90, length=0.05, col=colz[2])
  legend(legend_loc, title = "Luck component",
         legend=c("Survival", "Growth", "Fecundity", "Birth state"), 
         col=colz[c(6,4,2,8)], lty=1, lwd=2)
}
