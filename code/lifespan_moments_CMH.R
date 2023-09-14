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
