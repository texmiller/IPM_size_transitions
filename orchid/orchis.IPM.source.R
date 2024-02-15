##############################################################################################################
########## define demographic functions 
#Transforms all values below/above limits in min/max size to prevent eviction
##############################################################################################################
invlogit<-function(x){exp(x)/(1+exp(x))}

#PRODUCTION OF PROTOCORMS BY Y-SIZED INDIVIDUALS
fx<-function(x,params){
  xb<-pmin(pmax(x,params$minsize),params$maxsize) 
  pflow<-invlogit(params$flow.int+params$flow.size*xb) 
	nflow<-exp(params$flowers.int+params$flowers.size*xb)
	protocorms<-pflow*nflow*invlogit(params$fruits)*params$seeds*params$germ
	return(protocorms)
	}

#GROWTH FROM SIZE X TO Y
gxy_GAU<-function(x,y,params){
  xb<-pmin(pmax(x,params$minsize),params$maxsize) 
  pflow<-invlogit(params$flow.int+params$flow.size*xb) 
  gveg<-params$grow.int+params$grow.size*xb
  gflow<-params$grow.int+params$grow.flow+(params$grow.size+params$grow.size.flow)*xb
  mu<-gveg*(1-pflow)+gflow*pflow
	return(dnorm(x=y,mean=mu,sd=exp(params$growsd.int+params$growsd.fit*mu+params$growsd.fit2*mu^2)))
}

gxy_SST<-function(x,y,params){
  xb<-pmin(pmax(x,params$minsize),params$maxsize) 
  pflow<-invlogit(params$flow.int+params$flow.size*xb) 
  gveg<-params$grow.int+params$grow.size*xb
  gflow<-params$grow.int+params$grow.flow+(params$grow.size+params$grow.size.flow)*xb
  mu<-gveg*(1-pflow)+gflow*pflow
  return(dSST(x=y,mu=mu,
              sigma=exp(params$growsd.int+params$growsd.fit*mu+params$growsd.fit2*mu^2),
              nu=exp(params$growSST.nu.int + params$growSST.nu.fit*mu),
              tau=exp(params$growSST.tau.int + params$growSST.tau.fit*mu)+2))
}

## size distribution of plants that emerge from tubers
recruits<-function(y,params){
	dnorm(y,mean=params$kidsize,sd=params$kidsize.sd)
}

#recruits<-function(y,params){
#  fac1<-sqrt(2*pi)*params$kidsize.sd
#  fac2<-((y-params$kidsize)^2)/(2*params$kidsize.sd)
#  kids<-(exp(-fac2)/fac1)
#  return(kids)
#}


## size distribution of plants that emerge from dormancy
wakeup<-function(y,params){
	dnorm(y,mean=params$Dsize,sd=params$Dsizesd)
}
#wakeup<-function(y,params){
#  fac1<-sqrt(2*pi)*params$Dsizesd
#  fac2<-((y-params$Dsize)^2)/(2*params$Dsizesd)
#  plants<-(exp(-fac2)/fac1)
#  return(plants)
#}


#SURVIVAL AT SIZE X
sx<-function(x,params){
  xb<-pmin(pmax(x,params$minsize),params$maxsize) 
  invlogit(params$surv.int+params$surv.size*xb)
  #return(1)
	}

#Probability of dormancy at size x
dx<-function(x,params){
  xb<-pmin(pmax(x,params$minsize),params$maxsize) 
  invlogit(params$dorm.int+params$dorm.size*xb)
	}

#Survival * growth * not going dormant
pxy<-function(x,y,params,dist){
	sx(x,params)*(1-dx(x,params))*do.call(paste0("gxy_",dist),list(x,y,params))
	}

###R0 functions, one for each habitat
returnR0<-function(params,dist,lower.extend=0,upper.extend=0){

	L<-params$minsize - lower.extend
	U<-params$maxsize + upper.extend
	n<-params$matsize#matrix size
	h<-(U-L)/n
	b<-L+c(0:n)*h;
	y<-0.5*(b[1:n]+b[2:(n+1)])

	T<-matrix(0,(n+3),(n+3))## matrix size reflects length of size distribution + protocorms, tubers, and dormant plants
	T[3,4:(n+3)]<-sx(y,params)*dx(y,params) ## probabilities of dormancy go in this row, conditioned on survival
	T[4:(n+3),3]<-wakeup(y,params)*h	##distribution of plants emerging from dormancy -- we need a probability that they survival dormancy
	T[4:(n+3),4:(n+3)]<-h*t(outer(y,y,pxy,params=params,dist=dist)) ##growth function applied to the inner matrix					   
	T[4:(n+3),2]<-recruits(y,params)*params$sigmat*h #put recruit sizes in the matrix & reduce by probability tuber becomes a recruit
	T[2,1]<-params$sigmap##survival of protocorms to tubers
	
	F<-matrix(0,(n+3),(n+3))
	F[1,4:(n+3)]<-fx(y,params) ## production of new protocorms goes in the top row; no integration here bc this returns a count not a probability

	M<-T+F #sum the T & F matrix to get the whole matrix
		 #this is the discrete approximation to the continuous kernel
	
	## compute R0
	s <- length(diag(T))
	N <- try(solve(diag(s) - T), silent = TRUE)
	if (class(N)[1] == "try-error") {
	  r <- NA
	}
	else {
	  R <- F %*% N
	  r <- Re(eigen(R)$values[1])
	}
    
	#return(r)
  #use this return to get lambda      
	return(list(R0=r,matrix=M,Tmatrix=T,Fmatrix=F,meshpts=y))
	}

