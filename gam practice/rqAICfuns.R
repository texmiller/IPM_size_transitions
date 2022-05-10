## AIC model averaging of df=(1,2,3,4) quantile regression
## df=3 or 4 are natural B-spline bases including intercept.  
## Model df are over-weighted by default, as recommended for 
## spline smoothing in ordinary regression. Here it is just
## a weak hedge against overfitting.
rqAIC_avg = function(x,y,tau, L=NULL, U=NULL, gamma=2) {
    e = order(x); x=x[e]; y=y[e]; 
    if(is.null(U)) U = max(x); if(is.null(L)) L = min(x) 

  # fit constant (df=1)
    fit.1 = rq(y~1, tau = tau); 

  # fit linear (df=2)b
     fit.2 = rq(y~x, tau = tau); 

  # fit with df=3
    B3=create.bspline.basis(rangeval=c(L,U), nbasis=3, norder=2);
    X3 = eval.basis(B3,x);
    # X3 = ns(x,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.3 = rq(y ~ X3-1, tau=tau) 

  # fit with df=4
    B4=create.bspline.basis(rangeval=c(L,U), nbasis=4, norder=3);
    X4 = eval.basis(B4,x);
    # X4 = ns(x,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.4 = rq(y ~ X4-1, tau=tau) 
    
  # fit with df=5
    B5=create.bspline.basis(rangeval=c(L,U), nbasis=5, norder=4);
    X5 = eval.basis(B5,x);
    # X4 = ns(x,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.5 = rq(y ~ X5-1, tau=tau)     
   
    AICs = c(AIC(fit.1)[1],AIC(fit.2)[1], AIC(fit.3)[1], AIC(fit.4)[1],AIC(fit.5)[1]); 
  # For AIC weights with over-weighted model df 
    z = AICs + 2*(gamma-1)*c(1:5)  
    w = exp(-0.5*(z-mean(z))); w = w/sum(w); 
    
    px = seq(L,U,length=200); 
    # X3 = ns(px,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    # X4 = ns(px,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    X3 = eval.basis(B3,px); 
    X4 = eval.basis(B4,px); 
    X5 = eval.basis(B5,px);
    
    qhat = w[1]*coef(fit.1)[1] + w[2]*(coef(fit.2)[1]+px*coef(fit.2)[2]) + w[3]*X3%*%coef(fit.3) + w[4]*X4%*%coef(fit.4) + w[5]*X5%*%coef(fit.5)
    qfun = approxfun(px,qhat,rule=1); 
    return(list(qfun=qfun,wts=w)); 
} 

## AIC selection of df=(1,2,3) quantile regression
## df=1 and 2 are constant and linear functions
## df=3 is a quadratic B-spline basis from fda  
## Model df are over-weighted by default, as recommended for 
## spline smoothing in ordinary regression. Here it is just
## a weak hedge against overfitting.
rqAIC3 = function(x,y,tau, L=NULL, U=NULL, gamma=1.4) {
    e = order(x); x=x[e]; y=y[e]; 
    if(is.null(U)) U = max(x); if(is.null(L)) L = min(x) 

  # fit constant (df=1)
    fit.1 = rq(y~1, tau = tau); 

  # fit linear (df=2)b
     fit.2 = rq(y~x, tau = tau); 

  # fit with df=3
    B3=create.bspline.basis(rangeval=c(L,U), nbasis=3, norder=3);
    X3 = eval.basis(B3,x);
    fit.3 = rq(y ~ X3-1, tau=tau) 

    AICs = c(AIC(fit.1)[1],AIC(fit.2)[1], AIC(fit.3)[1]) 
  # Over-weight model df 
    z = AICs + 2*(gamma-1)*c(1:3);
    df_min = min(which(z==min(z))); 
    if(df_min==1) fit=fit.1;
    if(df_min==2) fit=fit.2; 
    if(df_min==3) fit=fit.3 
  
    return(fit); 
} 











