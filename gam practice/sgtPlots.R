x=seq(-5,5,length=500); 
plot(x,dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=100),type="l"); 
points(x,dsgt(x,mu=0,sigma=1,lambda=0,p=100,q=2),type="l",col="red")
points(x,dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=2),type="l",col="blue")
points(x,dsgt(x,mu=0,sigma=1,lambda=0,p=100,q=100),type="l",col="purple"); 

x=seq(-5,5,length=500); 
py1 = dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=1.25);
py2 = dsgt(x,mu=0,sigma=1,lambda=0,p=4,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 

require(splines);
z=rnorm(1000); z=sort(z); z=z*abs(z); hist(z); 

B = create.bspline.basis(rangeval=range(z), norder=4, nbasis=4)
x = seq(min(z),max(z),length=200); 
out = eval.basis(x,B);   
matplot(x,out,type="l"); 