#################################################################
# Recruit size data
#################################################################
recruitSizes=c(16.7, 6.4, 10, 4.1, 4.7, 11.2, 6.3, 9.7, 5); #these are areas of known recruits  

#xmin=min(log(recruitSizes)); xmax=max(log(recruitSizes));
#mins=maxs=numeric(500); 
#for(k in 1:500) {
#	z=rnorm(length(recruitSizes),mean=mean(log(recruitSizes)),sd=sd(log(recruitSizes)));
#	mins[k]=min(z); maxs[k]=max(z); 
#}	
# hist(mins); abline(v=xmin); hist(maxs); abline(v=xmax); 

##########################################################
# Setup the main data frame: do this before any analysis
##########################################################
X=read.csv("e:\\pubs\\mardis\\bruno\\data\\lesion_final_condensed2.csv",as.is=TRUE); 
X$Area.T5=as.numeric(X$Area.T5)
#### NOTE: the l_f_c2 file has correct NA's for state=I and areas not measured. 

# Combine into one list of transitions with Year (=first year) as factor. 
x23=subset(X,select=c(Site,Fan.number,State.T2,Area.T2,Infected.T2,Lesion.T2,State.T3,Area.T3,Infected.T3,Lesion.T3)); x23$Year="2002"; 
x34=subset(X,select=c(Site,Fan.number,State.T3,Area.T3,Infected.T3,Lesion.T3,State.T4,Area.T4,Infected.T4,Lesion.T4)); x34$Year="2003";
x45=subset(X,select=c(Site,Fan.number,State.T4,Area.T4,Infected.T4,Lesion.T4,State.T5,Area.T5,Infected.T5,Lesion.T5)); x45$Year="2004";

newnames=c("Site","Fan.number","State1","Area1","Infected1","Lesion1","State2","Area2","Infected2","Lesion2","Year"); 
names(x23)=newnames; names(x34)=newnames; names(x45)=newnames; 
XC=rbind(x23,x34,x45);
e=!is.na(XC$State1); XC=XC[e,];  

XC$H1=XC$Area1-XC$Infected1; XC$H1[XC$H1<0]=0; 
XC$H2=XC$Area2-XC$Infected2; XC$H2[XC$H2<0]=0; 
XC$State1[XC$H1==0]="D"; 
XC$State2[XC$H2==0]="D";

XH = subset(XC, (State1=="H")&(State2=="H") ); 
XH = subset(XH, select=c(Site,Fan.number,Year,Area1,Area2))
XH = na.omit(XH); 

XH$Site=factor(XH$Site); XH$Year=factor(XH$Year); 
XH$logarea.t0 = log(XH$Area1); 
XH$logarea.t1 = log(XH$Area2); 