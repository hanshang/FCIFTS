## Description: Simulation codes for eigenvalue ratio estimators and variance ratio tests.
## The code can produce the results given in Tables 1 and 2 of the paper.

## Loading packages

library(fda);library(tseries);library(sandwich);library(sde);
library(variables);library(basefun);library(polynom);
library(fracdiff);library(LongMemoryTS);library(arfima)

## Loading precalculated critical values of the variance-ratio test
setwd("../code/auxiliary")
load("CV_d1_0.5_mean_long_500_new.RData")

CVV=rbind(seq(0.7,1,0.01),MAXCV)
dalpha=0.5

####################################################
## Key simulation paramters and related settings ###
####################################################
SDIM=5    # K or q_{max}

d=0.95    # 1st memory paramter
d2=0.3    # 2nd memory paramter

LRS=2     # LRS=1 : LRS-type method
          # other values : the proposed method
####################################################
####################################################

################
## Simulation ##
################
inner = source("/auxiliary/inprod.R")
lrvar = source("/auxiliary/lr_var_v2_for_fractional.R")

## Setting for paratmeres and basis functions ##
nnbasis=31 ;  nt = 200 ; t = (0:(nt-1))/(nt-1); lbnumber=40
lbb<- Legendre_basis(numeric_var("x", support = c(0, 1)),order = lbnumber)
lb=lbb(t)
LB=matrix(0,nrow=length(t),ncol=lbnumber)
for(i in 2:lbnumber){
  for(j in 1:lbnumber)  {
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber){
  LB[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
}

lbnumber2=24
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}
LBF=cbind(rep(1,length(t)),LBF)
LBF0=LBF

PI0=NULL
for (j in 0:T)
{
if (j==0){aa=1}
if (j==1){aa=1*(d/1)}
if (j>=2){aa=PI0[j]*(j-1+d)/j}
PI0=append(PI0,aa)
}

PI=NULL
for(d in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,-0.5))
{
PItest2=NULL
for (j in 0:T)
{
if (j==0){aa=1}
if (j==1){aa=1*(d/1)}
if (j>=2){aa=PItest2[j]*(j-1+d)/j}
PItest2=append(PItest2,aa)
}
PI=cbind(PI,PItest2)
}

finteg=function(x,d)
{
fdtmp1=NULL
dstar=d*10*(d>=0)+ 11*(d<0)
 for(jji in 1:nrow(x))
   {
    if (jji==1){fdtmp1=rbind(fdtmp1,(PI[1:jji,dstar]*x[jji:1,]))}else{
    fdtmp1=rbind(fdtmp1,colSums(PI[1:jji,dstar]*x[jji:1,]))}
   }
   t(fdtmp1)
}
bdd=0.15;bdd2=0.15;nobs=c(200,350,500,1000)
DD1a=NULL;DD1b=NULL;DD1c=NULL;DD1d=NULL;DD1e=NULL;DD2a=NULL;DD2b=NULL;DD2c=NULL;DD2d=NULL;DD2e=NULL;DD3a=NULL;DD3b=NULL;DD3c=NULL;DD3d=NULL;DD2e=NULL;DD4a=NULL;DD4b=NULL;DD4c=NULL;DD4d=NULL;DD4e=NULL
set.seed(1) ; kernel=2
for (jjjj in 1:2000)
{
yy1=arima.sim(n=T,list(ar=runif(1,-bdd,bdd),ma=runif(1,-bdd,bdd)),sd =1, n.start=100)
yy2=arima.sim(n=T,list(ar=runif(1,-bdd,bdd),ma=runif(1,-bdd,bdd)),sd =1, n.start=100)
yy3=arima.sim(n=T,list(ar=runif(1,-bdd,bdd),ma=runif(1,-bdd,bdd)),sd =1, n.start=100)
yy4=arfima.sim(T, model = list(phi = runif(1,-bdd,bdd), theta = runif(1,-bdd,bdd) , dfrac = d2,  dint = 0))
yy5=arfima.sim(T, model = list(phi = runif(1,-bdd,bdd), theta = runif(1,-bdd,bdd) , dfrac = d2, dint = 0))
y1=NULL
y2=NULL
y3=NULL
y4=NULL
y5=NULL
for(jji in 1:T)
{
y1=append(y1,sum(PI0[1:jji]*yy1[jji:1]))
y2=append(y2,sum(PI0[1:jji]*yy2[jji:1]))
y3=append(y3,sum(PI0[1:jji]*yy3[jji:1]))
}
y4=yy4
y5=yy5

A=diag(runif(length(6:25),-bdd2,bdd2))
B=diag(runif(length(6:25),-bdd2,bdd2))
for (i in 1:length(6:25)){
for (j in 1:length(6:25)){
if (abs(i-j)<=2 & i !=j){A[i,j] = runif(1,-bdd2,bdd2);B[i,j] = runif(1,-bdd2,bdd2)}
}
}

vecy = cbind(rnorm(length(6:25),0,0.97^(0:19)))
epast= cbind(rnorm(length(6:25),0,0.97^(0:19)))
for(jji in 2:(T+100))
{
enow=cbind(rnorm(length(6:25),0,0.97^(0:19)))
vecy=cbind(vecy,A%*%vecy[,jji-1] + enow + B%*%epast) ##d2=0.3
epast=enow
}
vecy=vecy[,101:(T+100)]

set1=sample(1:5,size=5)
rpert=append(set1,sample(setdiff(1:ncol(LBF),set1),size=length(6:ncol(LBF))))
LBF=LBF0[,rpert]

eta = matrix(NA, nrow = nt, ncol = T)
x_mat = eta
mu = -2*(t-1/2)^2 + 0.5
for(i in 1:T){
eta = LBF[,6:ncol(LBF)]%*%vecy[,i]
x_mat[,i]=mu + y1[i]*LBF[,1]+y2[i]*LBF[,2]+y3[i]*LBF[,3]+y4[i]*LBF[,4]+y5[i]*LBF[,5]+ eta
}
	hh2=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
	xcoef=t(hh2)
	xcoef=t(xcoef)-rowMeans(t(xcoef))
    xcoef=t(xcoef)

	hh2=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
	ycoef=t(hh2)
	ycoef=t(ycoef-ycoef[,1])
	ycoef=t(ycoef)

	for (TTT in nobs)
	{

    ddd1=NULL

	for (jijj in 1:20)
	{
	if(LRS==1){input=t(ycoef[1:TTT,]%*%eig_col[,1])} else{input=rbind(ycoef[1:TTT,1]*rnorm(1,1,1)+ycoef[1:TTT,2]*rnorm(1,1,1)+ycoef[1:TTT,3]*rnorm(1,1,1)+ycoef[1:TTT,4]*rnorm(1,1,1)+ycoef[1:TTT,5]*rnorm(1,1,1))}
    dd1=local.W(input,m=floor(1+TTT^0.65),int=c(0.5,1))$d
	ddd1=append(ddd1,dd1)
   	}
    index=which(round(CVV[1,],digits=2)==round(max(max(ddd1),0.7),digits=2))

    sdim=SDIM
	while(sdim>=0)
	{
	if(sdim==1){sss=1;break}
	kdim=sdim+2
	lrx0 = crossprod(xcoef[1:TTT,])
    eig.lrx0 = eigen(lrx0)$vectors
    fdtmp0 = t(xcoef[1:TTT,]%*%eig.lrx0[,1:kdim])   # true=3
	fdtmp00=t(fdtmp0)

    fdtmp1 = finteg(fdtmp00,dalpha)
    lrx0 = crossprod(t(fdtmp0[1:kdim,]))
	lrx1 = crossprod(t(fdtmp1[1:kdim,]))
	EV=sort(eigen(lrx0%*%solve(lrx1),only.values = TRUE)$values)
    test=(TTT^(2*dalpha))*(EV) ;    CV=CVV[2:(sdim+1),index]
    testresult=(test[sdim]< CV[sdim])
	if(testresult==1){sss=sdim;break}else{sdim=sdim-1}
	}



if (TTT==nobs[1]){DD1a = append(DD1a,sss)}
if (TTT==nobs[2]){DD1b = append(DD1b,sss)}
if (TTT==nobs[3]){DD1c = append(DD1c,sss)}
if (TTT==nobs[4]){DD1d = append(DD1d,sss)}


	sdim=SDIM
    lrx0 = crossprod(xcoef[1:TTT,]) ; eig.lrx0 = eigen(lrx0)$vectors ; fdtmp0 = t(xcoef[1:TTT,]%*%eig.lrx0[,1:sdim]);EV=eigen(lrx0)$values[1:sdim]   # true=3
	lEV=length(EV)
	sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])

if (TTT==nobs[1]){DD2a = append(DD2a,sss)}
if (TTT==nobs[2]){DD2b = append(DD2b,sss)}
if (TTT==nobs[3]){DD2c = append(DD2c,sss)}
if (TTT==nobs[4]){DD2d = append(DD2d,sss)}



	sdim=SDIM
	bandw=trunc(1+TTT^(0.3))
	lrx0 = crossprod(xcoef[1:TTT,])
    eig.lrx0 = eigen(lrx0)$vectors
    fdtmp0 = t(xcoef[1:TTT,]%*%eig.lrx0[,4:(4+sdim-1)])
	EV=eigen(lrvar(t(fdtmp0))$omega)$values
	lEV=length(EV)
	sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])
if (TTT==nobs[1]){DD3a = append(DD3a,sss)}
if (TTT==nobs[2]){DD3b = append(DD3b,sss)}
if (TTT==nobs[3]){DD3c = append(DD3c,sss)}
if (TTT==nobs[4]){DD3d = append(DD3d,sss)}

	sdim=SDIM
	bandw=trunc(1+TTT^(2/5))
	lrx0 = crossprod(xcoef[1:TTT,])
    eig.lrx0 = eigen(lrx0)$vectors
    fdtmp0 = t(xcoef[1:TTT,]%*%eig.lrx0[,4:(4+sdim-1)])
	EV=eigen(lrvar(t(fdtmp0))$omega)$values
	lEV=length(EV)
	sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])

if (TTT==nobs[1]){DD4a = append(DD4a,sss)}
if (TTT==nobs[2]){DD4b = append(DD4b,sss)}
if (TTT==nobs[3]){DD4c = append(DD4c,sss)}
if (TTT==nobs[4]){DD4d = append(DD4d,sss)}
	}

	if(jjjj%%10==0){print(c(jjjj,mean(DD1a==3),mean(DD1b==3),mean(DD1c==3),mean(DD1d==3),mean(DD2a==3),mean(DD2b==3),mean(DD2c==3),mean(DD2d==3),mean(DD3a==2),mean(DD3b==2),mean(DD3c==2),mean(DD3d==2),mean(DD4a==2),mean(DD4b==2),mean(DD4c==2),mean(DD4d==2)))}
}

##DD1 :VRratio, DD2:eigenvalue ratio estimator of d, DD3: eigen ratio estimator of d-b with bandwidth T^0.3, DD4: eigen ratio estimator of d-b with bandwidth T^0.4

## Relative frequency of correct estimation ##
round(c(mean(DD1a==3), mean(DD1b==3),mean(DD1c==3),mean(DD1d==3)),digit=3)
round(c(mean(DD2a==3), mean(DD2b==3),mean(DD2c==3),mean(DD2d==3)),digit=3)
round(c(mean(DD3a==2), mean(DD3b==2),mean(DD3c==2),mean(DD3d==2)),digit=3)
round(c(mean(DD4a==2), mean(DD4b==2),mean(DD4c==2),mean(DD4d==2)),digit=3)
