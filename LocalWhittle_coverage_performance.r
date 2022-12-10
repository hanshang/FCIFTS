## Description: Simulation codes for local Whittle estimation of d and d-b.
## The code can produce the results given in Tables 5 and 10 of the paper (in order to produce Table 10 from this code, only a slight modification is required)

## Loading packages
library(fda);library(tseries);library(sandwich);library(sde);
library(variables);library(basefun);library(polynom);
library(fracdiff);library(LongMemoryTS);library(arfima)

setwd("../code/auxiliary")
source("inprod.R")

####################################################
## Key simulation paramters and related settings ###
####################################################
d=0.95    # 1st memory paramter
d2=0.3    # 2nd memory paramter

bw1=0.65  # Bandwidth for the local Whittle estimation

LRS=2     # LRS=1 : LRS-type method
          # other values : the proposed method
####################################################
####################################################

interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    return(c(mean(score), cpd))
}

nnbasis=31;nt = 200;t = (0:(nt-1))/(nt-1);neigen = 20;lbnumber=40
lbb<- Legendre_basis(numeric_var("x", support = c(0, 1)),order = lbnumber)
lb=lbb(t)
LB=matrix(0,nrow=length(t),ncol=lbnumber)
for(i in 2:lbnumber)
{
    for(j in 1:lbnumber)
    {
        if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

        for(i in 1:lbnumber)
        {
          LB[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
        }
        lbnumber2=24
        LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
        for (i in 1:(lbnumber2/2))
        {
            LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
            LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
        }
        LBF=cbind(rep(1,length(t)),LBF)
        LBF0=LBF
        T=1000
    PI0=NULL
    for (j in 0:T)
    {
    if (j==0){aa=1}
    if (j==1){aa=1*(d/1)}
    if (j>=2){aa=PI0[j]*(j-1+d)/j}
    PI0=append(PI0,aa)
    }

    PI=NULL
    for (j in 0:T)
    {
    if (j==0){aa=1}
    if (j==1){aa=1*(d2/1)}
    if (j>=2){aa=PI[j]*(j-1+d2)/j}
    PI=append(PI,aa)
    }

    bdd=0.15;bdd2=0.15;nobs=c(200,350,500,1000)
    DD1a=NULL;DD1b=NULL;DD1c=NULL;DD1d=NULL;DD1e=NULL;DD1aa=NULL ;DD1bb=NULL;DD1cc=NULL;DD1dd=NULL;DD4a=NULL;DD4b=NULL;DD4c=NULL;DD4d=NULL;DD4e=NULL;DD4aa=NULL;DD4bb=NULL;DD4cc=NULL;DD4dd=NULL

    set.seed(999); kernel=2
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
    if (abs(i-j)<=1 & i !=j){A[i,j] = runif(1,-bdd2,bdd2);B[i,j] = runif(1,-bdd2,bdd2)}
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
    	xcoef=t(xcoef-xcoef[,1])
    	xcoef=t(xcoef)

    	hh2=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
    	ycoef=t(hh2)
    	ycoef=t(ycoef)-rowMeans(t(ycoef))
        ycoef=t(ycoef)

    	for (TTT in nobs)
    	{
    	ddd1=NULL;	ddd11=NULL;	ddd4=NULL;	ddd44=NULL;

    	for (jijj in 1:20)
    	{
    	if(LRS==1){eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors; input=t(xcoef[1:TTT,]%*%eig_col[,1])} else{eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=rbind(xcoef[1:TTT,1]*rnorm(1,1,1)+xcoef[1:TTT,2]*rnorm(1,1,1)+xcoef[1:TTT,3]*rnorm(1,1,1)+xcoef[1:TTT,4]*rnorm(1,1,1)+xcoef[1:TTT,5]*rnorm(1,1,1))}
        dd1=local.W(input[2:TTT]-input[1:(TTT-1)],m=floor(1+TTT^bw1),int=c(-0.5,0.5))$d
    	ddd1=append(ddd1,dd1)
    	if(LRS==1){eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,4])} else{eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,4]*rnorm(1,1,0)+ycoef[1:TTT,]%*%eig_col[,5]*rnorm(1,0,1))}
    	dd4=local.W(input,m=floor(1+TTT^bw1),int=c(0,0.5))$d
    	ddd4=append(ddd4,dd4)
    	}
    	ddd1=quantile(ddd1,1)
    	ddd4=quantile(ddd4,1)

    	mm=floor(1+TTT^bw1)
    	cover1=(-0.05 < ddd1+1.96*(1/sqrt(4*mm)))*(-0.05 > ddd1-1.96*(1/sqrt(4*mm)))
    	cover2=(0.3< ddd4+1.96*(1/sqrt(4*mm)))*(0.3 > ddd4-1.96*(1/sqrt(4*mm)))

    if (TTT==nobs[1]){DD1aa= append(DD1aa,ddd1); DD4aa = append(DD4aa,ddd4)}
    if (TTT==nobs[2]){DD1bb = append(DD1bb,ddd1); DD4bb = append(DD4bb,ddd4)}
    if (TTT==nobs[3]){DD1cc = append(DD1cc,ddd1); DD4cc = append(DD4cc,ddd4)}
    if (TTT==nobs[4]){DD1dd = append(DD1dd,ddd1); DD4dd = append(DD4dd,ddd4)}


    if (TTT==nobs[1]){DD1a = append(DD1a,cover1); DD4a = append(DD4a,cover2)}
    if (TTT==nobs[2]){DD1b = append(DD1b,cover1); DD4b = append(DD4b,cover2)}
    if (TTT==nobs[3]){DD1c = append(DD1c,cover1); DD4c = append(DD4c,cover2)}
    if (TTT==nobs[4]){DD1d = append(DD1d,cover1); DD4d = append(DD4d,cover2)}
    	}
	print(c(jjjj,mean(DD1a),mean(DD1b),mean(DD1c),mean(DD1d),mean(DD4a),mean(DD4b),mean(DD4c),mean(DD4d)))
}

## Coverage rates

round(c(mean(DD1a),mean(DD1b),mean(DD1c),mean(DD1d)),digits=4)
round(c(mean(DD4a),mean(DD4b),mean(DD4c),mean(DD4d)),digits=4)

## Interval scores

MM=floor(1+nobs^bw1)
round(interval_score(-DD1aa-0.05,-1.96*(1/sqrt(4*MM[1])),1.96*(1/sqrt(4*MM[1])),0.05),digits=4)
round(interval_score(-DD1bb-0.05,-1.96*(1/sqrt(4*MM[2])),1.96*(1/sqrt(4*MM[2])),0.05),digits=4)
round(interval_score(-DD1cc-0.05,-1.96*(1/sqrt(4*MM[3])),1.96*(1/sqrt(4*MM[3])),0.05),digits=4)
round(interval_score(-DD1dd-0.05,-1.96*(1/sqrt(4*MM[4])),1.96*(1/sqrt(4*MM[4])),0.05),digits=4)

round(interval_score(-DD4aa+0.3,-1.96*(1/sqrt(4*MM[1])),1.96*(1/sqrt(4*MM[1])),0.05),digits=4)
round(interval_score(-DD4bb+0.3,-1.96*(1/sqrt(4*MM[2])),1.96*(1/sqrt(4*MM[2])),0.05),digits=4)
round(interval_score(-DD4cc+0.3,-1.96*(1/sqrt(4*MM[3])),1.96*(1/sqrt(4*MM[3])),0.05),digits=4)
round(interval_score(-DD4dd+0.3,-1.96*(1/sqrt(4*MM[4])),1.96*(1/sqrt(4*MM[4])),0.05),digits=4)
