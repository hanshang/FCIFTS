## Description: code for the empirical results reported in the paper
## The results of our proposed methods for the male data are reported in "male" and "female"
## The results of our proposed methods for the male data are reported in "malelrs" and "femalelrs"

setwd("code/auxiliary")
source("inprod.r")
source("lr_var_v2_for_fractional.r")
library(fda);library(tseries);library(sandwich);library(sde);
library(variables);library(basefun);library(polynom);library(geigen);
library(fracdiff);library(LongMemoryTS);library(arfima)

kernel=2

PI=NULL
for(d in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,-0.5))
{
PItest2=NULL
for (j in 0:500)
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
dalpha=0.5

setwd("../")
female_qx = t(matrix(read.table("SWE_lt_female_death.txt", header = TRUE)[,3], 111, 271))
male_qx   = t(matrix(read.table("SWE_lt_male_death.txt", header = TRUE)[,3], 111, 271))

n_col = ncol(female_qx)
n_row = nrow(female_qx)
female_pop = male_pop = matrix(NA, n_row, n_col)
for(ij in 1:n_row)
{
    start_pop_female = start_pop_male = 1
    for(ik in 1:n_col)
    {
        female_pop[ij,ik] = female_qx[ij,ik] * start_pop_female
        start_pop_female = start_pop_female - female_pop[ij,ik]

        male_pop[ij,ik] = male_qx[ij,ik] * start_pop_male
        start_pop_male = start_pop_male - male_pop[ij,ik]
        rm(ik)
    }
    print(ij); rm(ij)
}
colnames(female_pop) = colnames(male_pop) = 0:110
rownames(female_pop) = rownames(male_pop) = 1751:2021
SWE_female_pop = replace(female_qx, which(female_qx == 0), 10^-5)
SWE_male_pop   = replace(male_qx,   which(male_qx == 0), 10^-5)

# Basis approximation
nnbasis=31;nt = length(0:110);t = (0:(nt-1))/(nt-1);neigen = 20 ;lbnumber=40
lbb<- Legendre_basis(numeric_var("x", support = c(0, 1)),order = lbnumber)
lb=lbb(t)
LB=matrix(0,nrow=length(t),ncol=lbnumber)
for(i in 2:lbnumber){
  for(j in 1:lbnumber)  {
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber){
  LB[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
}

#######################################################################
#######################################################################
# MALE DATA
#######################################################################
#######################################################################
############
set.seed(0)
x_mat=t(log(SWE_male_pop))
############
x0_mat=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
xcoef=t(x0_mat)
xcoef=t(xcoef-xcoef[,1])
xcoef=t(xcoef)
x0_mat=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
ycoef=t(x0_mat)
ycoef=t(ycoef)-rowMeans(t(ycoef))
ycoef=t(ycoef)

## Highest memory estimation
TTT=nrow(xcoef)
bw1=0.65
bw2=0.3

LRS=2
ddd1=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])} else{input=rbind(xcoef[1:TTT,1]*rnorm(1,1,1)+xcoef[1:TTT,2]*rnorm(1,1,1)+xcoef[1:TTT,3]*rnorm(1,1,1)+xcoef[1:TTT,4]*rnorm(1,1,1)+xcoef[1:TTT,5]*rnorm(1,1,1))}
dd1=local.W(input,m=floor(1+TTT^bw1),int=c(0.5,1.5))$d
ddd1=append(ddd1,dd1)
}
ddd1=quantile(ddd1,1)
dest=ddd1

LRS=1
ddd1=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])} else{input=rbind(xcoef[1:TTT,1]*rnorm(1,1,1)+xcoef[1:TTT,2]*rnorm(1,1,1)+xcoef[1:TTT,3]*rnorm(1,1,1)+xcoef[1:TTT,4]*rnorm(1,1,1)+xcoef[1:TTT,5]*rnorm(1,1,1))}
dd1=local.W(input,m=floor(1+TTT^bw1),int=c(0.5,1.5))$d
ddd1=append(ddd1,dd1)
}
ddd1=quantile(ddd1,1)
dlrs=ddd1

## VR test
sdim=7
dalpha=0.5
CV=c(15.56287, 25.15453, 34.5383, 43.60106, 52.51703, 61.49654 ,69.60918 ,78.28453, 87.07905, 95.75672, 103.8943, 112.5736, 120.9732 ,129.5989 ,137.7424)  ## precalculated CV for dest

while(sdim>=0)
{
if(sdim==1){sss=1;break}
kdim=sdim+2
lrx0 = crossprod(ycoef[1:TTT,])
eig.lrx0 = eigen(lrx0)$vectors
fdtmp0 = t(ycoef[1:TTT,]%*%eig.lrx0[,1:kdim])   # true=3
fdtmp00=t(fdtmp0)
dalpha=0.5
fdtmp1 = finteg(fdtmp00,dalpha)
lrx0 = crossprod(t(fdtmp0[1:kdim,]))
lrx1 = crossprod(t(fdtmp1[1:kdim,]))
EV=sort(eigen(lrx0%*%solve(lrx1),only.values = TRUE)$values)
test=(TTT^(2*dalpha))*(EV) ;
testresult=(test[sdim]< CV[sdim])
if(testresult==1){sss=sdim;break}else{sdim=sdim-1}
}
dimest=sss

sdim=7
lrx0 = crossprod(ycoef[1:TTT,]) ; eig.lrx0 = eigen(lrx0)$vectors ; fdtmp0 = t(ycoef[1:TTT,]%*%eig.lrx0[,1:sdim]);EV=eigen(lrx0)$values[1:sdim]
lEV=length(EV)
sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])
dimlrs=sss


###
sdim=dimest
set.seed(0)
###
LRS=2
ddd4=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1])} else{eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1]*rnorm(1,1,0)+ycoef[1:TTT,]%*%eig_col[,sdim+2]*rnorm(1,0,1))}
dd4=local.W(input,m=floor(1+TTT^bw1),int=c(0,0.5))$d
ddd4=append(ddd4,dd4)
}
ddd4=quantile(ddd4,1)
dest2=ddd4

set.seed(0)
LRS=1
ddd4=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1])} else{input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1]*rnorm(1,1,0)+ycoef[1:TTT,]%*%eig_col[,sdim+2]*rnorm(1,0,1))}
dd4=local.W(input,m=floor(1+TTT^bw1),int=c(0,0.5))$d
ddd4=append(ddd4,dd4)
}
ddd4=quantile(ddd4,1)
dlrs2=ddd4

bandw=trunc(1+TTT^(bw2))
lrx0 = crossprod(ycoef[1:TTT,])
eig.lrx0 = eigen(lrx0)$vectors
fdtmp0 = t(ycoef[1:TTT,]%*%eig.lrx0[,(sdim+1):(sdim+7)])   # d true=3
# EV=eigen(lrvar(t(fdtmp0))$omega)$values
EV=eigen(lrvar(t(fdtmp0)))$values
lEV=length(EV)
sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])
dimest2=sss

### Results ###
male = c(dest,dimest,dest2,dimest2)
malelrs=c(dlrs,dimlrs,dlrs2)

nnn=ncol(log(SWE_male_pop))
color.gradient <- function(x, colors=c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#8000FF"), colsteps=nnn) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

##Graphs

#\commHS
require(rainbow)
plot(fts(0:110, log(t(SWE_male_pop))), xlab = "Age", ylab = "Log mortality rate", main = "Swedish males")

dpi=600    #pixels per square inch
png("data1.png", width=6*dpi, height=5*dpi, res=dpi)
CLR=color.gradient(1:nnn)
plot(log(SWE_male_pop[1,]),col=CLR[1],type="l",lwd=1,ylim=c(-6,0), xlab = "", ylab = "")
for(i in 2:nnn)
{
lines(log(SWE_male_pop[i,]),col=CLR[i],lwd=1)
}

title(main = "Swedish male central death rates (1751-2021)", sub = NULL, xlab = "Age", ylab = "Log-death rates",cex.lab=1.25)
   dev.off()
   eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])
   input=t(xcoef[1:TTT,]%*%eig_col[,3])
   dpi=600    #pixels per square inch
  png("score1.png", width=6*dpi, height=5*dpi, res=dpi)
  plot(input[1:TTT],type="l", ylab="",xlab="")
 title(#main="main title", sub="sub-title",
 xlab="time", ylab="",cex.lab=1.25)
   dev.off()

     eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])
   input=t(xcoef[1:TTT,]%*%eig_col[,7])
   dpi=600    #pixels per square inch
  png("score2.png", width=6*dpi, height=5*dpi, res=dpi)
  plot(input[1:TTT],type="l", ylab="",xlab="")
 title(#main="main title", sub="sub-title",
 xlab="time", ylab="",cex.lab=1.25)
   dev.off()

      eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])
   input=t(xcoef[1:TTT,]%*%eig_col[,12]  )
   dpi=600    #pixels per square inch
  png("score3.png", width=6*dpi, height=5*dpi, res=dpi)
  plot(input[1:TTT],type="l", ylab="",xlab="")
 title(#main="main title", sub="sub-title",
 xlab="time", ylab="",cex.lab=1.25)
   dev.off()

#######################################################################
#######################################################################
# FEMALE DATA
#######################################################################
#######################################################################
set.seed(0)
x_mat=t(log(SWE_female_pop))
############

x0_mat=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
xcoef=t(x0_mat)
xcoef=t(xcoef-xcoef[,1])
xcoef=t(xcoef)

x0_mat=t(LB[2:(nt),])%*%x_mat[2:(nt),]*(t[2]-t[1])
ycoef=t(x0_mat)
ycoef=t(ycoef)-rowMeans(t(ycoef))
ycoef=t(ycoef)


## Highest memory estimation
TTT=nrow(xcoef)

LRS=2
ddd1=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])} else{input=rbind(xcoef[1:TTT,1]*rnorm(1,1,1)+xcoef[1:TTT,2]*rnorm(1,1,1)+xcoef[1:TTT,3]*rnorm(1,1,1)+xcoef[1:TTT,4]*rnorm(1,1,1)+xcoef[1:TTT,5]*rnorm(1,1,1))}
#input=input[2:TTT]-input[1:(TTT-1)]
dd1=local.W(input,m=floor(1+TTT^bw1),int=c(0.5,1.5))$d
ddd1=append(ddd1,dd1)
}
ddd1=quantile(ddd1,1)
dest=ddd1

LRS=1
ddd1=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])} else{input=rbind(xcoef[1:TTT,1]*rnorm(1,1,1)+xcoef[1:TTT,2]*rnorm(1,1,1)+xcoef[1:TTT,3]*rnorm(1,1,1)+xcoef[1:TTT,4]*rnorm(1,1,1)+xcoef[1:TTT,5]*rnorm(1,1,1))}
#input=input[2:TTT]-input[1:(TTT-1)]
dd1=local.W(input,m=floor(1+TTT^bw1),int=c(0.5,1.5))$d
ddd1=append(ddd1,dd1)
}
ddd1=quantile(ddd1,1)
dlrs=ddd1

## VR test
sdim=7

dalpha=0.5
CV= c(14.59652, 24.05726 ,33.37007 ,41.91549 ,50.79171, 59.07409, 67.68511, 75.69294, 84.15756, 93.04522 ,101.1344, 109.1647, 117.4161, 125.6964, 133.8105) ## precalculated CV for dest

while(sdim>=0)
{
if(sdim==1){sss=1;break}
kdim=sdim+2
lrx0 = crossprod(ycoef[1:TTT,])
eig.lrx0 = eigen(lrx0)$vectors
fdtmp0 = t(ycoef[1:TTT,]%*%eig.lrx0[,1:kdim])   # true=3
fdtmp00=t(fdtmp0)
dalpha=0.5
fdtmp1 = finteg(fdtmp00,dalpha)
lrx0 = crossprod(t(fdtmp0[1:kdim,]))
lrx1 = crossprod(t(fdtmp1[1:kdim,]))
EV=sort(eigen(lrx0%*%solve(lrx1),only.values = TRUE)$values)
test=(TTT^(2*dalpha))*(EV) ;
testresult=(test[sdim]< CV[sdim])
if(testresult==1){sss=sdim;break}else{sdim=sdim-1}
}
dimest=sss

sdim=7
lrx0 = crossprod(ycoef[1:TTT,]) ; eig.lrx0 = eigen(lrx0)$vectors ; fdtmp0 = t(ycoef[1:TTT,]%*%eig.lrx0[,1:sdim]);EV=eigen(lrx0)$values[1:sdim]
lEV=length(EV)
sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])
dimlrs=sss



###
sdim=dimest
set.seed(0)
###
LRS=2
ddd4=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1])} else{eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1]*rnorm(1,1,0)+ycoef[1:TTT,]%*%eig_col[,sdim+2]*rnorm(1,0,1))}
dd4=local.W(input,m=floor(1+TTT^bw1),int=c(0,0.5))$d
ddd4=append(ddd4,dd4)
}
ddd4=quantile(ddd4,1)
dest2=ddd4

set.seed(0)
LRS=1
ddd4=NULL
for (jijj in 1:20)
{
if(LRS==1){eig_col = eigen(crossprod(ycoef[1:TTT,]))$vectors;input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1])} else{input=t(ycoef[1:TTT,]%*%eig_col[,sdim+1]*rnorm(1,1,0)+ycoef[1:TTT,]%*%eig_col[,sdim+2]*rnorm(1,0,1))}
dd4=local.W(input,m=floor(1+TTT^bw1),int=c(0,0.5))$d
ddd4=append(ddd4,dd4)
}
ddd4=quantile(ddd4,1)
dlrs2=ddd4


bandw=trunc(1+TTT^(bw2))
lrx0 = crossprod(ycoef[1:TTT,])
eig.lrx0 = eigen(lrx0)$vectors
fdtmp0 = t(ycoef[1:TTT,]%*%eig.lrx0[,(sdim+1):(sdim+7)])   # d true=3
EV=eigen(lrvar(t(fdtmp0))$omega)$values
lEV=length(EV)
sss=which(max(EV[1:(lEV-1)]/EV[2:(lEV)])==EV[1:(lEV-1)]/EV[2:(lEV)])
dimest2=sss

female = c(dest,dimest,dest2,dimest2)
femalelrs=c(dlrs,dimlrs,dlrs2)



##Graphs
 dpi=600    #pixels per square inch
  png("data2.png", width=6*dpi, height=5*dpi, res=dpi)
 CLR=color.gradient(1:nnn)
plot(log(SWE_female_pop[1,]),col=CLR[1],type="l",lwd=1,ylim=c(-6,0), xlab = "", ylab = "")
for(i in 2:nnn)
{
lines(log(SWE_female_pop[i,]),col=CLR[i],lwd=1)
}
title(main = "Swedish female central death rates (1751-2021)", sub = NULL, xlab = "Age", ylab = "Log-death rates",cex.lab=1.25)
   dev.off()

   eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])
   input=t(xcoef[1:TTT,]%*%eig_col[,3])
   dpi=600    #pixels per square inch
  png("score1a.png", width=6*dpi, height=5*dpi, res=dpi)
  plot(input[1:TTT],type="l", ylab="",xlab="")
 # axis(1,cex.axis=1.35)
 # axis(2,cex.axis=1.35)
 title(#main="main title", sub="sub-title",
 xlab="time", ylab="",cex.lab=1.25)
   dev.off()

     eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])
   input=t(xcoef[1:TTT,]%*%eig_col[,7])
   dpi=600    #pixels per square inch
  png("score2a.png", width=6*dpi, height=5*dpi, res=dpi)
  plot(input[1:TTT],type="l", ylab="",xlab="")
 # axis(1,cex.axis=1.35)
 # axis(2,cex.axis=1.35)
 title(#main="main title", sub="sub-title",
 xlab="time", ylab="",cex.lab=1.25)
   dev.off()

      eig_col = eigen(crossprod(xcoef[1:TTT,]))$vectors;input=t(xcoef[1:TTT,]%*%eig_col[,1])
   input=t(xcoef[1:TTT,]%*%eig_col[,12]  )
   dpi=600    #pixels per square inch
  png("score3a.png", width=6*dpi, height=5*dpi, res=dpi)
  plot(input[1:TTT],type="l", ylab="",xlab="")
 # axis(1,cex.axis=1.35)
 # axis(2,cex.axis=1.35)
 title(#main="main title", sub="sub-title",
 xlab="time", ylab="",cex.lab=1.25)
   dev.off()


# Estimators
male
female
malelrs
femalelrs
