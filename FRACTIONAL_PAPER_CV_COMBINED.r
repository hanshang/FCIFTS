## Description: Critical values of the VR test can be computed from this code
## Specified quantiles for the max-test are stores in MAXCV
## Specified quantiles for the trace-test are stores in TRACECV

####################################################
## Key simulation paramters and related settings ###
####################################################
outeriter=10000      ## number of repetitions
SET=seq(0.95)        ## sets of the values of d (the integration order)
etaquantile=0.95     ## quantile
leng=500             ## sample points of the fractional Brownian motions
dim = 1:10           ## the set of value of q
####################################################
####################################################

MEAN=1

MAXCV=NULL
TRACECV=NULL

#SET=c(0.95)###################################################################


T=leng
DD=cbind(rep(1,T))
d=0.5                     #####################################################
PItest=NULL
for (j in 0:T)
{
if (j==0){aa=1}
if (j==1){aa=1*(d/1)}
if (j>=2){aa=PItest[j]*(j-1+d)/j}
PItest=append(PItest,aa)
}
d000=d

MAXCV=matrix(0,nrow=length(dim),ncol=length(SET))
TRACECV=matrix(0,nrow=length(dim),ncol=length(SET))

for(nnn in 1:length(SET))
{
    d0=SET[nnn]

    DDD=NULL
    DDD2=NULL
    DDDD=NULL
    DDDD2=NULL

    d=d0-1
    fdcoeff=NULL
    for (j in 0:T)
    {
        if (j==0){aa=1}
        if (j==1){aa=1*(d/1)}
        if (j>=2){aa=fdcoeff[j]*(j-1+d)/j}
        fdcoeff=append(fdcoeff,aa)
    }

    CCC=matrix(0,nrow=outeriter,ncol=max(dim))
    CCC2=matrix(0,nrow=outeriter,ncol=max(dim))
    for(jjjij in 1:outeriter)
    {
        aaa=NULL
        for (i in 1:max(dim))
        {
            atem=as.vector(BM(x=0, t0=0, T=1, N=leng-1))
            aaa=rbind(aaa,atem)
        }
        aaa=t(aaa)
        aa=NULL

        for(jji in 1:nrow(aaa))
        {
	         if(max(dim)==1){ if (jji==1){ aa=rbind(aa,(fdcoeff[1:jji]*aaa[jji:1])) }else{
           aa=rbind(aa,sum(fdcoeff[1:jji]*aaa[jji:1])) }} else{
           if (jji==1){ aa=rbind(aa,(fdcoeff[1:jji]*aaa[jji:1,])) }else{
           aa=rbind(aa,colSums(fdcoeff[1:jji]*aaa[jji:1,])) }}
        }
        if (MEAN==1 &max(dim)==1){aa=aa-mean(aa)}
        if (MEAN==1 &max(dim)>1){aa= t(t(aa)-rowMeans(t(aa)))}
        aaa=(aa)
        bb=NULL
        for(jji in 1:nrow(aaa))
        {
	         if(max(dim)==1){ if (jji==1){ bb=rbind(bb,(PItest[1:jji]*aaa[jji:1])) }else{
           bb=rbind(bb,sum(PItest[1:jji]*aaa[jji:1])) }} else{
           if (jji==1){ bb=rbind(bb,(PItest[1:jji]*aaa[jji:1,])) }else{
           bb=rbind(bb,colSums(PItest[1:jji]*aaa[jji:1,])) }}
        }

        for(jdim in dim)
        {
            AA=t(aa[,1:jdim])%*%aa[,1:jdim]
            BB=t(bb[,1:jdim])%*%bb[,1:jdim]
            CC=AA%*%solve(BB)
            EV=eigen(CC)$values
            CCC[jjjij,jdim]=leng^{d000*2}*max(EV)
            CCC2[jjjij,jdim]=leng^{d000*2}*sum(EV)
        }
        if(jjjij%%100==0)
        {
            print(jjjij)
        }
    }
    for (ikj in dim)
    {
        MAXCV[ikj,nnn]=quantile(CCC[,ikj],etaquantile)
        TRACECV[ikj,nnn]=quantile(CCC2[,ikj],etaquantile)
    }
    print(nnn)
}
