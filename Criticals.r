#########################################################################################################
# Critical Value simulation Codes of New Combined Threshold Unit Root Test introduced in
# 
# "Does the shape of economic recovery matter? An alternative unit root test with new smooth transition model"
# 
# written by ...
#
# 1) Results are printed to the screen.
#
# 2) All steps of application process are explained below. Please check them before run the codes.
#
# *Instructions*
#
# 1) First, please install required R packcages stated below.
# 2) Set required inputs in Section -C "Simulation Settings".
# 3) Run all codes
#
#########################################################################################################
#########################################################################################################
# Section - A
# Required Packcages 
library(nloptr)
#########################################################################################################
# Section - B
# Required Functions
# Data Generation Process under Null 
nulldgp<-function(n){
nn<-n
nt<-as.matrix(rnorm(nn))
nt<-as.matrix(nt-mean(nt))
yt<-nt
for(t in 2:nn) {yt[t]<-yt[t-1]+nt[t]}
yt<-yt
}

# Starting value estimation for Models (A, B & C)
start0 <- function(yt){
tt<-length(yt)
tr<-as.matrix(c(1:tt))
cc<-as.matrix(rep(1,tt))

tauGs<-c()
for(kk in 1:tt){ tauGs[kk]<-(kk/tt) }
gamaGs<-seq(0,1,0.1)

SSRA<-matrix(0,length(tauGs),length(gamaGs))
SSRB<-matrix(0,length(tauGs),length(gamaGs))
SSRC<-matrix(0,length(tauGs),length(gamaGs))

for(vk in 1:3){
for(vi in 1:(length(gamaGs))){
for(vj in 1:(length(tauGs))){

    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>tauGs[vj]*tt) It[i] = 1

st <- (1-exp(-(gamaGs[vi])*(tr-(tauGs[vj]*tt))^2))

if(vk==1){X<-as.matrix(cbind(cc,st,(It*st)))}
if(vk==2){X<-as.matrix(cbind(cc,tr,st,(It*st)))}
if(vk==3){X<-as.matrix(cbind(cc,tr,st,(It*st),(tr*It*st)))}

XX<-(MASS:::ginv(t(X)%*%X, tol=1e-25))
betas<-XX%*%(t(X)%*%yt)
Yhat<-t((t(betas)%*%t(X)))
ehat<-(yt-Yhat)
if(vk==1){SSRA[vj,vi]<-sum((ehat)^2)}
if(vk==2){SSRB[vj,vi]<-sum((ehat)^2)}
if(vk==3){SSRC[vj,vi]<-sum((ehat)^2)} 
}}}

for(vs in 1:3){
if(vs==1){posA<-which(SSRA == min(SSRA), arr.ind = TRUE);
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>tauGs[posA[1]]*tt) It[i] = 1
st1 <- (1-exp(-(gamaGs[posA[2]])*(tr-(tauGs[posA[1]]*tt))^2));
X<-as.matrix(cbind(cc,st1,(It*st1)))}
if(vs==2){posB<-which(SSRB == min(SSRB), arr.ind = TRUE);
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>tauGs[posB[1]]*tt) It[i] = 1
st1 <- (1-exp(-(gamaGs[posB[2]])*(tr-(tauGs[posB[1]]*tt))^2));
X<-as.matrix(cbind(cc,tr,st1,(It*st1)))}
if(vs==3){posC<-which(SSRC == min(SSRC), arr.ind = TRUE);
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>tauGs[posC[1]]*tt) It[i] = 1
st1 <- (1-exp(-(gamaGs[posC[2]])*(tr-(tauGs[posC[1]]*tt))^2));
X<-as.matrix(cbind(cc,tr,st1,(It*st1),(tr*It*st1)))}

XX<-(MASS:::ginv(t(X)%*%X, tol=1e-25))

if(vs==1){betasA<-XX%*%(t(X)%*%yt);
betasA<-rbind(betasA,gamaGs[posA[2]],tauGs[posA[1]])}
if(vs==2){betasB<-XX%*%(t(X)%*%yt);
betasB<-rbind(betasB,gamaGs[posB[2]],tauGs[posB[1]])}
if(vs==3){betasC<-XX%*%(t(X)%*%yt);
betasC<-rbind(betasC,gamaGs[posC[2]],tauGs[posC[1]])}
}
list(betasA=betasA, betasB=betasB, betasC=betasC)
}

# Models (A, B & C) Estimation with Sequential Quadratic Programing
est<-function(yt, b0a, b0b, b0c){
tt<-length(yt)
tr<-as.matrix(c(1:tt))
cc<-as.matrix(rep(1,tt))
X   <- as.matrix(cbind(cc,tr))

stabA<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[5]*tt) It[i] = 1

    s1t <- (1-exp(-(b[4])*(tr-(b[5]*tt))^2))
    yth <- (b[1]*X[,1])+(b[2]*s1t)+(b[3]*It*s1t)

    sum((yt-yth)^2)
}

stabB<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[6]*tt) It[i] = 1

    s1t <- (1-exp(-(b[5])*(tr-(b[6]*tt))^2))
    yth <- (b[1]*X[,1])+(b[2]*X[,2])+(b[3]*s1t)+(b[4]*It*s1t)

    sum((yt-yth)^2)
}

stabC<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[7]*tt) It[i] = 1

    s1t <- (1-exp(-(b[6])*(tr-(b[7]*tt))^2))
    yth <- (b[1]*X[,1])+(b[2]*X[,2])+(b[3]*s1t)+(b[4]*It*s1t)+(b[5]*X[,2]*It*s1t)

    sum((yt-yth)^2)
}

lba<-c(-Inf,-Inf,-Inf,0,0)
uba<-c(Inf,Inf,Inf,Inf,1)

lbb<-c(-Inf,-Inf,-Inf,-Inf,0,0)
ubb<-c(Inf,Inf,Inf,Inf,Inf,1)

lbc<-c(-Inf,-Inf,-Inf,-Inf,-Inf,0,0)
ubc<-c(Inf,Inf,Inf,Inf,Inf,Inf,1)

SA<-nloptr:::slsqp(b0a,fn=stabA,lower=lba,upper=uba, control=list(xtol_rel=1e-25))
SB<-nloptr:::slsqp(b0b,fn=stabB,lower=lbb,upper=ubb, control=list(xtol_rel=1e-25))
SC<-nloptr:::slsqp(b0c,fn=stabC,lower=lbc,upper=ubc, control=list(xtol_rel=1e-25))

staba<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[5]*tt) It[i] = 1

    s1t <- (1-exp(-(b[4])*(tr-(b[5]*tt))^2))
    (b[1]*X[,1])+(b[2]*s1t)+(b[3]*It*s1t)
}
yhatA<-staba(SA$par)

stabb<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[6]*tt) It[i] = 1

    s1t <- (1-exp(-(b[5])*(tr-(b[6]*tt))^2))
    (b[1]*X[,1])+(b[2]*X[,2])+(b[3]*s1t)+(b[4]*It*s1t)
}
yhatB<-stabb(SB$par)

stabc<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[7]*tt) It[i] = 1

    s1t <- (1-exp(-(b[6])*(tr-(b[7]*tt))^2))
    (b[1]*X[,1])+(b[2]*X[,2])+(b[3]*s1t)+(b[4]*It*s1t)+(b[5]*X[,2]*It*s1t)
}
yhatC<-stabc(SC$par)

list(ehat=as.matrix(cbind((yt-yhatA),(yt-yhatB),(yt-yhatC))), yhat=as.matrix(cbind(yhatA,yhatB,yhatC)),BetaA=as.matrix(SA$par),BetaB=as.matrix(SB$par),BetaC=as.matrix(SC$par))
}

# ADF Estimator
linear<-function(dat,p,cons,trend){
    t <- nrow(dat)
    n <- t - 1 - p
    dy <- as.matrix(dat[2:t]-dat[1:(t-1)])
    y  <- as.matrix(dy[(1+p):(t-1)])
    if(p!=0) {x  <- cbind(dat[(1+p):(t-1)], dy[p:(t-2)])} else {x  <- cbind(dat[(1+p):(t-1)])}

    if(p>=2){for (k in 2:p) x <- cbind(x,dy[(p+1-k):(t-1-k)])}
    if(cons==1 & trend==0)  x <- cbind(rep(1,n), x)
    if(cons==1 & trend==1)  x <- cbind(rep(1,n),seq(1,n,1),x)

    kx     <- ncol(x)
    xx     <- solve(t(x)%*%x, tol=1e-25)
    bols   <- xx%*%(t(x)%*%y) 
    eols   <- y - x%*%bols
    sols   <- t(eols)%*%eols
    sigols <- sols/(n-kx)
    seols  <- sqrt(diag(xx)%*%sigols)
    rho    <- bols[1+cons+trend]
    tadf   <- rho/seols[1+cons+trend]
    ar0    <- as.matrix(bols[(cons+trend+2):(cons+trend+1+p)])
    bic<-log((sum(eols^2))/t)+(length(bols)*(log(t)/t))
    aic<-log((sum(eols^2))/t)+(length(bols)*(2/t))
    list(ar0=ar0,eols=eols,aic=aic,bic=bic,tadf=tadf,rho=rho)
}
#########################################################################################################
# Section - C
# Simulation Settings
rep      <- 10000
n        <- 500
p        <- 0

ptm <- proc.time()
set.seed(01210)

TDist <- matrix(0,rep,3)

for(i in 1:rep){
yt   <- as.matrix(nulldgp(n))
TF   <- start0(yt)
b0a  <- TF$betasA
b0b  <- TF$betasB
b0c  <- TF$betasC
EST  <- est(yt, b0a, b0b,b0c)
D    <- as.matrix(EST$ehat)

tA <- linear(as.matrix(D[,1]),p,0,0);
tB <- linear(as.matrix(D[,2]),p,0,0);
tC <- linear(as.matrix(D[,3]),p,0,0);

TDist[i,1] <- tA$tadf
TDist[i,2] <- tB$tadf
TDist[i,3] <- tC$tadf
}

TDist <- as.matrix(cbind(sort(TDist[,1],decreasing=TRUE),
sort(TDist[,2],decreasing=TRUE),sort(TDist[,3],decreasing=TRUE)))

CR      <- matrix(0,3,3)
for(i in 1:3){
if(i==1) crt <- 0.90
if(i==2) crt <- 0.95
if(i==3) crt <- 0.99
CR[i,1] <- TDist[,1][round(rep*crt)]
CR[i,2] <- TDist[,2][round(rep*crt)]
CR[i,3] <- TDist[,3][round(rep*crt)]
}

proc.time() - ptm
printres <- function(CR){
cat ("\n")
cat ("Critical Values for Selected Sample Size"                                           , "\n")
cat ("Number of Observations (T):", n                                                     , "\n")
cat ("\n")
cat ("   ","Model A"," ","Model B"," ","Model C"                                          , "\n")
cat ("   ","t stat","  ","t stat","  ","t stat"                                           , "\n")
cat ("10%", round(CR[1,1],3)," ",round(CR[1,2],3),"  ",round(CR[1,3],3)                   , "\n")
cat (" 5%", round(CR[2,1],3)," ",round(CR[2,2],3),"  ",round(CR[2,3],3)                   , "\n")
cat (" 1%", round(CR[3,1],3)," ",round(CR[3,2],3),"  ",round(CR[3,3],3)                   , "\n")
cat ("\n")
}

printres(CR)
# END #
