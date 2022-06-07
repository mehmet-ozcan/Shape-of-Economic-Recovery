#########################################################################################################
# Empirical Power simulation Codes of New Combined Threshold Unit Root Test introduced in
# 
# "Does the shape of economic recovery matter? An alternative unit root test with new smooth transition model"
# 
# written by ...
#
# 1) Results are printed to the screen.
#
# 2) All steps of simulation process are explained below. Please check them before run the codes.
#
# *Instructions*
#
# 1) First, please install required R packcages stated below.
# 2) You could simulate different parameter settings. Please check "simulation settings" part blow.
# 3) This codes includes critical values for T=100 and T=200.
# 4) Run all codes
#
#########################################################################################################
# Section - A
# Required Packcages 
library(nloptr)
#########################################################################################################
# Section - B
# Required Functions
#########################################################################################################

# Data Generation Process for Power Simulations  
idgp<-function(n,ph,a2,gm,tu){
et<-as.matrix(rnorm(n))
et<-as.matrix(et-mean(et))
wt<-matrix(0,n,1)
for(t in 2:n) {wt[t]<-(ph*wt[t-1]+et[t])}
It  <- matrix(0,n,1)
for(i in 1:n) if(i>tu*n) It[i] = 1
t<-c(1:n)
z0 <- (1-exp(-gm*(t-tu*n)^2))
yt <- as.matrix(1+(1*z0)+(a2*It*z0)+wt)
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

for(vi in 1:(length(gamaGs))){
for(vj in 1:(length(tauGs))){

    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>tauGs[vj]*tt) It[i] = 1

st <- (1-exp(-(gamaGs[vi])*(tr-(tauGs[vj]*tt))^2))

X<-as.matrix(cbind(cc,st,(It*st)))

XX<-(MASS:::ginv(t(X)%*%X, tol=1e-25))
betas<-XX%*%(t(X)%*%yt)
Yhat<-t((t(betas)%*%t(X)))
ehat<-(yt-Yhat)
SSRA[vj,vi]<-sum((ehat)^2)
}}

posA<-which(SSRA == min(SSRA), arr.ind = TRUE);
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>tauGs[posA[1]]*tt) It[i] = 1
st1 <- (1-exp(-(gamaGs[posA[2]])*(tr-(tauGs[posA[1]]*tt))^2));
X<-as.matrix(cbind(cc,st1,(It*st1)))

XX<-(MASS:::ginv(t(X)%*%X, tol=1e-25))

betasA<-XX%*%(t(X)%*%yt);
betasA<-rbind(betasA,gamaGs[posA[2]],tauGs[posA[1]])

list(betasA=betasA)
}

#########################################################################################################
# Models (A, B & C) Estimation with Sequential Quadratic Programing
est<-function(yt, b0a){
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

lba<-c(-Inf,-Inf,-Inf,0,0)
uba<-c(Inf,Inf,Inf,Inf,1)

SA<-nloptr:::slsqp(b0a,fn=stabA,lower=lba,upper=uba, control=list(xtol_rel=1e-25))

staba<-function(b){
    It  <- matrix(0,tt,1)
    for(i in 1:tt) if(i>b[5]*tt) It[i] = 1

    s1t <- (1-exp(-(b[4])*(tr-(b[5]*tt))^2))
    (b[1]*X[,1])+(b[2]*s1t)+(b[3]*It*s1t)
}
yhatA<-staba(SA$par)

list(ehat=as.matrix(cbind((yt-yhatA))), yhat=as.matrix(cbind(yhatA)),BetaA=as.matrix(SA$par))
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
#########################################################################################################
rep      <- 10000
n        <- 100
p        <- 4

ph       <- 0.5
a2       <- 10
gm       <- 0.1
tu       <- 0.8
###########################################################
# 5% Criticals T=100 ######################################
crA <- -4.582
adf <- -3.456
# 5% Criticals T=200 ######################################
#crA <- -4.258
#adf <- -3.433
###########################################################
# Simulation Process

ptm <- proc.time()
set.seed(01210)

TDist <- matrix(0,rep,2)

for(i in 1:rep){
yt <- as.matrix(idgp(n,ph,a2,gm,tu))
TF   <- start0(yt)
b0a  <- TF$betasA
EST  <- est(yt, b0a)
tp <- linear(as.matrix(EST$ehat),p,0,0);
tf <- linear(yt,p,1,1);

TDist[i,1] <- sum(tp$tadf>crA)
TDist[i,2] <- sum(tf$tadf>adf)
}

sA <- (100-(sum(TDist[,1])*100)/rep)
sB <- (100-(sum(TDist[,2])*100)/rep)

SZ <- rbind(sA,sB)
proc.time() - ptm


printsz <- function(SZ){
cat ("\n")
cat ("Power values for Selected Sample Size"                                              , "\n")
cat ("\n")
cat ("              Nominal Size: 5%"                                                     , "\n")
cat ("Number of Observations (T):", n                                                     , "\n")
cat ("               Replication:", rep                                                   , "\n")
cat ("    Autoregressive Lag (k):", p                                                     , "\n")
cat ("                       Tau:", tu                                                    , "\n")
cat ("                   Alpha 2:", a2                                                    , "\n")
cat ("                     Gamma:", gm                                                    , "\n")
cat ("                       Phi:", ph                                                    , "\n")
cat ("\n")
cat ("   ","Model A Pa"," ","ADF Tt"                                                      , "\n")
cat ("   ", round(SZ[1],3),"       ",round(SZ[2],3)                                       , "\n")
cat ("\n")
}

printsz(SZ)
# END #
