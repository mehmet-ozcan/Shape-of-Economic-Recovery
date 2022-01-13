#########################################################################################################
# Empirical Application Codes of
# 
# "Does the shape of economic recovery matter? An alternative unit root test with new smooth transition model"
# 
# written by ...
# 
# 
# 1) Results are printed to the screen.
#
# 2) All steps of application process are explained below. Please check them before run the codes.
#
# *Instructions*
#
# 1) First, please install required R packcages stated in PART-1
# 2) Run all codes in PART-1
# 3) Before run PART-2, please chose working directory folder that contains your "data.txt" file.
# 4) PART-2 codes analyse both the USA's and Turkey's data. Threfore, before analyse your data try
#    data of this study. It could be found at "data.txt" file in same repository of Dr. Özcan's GitHub page.
#
#########################################################################################################
# PART -1 
#########################################################################################################
# Required Packcages 
library(nloptr)
# Required Functions
#########################################################################################################
#########################################################################################################
# Required Functions
#########################################################################################################
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
#########################################################################################################
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
# PART - 2
#########################################################################################################
# Aplication
# select working directory file.
# setwd("your folder directory") 
# read data

setwd("")

D   <- read.table("data.txt")
cai <- D[36:135,1]
fri <- D[36:135,2]
dei <- D[36:135,3]
iti <- D[36:135,4]
jpi <- D[36:135,5]
uki <- D[36:135,6]
usi <- D[36:135,7]
rui <- D[36:135,8]
# logarithms
cai <-as.matrix(log(cai))
fri <-as.matrix(log(fri))
dei <-as.matrix(log(dei))
iti <-as.matrix(log(iti))
jpi <-as.matrix(log(jpi))
uki <-as.matrix(log(uki))
usi <-as.matrix(log(usi))
rui <-as.matrix(log(rui))

# Mm  <-as.matrix(cbind(cai,fri,dei,iti,jpi,uki,usi,rui))
# estimation parameter of transition functions & Unit Root Test Statistics
vrb <- rui

TF   <- start0(vrb)
b0a  <- TF$betasA
b0b  <- TF$betasB
b0c  <- TF$betasC
EST  <- est(vrb, b0a, b0b,b0c)
Ytr   <- as.matrix(EST$yhat)

Crit <- matrix(0,10,2)
for(k in 1:10){
tp <- linear(as.matrix(EST$ehat[,3]),k,0,0);
tf <- linear(vrb,k,1,1);
Crit[k,1] <- tp$bic
Crit[k,2] <- tf$bic
}

tP <- linear(as.matrix(EST$ehat[,3]),which.min(Crit[,1]),0,0)
tF <- linear(vrb,which.min(Crit[,2]),1,1)

# Plots
par(mfrow=c(3,1))
plot(ts(Ytr[,1], start=c(2012,12),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",main="Model A",ylim=c(min(vrb),max(vrb)))
mtext(side=3,text=paste("Industrial Production and Fitted Transition"),line=2.5,lwd=2,cex=1)
lines(ts(vrb, start=c(2012,12), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2012:2021), labels=c(2012:2021),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb),1),round(max(vrb),1),0.1),lwd=1.5)
abline(v=(seq(2012,2021,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
plot(ts(Ytr[,2], start=c(2012,12),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",main="Model B",ylim=c(min(vrb),max(vrb)))
lines(ts(vrb, start=c(2012,12), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2012:2021), labels=c(2012:2021),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb),1),round(max(vrb),1),0.1),lwd=1.5)
abline(v=(seq(2012,2021,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
plot(ts(Ytr[,3], start=c(2012,12),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",main="Model C",ylim=c(min(vrb),max(vrb)))
lines(ts(vrb, start=c(2012,12), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2012:2021), labels=c(2012:2021),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb),1),round(max(vrb),1),0.1),lwd=1.5)
abline(v=(seq(2012,2021,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")

# Table
printsz <- function(SZ){
cat ("\n")
cat ("Unit Root Test Results for Selected I.Production Index Series"                      , "\n")
cat ("\n")
cat ("Sample Period: Sep-2012 - March-2020"                                                 , "\n")
cat ("Sample Size  : 100"                                                                 , "\n")
cat ("\n")
cat ("ADF Test Results:"                                                                  , "\n")
cat ("           Test Statistics:", round(tF$tadf,3)                                      , "\n")
cat ("    Autoregressive Lag (k):", which.min(Crit[,2])                                   , "\n")
cat ("           Critical Values:"                                                        , "\n")
cat ("                      1%: -4.078"                                                   , "\n")
cat ("                      5%: -3.456"                                                   , "\n")
cat ("                     10%: -3.154"                                                   , "\n")
cat ("\n")
cat ("P(alphabeta) Test Results:"                                                         , "\n")
cat ("           Test Statistics:", round(tP$tadf,3)                                      , "\n")
cat ("    Autoregressive Lag (k):", which.min(Crit[,1])                                   , "\n")
cat ("           Critical Values:"                                                        , "\n")
cat ("                      1%: -6.060"                                                   , "\n")
cat ("                      5%: -5.431"                                                   , "\n")
cat ("                     10%: -5.078"                                                   , "\n")
cat ("\n")
}
SZ<-1
printsz(SZ)
EST$BetaC

# end
