#########################################################################################################
# Plot Codes of
# 
# "Does the shape of economic recovery matter? An alternative unit root test with new smooth transition model"
# 
# written by Dr. Mehmet ÖZCAN / www.ozcanmehmet.com / mehmetozcan@kmu.edu.tr 
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
# 4) PART-2 codes draw plots for selected G8 country's data. Threfore, before analyse your data try
#    data of this study. It could be found at "cpi.txt", "ind.txt", and "unemp.txt" files in same repository of Dr. Özcan's GitHub page.
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
#########################################################################################################
# PART - 2
#########################################################################################################
# Aplication
# select working directory folder
# setwd("your folder directory") 
# read data
# the plot codes automatically draw and save .png file into working directory folder

setwd("C:/Users/Ozc/Desktop/Destek")

D1   <- as.matrix(log(read.table("ipi.txt")))
D2   <- as.matrix(read.table("unp.txt"))
D3   <- as.matrix(log(read.table("cpi.txt")))
#########################################################################################################
# IPI Plots
YTR <- matrix(0,nrow(D1),8)
for(i in 1:8){
vrb  <- D1[,i]
TF   <- start0(vrb)
b0a  <- TF$betasA
b0b  <- TF$betasB
b0c  <- TF$betasC
EST  <- est(vrb, b0a, b0b,b0c)
Ytr  <- as.matrix(EST$yhat)
YTR[,i]  <- Ytr[,3]}

graphics.off()
png(filename="IPIs.png",width=1380,height=1920,pointsize=10,res=250)
par(mfrow=c(4,2))
par(cex=0.8, mar=c(6.15,2,1,1))
# Canada
vrb<-D1[,1]
plot(ts(YTR[,1], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Canada IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd =TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# France
vrb<-D1[,2]
plot(ts(YTR[,2], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("France IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# Germany
vrb<-D1[,3]
plot(ts(YTR[,3], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Germany IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# Italy
vrb<-D1[,4]
plot(ts(YTR[,4], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Italy IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# Japan
vrb<-D1[,5]
plot(ts(YTR[,5], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Japan IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# United Kingdom
vrb<-D1[,6]
plot(ts(YTR[,6], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("U. Kingdom IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# United States
vrb<-D1[,7]
plot(ts(YTR[,7], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("U. States IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
# Russia
vrb<-D1[,8]
plot(ts(YTR[,8], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Russia IPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb),1),round(max(vrb),1),0.1), col="lightgray", lty="dotted")
dev.off()
graphics.off()

#########################################################################################################
# UNP Plots
YTR <- matrix(0,nrow(D2),7)
for(i in 1:7){
vrb  <- D2[,i]
TF   <- start0(vrb)
b0a  <- TF$betasA
b0b  <- TF$betasB
b0c  <- TF$betasC
EST  <- est(vrb, b0a, b0b,b0c)
Ytr  <- as.matrix(EST$yhat)
YTR[,i]  <- Ytr[,3]}

# Plots
graphics.off()
png(filename="UNPs.png",width=1380,height=1920,pointsize=10,res=250)
par(mfrow=c(4,2))
par(cex=0.8, mar=c(6.15,2,1,1))
# Canada
vrb<-D2[,1]
plot(ts(YTR[,1], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Canada UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd =TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
# France
vrb<-D2[,2]
plot(ts(YTR[,2], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("France UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
# Germany
vrb<-D2[,3]
plot(ts(YTR[,3], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Germany UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
# Italy
vrb<-D2[,4]
plot(ts(YTR[,4], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Italy UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
# Japan
vrb<-D2[,5]
plot(ts(YTR[,5], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Japan UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
# United Kingdom
vrb<-D2[,6]
plot(ts(YTR[,6], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("U. Kingdom UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
# United States
vrb<-D2[,7]
plot(ts(YTR[,7], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("U. States UNP", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),1), col="lightgray", lty="dotted")
dev.off()
graphics.off()

#########################################################################################################
# CPI Plots
YTR <- matrix(0,nrow(D3),8)
for(i in 1:8){
vrb  <- D3[,i]
TF   <- start0(vrb)
b0a  <- TF$betasA
b0b  <- TF$betasB
b0c  <- TF$betasC
EST  <- est(vrb, b0a, b0b,b0c)
Ytr  <- as.matrix(EST$yhat)
YTR[,i]  <- Ytr[,3]}

# Plots

graphics.off()
png(filename="CPIs.png",width=1380,height=1920,pointsize=10,res=250)
par(mfrow=c(4,2))
par(cex=0.8, mar=c(6.15,2,1,1))
# Canada
vrb<-D3[,1]
plot(ts(YTR[,1], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Canada CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd =TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# France
vrb<-D3[,2]
plot(ts(YTR[,2], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("France CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# Germany
vrb<-D3[,3]
plot(ts(YTR[,3], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Germany CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# Italy
vrb<-D3[,4]
plot(ts(YTR[,4], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Italy CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# Japan
vrb<-D3[,5]
plot(ts(YTR[,5], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Japan CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# United Kingdom
vrb<-D3[,6]
plot(ts(YTR[,6], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("U. Kingdom CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# United States
vrb<-D3[,7]
plot(ts(YTR[,7], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("U. States CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
# Russia
vrb<-D3[,8]
plot(ts(YTR[,8], start=c(2013,9),frequency=12),lty=2,axes=FALSE,lwd=2,col="black",ylab=NA,xlab="Years",ylim=c(min(vrb),max(vrb)))
legend(x = "bottom", inset = c(0, -0.8), legend = c("Russia CPI", "Fitted Transition"), lty = c(1, 2),col = c("gray42", "black"),lwd=2,xpd=TRUE,horiz=TRUE,text.font=2,box.lty=0) 
lines(ts(vrb, start=c(2013,9), frequency=12),lty=1, lwd=2, col="gray42")
axis(1, at =c(2013:2022), labels=c(2013:2022),las=3, gap.axis=0.1)
axis(2, at = seq(round(min(vrb)-0.5,1),round(max(vrb)+0.5,1),0.05),lwd=1.5)
abline(v=(seq(2013,2022,1)), col="lightgray", lty="dotted")
abline(h=seq(round(min(vrb)-0.5,1),round(max(vrb),1)+0.5,0.05), col="lightgray", lty="dotted")
dev.off()
graphics.off()

# END 