###
# Functions for running the analyses
# Zane Billings
# 2022-12-12
###

# Log Likelihood function for censored trivariate normal with mean function Az
# (where A is 3 x 4 matrix of coefficients on intercept, sex, race, age)
Maxing <- function(X,mn,m1,m2,m3,m12,m13,m23,ma,dmatrix) {
	like<-function(t) {
		fc<-rep(0,n)
		sig123 <- matrix(c(t[13],t[14],t[15],t[14],t[16],t[17],t[15],t[17],t[18]),3,3)
		sig12 <- matrix(c(t[13],t[14],t[14],t[16]),2,2)
		sig13 <- matrix(c(t[13],t[15],t[15],t[18]),2,2)
		sig23 <- matrix(c(t[16],t[17],t[17],t[18]),2,2)
		if(mn>0) fc[1:mn] <- -log(dnorm(X[1:mn,2], mean=(t[1]+t[2]*X[1:mn,5]+t[3]*X[1:mn,6]+t[4]*X[1:mn,7]+as.double(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%t((X[1:mn,3:4]-matrix(c(t[5]+t[6]*X[1:mn,5]+t[7]*X[1:mn,6]+t[8]*X[1:mn,7],t[9]+t[10]*X[1:mn,5]+t[11]*X[1:mn,6]+t[12]*X[1:mn,7]),mn,2))))),sd=sqrt((t[13]-as.double(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%matrix(c(t[14],t[15]),2,1)))))*dnorm(X[1:mn,3],mean=(t[5]+t[6]*X[1:mn,5]+t[7]*X[1:mn,6]+t[8]*X[1:mn,7]+t[17]/t[18]*(X[1:mn,4]-t[9]-t[10]*X[1:mn,5]-t[11]*X[1:mn,6]-t[12]*X[1:mn,7])),sd=sqrt(t[16]-t[17]^2/t[18]))*dnorm(X[1:mn,4],mean=(t[9]+t[10]*X[1:mn,5]+t[11]*X[1:mn,6]+t[12]*X[1:mn,7]),sd=sqrt(t[18])))	
		if(m1>0) fc[(mn+1):(mn+m1)] <- -log(dnorm(X[(mn+1):(mn+m1),3],mean=(t[5]+t[6]*X[(mn+1):(mn+m1),5]+t[7]*X[(mn+1):(mn+m1),6]+t[8]*X[(mn+1):(mn+m1),7]+t[17]/t[18]*(X[(mn+1):(mn+m1),4]-t[9]-t[10]*X[(mn+1):(mn+m1),5]-t[11]*X[(mn+1):(mn+m1),6]-t[12]*X[(mn+1):(mn+m1),7])),sd=sqrt(t[16]-t[17]^2/t[18]))*dnorm(X[(mn+1):(mn+m1),4],mean=(t[9]+t[10]*X[(mn+1):(mn+m1),5]+t[11]*X[(mn+1):(mn+m1),6]+t[12]*X[(mn+1):(mn+m1),7]),sd=sqrt(t[18]))*pnorm(d[1],mean=(t[1]+t[2]*X[(mn+1):(mn+m1),5]+t[3]*X[(mn+1):(mn+m1),6]+t[4]*X[(mn+1):(mn+m1),7]+as.double(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%t((X[(mn+1):(mn+m1),3:4]-matrix(c(t[5]+t[6]*X[(mn+1):(mn+m1),5]+t[7]*X[(mn+1):(mn+m1),6]+t[8]*X[(mn+1):(mn+m1),7],t[9]+t[10]*X[(mn+1):(mn+m1),5]+t[11]*X[(mn+1):(mn+m1),6]+t[12]*X[(mn+1):(mn+m1),7]),(m1),2))))),sd=sqrt((t[13]-as.double(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%matrix(c(t[14],t[15]),2,1))))))
		if(m2>0) fc[(mn+m1+1):(mn+m1+m2)] <- -log(dnorm(X[(mn+m1+1):(mn+m1+m2),2],mean=(t[1]+t[2]*X[(mn+m1+1):(mn+m1+m2),5]+t[3]*X[(mn+m1+1):(mn+m1+m2),6]+t[4]*X[(mn+m1+1):(mn+m1+m2),7]+t[15]/t[18]*(X[(mn+m1+1):(mn+m1+m2),4]-t[9]-t[10]*X[(mn+m1+1):(mn+m1+m2),5]-t[11]*X[(mn+m1+1):(mn+m1+m2),6]-t[12]*X[(mn+m1+1):(mn+m1+m2),7])),sd=sqrt(t[13]-t[15]^2/t[18]))*dnorm(X[(mn+m1+1):(mn+m1+m2),4],mean=(t[9]+t[10]*X[(mn+m1+1):(mn+m1+m2),5]+t[11]*X[(mn+m1+1):(mn+m1+m2),6]+t[12]*X[(mn+m1+1):(mn+m1+m2),7]),sd=sqrt(t[18]))*pnorm(d[2],mean=(t[5]+t[6]*X[(mn+m1+1):(mn+m1+m2),5]+t[7]*X[(mn+m1+1):(mn+m1+m2),6]+t[8]*X[(mn+m1+1):(mn+m1+m2),7]+as.double(matrix(c(t[14],t[17]),1,2)%*%solve(sig13)%*%t((X[(mn+m1+1):(mn+m1+m2),c(2,4)]-matrix(c(t[1]+t[2]*X[(mn+m1+1):(mn+m1+m2),5]+t[3]*X[(mn+m1+1):(mn+m1+m2),6]+t[4]*X[(mn+m1+1):(mn+m1+m2),7],t[9]+t[10]*X[(mn+m1+1):(mn+m1+m2),5]+t[11]*X[(mn+m1+1):(mn+m1+m2),6]+t[12]*X[(mn+m1+1):(mn+m1+m2),7]),m2,2))))),sd=sqrt((t[16]-as.double(matrix(c(t[14],t[17]),1,2)%*%solve(sig13)%*%matrix(c(t[14],t[17]),2,1))))))
		if(m3>0) fc[(mn+m1+m2+1):(mn+m1+m2+m3)] <- -log(dnorm(X[(mn+m1+m2+1):(mn+m1+m2+m3),2],mean=(t[1]+t[2]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[3]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[4]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]+t[14]/t[16]*(X[(mn+m1+m2+1):(mn+m1+m2+m3),3]-t[5]-t[6]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]-t[7]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]-t[8]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7])),sd=sqrt(t[13]-t[14]^2/t[16]))*dnorm(X[(mn+m1+m2+1):(mn+m1+m2+m3),3],mean=(t[5]+t[6]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[7]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[8]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]),sd=sqrt(t[16]))*pnorm(d[3],mean=(t[9]+t[10]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[11]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[12]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]+as.double(matrix(c(t[15],t[17]),1,2)%*%solve(sig12)%*%t((X[(mn+m1+m2+1):(mn+m1+m2+m3),2:3]-matrix(c(t[1]+t[2]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[3]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[4]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7],t[5]+t[6]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[7]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[8]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]),m3,2))))),sd=sqrt((t[18]-as.double(matrix(c(t[15],t[17]),1,2)%*%solve(sig12)%*%matrix(c(t[15],t[17]),2,1))))))
		if(m12>0) fc[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12)] <- -log(dnorm(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),4],mean=(t[9]+t[10]*X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),5]+t[11]*X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),6]+t[12]*X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),7]),sd=sqrt(t[18]))*apply(cbind(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),],dmatrix[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),1:2]),1,Int1,Lower=c(-Inf,-Inf),t=t,sig12=sig12))
		if(m13>0) fc[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13)] <- -log(dnorm(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),3],mean=(t[5]+t[6]*X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),5]+t[7]*X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),6]+t[8]*X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),7]),sd=sqrt(20))*apply(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),],1,Int2,Lower=c(-Inf,-Inf),Upper=c(d[1],d[3]),t=t,sig13=sig13))
		if(m23>0) fc[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23)] <- -log(dnorm(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),2],mean=(t[1]+t[2]*X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),5]+t[3]*X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),6]+t[4]*X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),7]),sd=sqrt(t[13]))*apply(cbind(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),],dmatrix[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),2:3]),1,Int3,Lower=c(-Inf,-Inf),t=t,sig23=sig23))
		if(ma>0) fc[(mn+m1+m2+m3+m12+m13+m23+1):n] <- -log(apply(cbind(X[(mn+m1+m2+m3+m12+m13+m23+1):n,],dmatrix[(mn+m1+m2+m3+m12+m13+m23+1):n,]),1,Int4,t=t,sig123=sig123,Lower=c(-Inf,-Inf,-Inf)))
		like <- sum(fc)
	}
	like
}

#functions to calculate cdf in order to be able to use *apply* function, eliminating the need for a loop which slows maximization
Int1 <- function(X,t,sig12,Lower) sadmvn(lower=Lower,upper=X[8:9],mean=(c(t[1]+t[2]*X[5]+t[3]*X[6]+t[4]*X[7],t[5]+t[6]*X[5]+t[7]*X[6]+t[8]*X[7])+as.vector(matrix(c(t[15],t[17]),2,1)*1/t[18]*(X[4]-t[9]-t[10]*X[5]-t[11]*X[6]-t[12]*X[7]))),varcov=make.positive.definite(sig12-1/t[18]*matrix(c(t[15],t[17]),2,1)%*%t(matrix(c(t[15],t[17]),2,1))),abseps=0.000000001)
Int2 <- function(X,t,sig13,Lower,Upper) sadmvn(lower=Lower,upper=Upper,mean=(c(t[1]+t[2]*X[5]+t[3]*X[6]+t[4]*X[7],t[9]+t[10]*X[5]+t[11]*X[6]+t[12]*X[7])+as.vector(matrix(c(t[14],t[17]),2,1)*1/t[16]*(X[3]-t[5]-t[6]*X[5]-t[7]*X[6]-t[8]*X[7]))),varcov=make.positive.definite(sig13-1/t[16]*matrix(c(t[14],t[17]),2,1)%*%t(matrix(c(t[14],t[17]),2,1))),abseps=0.000000001)
Int3 <- function(X,t,sig23,Lower) sadmvn(lower=Lower,upper=X[8:9],mean=(c(t[5]+t[6]*X[5]+t[7]*X[6]+t[8]*X[7],t[9]+t[10]*X[5]+t[11]*X[6]+t[12]*X[7])+as.vector(matrix(c(t[14],t[15]),2,1)*1/t[13]*(X[2]-t[1]-t[2]*X[5]-t[3]*X[6]-t[4]*X[7]))),varcov=make.positive.definite(sig23-1/t[13]*matrix(c(t[14],t[15]),2,1)%*%t(matrix(c(t[14],t[15]),2,1))),abseps=0.000000001)
Int4 <- function(X,t,sig123,Lower) sadmvn(lower=Lower,upper=X[8:10],mean=c(t[1]+t[2]*X[5]+t[3]*X[6]+t[4]*X[7],t[5]+t[6]*X[5]+t[7]*X[6]+t[8]*X[7],t[9]+t[10]*X[5]+t[11]*X[6]+t[12]*X[7]),varcov=make.positive.definite(sig123),abseps=0.000000001)


#Function to compute weighted maximum likelihood estimates for gamma parameters
EMmaxing <- function(X,mn,EM) {
	like<-function(t) {
		sig23 <- matrix(c(t[16],t[17],t[17],t[18]),2,2)
		fc<-rep(0,(mn+EM*(n-mn)))
		if(mn>0) fc[1:mn] <-  -log(dnorm(X[1:mn,2], mean=(t[1]+t[2]*X[1:mn,5]+t[3]*X[1:mn,6]+t[4]*X[1:mn,7]+as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%t((X[1:mn,3:4]-matrix(c(t[5]+t[6]*X[1:mn,5]+t[7]*X[1:mn,6]+t[8]*X[1:mn,7],t[9]+t[10]*X[1:mn,5]+t[11]*X[1:mn,6]+t[12]*X[1:mn,7]),mn,2))))),sd=sqrt((t[13]-as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%matrix(c(t[14],t[15]),2,1)))))*dnorm(X[1:mn,3],mean=(t[5]+t[6]*X[1:mn,5]+t[7]*X[1:mn,6]+t[8]*X[1:mn,7]+t[17]/t[18]*(X[1:mn,4]-t[9]-t[10]*X[1:mn,5]-t[11]*X[1:mn,6]-t[12]*X[1:mn,7])),sd=sqrt(t[16]-t[17]^2/t[18]))*dnorm(X[1:mn,4],mean=(t[9]+t[10]*X[1:mn,5]+t[11]*X[1:mn,6]+t[12]*X[1:mn,7]),sd=sqrt(t[18])))	
		if((n-mn)>0) fc[(mn+1):(mn+EM*(n-mn))] <-  -1/EM*log(dnorm(X[(mn+1):(mn+EM*(n-mn)),2], mean=(t[1]+t[2]*X[(mn+1):(mn+EM*(n-mn)),5]+t[3]*X[(mn+1):(mn+EM*(n-mn)),6]+t[4]*X[(mn+1):(mn+EM*(n-mn)),7]+as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%t((X[(mn+1):(mn+EM*(n-mn)),3:4]-matrix(c(t[5]+t[6]*X[(mn+1):(mn+EM*(n-mn)),5]+t[7]*X[(mn+1):(mn+EM*(n-mn)),6]+t[8]*X[(mn+1):(mn+EM*(n-mn)),7],t[9]+t[10]*X[(mn+1):(mn+EM*(n-mn)),5]+t[11]*X[(mn+1):(mn+EM*(n-mn)),6]+t[12]*X[(mn+1):(mn+EM*(n-mn)),7]),EM*(n-mn),2))))),sd=sqrt((t[13]-as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%matrix(c(t[14],t[15]),2,1)))))*dnorm(X[(mn+1):(mn+EM*(n-mn)),3],mean=(t[5]+t[6]*X[(mn+1):(mn+EM*(n-mn)),5]+t[7]*X[(mn+1):(mn+EM*(n-mn)),6]+t[8]*X[(mn+1):(mn+EM*(n-mn)),7]+t[17]/t[18]*(X[(mn+1):(mn+EM*(n-mn)),4]-t[9]-t[10]*X[(mn+1):(mn+EM*(n-mn)),5]-t[11]*X[(mn+1):(mn+EM*(n-mn)),6]-t[12]*X[(mn+1):(mn+EM*(n-mn)),7])),sd=sqrt(t[16]-t[17]^2/t[18]))*dnorm(X[(mn+1):(mn+EM*(n-mn)),4],mean=(t[9]+t[10]*X[(mn+1):(mn+EM*(n-mn)),5]+t[11]*X[(mn+1):(mn+EM*(n-mn)),6]+t[12]*X[(mn+1):(mn+EM*(n-mn)),7]),sd=sqrt(t[18])))	
		like <- sum(fc)
	}
	like
}

#Acceptance function for acceptance-rejection sampler (probability density of Y|X
AR <- function(x1,x2,x3,x4,y,beta) (1/(1+exp(-(beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x3+beta[5]*x4[1]+beta[6]*x4[2]+beta[7]*x4[3]))))^(y)*(1-1/(1+exp(-(beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x3+beta[5]*x4[1]+beta[6]*x4[2]+beta[7]*x4[3]))))^(1-y)

#Creating ordered data matrix (no censoring - all three covariates censored); uses 0 as cut-off since 
#censored values are coded as -4,-2,-5 and logs of all covariates here are positive
cendata <-function(Y1,X1){
	Xcc <- X1[(X1[,2]>0 & X1[,3]>0 & X1[,4]>0),]==
	Xmiss1 <- X1[(X1[,2]<0 & X1[,3]>0 & X1[,4]>0),]
	Xmiss2 <- X1[(X1[,2]>0 & X1[,3]<0 & X1[,4]>0),]
	Xmiss3 <- X1[(X1[,2]>0 & X1[,3]>0 & X1[,4]<0),]
	Xmiss12 <- X1[(X1[,2]<0 & X1[,3]<0 & X1[,4]>0),]
	Xmiss13 <- X1[(X1[,2]<0 & X1[,3]>0 & X1[,4]<0),]
	Xmiss23 <- X1[(X1[,2]>0 & X1[,3]<0 & X1[,4]<0),]
	Xmissa <- X1[(X1[,2]<0 & X1[,3]<0 & X1[,4]<0),]
	X <- rbind(Xcc,Xmiss1,Xmiss2,Xmiss3,Xmiss12,Xmiss13,Xmiss23,Xmissa)
	Ycc <- Y1[(X1[,2]>0 & X1[,3]>0 & X1[,4]>0)]
	Ymiss1 <- Y1[(X1[,2]<0 & X1[,3]>0 & X1[,4]>0)]
	Ymiss2 <- Y1[(X1[,2]>0 & X1[,3]<0 & X1[,4]>0)]
	Ymiss3 <- Y1[(X1[,2]>0 & X1[,3]>0 & X1[,4]<0)]
	Ymiss12 <- Y1[(X1[,2]<0 & X1[,3]<0 & X1[,4]>0)]
	Ymiss13 <- Y1[(X1[,2]<0 & X1[,3]>0 & X1[,4]<0)]
	Ymiss23 <- Y1[(X1[,2]>0 & X1[,3]<0 & X1[,4]<0)]
	Ymissa <- Y1[(X1[,2]<0 & X1[,3]<0 & X1[,4]<0)]
	Y <- c(Ycc,Ymiss1,Ymiss2,Ymiss3,Ymiss12,Ymiss13,Ymiss23,Ymissa)
	m <- c(length(Ycc),length(Ymiss1),length(Ymiss2),length(Ymiss3),length(Ymiss12),length(Ymiss13),length(Ymiss23),length(Ymissa))
	list(Y,X,m)
}

#Function to find imputations for censored values
Ximps <- function(Y,X,ests,Bcc,mn,m1,m2,m3,m12,m13,m23,ma,M,d,twocens) {
	
	Xmissing1 <- matrix(0,M,m1)
	Xmissing2 <- matrix(0,M,m2)
	Xmissing3 <- matrix(0,M,m3)
	Xmissing12 <- matrix(0,2*M,m12)
	Xmissing13 <- matrix(0,2*M,m13)
	Xmissing23 <- matrix(0,2*M,m23)
	Xmissinga <- matrix(0,3*M,ma)
	sig123 <- matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3)
	sig12 <- matrix(c(ests[13],ests[14],ests[14],ests[16]),2,2)
	sig13 <- matrix(c(ests[13],ests[15],ests[15],ests[18]),2,2)
	sig23 <- matrix(c(ests[16],ests[17],ests[17],ests[18]),2,2)
	
	for(i in (mn+1):(mn+m1)){
		r<-0
		while(r<M){
			U <- runif(1,0,1)
			x <- rtnorm(1,mean=(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7]+as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix((X[i,3:4]-c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[13]-as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix(c(ests[14],ests[15]),2,1)))),lower=-Inf,upper=d[1])
			if(U<AR(x,X[i,3],X[i,4],X[i,5:7],Y[i],Bcc)) {Xmissing1[(r+1),(i-mn)]<-x 
			r=r+1}
		}
	}
	
	for(i in (mn+m1+1):(mn+m1+m2)){
		r<-0
		while(r<M){
			U <- runif(1,0,1)
			x <- rtnorm(1,mean=(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]+as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix((X[i,c(2,4)]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[16]-as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix(c(ests[14],ests[17]),2,1)))),lower=-Inf,upper=d[2])
			if(U<AR(X[i,2],x,X[i,4],X[i,5:7],Y[i],Bcc)) {Xmissing2[(r+1),(i-mn-m1)]<-x 
			r=r+1}
		}
	}
	
	for(i in (mn+m1+m2+1):(mn+m1+m2+m3)){
		r<-0
		while(r<M){
			U <- runif(1,0,1)
			x <- rtnorm(1,mean=(ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]+as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix((X[i,2:3]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])),2,1))),sd=sqrt((ests[18]-as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix(c(ests[15],ests[17]),2,1)))),lower=-Inf,upper=d[3])
			if(U<AR(X[i,2],X[i,3],x,X[i,5:7],Y[i],Bcc)) {Xmissing3[(r+1),(i-mn-m1-m2)]<-x 
			r=r+1}
		}
	}
	
	for(i in (mn+m1+m2+m3+1):(mn+m1+m2+m3+m12)){
		if(any(twocens==i)) Upper=c(d[1],log(2)) else Upper=c(d[1],d[2])
		x <- rtmvnorm(3*M,algorithm="gibbs",mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])+as.vector(matrix(c(ests[15],ests[17]),2,1)*1/ests[18]*(X[i,4]-ests[9]-ests[10]*X[i,5]-ests[11]*X[i,6]-ests[12]*X[i,7]))),sigma=sig12-1/ests[18]*matrix(c(ests[15],ests[17]),2,1)%*%t(matrix(c(ests[15],ests[17]),2,1)),lower=c(-Inf,-Inf),upper=Upper)
		r <- t <- 0
		while(r<M){
			t <- t +1
			if(t==3*M) 	{x <- rtmvnorm(3*M,algorithm="gibbs",mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])+as.vector(matrix(c(ests[15],ests[17]),2,1)*1/ests[18]*(X[i,4]-ests[9]-ests[10]*X[i,5]-ests[11]*X[i,6]-ests[12]*X[i,7]))),sigma=sig12-1/ests[18]*matrix(c(ests[15],ests[17]),2,1)%*%t(matrix(c(ests[15],ests[17]),2,1)),lower=c(-Inf,-Inf),upper=Upper)
			t <- 1}
			U <- runif(1,0,1)
			if(U<AR(x[t,1],x[t,2],X[i,4],X[i,5:7],Y[i],Bcc)) {Xmissing12[(2*r+1):(2*r+2),(i-mn-m1-m2-m3)]<-x[t,] 
			r=r+1}
		}
	}
	
	for(i in (mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13)){
		x <- rtmvnorm(3*M,algorithm="gibbs",mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[17]),2,1)*1/ests[16]*(X[i,3]-ests[5]-ests[6]*X[i,5]-ests[7]*X[i,6]-ests[8]*X[i,7]))),sigma=sig13-1/ests[16]*matrix(c(ests[14],ests[17]),2,1)%*%t(matrix(c(ests[14],ests[17]),2,1)),lower=c(-Inf,-Inf),upper=c(d[1],d[3]))
		r <- t <- 0
		while(r<M){
			t <- t +1
			if(t==3*M) 	{x <- rtmvnorm(3*M,algorithm="gibbs",mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[17]),2,1)*1/ests[16]*(X[i,3]-ests[5]-ests[6]*X[i,5]-ests[7]*X[i,6]-ests[8]*X[i,7]))),sigma=sig13-1/ests[16]*matrix(c(ests[14],ests[17]),2,1)%*%t(matrix(c(ests[14],ests[17]),2,1)),lower=c(-Inf,-Inf),upper=c(d[1],d[3]))
			t <- 1}
			U <- runif(1,0,1)
			if(U<AR(x[t,1],X[i,3],x[t,2],X[i,5:7],Y[i],Bcc)) {Xmissing13[(2*r+1):(2*r+2),(i-mn-m1-m2-m3-m12)]<-x[t,] 
			r=r+1}
		}
	}
	
	for(i in (mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23)){
		if(any(twocens==i)) Upper=c(log(2),d[3]) else Upper=c(d[2],d[3])
		x <- rtmvnorm(3*M,algorithm="gibbs",mean=(c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[15]),2,1)*1/ests[13]*(X[i,2]-ests[1]-ests[2]*X[i,5]-ests[3]*X[i,6]-ests[4]*X[i,7]))),sigma=(sig23-1/ests[13]*matrix(c(ests[14],ests[15]),2,1)%*%t(matrix(c(ests[14],ests[15]),2,1))),lower=c(-Inf,-Inf),upper=Upper)
		r <- t <- 0
		while(r<M){
			t <- t +1
			if(t==3*M) 	{x <- rtmvnorm(3*M,algorithm="gibbs",mean=(c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[15]),2,1)*1/ests[13]*(X[i,2]-ests[1]-ests[2]*X[i,5]-ests[3]*X[i,6]-ests[4]*X[i,7]))),sigma=(sig23-1/ests[13]*matrix(c(ests[14],ests[15]),2,1)%*%t(matrix(c(ests[14],ests[15]),2,1))),lower=c(-Inf,-Inf),upper=Upper)
			t <- 1}
			U <- runif(1,0,1)
			if(U<AR(X[i,2],x[t,1],x[t,2],X[i,5:7],Y[i],Bcc)) {Xmissing23[(2*r+1):(2*r+2),(i-mn-m1-m2-m3-m12-m13)]<-x[t,]
			r=r+1}
		}
	}
	
	for(i in (mn+m1+m2+m3+m12+m13+m23+1):n){
		if(any(twocens==i)) Upper=c(d[1],log(2),d[3]) else Upper=d
		x <- rtmvnorm(2*M,algorithm="gibbs",mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sigma=matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3),lower=c(-Inf,-Inf,-Inf),upper=Upper)
		r <- t <- 0
		while(r<M){
			t <- t +1
			if(t==2*M) {x <- rtmvnorm(2*M,algorithm="gibbs",mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sigma=matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3),lower=c(-Inf,-Inf,-Inf),upper=Upper)
			t <- 1}
			U <- runif(1,0,1)
			if(U<AR(x[t,1],x[t,2],x[t,3],X[i,5:7],Y[i],Bcc))  {Xmissinga[(3*r+1):(3*r+3),(i-mn-m1-m2-m3-m12-m13-m23)]<-x[t,]
			r=r+1}
		}
	}
	
	
	Xmi<- NULL
	for(i in 1:M){
		Xmi <- cbind(Xmi,rbind(X[1:mn,],cbind(1,Xmissing1[i,],X[(mn+1):(mn+m1),3:7]),cbind(X[(mn+m1+1):(mn+m1+m2),1:2],Xmissing2[i,],X[(mn+m1+1):(mn+m1+m2),4:7]),cbind(X[(mn+m1+m2+1):(mn+m1+m2+m3),1:3],Xmissing3[i,],X[(mn+m1+m2+1):(mn+m1+m2+m3),5:7]),cbind(1,t(Xmissing12[(2*i-1):(2*i),]),X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),4:7]),cbind(1,Xmissing13[(2*i-1),],X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),3],Xmissing13[(2*i),],X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),5:7]),cbind(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),1:2],t(Xmissing23[(2*i-1):(2*i),]),X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),5:7]),cbind(1,t(Xmissinga[(3*i-2):(3*i),]),X[(mn+m1+m2+m3+m12+m13+m23+1):n,5:7])))
	}
	
	Xem <- rbind(X[1:mn,], cbind(1,as.vector(Xmissing1),rep(X[(mn+1):(mn+m1),3],each=M),rep(X[(mn+1):(mn+m1),4],each=M),rep(X[(mn+1):(mn+m1),5],each=M),rep(X[(mn+1):(mn+m1),6],each=M),rep(X[(mn+1):(mn+m1),7],each=M)),
							 cbind(1,rep(X[(mn+m1+1):(mn+m1+m2),2],each=M),as.vector(Xmissing2),rep(X[(mn+m1+1):(mn+m1+m2),4],each=M),rep(X[(mn+m1+1):(mn+m1+m2),5],each=M),rep(X[(mn+m1+1):(mn+m1+m2),6],each=M),rep(X[(mn+m1+1):(mn+m1+m2),7],each=M)),
							 cbind(1,rep(X[(mn+m1+m2+1):(mn+m1+m2+m3),2],each=M),rep(X[(mn+m1+m2+1):(mn+m1+m2+m3),3],each=M),as.vector(Xmissing3),rep(X[(mn+m1+m2+1):(mn+m1+m2+m3),5],each=M),rep(X[(mn+m1+m2+1):(mn+m1+m2+m3),6],each=M),rep(X[(mn+m1+m2+1):(mn+m1+m2+m3),7],each=M)),
							 cbind(1,as.vector(Xmissing12[seq(1,2*M,2),]),as.vector(Xmissing12[seq(2,2*M,2),]),rep(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),4],each=M),rep(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),5],each=M),rep(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),6],each=M),rep(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),7],each=M)),
							 cbind(1,as.vector(Xmissing13[seq(1,2*M,2),]),rep(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),3],each=M),as.vector(Xmissing13[seq(2,2*M,2),]),rep(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),5],each=M),rep(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),6],each=M),rep(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),7],each=M)),
							 cbind(1,rep(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),2],each=M),as.vector(Xmissing23[seq(1,2*M,2),]),as.vector(Xmissing23[seq(2,2*M,2),]),rep(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),5],each=M),rep(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),6],each=M),rep(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),7],each=M)),
							 cbind(1,as.vector(Xmissinga[seq(1,3*M,3),]),as.vector(Xmissinga[seq(2,3*M,3),]),as.vector(Xmissinga[seq(3,3*M,3),]),rep(X[(mn+m1+m2+m3+m12+m13+m23+1):(mn+m1+m2+m3+m12+m13+m23+ma),5],each=M),rep(X[(mn+m1+m2+m3+m12+m13+m23+1):(mn+m1+m2+m3+m12+m13+m23+ma),6],each=M),rep(X[(mn+m1+m2+m3+m12+m13+m23+1):(mn+m1+m2+m3+m12+m13+m23+ma),7],each=M)))
	
	return(list(Xmi,Xem))
	
}  

#Calculates Improper multiple imputation estimates (and standard errors)
IMI <- function(Y,X,InitBN,Bcc,InitVar,mn,m1,m2,m3,m12,m13,m23,ma,M,d,twocens){
	Xmimp <- Ximps(Y,X,InitBN,Bcc,mn,m1,m2,m3,m12,m13,m23,ma,M,d,twocens)[[1]]
	Bmis <- matrix(0,7,M)
	FinalVar <- matrix(0,25,25)
	for(g in 1:M){
		ImProp <- brglm(Y ~ Xmimp[,(7*g-5)] + Xmimp[,(7*g-4)]+ Xmimp[,(7*g-3)]+ Xmimp[,(7*g-2)]+ Xmimp[,(7*g-1)]+ Xmimp[,(7*g)],family = binomial(logit))
		Bmis[,g] <- as.vector(coef(ImProp))
		FinalVar[1:7,1:7] <- FinalVar[1:7,1:7] + as.matrix(summary(ImProp)$cov.unscaled)/M
		Ests <- optim(InitBN,Maxing(Xmimp[,(7*g-6):(7*g)],n,0,0,0,0,0,0,0,dmatrix),hessian=TRUE,method="BFGS")
		FinalBN[,g] <-as.vector(Ests$par)
		FinalVar[8:25,8:25] <- FinalVar[8:25,8:25] + solve(Ests$hessian)/M
	}
	
	return(list(rowMeans(Bmis),sqrt(diag(ApproxVar(Bmis,FinalBN,InitVar,FinalVar)[1:7,1:7]))))
}


#Calculates approximate Standard Error for IMI assuming moderate sized M
ApproxVar <- function(beta,ests,I,F) {
	par <-rbind(beta,ests)
	BW <- matrix(0,25,25)
	for(i in 1:M) BW <- BW + (par[,i]-rowMeans(par))%*%t(par[,i]-rowMeans(par))/(M-1)
	(F + (1+1/M)*BW+BW%*%solve(F)%*%I%*%solve(F)%*%BW)
}

#Fully bayesian method; obtains draws of beta and gamma parameters; also samples X 
#(treating as parameter for Bayesian method; posterior predictive for PMI)
Bayes<-function(Y,X,Bcc,ests,mn,m1,m2,m3,m12,m13,m23,ma,d,n.samples,twocens){
	
	#Initial values/defs:
	p<-ncol(X)
	prior.mn2 <- ests[c(1,5,9)]
	prior.mns <- 0
	prior.sds <- 10
	canbeta.sd <- c(0.2,.05,.02,.05,.1,.1,.002)	#Fiddled with these a little to get acceptance rates in recommended ranges (~20-40%)
	cangamma.sd <- c(0.5,.1,.1,.001,1,.1,.1,.001,.6,.1,.1,.001)
	beta <- Bcc
	gamma <- ests
	tau <- 0.01
	keep.beta <-matrix(0,n.samples,p)
	keep.gamma <- matrix(0,n.samples,18)
	keep.X <- matrix(0,n,n.samples*7)
	
	
	for(i in 1:n.samples){
		
		#Update Xcensored
		X <- Ximps(Y,X,gamma,Bcc,mn,m1,m2,m3,m12,m13,m23,ma,1,d,twocens)[[1]]
		keep.X[,(i*7-6):(i*7)] <- X
		
		#Update Beta
		for(l in 1:10){
			for(j in 1:p){
				canbeta<-beta
				canbeta[j]<-rnorm(1,beta[j],canbeta.sd[j])      #Draw candidate:
				R<-prod(dbinom(Y,1,exp(X%*%canbeta)/(1+exp(X%*%canbeta)))/(dbinom(Y,1,exp(X%*%beta)/(1+exp(X%*%beta)))))*(det(t(X)%*%diag(as.vector(exp(X%*%canbeta)/(1+exp(X%*%canbeta))^2))%*%X))/  #Compute acceptance ratio:
					((det(t(X)%*%diag(as.vector(exp(X%*%beta)/(1+exp(X%*%beta))^2))%*%X)))
				U<-runif(1)                          
				if(U<R){  beta<-canbeta }  #Accept the candidate w/ prob min(R,1)
			}}   
		keep.beta[i,]<-beta
		
		#Update mean parameters in gamma (except intercept terms
		for(l in 1:10){
			for(j in c(2,3,4,6,7,8,10,11,12)){
				cangamma<-gamma
				cangamma[j]<-rnorm(1,gamma[j],cangamma.sd[j])      #Draw candidate:
				R<-prod((dnorm(X[,2], mean=(cangamma[1]+cangamma[2]*X[,5]+cangamma[3]*X[,6]+cangamma[4]*X[,7]+as.real(matrix(c(cangamma[14],cangamma[15]),1,2)%*%solve(matrix(gamma[c(5,6,6,9)],2,2))%*%t((X[,3:4]-matrix(c(cangamma[5]+cangamma[6]*X[,5]+cangamma[7]*X[,6]+cangamma[8]*X[,7],cangamma[9]+cangamma[10]*X[,5]+cangamma[11]*X[,6]+cangamma[12]*X[,7]),n,2))))),sd=sqrt((cangamma[13]-as.real(matrix(c(cangamma[14],cangamma[15]),1,2)%*%solve(matrix(gamma[c(5,6,6,9)],2,2))%*%matrix(c(cangamma[14],cangamma[15]),2,1)))))*dnorm(X[,3],mean=(cangamma[5]+cangamma[6]*X[,5]+cangamma[7]*X[,6]+cangamma[8]*X[,7]+cangamma[17]/cangamma[18]*(X[,4]-cangamma[9]-cangamma[10]*X[,5]-cangamma[11]*X[,6]-cangamma[12]*X[,7])),sd=sqrt(cangamma[16]-cangamma[17]^2/cangamma[18]))*dnorm(X[,4],mean=(cangamma[9]+cangamma[10]*X[,5]+cangamma[11]*X[,6]+cangamma[12]*X[,7]),sd=sqrt(cangamma[18])))/(dnorm(X[,2], mean=(gamma[1]+gamma[2]*X[,5]+gamma[3]*X[,6]+gamma[4]*X[,7]+as.real(matrix(c(gamma[14],gamma[15]),1,2)%*%solve(matrix(gamma[c(5,6,6,9)],2,2))%*%t((X[,3:4]-matrix(c(gamma[5]+gamma[6]*X[,5]+gamma[7]*X[,6]+gamma[8]*X[,7],gamma[9]+gamma[10]*X[,5]+gamma[11]*X[,6]+gamma[12]*X[,7]),n,2))))),sd=sqrt((gamma[13]-as.real(matrix(c(gamma[14],gamma[15]),1,2)%*%solve(matrix(gamma[c(5,6,6,9)],2,2))%*%matrix(c(gamma[14],gamma[15]),2,1)))))*dnorm(X[,3],mean=(gamma[5]+gamma[6]*X[,5]+gamma[7]*X[,6]+gamma[8]*X[,7]+gamma[17]/gamma[18]*(X[,4]-gamma[9]-gamma[10]*X[,5]-gamma[11]*X[,6]-gamma[12]*X[,7])),sd=sqrt(gamma[16]-gamma[17]^2/gamma[18]))*dnorm(X[,4],mean=(gamma[9]+gamma[10]*X[,5]+gamma[11]*X[,6]+gamma[12]*X[,7]),sd=sqrt(gamma[18]))))*prod(dnorm(cangamma,prior.mns,prior.sds))/
					(prod(dnorm(cangamma,prior.mns,prior.sds)))
				U<-runif(1)                          
				if(U<R){  gamma<-cangamma }  #Accept the candidate w/ prob min(R,1)
			}}   
		
		#Update covariance parameters in gamma and intercept terms given mean paramters found above
		Xcent <- X[,2:4]-cbind(X[,5:7]%*%matrix(c(gamma[2],gamma[3],gamma[4]),3,1),X[,5:7]%*%matrix(c(gamma[6],gamma[7],gamma[8]),3,1),X[,5:7]%*%matrix(c(gamma[10],gamma[11],gamma[12]),3,1))
		meanx <- matrix(colMeans(Xcent),n,3,byrow=TRUE)
		CCSP <-t((Xcent-meanx))%*%(Xcent-meanx)
		Sig <- riwish((n-1+2.1),(CCSP+(0.1*diag(3))+n*tau/(tau+n)*(colMeans(Xcent)-prior.mn2)%*%t((colMeans(Xcent)-prior.mn2))))
		gamma[13:18] <- as.vector(Sig)[c(1,2,3,5,6,9)]
		gamma[c(1,5,9)] <- rmvnorm(1,mean=(1/(n+tau)*(n*colMeans(Xcent)+tau*prior.mn2)),sigma=(1/(n+tau)*Sig))
		keep.gamma[i,]<-gamma
	}
	
	#Return the posterior sample beta and the MC acceptance rates acc.rate:
	list(beta=keep.beta,gamma=keep.gamma,Ximp=keep.X)
}

#Compute PMI given 
PMI <- function(Y,XbayesMI,M){
	draws <- sort(sample(50:(dim(XbayesMI)[[2]]/7),M,replace=FALSE))	#sampling posterior predictive datasets
	VAR <- matrix(0,7,7)
	Bmis <- matrix(0,7,M)
	
	#Calculating estimates for each imputed dataset
	for(g in 1:M){
		PropMI<- brglm(Y ~ XbayesMI[,(7*draws[g]-5)]+ XbayesMI[,(7*draws[g]-4)]+XbayesMI[,(7*draws[g]-3)]+XbayesMI[,(7*draws[g]-2)]+XbayesMI[,(7*draws[g]-1)]+XbayesMI[,(7*draws[g])],family = binomial(logit))
		Bmis[,g] <- as.vector(coef(PropMI))
		VAR <- VAR + as.matrix(summary(PropMI)$cov.unscaled)/M
	}
	
	return(list(rowMeans(Bmis),as.vector(sqrt(diag(VAR[1:7,1:7])+(1+1/M)*apply(Bmis,1,var)))))
}

#Implements Monte Carlo EM algorithm
MCEM <- function(Y,X,InitBN,Bcc,mn,m1,m2,m3,m12,m13,m23,ma,n,d,EM,twocens){
	#Initial estimates of parameters in dist of X (1st 18) and Y|X (last 7)
	Parest <- matrix(c(rep(0,25),as.vector(c(InitBN,Bcc))),2,25,byrow=TRUE)
	g <-2  #counter to update estimates and decide convergence
	
	while(abs((mean(abs(Parest[max(1,(g-4)):g,]))-mean(abs(Parest[max(1,(g-9)):max(1,g-4),]))))>0.001){
		Ximp <- Ximps(Y,X,Parest[g,1:18],Parest[g,19:25],mn,m1,m2,m3,m12,m13,m23,ma,EM,d,twocens)[[2]] #obtain EM imputations given parameters
		Yimp <- c(Y[1:mn],rep(Y[(mn+1):n],each=EM))
		
		#Find weighted EM estimates at step g-1
		mylogitmimpc <- brglm(Yimp ~ Ximp[,2]+Ximp[,3]+Ximp[,4]+Ximp[,5]+Ximp[,6]+Ximp[,7],weights=c(rep(1,mn),rep(1/EM,EM*(n-mn))),family = binomial(logit))
		Gam <- optim(Parest[g,1:18], EMmaxing(Ximp,mn,EM), method="BFGS", hessian=FALSE)
		meanX <- Gam$par
		Parest <- rbind(Parest,c(as.vector(meanX), as.vector(coefficients(mylogitmimpc))))
		g <- g+1
		print(g)
		print(abs((mean(abs(Parest[max(1,(g-4)):g,]))-mean(abs(Parest[max(1,(g-9)):max(1,g-4),])))))
	}
	
	return(list(colMeans(Parest[(g-4):g,19:25]),colMeans(Parest[(g-4):g,1:18]),g))
}

#Finds standard error of Monte Carlo EM agorithm
sdMCEM <- function(Y,X,ests,Bcc,n,d,B,EM,twocens){
	Bemi <- matrix(0,7,B)
	for(i in 1:B){
		samp <- sample(1:n,n,replace=TRUE) #Bootstrap sample points
		CD <- cendata(Y[samp],X[samp,]) #Bootstrap sample
		Y1 <- CD[[1]]
		X1 <- CD[[2]]
		mn <- CD[[3]][1] 
		m1 <- CD[[3]][2] 
		m2 <- CD[[3]][3]
		m3 <- CD[[3]][4] 
		m12 <- CD[[3]][5] 
		m13 <- CD[[3]][6] 
		m23 <- CD[[3]][7] 
		ma <- CD[[3]][8] 
		Bemi[,i] <- MCEM(Y1,X1,ests,Bcc,m1,m2,m3,m12,m13,m23,ma,n,d,EM,twocens)
	}
	apply(Bemi,1,sd)
}
