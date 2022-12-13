###################Necessary Packages########################
library(mvtnorm)
library(tmvtnorm)
library(msm)
library(brglm)
library(cubature)
library(Matrix)
library(MCMCpack)

##################Simulation Parameters######################
N <- 2000	#number of simulated datasets
n <- 500	#number of obsevations per dataset
M <- 15	#number of imputations for multiple imputation procedures
B1 <- 50	#number of bootstrap samples for finding variance of IMI
B2 <- 25	#number of bootstrap samples for finding variance of ML estimates 
EM <- 250	#number of draws for MC EM algorithm procedure
n.samples <- 500	#number of samples from posterior fully Bayesian method
Beta <- c(3,1.5,-3) #logistic regression coefficients
mu <- c(0,0)	#mean of covariates
var <- c(1,1)	#variance of covariates
rho <- 0.5		#covariance (correlation since var is 1) of covariates
d <-c(-1.1685,-1.1685) #detections limits for 20% censoring
#d <-c(-0.6259,-0.6259) #detections limits for 40% censoring
#d <-c(-0.1619,-0.1619) #detections limits for 60% censoring

#################Data Generation Function####################

#Generating data for simulations (X first, then Y give X)
Datagen <- function(N,n,Beta,mu,var,rho) {
	Xgen <- matrix(0,n,3*N)
	Ygen <- matrix(0,n,N)
		for(k in 1:N){
		Xgen[,((3*k-2):(3*k))]<-cbind(rep(1,n),rmvnorm(n,c(mu[1],mu[2]),matrix(c(var[1],rho,rho,var[2]),2,2)))
		Ygen[,k]<-rbinom(n,1,1/(1+exp(-(Xgen[,((3*k-2):(3*k))]%*%Beta))))
	}
	return(list(Ygen,Xgen))
}

####################Computational Functions###################

#Puts data in order CCs, censored X1, censored X2, censored both
cendata <-function(Y1,X1,d){
	Xcc <- X1[(X1[,2]>d[1] & X1[,3]>d[2]),]
	Xmiss1 <- X1[(X1[,2]<d[1] & X1[,3]>d[2]),]
	Xmiss2 <- X1[(X1[,2]>d[1] & X1[,3]<d[2]),]
	Xmissb <- X1[(X1[,2]<d[1] & X1[,3]<d[2]),]
	X <- rbind(Xcc,Xmiss1,Xmiss2,Xmissb)
	Ycc <- Y1[(X1[,2]>d[1] & X1[,3]>d[2])]
	Ymiss1 <- Y1[(X1[,2]<d[1] & X1[,3]>d[2])]
	Ymiss2 <- Y1[(X1[,2]>d[1] & X1[,3]<d[2])]
	Ymissb <- Y1[(X1[,2]<d[1] & X1[,3]<d[2])]
	Y <- c(Ycc,Ymiss1,Ymiss2,Ymissb)
	m <- c(length(Ycc),length(Ymiss1),length(Ymiss2),length(Ymissb))
	list(Y,X,m)
}

#Function to get substitution estimates (D/sqrt(2)), assumes simulated data was log-transformed from a positive r.v.
substitution <- function(Y,X,mn,m1,m2,mb,d) {
	Xsub <- X
	plug1 <- log(exp(d[1])/sqrt(2)) 
	plug2 <- log(exp(d[2])/sqrt(2)) 
	Xsub[(mn+1):(mn+m1),2] <- rep(plug1,m1)
	Xsub[(mn+m1+1):(mn+m1+m2),3] <- rep(plug2,m2)
	if(mb>0) Xsub[(mn+m1+m2+1):(mn+m1+m2+mb),2:3] <- matrix(c(rep(plug1,mb),rep(plug2,mb)),mb,2)
	Sub <- brglm(Y ~ Xsub[,2] + Xsub[,3],family = binomial(logit))
	return(Sub)
}


#negative Log-Likelihood function for censored bivariate normal 
Maxing <- function(X,mn,m1,m2,mb) {
	like<-function(t) {
	fc<-rep(0,n)
	if(mn>0) fc[1:mn] <- -log(dnorm(X[1:mn,3],(t[2]),sqrt(max(t[5],0.001)))*dnorm(X[1:mn,2],(t[1]+t[4]/t[5]*(X[1:mn,3]-t[2])),sqrt(max((t[3]-t[4]^2/t[5]),0.001))))
	if(m1>0) fc[(mn+1):(mn+m1)] <- -log(dnorm(X[(mn+1):(mn+m1),3],(t[2]),sqrt(max(t[5],0.001)))*pnorm(d[1],(t[1]+t[4]/t[5]*(X[(mn+1):(mn+m1),3]-t[2])),sqrt(max((t[3]-t[4]^2/t[5]),0.001))))
      if(m2>0) fc[(mn+m1+1):(mn+m1+m2)] <- -log(dnorm(X[(mn+m1+1):(mn+m1+m2),2],(t[1]),sqrt(max(t[3],0.001)))*pnorm(d[2],(t[2]+t[4]/t[3]*(X[(mn+m1+1):(mn+m1+m2),2]-t[1])),sqrt(max((t[5]-t[4]^2/t[3]),0.001))))	
	if(mb>0) fc[(mn+m1+m2+1):n] <- rep(-log(pmvnorm(lower=c(-Inf,-Inf),upper=c(d[1],d[2]),mean=c((t[1]),(t[2])),sigma=matrix(c(max(t[3],0.001),t[4],t[4],max(t[5],0.001)),2,2))),mb)
	like <- sum(fc)
	}
	like
}

#negative Log-Likelihood function for weighted bivariate normal
EMmaxing <- function(X,mn,EM) {
	like<-function(t) {
	fc<-rep(0,(mn+EM*(n-mn)))
	if(mn>0) fc[1:mn] <- -log(dnorm(X[1:mn,3],(t[2]),sqrt(max(t[5],0.001)))*dnorm(X[1:mn,2],(t[1]+t[4]/t[5]*(X[1:mn,3]-t[2])),sqrt(max((t[3]-t[4]^2/t[5]),0.001))))
	if((n-mn)>0) fc[(mn+1):(mn+EM*(n-mn))] <-  -1/EM*log(dnorm(X[(mn+1):(mn+EM*(n-mn)),3],(t[2]),sqrt(max(t[5],0.001)))*dnorm(X[(mn+1):(mn+EM*(n-mn)),2],(t[1]+t[4]/t[5]*(X[(mn+1):(mn+EM*(n-mn)),3]-t[2])),sqrt(max((t[3]-t[4]^2/t[5]),0.001))))
	like <- sum(fc)
	}
	like
}

#Function to get mean imputations for X given other covariates
Xmeanimp <- function(X,ests,mn,m1,m2,mb,d) {
	Xmiss <- Xnew <- X[(mn+1):n,]
	distnuma <- function(x) x[1]*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]), sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	distnumb <- function(x) x[2]*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]), sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	int1 <- adaptIntegrate(distnuma,lowerLimit=c(-10,-10),upperLimit=c(d[1],d[2]))$integral
	int2 <- adaptIntegrate(distnumb,lowerLimit=c(-10,-10),upperLimit=c(d[1],d[2]))$integral

	for(i in 1:m1) Xnew[i,2] <- (ests[1]+ests[4]/ests[5]*(Xmiss[i,3]-ests[2]))-dnorm((d[1]-(ests[1]+ests[4]/ests[5]*(X[i,3]-ests[2])))/sqrt(ests[3]-ests[4]^2/ests[5]))/pnorm(d[1],mean=(ests[1]+ests[4]/ests[5]*(X[i,3]-ests[2])),sd=sqrt(ests[3]-ests[4]^2/ests[5]))*sqrt(ests[3]-ests[4]^2/ests[5])
	for(i in (m1+1):(m1+m2)) Xnew[i,3] <- (ests[2]+ests[4]/ests[3]*(Xmiss[i,2]-ests[1]))-dnorm((d[2]-(ests[2]+ests[4]/ests[3]*(X[i,2]-ests[1])))/sqrt(ests[5]-ests[4]^2/ests[3]))/pnorm(d[2],mean=(ests[2]+ests[4]/ests[3]*(X[i,2]-ests[1])),sd=sqrt(ests[5]-ests[4]^2/ests[3]))*sqrt(ests[5]-ests[4]^2/ests[3])
	if(mb>0) Xnew[(m1+m2+1):(m1+m2+mb),2:3] <- matrix(c(int1/pmvnorm(lower=c(-Inf,-Inf),upper=c(d[1],d[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2)),int2/pmvnorm(lower=c(-Inf,-Inf),upper=c(d[1],d[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))),mb,2,byrow=TRUE)	

	X <- rbind(X[1:mn,],Xnew)
}

#Function to get mean imputations for X given other covariates and response
Xmeanimp2 <- function(Y,X,ests,Bcc,mn,m1,m2,mb,d) {
	Ymiss <- Y[(mn+1):n]
	Xnew <- Xmiss <- X[(mn+1):n,]
	distnum1 <- function(x,xmiss,ymiss) x*(1/(1+exp(-(Bcc[1]+Bcc[2]*x+Bcc[3]*xmiss))))^(ymiss)*(1-(1/(1+exp(-(Bcc[1]+Bcc[2]*x+Bcc[3]*xmiss)))))^(1-ymiss)*dnorm(x,mean=(ests[1]+ests[4]/ests[5]*(xmiss-ests[2])),sd=sqrt(ests[3]-ests[4]^2/ests[5]))
	distden1 <- function(x,xmiss,ymiss) (1/(1+exp(-(Bcc[1]+Bcc[2]*x+Bcc[3]*xmiss))))^(ymiss)*(1-(1/(1+exp(-(Bcc[1]+Bcc[2]*x+Bcc[3]*xmiss)))))^(1-ymiss)*dnorm(x,mean=(ests[1]+ests[4]/ests[5]*(xmiss-ests[2])),sd=sqrt(ests[3]-ests[4]^2/ests[5]))
	distnum2 <- function(x,xmiss,ymiss) x*(1/(1+exp(-(Bcc[1]+Bcc[2]*xmiss+Bcc[3]*x))))^(ymiss)*(1-(1/(1+exp(-(Bcc[1]+Bcc[2]*xmiss+Bcc[3]*x)))))^(1-ymiss)*dnorm(x,mean=(ests[2]+ests[4]/ests[3]*(xmiss-ests[1])),sd=sqrt(ests[5]-ests[4]^2/ests[3]))
	distden2 <- function(x,xmiss,ymiss) (1/(1+exp(-(Bcc[1]+Bcc[2]*xmiss+Bcc[3]*x))))^(ymiss)*(1-(1/(1+exp(-(Bcc[1]+Bcc[2]*xmiss+Bcc[3]*x)))))^(1-ymiss)*dnorm(x,mean=(ests[2]+ests[4]/ests[3]*(xmiss-ests[1])),sd=sqrt(ests[5]-ests[4]^2/ests[3]))
	distnuma1 <- function(x) x[1]*(1/(1+exp(-(Bcc[1]+Bcc[2]*x[1]+Bcc[3]*x[2]))))*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	distnuma0 <- function(x) x[1]*(1-(1/(1+exp(-(Bcc[1]+Bcc[2]*x[1]+Bcc[3]*x[2])))))*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	distnumb1 <- function(x) x[2]*(1/(1+exp(-(Bcc[1]+Bcc[2]*x[1]+Bcc[3]*x[2]))))*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	distnumb0 <- function(x) x[2]*(1-(1/(1+exp(-(Bcc[1]+Bcc[2]*x[1]+Bcc[3]*x[2])))))*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	distdenb1 <- function(x) (1/(1+exp(-(Bcc[1]+Bcc[2]*x[1]+Bcc[3]*x[2]))))*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	distdenb0 <- function(x) (1-(1/(1+exp(-(Bcc[1]+Bcc[2]*x[1]+Bcc[3]*x[2])))))*dmvnorm(c(x[1],x[2]),mean=c(ests[1],ests[2]),sigma=matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
	int1 <- adaptIntegrate(distdenb1,lower=c(-10,-10),upper=c(d[1],d[2]))$integral
	int0 <- adaptIntegrate(distdenb0,lower=c(-10,-10),upper=c(d[1],d[2]))$integral

	for(i in 1:m1) Xnew[i,2] <- integrate(distnum1,lower=-Inf,upper=d[1],xmiss=Xmiss[i,3],ymiss=Ymiss[i])$value/integrate(distden1,lower=-Inf,upper=d[1],xmiss=Xmiss[i,3],ymiss=Ymiss[i])$value
	for(i in (m1+1):(m1+m2)) Xnew[i,3] <- integrate(distnum2,lower=-Inf,upper=d[2],xmiss=Xmiss[i,2],ymiss=Ymiss[i])$value/integrate(distden2,lower=-Inf,upper=d[2],xmiss=Xmiss[i,2],ymiss=Ymiss[i])$value
	Xnew[((m1+m2+1):(m1+m2+mb))[Ymiss[(m1+m2+1):(m1+m2+mb)]==1],2:3] <- c(adaptIntegrate(distnuma1,lowerLimit=c(-10,-10),upperLimit=c(d[1],d[2]))$integral/int1,adaptIntegrate(distnumb1,lowerLimit=c(-10,-10),upperLimit=c(d[1],d[2]))$integral/int1)
	Xnew[((m1+m2+1):(m1+m2+mb))[Ymiss[(m1+m2+1):(m1+m2+mb)]==0],2:3] <- c(adaptIntegrate(distnuma0,lowerLimit=c(-10,-10),upperLimit=c(d[1],d[2]))$integral/int0,adaptIntegrate(distnumb0,lowerLimit=c(-10,-10),upperLimit=c(d[1],d[2]))$integral/int0)

	rbind(X[1:mn,],Xnew)
}

#Employs Monte Carlo EM algorithm for maximizing likelihood 
MCEM <- function(Y,X,ests,Bcc,mn,m1,m2,mb,d,EM){
	Parest <- matrix(c(rep(0,8),as.vector(c(ests,Bcc))),2,8,byrow=TRUE)
	g <-2  #counter to update estimates and decide convergence
	while(abs((mean(abs(Parest[max(1,(g-4)):g,]))-mean(abs(Parest[max(1,(g-9)):max(1,g-4),]))))>0.001){
		Ximp <- Ximps(Y,X,Parest[g,1:5],Parest[g,6:8],mn,m1,m2,mb,d,EM)[[2]] #obtain EM imputations given parameters
		Yimp <- c(Y[1:mn],rep(Y[(mn+1):n],each=EM))

		#Find weighted EM estimates at step g-1
		mylogitmimpc <- brglm(Yimp ~ Ximp[,2]+Ximp[,3],weights=c(rep(1,mn),rep(1/EM,EM*(n-mn))),family = binomial(logit))
		meanX <- nlm(EMmaxing(Ximp,mn,EM), Parest[g,1:5], hessian=TRUE)$estimate
		Parest <- rbind(Parest,c(as.vector(meanX), as.vector(coefficients(mylogitmimpc))))
		g <- g+1
	}
	return(colMeans(Parest[(g-4):g,6:8]) )
}

#Computes variance of MCEM ML method via bootstrap
varMCEM <- function(Y,X,ests,Bcc,n,d,B2,EM){
	Bemi <- matrix(0,3,B2)
	for(i in 1:B2){
		samp <- sample(1:n,n,replace=TRUE) #Bootstrap sample points
		CD <- cendata(Y[samp],X[samp,],d) #Bootstrap sample
		Y1 <- CD[[1]]
		X1 <- CD[[2]]
		mn <- CD[[3]][1]
		m1 <- CD[[3]][2]
		m2 <- CD[[3]][3]
		mb <- CD[[3]][4]
		Bemi[,i] <- MCEM(Y1,X1,ests,Bcc,mn,m1,m2,mb,d,EM)
	}
	apply(Bemi,1,sd)
}

#Function to produce vector of imputations using acceptance-rejection with Metropolis-Hastings back-up
#Outputs both matrix convenient for M-algorithm and general imputations	
Ximps <- function(Y,X,ests,Bcc,mn,m1,m2,mb,d,M) {
	Xmiss <- X[(n-m1-m2-mb+1):n,]
	Ymiss <- Y[(n-m1-m2-mb+1):n]

	Xmissing1 <- rep(0,m1*M)
	for(i in 1:m1){
		r <- t <- 0   #t is a counter that when too high calls an MH algorithm
		x <- rtnorm(M*10,(ests[1]+ests[4]/ests[5]*(Xmiss[i,3]-ests[2])),sqrt(ests[3]-ests[4]^2/ests[5]),lower=-Inf,upper=d[1])
		while(r<M){
		t <- t+1
		if(runif(1,0,1)<AR(x[t],Xmiss[i,3],Ymiss[i],Bcc)) {Xmissing1[(i*M-r)]<- x[t]
											r=r+1}
		if(t>=M*10){Xmissing1[((i-1)*M+r+1):(i*M)] <- Met1(Xmiss[i,3],Ymiss[i],Bcc,ests,d,(M-r))
				r=M}
		}
	}

	Xmissing2 <- rep(0,m2*M)
	for(i in (m1+1):(m1+m2)){
		r <- t <- 0
		x <- rtnorm(M*10,(ests[2]+ests[4]/ests[3]*(Xmiss[i,2]-ests[1])),sqrt(ests[5]-ests[4]^2/ests[3]),lower=-Inf,upper=d[2])
		while(r<M){
		t <- t+1
		if(runif(1,0,1)<AR(Xmiss[i,2],x[t],Ymiss[i],Bcc)) {Xmissing2[((i-m1)*M-r)]<-x[t]										
		r=r+1}
		if(t>=M*10){Xmissing2[((i-m1-1)*M+r+1):((i-m1)*M)] <- Met2(Xmiss[i,2],Ymiss[i],Bcc,ests,d,(M-r))
				r=M}
		}
	}

	if(mb>0){
	quant <- 5
	x <- rtmvnorm(mb*quant*M,c(ests[1],ests[2]),matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2),lower=c(-Inf,-Inf),upper=c(d[1],d[2]))
	Xmissingb <- matrix(0,mb*M,2)
	count <- 0
	if(mb>0){
	for(i in (m1+m2+1):(m1+m2+mb)){
		r <-  t <-0
		while(r<M){
		t <- t + 1
		if(runif(1,0,1)<AR(x[((i-m1-m2)*quant*M+1-t),1],x[((i-m1-m2)*quant*M+1-t),2],Ymiss[i],Bcc))  {Xmissingb[((i-m1-m2)*M-r),]<-x[((i-m1-m2)*quant*M+1-t),]
											   r=r+1}
		if(t>=M*quant & r<M){Xmissingb[((i-m1-m2-1)*M+r+1):((i-m1-m2)*M),] <- Met12(Ymiss[i],Bcc,ests,d,(M-r))
				    r=M}
		} 
	}}}
	if(mb==0) Xmissingb <- NULL

	XEM <- cbind(1,rbind(X[1:mn,2:3],cbind(Xmissing1,rep(X[(mn+1):(mn+m1),3],each=M)),cbind(rep(X[(mn+m1+1):(mn+m1+m2),2],each=M),Xmissing2),Xmissingb))
	
	Ximputed <- NULL
	if(mb>1){
	for(i in 1:M){
		Ximputed <- cbind(Ximputed,rbind(X[1:(n-m1-m2-mb),],cbind(X[(n-m1-m2-mb+1):(n-m2-mb),1],Xmissing1[seq(i,(M*m1),M)],X[(n-m1-m2-mb+1):(n-m2-mb),3]),cbind(X[(n-m2-mb+1):(n-mb),1:2],Xmissing2[seq(i,(M*m2),M)]),cbind(1,Xmissingb[seq(i,(M*mb),M),])))
	}}
	if(mb==1){
	for(i in 1:M){
		Ximputed <- cbind(Ximputed,rbind(X[1:(n-m1-m2-mb),],cbind(X[(n-m1-m2-mb+1):(n-m2-mb),1],Xmissing1[seq(i,(M*m1),M)],X[(n-m1-m2-mb+1):(n-m2-mb),3]),cbind(X[(n-m2-mb+1):(n-mb),1:2],Xmissing2[seq(i,(M*m2),M)]),c(1,Xmissingb[seq(i,(M*mb),M),])))
	}}
	if(mb==0){
	for(i in 1:M){
		Ximputed <- cbind(Ximputed,rbind(X[1:(n-m1-m2-mb),],cbind(X[(n-m1-m2-mb+1):(n-m2-mb),1],Xmissing1[seq(i,(M*m1),M)],X[(n-m1-m2-mb+1):(n-m2-mb),3]),cbind(X[(n-m2-mb+1):(n-mb),1:2],Xmissing2[seq(i,(M*m2),M)])))
	}}

	return(list(Ximputed,XEM))
}

#Back-up sampler for censored X1's if Acceptance-rejection is slow
Met1 <- function(X,y,Bcc,ests,d,M){
	X1imp <- rep(0,M)
	if(m1>0){
	for(i in 1:M){
		Xi <- rep(0,20)
		Xi[1] <- rtnorm(1,ests[1],ests[3],upper=d[1])
		rden <- AR(Xi[1],X,y,Bcc)*dnorm(Xi[1],(ests[1]+ests[4]/ests[5]*(X-ests[2])),sqrt(ests[3]-ests[4]^2/ests[5]))
		for(t in 2:20) {
				Xstar <- rnorm(1,Xi[t-1],1)
				rnum <- AR(Xstar,X,y,Bcc)*dnorm(Xstar,(ests[1]+ests[4]/ests[5]*(X-ests[2])),sqrt(ests[3]-ests[4]^2/ests[5]))
				r <- rnum/max(rden,1e-300)
				u <- runif(1)
				if(r>u & Xstar<d[1]) {rden <- rnum
					   Xi[t] <- Xstar}  else Xi[t] <- Xi[(t-1)]		 
		}
	X1imp[i] <- Xi[20]
	}}
	return(X1imp)
}

#Back-up sampler for censored X2's if Acceptance-rejection is slow
Met2 <- function(X,y,Bcc,ests,d,M){
	X2imp <- rep(0,M)
	if(m2>0){
	for(i in 1:M){
		Xi <- rep(0,20)
		Xi[1] <- rtnorm(1,ests[2],ests[5],upper=d[2])	
		rden <- AR(X,Xi[1],y,Bcc)*dnorm(Xi[1],(ests[2]+ests[4]/ests[3]*(X-ests[1])),sqrt(ests[5]-ests[4]^2/ests[3]))
		for(t in 2:20) {
				Xstar <- rnorm(1,Xi[t-1],1)
				rnum <- AR(Xstar,X,y,Bcc)*dnorm(Xstar,(ests[2]+ests[4]/ests[3]*(X-ests[1])),sqrt(ests[5]-ests[4]^2/ests[3]))
				r <- rnum/max(rden,1e-300)
				u <- runif(1)
				if(r>u & Xstar<d[2]) {rden <- rnum
					   Xi[t] <- Xstar}  else Xi[t] <- Xi[(t-1)]			 
		}
	X2imp[i] <- Xi[20]
	}}
	return(X2imp)
}

#Back-up sampler for censored X1 & X2's if Acceptance-rejection is slow
Met12 <- function(y,Bcc,ests,d,M){
	start1 <- ests[1]-dnorm((d[1]- ests[1])/sqrt(ests[3]))/pnorm((d[1]- ests[1])/sqrt(ests[3]))*sqrt(ests[3])
	start2 <- ests[2]-dnorm((d[2]- ests[2])/sqrt(ests[5]))/pnorm((d[2]- ests[2])/sqrt(ests[5]))*sqrt(ests[5])

	Genvals <- rmvnorm(20*M,c(0,0),matrix(c(ests[3]/2,ests[4]/2,ests[4]/2,ests[5]/2),2,2))
	l <- 1
	X12imp <- matrix(0,M,2)
	if(mb>0){
	for(i in 1:M){
		Xi <- matrix(0,20,2)
		Xi[1,] <- rtmvnorm(1,c(ests[1],ests[2]),matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2),upper=d)
		rden <- AR(Xi[1,1],Xi[1,2],y,Bcc)*dmvnorm(Xi[1,],c((ests[1]),(ests[2])),matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))	
		for(t in 2:20) {
				Xstar <- Xi[(t-1),]+Genvals[l,]
				l <- l+1
				rnum <- AR(Xstar[1],Xstar[2],y,Bcc)*dmvnorm(Xstar,c((ests[1]),(ests[2])),matrix(c(ests[3],ests[4],ests[4],ests[5]),2,2))
				r <- rnum/max(rden,1e-300)
				u <- runif(1)
				if(r>u & Xstar[1]<d[1] & Xstar[2]<d[2]) {rden <- rnum
					   Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]	 
		}
	X12imp[i,] <- Xi[20,]
	}}
	return(X12imp)
}

#Acceptance Function for Acceptance-Rejection Sampler (probability density of Y|X)
AR <- function(x1,x2,y,beta) (1/(1+exp(-(beta[1]+beta[2]*x1+beta[3]*x2))))^(y)*(1-1/(1+exp(-(beta[1]+beta[2]*x1+beta[3]*x2))))^(1-y)

#Fully bayesian method; obtains draws of beta and gamma parameters; also samples X (treating as parameter for Bayesian method; posterior predictive for PMI)
Bayes<-function(Y,X,ests,Bcc,mn,m1,m2,mb,d,n.samples){
    
	#Initial values/defs:
	p<-ncol(X)
	prior.mn <- ests[1:2]
	beta <- Bcc
	gamma <- ests
	tau <- 0.01
	keep.beta <-matrix(0,n.samples,3)
 	keep.gamma <- matrix(0,n.samples,5)  #Used to monitor convergence of gamma parameters
	keep.x <- matrix(0,n,n.samples*3)

	for(i in 1:n.samples){

		#Update X (censored values)
		X <- Ximps(Y,X,gamma,beta,mn,m1,m2,mb,d,1)[[1]]
	    	keep.x[,(3*i-2):(3*i)] <- X

     		#Update Beta
     		for(l in 1:10){
    		for(j in 1:p){
      		canbeta<-beta
      		canbeta[j]<-rnorm(1,beta[j],1)      #Draw candidate:
       		R<-prod(dbinom(Y,1,exp(X%*%canbeta)/(1+exp(X%*%canbeta))))*sqrt(det(t(X)%*%diag(as.vector(exp(X%*%canbeta)/(1+exp(X%*%canbeta))^2))%*%X))/  #Compute acceptance ratio:
          		(prod(dbinom(Y,1,exp(X%*%beta)/(1+exp(X%*%beta))))*sqrt(det(t(X)%*%diag(as.vector(exp(X%*%beta)/(1+exp(X%*%beta))^2))%*%X)))                        
       		if(runif(1)<R){  beta<-canbeta }  #Accept the candidate w/ prob min(R,1)
     		}}   
		keep.beta[i,]<-beta

		#Update gamma
    		meanx <- matrix(colMeans(X[,2:3]),n,2,byrow=TRUE)
     		CCSP <-t((X[,2:3]-meanx))%*%(X[,2:3]-meanx)
     		Sig <- riwish((n-1+2.1),(CCSP+(0.1*diag(2))+n*tau/(tau+n)*(colMeans(X[,2:3])-prior.mn)%*%t((colMeans(X[,2:3])-prior.mn))))
     		gamma[3:5] <- as.vector(Sig)[c(1,2,4)]
     		gamma[1:2] <- rmvnorm(1,mean=(1/(n+tau)*(n*colMeans(X[,2:3])+tau*prior.mn)),sigma=(1/(n+tau)*Sig))
    		keep.gamma[i,]<-gamma
}

#Return the posterior sample beta, gamma, and x
list(beta=keep.beta,gamma=keep.gamma,x=keep.x)
}

#Proper multiple imputation using posterior predictive datasets
Pmi <- function(Y,Xpred,ests,M){
	draws <- sort(sample(50:(dim(Xpred)[[2]]/3),M,replace=FALSE))	#sampling posterior predictive datasets
	VAR <- matrix(0,3,3)
	Bmis <- matrix(0,3,M)

	#Calculating estimates for each imputed dataset
	for(g in 1:M){
		propmi <- brglm(Y ~ Xpred[,(3*draws[g]-1)] + Xpred[,(3*draws[g])],family = binomial(logit))
		Bmis[,g] <- as.vector(coef(propmi))
		VAR <- VAR + as.matrix(summary(propmi)$cov.unscaled)/M
	}
	return(list(rowMeans(Bmis),sqrt(diag(VAR[1:3,1:3])+(1+1/M)*apply(Bmis,1,var))))
}

#Improper multiple imputation method
Imi <- function(Y,X,ests,Bcc,InitVar,mn,m1,m2,mb,n,d,M){
	Xmimp <- Ximps(Y,X,ests,Bcc,mn,m1,m2,mb,d,M)[[1]]
	FinalVar <- matrix(0,8,8)
	Bmis <- matrix(0,3,M)
	FinalBN <- matrix(0,5,M)

	#Calculating Parameter Estimates for each Imputation
	for(g in 1:M){
		impropmi <- brglm(Y ~ Xmimp[,(3*g-1)] + Xmimp[,(3*g)],family = binomial(logit))
		Bmis[,g] <- as.vector(coef(impropmi))
		FinalVar[1:3,1:3] <- FinalVar[1:3,1:3] + as.matrix(summary(impropmi)$cov.unscaled)/M
		Ests <- nlm(Maxing(Xmimp[,(3*g-2):(3*g)],n,0,0,0), ests,hessian=TRUE)
		FinalBN[,g] <-as.vector(Ests$estimate)
		FinalVar[4:8,4:8] <- FinalVar[4:8,4:8] + solve(Ests$hessian)/M
	}

	return(list(rowMeans(Bmis),sqrt(diag(ConsistVar(Y,Xmimp,Bmis,FinalBN,InitVar,FinalVar)[1:3,1:3])), sqrt(diag(ApproxVar(Bmis,FinalBN,InitVar,FinalVar)[1:3,1:3])), sqrt(diag(FinalVar[1:3,1:3])+(1+1/M)*apply(Bmis,1,var))))
}

#Calculates consistent standard error for IMI
ConsistVar <- function(Y,X,beta,ests,I,F){
	par <- rbind(beta,ests)
	IFmIN <-matrix(0,8,8)
	ScoreMatrix <- matrix(0,8,M*n)
	for(g in 1:M){
		ScoreMatrix[1,(g*n+1-n):(g*n)] <- (Y-1/(1+exp(-(par[1,g]+par[2,g]*X[,(3*g-1)]+par[3,g]*X[,(3*g)]))))
		ScoreMatrix[2,(g*n+1-n):(g*n)] <- X[,(3*g-1)]*(Y-1/(1+exp(-(par[1,g]+par[2,g]*X[,(3*g-1)]+par[3,g]*X[,(3*g)]))))
		ScoreMatrix[3,(g*n+1-n):(g*n)] <- X[,(3*g)]*(Y-1/(1+exp(-(par[1,g]+par[2,g]*X[,(3*g-1)]+par[3,g]*X[,(3*g)]))))
		ScoreMatrix[4,(g*n+1-n):(g*n)] <- -(-(2*(X[,(3*g-1)]-par[4,g]))/par[6,g]+2*par[7,g]*(X[,(3*g)]-par[5,g])/(par[6,g]*par[8,g]))/(2-2*par[7,g]^2/(par[6,g]*par[8,g]))
		ScoreMatrix[5,(g*n+1-n):(g*n)] <- -(-(2*(X[,(3*g)]-par[5,g]))/par[8,g]+2*par[7,g]*(X[,(3*g-1)]-par[4,g])/(par[6,g]*par[8,g]))/(2-2*par[7,g]^2/(par[6,g]*par[8,g]))
		ScoreMatrix[6,(g*n+1-n):(g*n)] <- -1/(2*par[6,g])-(1/2)*par[7,g]^2/(par[6,g]^2*par[8,g]*(1-par[7,g]^2/(par[6,g]*par[8,g])))+(2*((X[,(3*g-1)]-par[4,g])^2/par[6,g]+(X[,(3*g)]-par[5,g])^2/par[8,g]-2*par[7,g]*(X[,(3*g-1)]-par[4,g])*(X[,(3*g)]-par[5,g])/(par[6,g]*par[8,g])))*par[7,g]^2/((2-2*par[7,g]^2/(par[6,g]*par[8,g]))^2*par[6,g]^2*par[8,g])-(-(X[,(3*g-1)]-par[4,g])^2/par[6,g]^2+2*par[7,g]*(X[,(3*g-1)]-par[4,g])*(X[,(3*g)]-par[5,g])/(par[6,g]^2*par[8,g]))/(2-2*par[7,g]^2/(par[6,g]*par[8,g]))
		ScoreMatrix[7,(g*n+1-n):(g*n)] <- par[7,g]/(par[6,g]*par[8,g]*(1-par[7,g]^2/(par[6,g]*par[8,g])))-(4*((X[,(3*g-1)]-par[4,g])^2/par[6,g]+(X[,(3*g)]-par[5,g])^2/par[8,g]-2*par[7,g]*(X[,(3*g-1)]-par[4,g])*(X[,(3*g)]-par[5,g])/(par[6,g]*par[8,g])))*par[7,g]/((2-2*par[7,g]^2/(par[6,g]*par[8,g]))^2*par[6,g]*par[8,g])+(2*(X[,(3*g-1)]-par[4,g]))*(X[,(3*g)]-par[5,g])/((2-2*par[7,g]^2/(par[6,g]*par[8,g]))*par[6,g]*par[8,g])              
		ScoreMatrix[8,(g*n+1-n):(g*n)] <- -1/(2*par[8,g])-(1/2)*par[7,g]^2/(par[6,g]*par[8,g]^2*(1-par[7,g]^2/(par[6,g]*par[8,g])))+(2*((X[,(3*g-1)]-par[4,g])^2/par[6,g]+(X[,(3*g)]-par[5,g])^2/par[8,g]-2*par[7,g]*(X[,(3*g-1)]-par[4,g])*(X[,(3*g)]-par[5,g])/(par[6,g]*par[8,g])))*par[7,g]^2/((2-2*par[7,g]^2/(par[6,g]*par[8,g]))^2*par[6,g]*par[8,g]^2)-(-(X[,(3*g)]-par[5,g])^2/par[8,g]^2+2*par[7,g]*(X[,(3*g-1)]-par[4,g])*(X[,(3*g)]-par[5,g])/(par[6,g]*par[8,g]^2))/(2-2*par[7,g]^2/(par[6,g]*par[8,g]))
	}
	for(i in 1:n){
		Means <- rowMeans(ScoreMatrix[,seq(i,M*n,n)])
		for(g in 1:M){
			IFmIN <- IFmIN + (ScoreMatrix[,(n*(g-1)+i)]-Means)%*%t((ScoreMatrix[,(n*(g-1)+i)]-Means))/(M-1)
	}}
	(F + (1+1/M)*F%*%IFmIN%*%F+ F%*%IFmIN%*%I%*%IFmIN%*%F)
}

#Calculates approximate Standard Error for IMI assuming moderate sized M
ApproxVar <- function(beta,ests,I,F) {
	par <-rbind(beta,ests)
	BW <- matrix(0,8,8)
	for(i in 1:M) BW <- BW + (par[,i]-rowMeans(par))%*%t(par[,i]-rowMeans(par))/(M-1)
	(F + (1+1/M)*BW+BW%*%I%*%BW)
}

#Calculates a bootstrap standard error estimate for the improper multiple imputation method
BSse <- function(Y,X,d,B1,M){
	Bmis <- matrix(0,3,M)
	Bmisboot <- matrix(0,3,B1)
	for(g in 1:B1){
		oksamp <- 0 #OKsamp 0 if bad sample (all 1's/0's in CC's)
		while(oksamp<1){
			samp <- sample(1:n,n,replace=TRUE) #Bootstrap sample points
			tm <- cendata(Y[samp],X[samp,],d) #Bootstrap sample
			Y1 <- tm[[1]]
			X1 <- tm[[2]]
			mn <- tm[[3]][1]
			m1 <- tm[[3]][2]
			m2 <- tm[[3]][3]
			mb <- tm[[3]][4]

			#If weird sample (and weird parameter estimates because of this), throws out sample and starts over
			tryCatch(IMIm<-glm(Y1[1:mn] ~ X1[1:mn,2] + X1[1:mn,3],family = binomial(logit)), warning=function(...) oksamp <<- 0)
			if(max(abs(as.vector(coef(IMIm))))<8) oksamp <-1
		}

		#Complete Case Estimates for Bootstrap sample g
		CC <- brglm(Y1[1:mn] ~ X1[1:mn,2] + X1[1:mn,3],family = binomial(logit))
		Bcc <- as.vector(coef(CC))
		InitBN <- as.vector(nlm(Maxing(X1,mn,m1,m2,mb), c(0.001,0.001,1,0.5,1))$estimate)
		Xmimpb <- Ximps(Y1,X1,InitBN,Bcc,mn,m1,m2,mb,d,M)[[1]]
		for(i in 1:M){
			MIb <- brglm(Y1 ~ Xmimpb[,(3*i-1)] + Xmimpb[,(3*i)],family = binomial(logit))
			Bmis[,i] <- as.vector(coef(MIb))
		}

		Bmisboot[,g] <- rowMeans(Bmis)
	}
	return(apply(Bmisboot,1,sd))
}


#########################Simulation###########################
set.seed(8089205)
Data <- Datagen(N,n,Beta,mu,var,rho)  #Generates data for all N datasets
Ygen <- Data[[1]]
Xgen <- Data[[2]]

#Initialization of Method Estimate Vectors and SEs
Btrue <- Bcc <- Bsub <- Bmean1 <- Bmean2 <- Bml <- Bbayes <- Bpmi  <- Bimi <- matrix(0,3,N)
SEtrue <- SEcc <- SEsub <- SEmean1 <- SEmean2 <- SEml <- SEbayes <- SEpmi  <- SEimi_c <- SEimi_a <- SEimi_r <- SEimi_b <- matrix(0,3,N)
Ttrue <- Tcc <- Tsub <- Tcov <- Tmean1 <- Tmean2 <- Tml <- Tbayes <- Timi <- rep(0,N)
InitVar <- matrix(0,8,8)

#Loops through all N datasets and calculates Beta and SE estimates for each method
for(k in 1:N){
	CD <- cendata(Ygen[,k],Xgen[,((3*k-2):(3*k))],d)
	Y <- CD[[1]]
	X <- CD[[2]]
	mn <- CD[[3]][1]	#no censoring (missing)
	m1 <- CD[[3]][2]  #x1 censored
	m2 <- CD[[3]][3]	#x2 censored
	mb <- CD[[3]][4]	#x1 and x2 censored

	####Method 1: Omnicient Method####
	Ttrue[k] <- system.time({TRUEm <- brglm(Y ~ X[,2] + X[,3],family = binomial(logit)) #brglm is biased correction (for 1/nth order bias) version of glm
	Btrue[,k] <- as.vector(coef(TRUEm))
	SEtrue[,k] <- as.vector(summary(TRUEm)$coefficient[,2])})[3]


	####Method 2: Complete Case Method####
	Tcc[k] <- system.time({CC <- brglm(Y[1:mn] ~ X[1:mn,2] + X[1:mn,3],family = binomial(logit))
	Bcc[,k] <- as.vector(coef(CC))
	SEcc[,k] <- as.vector(summary(CC)$coefficient[,2])
	InitVar[1:3,1:3] <- as.matrix(summary(CC)$cov.unscaled) #Getting Inital var ests. for improper MI var estimation
	})[3]


	####Method 3: Substituation with D/sqrt(2)####
	Tsub[k] <- system.time({SUB <- substitution(Y,X,mn,m1,m2,mb,d)
	Bsub[,k] <- as.vector(coef(SUB))
	SEsub[,k] <- as.vector(summary(SUB)$coefficient[,2])
	})[3]
	
	#Maximum Likelihood Estimates for BVN distribution for X covariates (needed for Mean1, Mean2, and IMI methods and helpful for 
	#initial values in Bayesian and ML methods
	Tcov[k] <- system.time({Ests <- nlm(Maxing(X,mn,m1,m2,mb), c(0.001,0.001,1,rho,1),hessian=TRUE)
	InitBN <- as.vector(Ests$estimate)
	InitVar[4:8,4:8] <- solve(Ests$hessian) #Getting Inital var ests. for improper MI var estimation
	})[3]

	####Method 4: Conditional Mean method ignoring Y####
	Tmean1[k] <- system.time({Ximputed <- Xmeanimp(X,InitBN,mn,m1,m2,mb,d)
	Mean1 <- brglm(Y ~ Ximputed[,2] + Ximputed[,3],family = binomial(logit))
	Bmean1[,k] <- as.vector(coef(Mean1))
	SEmean1[,k] <- as.vector(summary(Mean1)$coefficient[,2])
	})[3]

	####Method 5: Conditional Mean method conditioning on Y and observed X####
	Tmean2[k] <- system.time({Ximputed2 <- Xmeanimp2(Y,X,InitBN,Bcc[,k],mn,m1,m2,mb,d)
	Mean2 <- brglm(Y ~ Ximputed2[,2] + Ximputed2[,3],family = binomial(logit))
	Bmean2[,k] <- as.vector(coef(Mean2))
	SEmean2[,k] <- as.vector(summary(Mean2)$coefficient[,2])
	})[3]

	####Method 6: Maximum Likelihood via Monte Carlo EM algorithm####
	Tml[k] <- system.time({Bml[,k] <- MCEM(Y,X,InitBN,Bcc[,k],mn,m1,m2,mb,d,EM)
	
	#In order to get variance estimate in general, use bootstrap below;  note that we can use Louis (1982) as well to obtain a bootstrap
	#variance estimate if we do not average across last 5 (or any number >1) iterations to get estimate
	#SEml[,k] <- varMCEM(Y,X,InitBN,Bcc[,k],n,d,B2,EM)	
	})[3]

	####Method 7: Fully Bayesian Method####
	Tbayes[k] <- system.time({BAYES <- Bayes(Y,X,InitBN,Bcc[,k],mn,m1,m2,mb,d,n.samples)
	Bbayes[,k] <- c(density(BAYES[[1]][50:n.samples,1])$x[density(BAYES[[1]][50:n.samples,1])$y==max(density(BAYES[[1]][50:n.samples,1])$y)],density(BAYES[[1]][50:n.samples,2])$x[density(BAYES[[1]][50:n.samples,2])$y==max(density(BAYES[[1]][50:n.samples,2])$y)],density(BAYES[[1]][50:n.samples,3])$x[density(BAYES[[1]][50:n.samples,3])$y==max(density(BAYES[[1]][50:n.samples,3])$y)])
	SEbayes[,k] <- apply(BAYES[[1]][50:n.samples,],2,sd)
	})[3]

	####Method 8: Proper Multiple Imputation####
	PMI <-Pmi(Y,BAYES[[3]],InitBN,M)
	Bpmi[,k] <- PMI[[1]]
	SEpmi[,k] <- PMI[[2]]

	
	####Method 9: Improper Multiple Imputation####
	Timi[k] <- system.time({IMI <- Imi(Y,X,InitBN,Bcc[,k],InitVar,mn,m1,m2,mb,n,d,M)
	Bimi[,k] <- IMI[[1]]
	SEimi_c[,k] <- IMI[[2]]		#Consistent variance estimator
	SEimi_a[,k] <- IMI[[3]]		#Approximate variance estimator (consistent as m goes to inf.)
	SEimi_r[,k] <- IMI[[4]]		#Rubin's variance estimator (works poorly for improper imputation
	})[3]
	SEimi_b[,k] <- BSse(Y,X,d,B1,M) #Bootstrap variance estimator
	#multiply time for method 9 by  5/6 to get final time

	print(k) #keep track of complete datasets
}




