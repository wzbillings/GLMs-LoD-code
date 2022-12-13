###################Necessary Packages########################
library(mnormt)
library(brglm)
library(tmvtnorm)
library(cubature)
library(msm)
library(corpcor)
library(MCMCpack)

set.seed(8089205)
####################Data Read and Clean-Up###################
## Reading in the only data that was provided with the code
dat <-
	readxl::read_excel(here::here("originals", "B2013SimData-2nfu6iq.xls")) |>
	as.data.frame()

# Defining Response, Y, the survival indicator (90 -> survive, else dead)
Y <- as.numeric(dat[ , "Alive?"])

# Defining Covariate matrix:
#  covariate 1 is TNF,
#  covariate 2 is IL-6,
#  covariate 3 is IL-10,
#  covariate 4 is sex (1=male),
#  covariate 5 is race,
#  covariate 6 is age
X <- dat |>
	dplyr::select(-`Observation..`, -`Alive?`) |>
	`colnames<-`(c("TNF", "IL-6", "IL-10", "Sex", "Race", "Age")) |>
	dplyr::mutate(dplyr::across(tidyselect::everything(), as.numeric))

# Eliminating those without first day observations (should not affect analysis
# greatly since lack of first day observation occured only when patient arrived
# on a weekend
ResCov <- cbind(Y, X)
ResCov <- na.omit(ResCov)
Y <- ResCov[,1]
X <- ResCov[,2:7]

# Defining n = number of individuals (1418) and censoring values after
# log-transformation; note that for 16 individuals, the censoring value for the
# IL-6 covariate was 2 instead of 5; this is taken into account in the methods
# below
n <- length(ResCov[,1])
d <- c(log(4), log(5), log(5))

# Creating Q-Q Plots for log-transormed data considering only available cases
# for each covariate; note that for IL-6, only data points with observations
# above (or equal to) 5 are considered since there are two DLs, depending on
# the patient
# Instead of calculating the quantiles starting from 0, we want to start from
# (n - # of censored data points + 1), because we are only plotting the
# non-censored values. The minimum non-censored value will estimate this order
# statistic instead of the minimum.
layout(matrix(c(1, 2, 3), nrow = 1))

qqnorm_cens <- function(x, cens, ...) {
	yplot <- log(sort(x[x >= cens]))
	n <- length(yplot)
	xvals <- seq(
		from = (n + 1) / length(x),
		to = 1,
		by = 1 / length(x)
	)
	xq <- qnorm(xvals)
	qqplot(xq, yplot, ...)
	qqline(yplot)
	
	invisible(1)
}

qqnorm_cens(
	X[, 1], 4,
	main="TNF", xlab="Quantiles of N(0,1)", ylab="log(TNF)",
	cex.main = 2, cex.lab = 1.5
)

qqnorm_cens(
	X[, 2], 5,
	main="IL-6", xlab="Quantiles of N(0,1)", ylab="log(IL-6)",
	cex.main = 2, cex.lab=1.5
)

qqnorm_cens(
	X[, 3], 5,
	main = "IL-10", xlab = "Quantiles of N(0,1)", ylab = "log(IL-10)",
	cex.main = 2, cex.lab=1.5
)

hist_cens <- function(x, cens, ...) {
	xplot <- log(ifelse(x >= cens, x, 0.5))
	myhist <- hist(xplot, xaxt = "n", ...)
	axis(
		1,
		at = c(myhist$mids[1], myhist$breaks[-c(1:2)]),
		labels = c(paste("<", cens), myhist$breaks[-c(1:2)])
	)
}

hist_cens(
	X[, 1], 4,
	main = "TNF", xlab = NULL,
	cex.main = 2, cex.lab = 1.4
)
hist_cens(
	X[, 2], 5,
	main = "IL-6", xlab = NULL,
	cex.main = 2, cex.lab = 1.4
)
hist_cens(
	X[, 3], 5,
	main = "IL-10", xlab = NULL,
	cex.main = 2, cex.lab = 1.4
)

layout(1)

# Taking log tranformation of covariates to make them approximately normal;
# defining race as a binary variable; 
# defining response as 1/0 (1 indicating event of survival)
for(i in 1:n){
  if(X[i, 1] > 4) X[i,1] <- log(X[i,1])
  if(X[i, 2] > 2) X[i,2] <- log(X[i,2])
  if(X[i, 3] > 5) X[i,3] <- log(X[i,3])
  if(X[i, 5] > 1) X[i,5] <- 0
  if(Y[i] == 90) Y[i] <-1 else Y[i] <- 0
}

#Adding in an intercept to X matrix
X <- cbind(1,X)

#creating an n x 3 matrix of the detection limits for convenience below
# (all rows are the same except for those individuals with censoring
# at 2 for IL-6
# Finds which observations for IL-6 are censored at 2 rather than 5
# (5 is much more common); none occur with just IL-6 censored.
dmatrix <- matrix(d, n, 3, byrow=TRUE)
twocens <- which(X[,3] == -2)
dmatrix[twocens,2] <- log(2)

###########Parameters and Matrix Initializations##############
M <- 15
EM <- 250
B <- 25
n.samples <- 1000
InitVar <- FinalVar <- matrix(0,25,25)
FinalBN <- FinalBN2 <- matrix(0,18,M)


######################Analysis of Data########################

CD <- cendata(Y,X)  #Ordering censored data according to censoring
Y <- CD[[1]]
X <- CD[[2]]
mn <- CD[[3]][1] #no censoring
m1 <- CD[[3]][2] #TNF censored
m2 <- CD[[3]][3] #IL-6 censored
m3 <- CD[[3]][4] #IL-10 censored
m12 <- CD[[3]][5] #TNF and IL-6 censored
m13 <- CD[[3]][6] #TNF and IL-10 censored
m23 <- CD[[3]][7] #IL-6 and IL-10 censored
ma <- CD[[3]][8] #TNF, IL-6, and IL-10 censored

source(here::here("funs.R"))

# maximizes censored trivariate normal

Ests <-
	optim(
		par = c(1.7, 0, 0, 0, 3.6, 0, 0, 0, 1.7, 0, 0, 0, 1.09, 0.96, 0.56, 5.28,
						1.80, 2.88),
		fn = Maxing(X, mn, m1, m2, m3, m12, m13, m23, ma, dmatrix),
		hessian = TRUE,
		method="BFGS",
		control = list(trace = 3, reltol = 1e-10, maxit = 1000)
	)

InitBN <- Ests$par  #Initial Parameter Estimates in trivariate normal
InitVar[8:25,8:25] <- solve(Ests$hessian) #Initial variance estimates

####1. Complete Case Estimates (ignoring data with Censoring)####
CC<- brglm(Y[1:mn] ~ X[1:mn,2] + X[1:mn,3]+X[1:mn,4] + X[1:mn,5]+X[1:mn,6]+X[1:mn,7],family = binomial(logit))
Bcc <- as.vector(coef(CC))
SEcc <- as.vector(summary(CC)$coefficient[,2])
InitVar[1:7,1:7] <- as.matrix(summary(CC)$cov.unscaled)  #Initial variance estimates


####2. Improper MI####
improp <- IMI(Y,X,InitBN,Bcc,InitVar,mn,m1,m2,m3,m12,m13,m23,ma,M,d,twocens)
Bmi <- improp[[1]]
SEmi <- improp[[2]]


####3. Fully Bayesian Method####
BAYES <- Bayes(Y,X,Bcc,InitBN,mn,m1,m2,m3,m12,m13,m23,ma,d,n.samples,twocens)
Bbayes <-  c(density(BAYES[[1]][50:n.samples,1])$x[density(BAYES[[1]][50:n.samples,1])$y==max(density(BAYES[[1]][50:n.samples,1])$y)],density(BAYES[[1]][50:n.samples,2])$x[density(BAYES[[1]][50:n.samples,2])$y==max(density(BAYES[[1]][50:n.samples,2])$y)],density(BAYES[[1]][50:n.samples,3])$x[density(BAYES[[1]][50:n.samples,3])$y==max(density(BAYES[[1]][50:n.samples,3])$y)],density(BAYES[[1]][50:n.samples,4])$x[density(BAYES[[1]][50:n.samples,4])$y==max(density(BAYES[[1]][50:n.samples,4])$y)],density(BAYES[[1]][50:n.samples,5])$x[density(BAYES[[1]][50:n.samples,5])$y==max(density(BAYES[[1]][50:n.samples,5])$y)],density(BAYES[[1]][50:n.samples,6])$x[density(BAYES[[1]][50:n.samples,6])$y==max(density(BAYES[[1]][50:n.samples,6])$y)],density(BAYES[[1]][50:n.samples,7])$x[density(BAYES[[1]][50:n.samples,7])$y==max(density(BAYES[[1]][50:n.samples,7])$y)])
SEbayes <- apply(BAYES[[1]][50:n.samples,],2,sd)
FinalBNbayes <- colMeans(BAYES[[2]][50:n.samples,]) # to monitor convergence

####4. Proper MI ####
propmi <- PMI(Y, BAYES[[3]],M)
BbayesMI <- propmi[[1]]
SEbayesMI <- propmi[[2]]


####5. MCEM Maximum Liklihood####
em <- MCEM(Y,X,InitBN,Bcc,mn,m1,m2,m3,m12,m13,m23,ma,n,d,EM,twocens)
Bem <- em[[1]]
SEem <- sdMCEM(Y,X,InitBN,Bcc,n,d,B,EM,twocens)
