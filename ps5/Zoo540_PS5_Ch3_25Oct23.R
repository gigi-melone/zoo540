# Zoology/Entomology 540: Problem Set 5

# Due 24 October (midnight)

	# For the homework, ONLY TURN IN THE R CODE THAT YOU USED. Start with the code from correlated_data which is pasted below. You can then add your code, remove any of the code below that you don't need for the homework questions, and save it to a separate file.

# 1. What questions do you have about the material in Ch 3? What needs more explanation? I'm serious about asking this question, because I want to improve the book. (NOTE: This requires NO NEW R CODE.)

# 2. Construct the right panel of figure 3.8 giving power curves for the detection of phylogenetic signal in simulated data using an OU transform. You can modify the code I provided for the left panel which simulates data using Pagel's λ transform which is below labeled "Fig. 3.8 left panel".

# 3. The parametric bootstrap of H0 for phylogenetic signal (subsection 3.5.3) uses the log likelihood ratio (LLR) as the test statistic. It is also possible to use the phylogenetic signal parameter (λ or α) as the test statistic. This would involve the following: (i) Fit the model to the data and calculate the value of the phylogenetic signal parameter (λ or α). (ii) Refit the model with no phylogenetic signal and use the resulting parameter values to simulate a large number of datasets. (iii) Refit the model including the phylogenetic signal parameter (λ or α) for each dataset and collect these values. (iv) The resulting distribution of λ or α estimated from the simulated datasets approximates the distribution of the estimator of λ or α, allowing P-values to be calculated. Perform this bootstrap and compare the results to the bootstrap of H0 using the LLR.

# 4. When investigating mixed models (chapters 1 and 2), we focused on testing hypotheses concerning the fixed effects (the slope). Knowing what you know about testing for phylogenetic signal (sections 3.4 and 3.5), what methods could you use to test hypotheses regarding the random effects in a mixed model? You don't need to do this in R: I'm only asking for ideas. But you will impress me if you can take the models from chapter 1 and perform a test of a null hypothesis concerning a random effect. 

#################################################################
# load libraries
#################################################################
#install.packages("ape")
library(ape)
#install.packages("phylolm")
library(phylolm)
library(lme4)
library(nlme)
library(mvtnorm)
#install.packages("phytools")
library(phytools)
library(phangorn)
#install.packages("rr2")
library(rr2)
library(phyr)
library(logistf)


#~~~~~~~~~~~~~~~~~~~~
n <- 30
lam<-0
alpha<-50
nboot <- 2000

phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
phy.ou <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha=alpha))$tree

boot0 <- data.frame(LLR.lam=rep(NA, nboot), LLR.OU=NA)
for(i in 1:nboot){
  Y.sim <- rTraitCont(phy.ou, model = "BM", sigma = 1)
  
  z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
  z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
  z.0.sim <- lm(Y.sim ~ 1)
  
  boot0$LLR.lam[i] <- 2*(z.lam.sim$logLik - logLik(z.0.sim))
  boot0$LLR.OU[i] <- 2*(z.OU.sim$logLik - logLik(z.0.sim))
}
LLR.lam.crit <- sort(boot0$LLR.lam)[floor(0.95 * nboot)]
LLR.OU.crit <- sort(boot0$LLR.OU)[floor(0.95 * nboot)]

# This plots the bootstrap distributions of LLR for Pagel's lambda and the OU transform.
par(mfrow=c(2,1))
hist(boot0$LLR.lam)
lines(.1*(1:100), nboot*dchisq(.1*(1:100), df=1), col="blue")
hist(boot0$LLR.OU)
lines(.1*(1:100), nboot*dchisq(.1*(1:100), df=1), col="blue")

# Compute the rejection rates for the 5 methods.
# Note that the same starter phylogeny is used for the simulations, since the critical values of lam and alpha used in the bootstraps around H0 could depend on the topology.
alpha.list <- c(1,2,5,10,20,50)
nsims <- 2000

reject <- data.frame(alpha=rep(NA, nsims*length(alpha.list)), lam.LRT=NA, OU.LRT=NA, lam.boot0=NA, OU.boot0=NA, K.perm=NA)
counter <- 0
for(alpha in alpha.list){
  phy.OU <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha=alpha))$tree
  for(i in 1:nsims){
    counter <- counter + 1
    reject$lam[counter] <- lam      
    Y.sim <- rTraitCont(phy.ou, model = "BM", sigma = 1)
    
    # LRTs
    z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
    z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
    z.0.sim <- lm(Y.sim ~ 1)
    LLR.lam <- 2*(z.lam.sim$logLik - logLik(z.0.sim)[1])
    LLR.OU <- 2*(z.OU.sim$logLik - logLik(z.0.sim)[1])
    reject$lam.LRT[counter] <- (pchisq(LLR.lam, df=1, lower.tail=F) < 0.05)
    reject$OU.LRT[counter] <- (pchisq(LLR.OU, df=1, lower.tail=F) < 0.05)
    
    # Bootstraps around H0
    reject$lam.boot0[counter] <- (LLR.lam > LLR.lam.crit)
    reject$OU.boot0[counter] <- (LLR.OU > LLR.OU.crit)
    
    # Permutation test with Blomberg's K
    z.phylosig.K <- phylosig(Y.sim, tree=phy, method="K", test=TRUE, nsim=2000)
    reject$K.perm[counter] <- (z.phylosig.K$P < 0.05)
  }
}
# This code took a long time to run, so I just saved the output and then reloaded it later. I've commented out these lines though.
# write.table(reject, file="Phylo signal power curves lam n=30.csv", sep=",", row.names=F)
# reject <- read.csv(file="Phylo signal power curves lam n=30.csv")

#Create power curves
power <- aggregate(reject[,c("lam.LRT","OU.LRT","lam.boot0","OU.boot0","K.perm")], by = list(reject$alpha), FUN = mean)
names(power)[1] <- "alpha"

# Fig. 3.8, right panel


par(mfrow=c(2,1))
plot(OU.LRT ~ alpha, data=power, typ="l", xlim=c(0,50), ylim=c(0, 1), xlab=expression(paste("Phylogenetic signal (", alpha, ")")), ylab="Fraction rejected")
lines(lam.LRT ~ alpha, data=power, col=2)
lines(lam.boot0 ~ alpha, data=power, col=3)
lines(OU.boot0 ~ alpha, data=power, col=4)
lines(K.perm ~ alpha, data=power, col=5)


lines(c(0,10), c(.05,.05), lty=2)
legend(.0,1,c("lam.LRT","OU.LRT", "lam.boot(H0)", "OU.boot(H0)", "K.perm"), col=1:5, lty=1)

--

#################################################################
# 3.3 Phylogenetic correlation
#################################################################

######
# 3.3.1 Phylogenetic covariances

# Generate a "phylogeny" for an LMM of 3 routes containing 4 stations
VCV.lmm <- kronecker(diag(3), matrix(1,nrow=4, ncol=4)) + 2*diag(12)
# to see what the kronecker product is, print out diag(3) and matrix(1,nrow=4, ncol=4) separately. It takes the second matrix and creates a large matrix by multiplying it to each element of the first matrix. As other examples, try kronecker(diag(c(1,2,3)), matrix(1,nrow=4, ncol=4)) and kronecker(matrix(1,nrow=3, ncol=3), diag(c(10,20,30)))
VCV.lmm <- matrix(VCV.lmm, nrow=12, dimnames=list(1:12,1:12))
phy.lmm <- vcv2phylo(VCV.lmm)
phy.lmm$tip.label <- rep(4:1, times=3)
phy.lmm$node.label <- c("","route 3","","","","route 2","","","route 1","","")
VCV.lmm 

# Fig. 3.1
par(mfrow=c(1,2))
plot(phy.lmm, show.node.label = T, node.pos=2, adj=-1.2, cex=1)

# Create a phylogeny using rtree(). compute.brlen() is used to make the tree ultrametric (having contemporaneous tips)
phy <- compute.brlen(rtree(n=12), method = "Grafen", power = 1)
phy$tip.label <- 12:1

round(vcv(phy)[12:1,12:1],2)
plot(phy)

######
# 3.3.2 The strength of phylogenetic signal 

# Crate a phylogeny using rtree(). compute.brlen() is used to make the tree ultrametric (having contemporaneous tips)
phy <- compute.brlen(rtree(n=12), method = "Grafen", power = 1)
phy$tip.label <- 12:1

# Fig. 3.2
layout(matrix(c(0,2,3,0,1,2,3,4,1,5,6,4,0,5,6,0),nrow=4, byrow=T))
plot(phy, no.margin=F)
mtext("Starter Tree")

# Create a set of trees changing the strength of phylogenetic signal.
lambda <- .75
plot(transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lambda))$tree)
mtext(expression(paste(lambda," = 0.75")))
lambda <- .25
plot(transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lambda))$tree)
mtext(expression(paste(lambda," = 0.25")))
lambda <- .001
plot(transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lambda))$tree)
mtext("Star Tree")

alpha <- 2
plot(transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree)
mtext(expression(paste(alpha," = 2")))
alpha <- 6
plot(transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree)
mtext(expression(paste(alpha," = 6")))

# This is an aside that you can skip: To see another version of the OU transform, this is the one from ape. Notice that the transform never recovers the
# starter tree and does not contain zeros. corMartins in ape should be the same as model="OUrandomRoot" in phylolm,
# although model="OUrandomRoot" and model="OUfixedRoot" give the same results. You can get the "OUfixedRoot"
# covariance matrix by subtracting off the lowest value from the corMartins covariance matrix.

alpha <- 1
# The OU transform in ape uses the corStruct "corMartins"
corMatrix(Initialize(corMartins(alpha, phy), data=as.data.frame(1:12)))
vcv(transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree)[12:1,12:1]

#################################################################
# 3.4 Estimating phylogenetic signal
#################################################################

######
# 3.4.1 Phylogenetic regression model

# Create a phylogeny.
n <- 12
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
phy$tip.label <- 12:1

# This function simulates data according to a given model of evolution (in this case BM). Try running this multiple times to see how much variation there is. Can you (mosdt of the time) see the correlatoins anticipated from the phylogeny?
Y <- rTraitCont(phy, model = "BM", sigma = 1)

# Fig. 3.3
par(mfrow=c(1,2))
contMap(phy, x=Y)
par(mai=c(.8,.1,.1,.1))
plot(Y,1:length(Y), xlab="Y", ylab="", yaxt="n")

######
# 3.4.2 Analyzing hierarchical data as a phylogeny

# Construct a dataset with 4 stations in each of 6 routes
d <- data.frame(route=rep(1:4, each=6), plot=rep(1:6, times=4), x=0)

# Parameter values for the model
b0 <- 0
sd.b <- 1
sd.e <- .5

# Simulate the dataset.
for(i in d$route){
	dd <- d[d$route == i,]
	nn <- nrow(dd)
	b1.route <- rnorm(n=1, mean=0, sd=sd.b)
	d$Y[d$route == i] <- b0 + b1.route + rnorm(n=nn, sd=sd.e)
}

# Construct a "phylo" object from the covariance matrix for the LMM. this requires the function vcv2phylo() which does what the name suggests. The function phy.lm <- multi2di() is needed to resolve polytomies (situations in which branches split into more than 2 new branches, so 3 or more branches come out of the same node).
vcv <- kronecker(sd.b^2*diag(nrow=4), matrix(1, nrow=6, ncol=6)) + sd.e^2*diag(dim(d)[1])
vcv <- vcv/max(vcv)
phy.lm <- vcv2phylo(vcv, tolerance = 1e-3)
phy.lm$tip.label <- 24:1
phy.lm <- multi2di(phy.lm)

# Fig. 3.4
par(mfrow=c(1,2), mai=c(.8,.1,.1,.1))
plot(phy.lm, node.pos=2)
plot(d$Y, 1:length(d$Y), xlab="Route", ylab="", yaxt="n")

# Fit the model as an LMM and as an PGLS. Note that these are not the same if the estimated lambda = 1, since this is the highest
# value allowed in phylolm(), even though it can go higher in lmer(). If they aren't the same, try re-simulating the dataset.

z.lmm <- lmer(Y ~ 1 + (1|route), REML=F, data=d)
summary(z.lmm)

# This uses phylolm from {phylolm}
z.pgls <- phylolm(Y ~ 1, phy=phy.lm, model = "lambda", data=d)
summary(z.pgls)

# For the PGLS, the route-level variance is lambda*vcv[2,1]*sigma2
cov.pgls <- z.pgls$optpar * vcv[2,1] * z.pgls$sigma2
cov.lmm <- summary(z.lmm)$varcor[[1]][1]
c(cov.pgls, cov.lmm)

# This is an aside that you can skip. This uses gls from {ape} and {nlme} to give the same results as phylom(), in case you are curious. This will give an error message, but comparing the results to those for phylolm and pgls, you can see that there is no problem.
summary(gls(Y ~ 1, correlation = corPagel(0.5 , phy = phy.lm) , data = d, method="ML"))

######
# 3.4.3 Choosing a branch-length transform 

# Create a phylogeny.
n <- 100
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

# Simulate data with Pagel's lambda.
lam <- .5
phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
Y <- rTraitCont(phy.lam, model = "BM", sigma = 1)

# Fit the data with both Pagel's lambda and the OU transform.
summary(phylolm(Y ~ 1, phy=phy, model = "lambda"))
summary(phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot"))

# Simulate data with the OU transform.
alpha <- 1
phy.alpha <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree
Y <- rTraitCont(phy.alpha, model = "BM", sigma = 1)

# Fit the data with both Pagel's lambda and the OU transform.
summary(phylolm(Y ~ 1, phy=phy, model = "lambda"))
summary(phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot"))

# How often does lambda best-fit data simulated with Pagel's lambda?
n <- 100
nsim <- 200
lam <- .5
sim <- data.frame(LL.lambda=rep(NA, nsim), LL.OU=NA)
for(i in 1:nsim){
	phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
	phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
	Y <- rTraitCont(phy.lam, model = "BM", sigma = 1)

	sim$LL.lambda[i] <- phylolm(Y ~ 1, phy=phy, model = "lambda")$logLik
	sim$LL.OU[i] <- phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot")$logLik
}
# This gives the proportion of simulations in which Pagel's lambda did better.
mean(sim$LL.lambda > sim$LL.OU)

# How often does OU best-fit data simulated with the OU transform?
alpha <- 1
sim <- data.frame(LL.lambda=rep(NA, nsim), LL.OU=NA)
for(i in 1:nsim){
	phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
	phy.alpha <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree
	Y <- rTraitCont(phy.alpha, model = "BM", sigma = 1)

	sim$LL.lambda[i] <- phylolm(Y ~ 1, phy=phy, model = "lambda")$logLik
	sim$LL.OU[i] <- phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot")$logLik
}
# This gives the proportion of simulations in which OU did better.
mean(sim$LL.lambda < sim$LL.OU)

#################################################################
# 3.5. Statistical tests for phylogenetic signal
#################################################################

######
# 3.5.1 Likelihood ratio test
n <- 30
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

# Simulate data with Pagel's lambda transform.
lam <- .5
phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
Y <- rTraitCont(phy.lam, model = "BM", sigma = 1)
# Note that all of the analyses for the subsections from here to subsection 3.5.5 all use this dataset Y. The results will depend on this dataset, so the figures will not be the same as those in the book. You can rerun the last line to give different datasets from the same phylogeny, or rerun the last two lines to get a new dataset from a new phylogeny.

# Estimate phylogenetic signal with Pagel's lambda and the OU transform.
z.lam <- phylolm(Y ~ 1, phy=phy, model = "lambda")
z.OU <- phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot")
summary(z.lam)
summary(z.OU)

# Perform a Likelihood Ratio Test (LRT).
z.0 <- lm(Y ~ 1)
logLik(z.0)

LLR.lam <- 2*(z.lam$logLik - logLik(z.0)[1])
LLR.OU <- 2*(z.OU$logLik - logLik(z.0)[1])

c(LLR.lam, pchisq(LLR.lam, df=1, lower.tail=F))
c(LLR.OU, pchisq(LLR.OU, df=1, lower.tail=F))
# NOTE: It is sometimes possible to get small negative LLR values. This occurs when there is a fitting problem in the full model -- the true maximum likelihood isn't found.

######
# 3.5.2 Parametric bootstrap of the parameter

# Use the built-in bootstrap capability of phylolm() to bootstrap around the estimate.
nboot <- 2000
z.lam.boot <- phylolm(Y ~ 1, phy=phy, model = "lambda", boot=nboot, full.matrix = TRUE)
z.OU.boot <- phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot", boot=nboot, full.matrix = TRUE)

P.boot.lam <- mean(z.lam.boot$bootstrap[,3] < 0.01)
P.boot.OU <- mean(z.OU.boot$bootstrap[,3] > 49)

# Fig 3.5
par(mfrow=c(1,2), mai=c(1,1,.4,.4))
hist(z.lam.boot$bootstrap[,3], main=expression(paste("Pagel's ", lambda, " boot(", lambda,")")), xlab=expression(lambda), ylim=c(0,1000), cex.main=1.5)
lines(z.lam.boot$optpar * c(1,1), c(0,1000), col="red")
lines(mean(z.lam.boot$bootstrap[,3]) * c(1,1), c(0,1000), col="green")
text(.8,600,paste("P =", round(P.boot.lam, digits=3)), cex=1.5)

hist(z.OU.boot$bootstrap[,3], main=expression(paste("OU boot(", alpha,")")), xlab=expression(alpha), ylim=c(0,1000), cex.main=1.5)
lines(z.OU.boot$optpar * c(1,1), c(0,1000), col="red")
lines(mean(z.OU.boot$bootstrap[,3]) * c(1,1), c(0,1000), col="green")
text(40,600,paste("P =", round(P.boot.OU, digits=3)), cex=1.5)

######
# 3.5.3 Parametric bootstrap of H0

# Use lm() to fit the null model (species are independent).
z.0 <- lm(Y ~ 1)

# Perform the bootstrap using the parameters estimated under H0. The function simulate() is applied to the model z.0. The simulated datasets are fit with both Pagel's lambda and the OU transform.
nboot <- 2000
Y.sim <- Y
boot <- data.frame(LLR.lam=rep(NA,nboot), LLR.OU=NA)
for(i in 1:nboot){
	Y.sim <- simulate(z.0)[[1]]
	names(Y.sim) <- names(Y)
	z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
	z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
	z.0.sim <- lm(Y.sim ~ 1)
	boot$LLR.lam[i] <- 2*(z.lam.sim$logLik - logLik(z.0.sim))
	boot$LLR.OU[i] <- 2*(z.OU.sim$logLik - logLik(z.0.sim))
}

P.boot.LLR.lam <- mean(boot$LLR.lam > LLR.lam)
P.boot.LLR.OU <- mean(boot$LLR.OU > LLR.OU)
c(P.boot.LLR.lam, P.boot.LLR.OU)

# Due to numerical issues, sometimes phylolm doesn't converge and gives negative LLR. When this happens, I set the LLR to zero.
boot$LLR.lam[boot$LLR.lam < 0] <- 0
boot$LLR.OU[boot$LLR.OU < 0] <- 0

# Fig 3.6
par(mfrow=c(1,2), mai=c(1,1,.4,.4))
hist(boot$LLR.lam, main=expression(paste("Pagel's ", lambda, " boot(H0)")), xlab="LLR", ylim=c(0,1800), cex.main=1.5)
lines(LLR.lam * c(1,1), c(0,nboot), col="red")
lines(.1*(1:100), nboot*dchisq(.1*(1:100), df=1), col="blue")
text(5,800,paste("P =", round(P.boot.LLR.lam, digits=3)), cex=1.5)

hist(boot$LLR.OU, main=expression("OU boot(H0)"), xlab="LLR", ylim=c(0,1800), cex.main=1.5)
lines(LLR.OU * c(1,1), c(0,nboot), col="red")
lines(.1*(1:100), nboot*dchisq(.1*(1:100), df=1), col="blue")
text(8,800,paste("P =", round(P.boot.LLR.OU, digits=3)), cex=1.5)

######
# 3.5.4 Permutation test
nperm <- 2000
perm <- data.frame(lam=rep(NA,nperm), alpha=NA)
for(i in 1:nperm){
	Y.perm <- Y[sample(1:length(Y))]
	names(Y.perm) <- names(Y)
	z.lam.perm <- phylolm(Y.perm ~ 1, phy=phy, model = "lambda")
	z.OU.perm <- phylolm(Y.perm ~ 1, phy=phy, model = "OUfixedRoot")
	perm$lam[i] <- z.lam.perm$optpar
	perm$alpha[i] <- z.OU.perm$optpar
}

P.perm.lam <- mean(perm$lam > z.lam$optpar)
P.perm.OU <- mean(perm$alpha < z.OU$optpar)

# Fig 3.7
par(mfrow=c(1,2), mai=c(1,1,.4,.4))
hist(perm$lam, main=expression(paste("Pagel's ", lambda, " perm")), xlab=expression(lambda), ylim=c(0,1800), cex.main=1.5)
lines(z.lam$optpar * c(1,1), c(0,nboot), col="red")
text(.6,800,paste("P =", round(P.perm.lam, digits=3)), cex=1.5)

hist(perm$alpha, main=expression("OU perm"), xlab=expression(alpha), ylim=c(0,1800), cex.main=1.5)
lines(z.OU$optpar * c(1,1), c(0,nboot), col="red")
text(30,800,paste("P =", round(P.perm.OU, digits=3)), cex=1.5)

######
# 3.5.5 Blomberg's K
z.phylosig.K <- phylosig(Y, tree=phy, method="K", test=TRUE, nsim=2000)
z.phylosig.K$P

z.phylosig.lam <- phylosig(Y, tree=phy, method="lambda", test=TRUE)
z.phylosig.lam$P

######
# 3.5.6 Type I errors and power for tests of phylogenetic signal

# For the parametric bootstrap of the LRT under H0, it is necessary to generate the critical values of the LLR. These are placed in LLR.crit.
# Note that these bootstraps to get the critical values of LLR.lam and LLR.OU are performed for the same phylogenetic tree as used for the power simulations later. This is because the critical LLR.lam and LLR.OU could depend on the topology of the phylogeny.
n <- 30
lam <- 0
nboot <- 2000

phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree

boot0 <- data.frame(LLR.lam=rep(NA, nboot), LLR.OU=NA)
for(i in 1:nboot){
	Y.sim <- rTraitCont(phy.lam, model = "BM", sigma = 1)
	
	z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
	z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
	z.0.sim <- lm(Y.sim ~ 1)
	
	boot0$LLR.lam[i] <- 2*(z.lam.sim$logLik - logLik(z.0.sim))
	boot0$LLR.OU[i] <- 2*(z.OU.sim$logLik - logLik(z.0.sim))
}
LLR.lam.crit <- sort(boot0$LLR.lam)[floor(0.95 * nboot)]
LLR.OU.crit <- sort(boot0$LLR.OU)[floor(0.95 * nboot)]

# This plots the bootstrap distributions of LLR for Pagel's lambda and the OU transform.
par(mfrow=c(2,1))
hist(boot0$LLR.lam)
lines(.1*(1:100), nboot*dchisq(.1*(1:100), df=1), col="blue")
hist(boot0$LLR.OU)
lines(.1*(1:100), nboot*dchisq(.1*(1:100), df=1), col="blue")

# Compute the rejection rates for the 5 methods.
# Note that the same starter phylogeny is used for the simulations, since the critical values of lam and alpha used in the bootstraps around H0 could depend on the topology.
lam.list <- c(0, .2, .4, .6, .8, 1)
nsims <- 2000

reject <- data.frame(lam=rep(NA, nsims*length(lam.list)), lam.LRT=NA, OU.LRT=NA, lam.boot0=NA, OU.boot0=NA, K.perm=NA)
counter <- 0
for(lam in lam.list){
	phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
	for(i in 1:nsims){
		counter <- counter + 1
		reject$lam[counter] <- lam	
		Y.sim <- rTraitCont(phy.lam, model = "BM", sigma = 1)
		
		# LRTs
		z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
		z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
		z.0.sim <- lm(Y.sim ~ 1)
		LLR.lam <- 2*(z.lam.sim$logLik - logLik(z.0.sim)[1])
		LLR.OU <- 2*(z.OU.sim$logLik - logLik(z.0.sim)[1])
		reject$lam.LRT[counter] <- (pchisq(LLR.lam, df=1, lower.tail=F) < 0.05)
		reject$OU.LRT[counter] <- (pchisq(LLR.OU, df=1, lower.tail=F) < 0.05)
		
		# Bootstraps around H0
		reject$lam.boot0[counter] <- (LLR.lam > LLR.lam.crit)
		reject$OU.boot0[counter] <- (LLR.OU > LLR.OU.crit)

		# Permutation test with Blomberg's K
		z.phylosig.K <- phylosig(Y.sim, tree=phy, method="K", test=TRUE, nsim=2000)
		reject$K.perm[counter] <- (z.phylosig.K$P < 0.05)
	}
}
# This code took a long time to run, so I just saved the output and then reloaded it later. I've commented out these lines though.
# write.table(reject, file="Phylo signal power curves lam n=30.csv", sep=",", row.names=F)
# reject <- read.csv(file="Phylo signal power curves lam n=30.csv")

#Create power curves
power <- aggregate(reject[,c("lam.LRT","OU.LRT","lam.boot0","OU.boot0","K.perm")], by = list(reject$lam), FUN = mean)
names(power)[1] <- "lam"

# Fig. 3.8, left panel. I've left the right panel as an exercise.
par(mfrow=c(1,2))
plot(lam.LRT ~ lam, data=power, typ="l", xlim=c(0,1), ylim=c(0, 1), xlab=expression(paste("Phylogenetic signal (", lambda, ")")), ylab="Fraction rejected")
lines(OU.LRT ~ lam, data=power, col=2)
lines(lam.boot0 ~ lam, data=power, col=3)
lines(OU.boot0 ~ lam, data=power, col=4)
lines(K.perm ~ lam, data=power, col=5)

lines(c(0,10), c(.05,.05), lty=2)
legend(.0,1,c("lam.LRT","OU.LRT", "lam.boot(H0)", "OU.boot(H0)", "K.perm"), col=1:5, lty=1)

#################################################################
# 3.6 Estimating regression coefficients
#################################################################

# Generate data for a phylogenetic regression.
n <- 20
b0 <- 0
b1 <- .5
lam.x <- 1
lam.e <- .5
s2.x <- 1
s2.e <- 1

phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree

# Plot the starter trees phy.x and phy.e
par(mfrow=c(1,2))
plot(phy.x, main="phy.x")
plot(phy.e, main="phy.e")

x <- rTraitCont(phy.x, model = "BM", sigma = s2.x)
e <- rTraitCont(phy.e, model = "BM", sigma = s2.e)
x <- x[match(names(e), names(x))]
Y <- b0 + b1 * x + e
Y <- array(Y)
rownames(Y) <- phy$tip.label
	
# Fit with Pagel's lambda.
summary(phylolm(Y ~ x, phy=phy, model = "lambda"))

# Fit without phylogenetic correlations.Note that this will be the same as when fit with Pagel's lambda in the case when the estimate of lambda = 0.
summary(lm(Y ~ x))

######
# 3.6.2 Type I errors and power

# Get critical values of the LLR for the parametric bootstrap around H0 foir the case when there is no phylogenetic signal in x. These are placed in LLR.crit.
n <- 30

phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
lam.x <- 0
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
lam.e <- 1
phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree

b0 <- 0
b1 <- 0
s2.x <- 1
s2.e <- 1

lam.x.list <- c(0,1)
nboot <- 2000
LLR.crit <- data.frame(lam.x=lam.x.list, pgls=NA)
lam.x.i <- 0
for(lam.x in lam.x.list){
	lam.x.i <- lam.x.i + 1
	phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
	LLRdeviance <- data.frame(pgls=rep(NA, nboot))
	for(i in 1:nboot){
		x <- rTraitCont(phy.x, model = "BM", sigma = s2.x)
		e <- rTraitCont(phy.e, model = "BM", sigma = s2.e)
		x <- x[match(names(e), names(x))]
		y <- b1 * x + e
		y <- array(y)
		rownames(y) <- phy$tip.label
			
		mod <- phylolm(y ~ x, phy=phy, model = "lambda")
		mod0 <- phylolm(y ~ 1, phy=phy, model = "lambda")
		
		LLRdeviance$pgls[i] <- 2*(mod$logLik - mod0$logLik)[1]
	}
	LLR.crit$pgls[lam.x.i] <- sort(LLRdeviance$pgls)[floor(.95*nboot)]
	show(LLR.crit)
}
round(LLR.crit, digits=2)

# Compute the power curves as b1 ranges from 0 to 1 for the case when there is no phylogenetic signal in x.

lam.x <- 0
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree

b1range <- c(0, .2, .4, .6, .8, 1)

nsims <- 2000
reject <- data.frame(b1.true=array(0,dim=nsims*length(b1range)), P.pgls=0, P.lm=0, pgls.boot0=0)

counter <- 0
for(b1 in b1range) for(i in 1:nsims){
	counter <- counter + 1
	
	x <- rTraitCont(phy.x, model = "BM", sigma = 1)
	e <- rTraitCont(phy.e, model = "BM", sigma = 1)
	x <- x[match(names(e), names(x))]
	y <- b1 * x + e
	y <- array(y)
	rownames(y) <- phy$tip.label
		
	mod.pgls <- phylolm(y ~ x, phy=phy, model = "lambda")
	mod.pgls0 <- phylolm(y ~ 1, phy=phy, model = "lambda")
	mod.lm <- lm(y ~ x)
	
	reject$b1.true[counter] <- b1
	reject$P.pgls[counter] <- summary(mod.pgls)$coef[2,4]
	reject$P.lm[counter] <- summary(mod.lm)$coef[2,4]
	deviance <- 2*(mod.pgls$logLik - mod.pgls0$logLik)[1]
	reject$pgls.boot0[counter] <- (deviance > LLR.crit[LLR.crit$lam.x == lam.x,2])
}
# This code took a long time to run, so I just saved the output and then reloaded it later. I've commented out these lines.
# write.table(reject, file="b1 power curves lam.x=0 n=30.csv", sep=",", row.names=F)
# reject <- read.csv(file="b1 power curves lam.x=0 n=30.csv")

# Aggregate the data to get the proportion of simulation datasets for which the null hypothesis of b1=0 is rejected.
reject$pgls <- reject$P.pgls < 0.05
power.pgls <- aggregate(reject$pgls, by = list(reject$b1.true), FUN = mean)
names(power.pgls) <- c("b1", "rejected")

reject$lm <- reject$P.lm < 0.05
power.lm <- aggregate(reject$lm, by = list(reject$b1.true), FUN = mean)
names(power.lm) <- c("b1", "rejected")

power.pgls.boot0 <- aggregate(reject$pgls.boot0, by = list(reject$b1.true), FUN = mean)
names(power.pgls.boot0) <- c("b1", "rejected")

# Figure 3.9 (left panel) with lam.x = 0. The right panel of Figure 3.9 is generated with the same code but setting lam.x = 1.
par(mfrow=c(1,2))
if(lam.x == 0){
  plot(rejected ~ b1, data=power.pgls, typ="l", main=expression(paste("No phylogenetic signal in ", italic(x))), ylab="Fraction rejected", ylim=c(0, max(power.pgls$rejected, power.lm$rejected)))
  legend(.42,.45,c("PGLS","LM","PGLS.boot0"), col=1:3, lty=1)
}else{
  plot(rejected ~ b1, data=power.pgls, typ="l", main=expression(paste("Phylogenetic signal in ", italic(x))), ylab="Fraction rejected", ylim=c(0, max(power.pgls$rejected, power.lm$rejected)))
}
lines(rejected ~ b1, data=power.lm, col=2)
lines(rejected ~ b1, data=power.pgls.boot0, col=3)
lines(c(0,10), c(.05,.05), lty=2)

######
# 3.6.3 Partial R2s

# Simulate a dataset with specified strengths of phylogenetic signal in x and e.
n <- 100
b0 <- 0
b1 <- .5
lam.x <- 0
lam.e <- .8

phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree

x <- rTraitCont(phy.x, model = "BM", sigma = 1)
e <- rTraitCont(phy.e, model = "BM", sigma = 1)
x <- x[match(names(e), names(x))]
Y <- b0 + b1 * x + e
Y <- array(Y)
rownames(Y) <- phy$tip.label
	
# Fit full model with Pagel's lambda
mod.f <- phylolm(Y ~ x, phy=phy, model = "lambda")
summary(mod.f)

# Fit reduced model without x
mod.x <- phylolm(Y ~ 1, phy=phy, model = "lambda")
summary(mod.x)

# Fit reduced model without phylogenetic correlations
mod.phy <- lm(Y ~ x)
summary(mod.phy)

# Fit reduced models without x or phylogenetic correlations
mod.0 <- lm(Y ~ 1)
summary(mod.0)

# Total R2 using the function R2.lik from the rr2 package. If only one model is specified, R2.lik assumes that the reduced model is the linear model with only an intercept.
round(R2_lik(mod.f), digits=2)

# Partial R2 for phy. The first model specified in R2.lik is the full model, and the second is the reduced model.
round(R2_lik(mod.f, mod.phy), digits=2)

# Partial R2 for x.
round(R2_lik(mod.f, mod.x), digits=2)

# Total R2 for phy without x.
round(R2_lik(mod.x), digits=2)

# Total R2 for x without phy.
round(R2_lik(mod.phy), digits=2)

############################################################
# 3.7 How good must the phylogeny be?
############################################################

# The code below gives an example of permuting nodes within a phylogeny, and then does a simple permutation illustration to see the effect of topological mistakes on regression analyses. The code after that gives all the procedures used to construct the power curve (Fig. 3.10). I've included it for completeness, but unless you have some specific need, I'd stop with the simple permutation illustration.

# The node above which species are permuted among tips.
node <- 35

n <- 20
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
phy.perm <- rotateNodes(phy, node=node)
phy.perm$tip.label <- phy$tip.label

layout(matrix(c(1,2,2,2,2), ncol=1))
tdist <- treedist(phy, phy.perm)[1]/2
nd <- node.depth(phy.perm, method=2)
nd <- round((nd[node]-1)/max(nd-1), digits=2)
plot(1 ~ 1, type="n", xaxt="n", )
text(1,1,paste("RF = ", tdist, "   nd = ", nd, "   node = ", node), cex=1.5)
plot(cophylo(phy,phy.perm), no.margin=F)

# This code gives an example of the analysis in which phylogenetic signal in x is determined by lam.x.
b0 <- 0
b1 <- .1
lam.x <- 1
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree

e <- rTrait(n=1, phy=phy)
x <- rTrait(n=1, phy=phy.x)

Y <- b0 + b1*x + e

phy.perm <- rotateNodes(phy, node=node)
phy.perm$tip.label <- phy$tip.label

summary(phylolm(Y ~ x, phy=phy, model = "lambda"))
summary(phylolm(Y ~ x, phy=phy.perm, model = "lambda"))

# Simple permutation illustration for how to compare regression results for differing degrees of degradation of the topology of the phylogeny. The function treedist in the package {phangorn} computes the Robinson-Foulds distance between two phylogenies, in this case between the original phylogeny and the permuted phylogeny. This simple example with only 20 species is designed to show you how the overall analysis to produce the power curves is performed.
nrep <- 50
sim <- data.frame(dist=rep(NA,nrep), est=NA, P=NA)
for(i in 1:nrep){
	phy.perm <- rotateNodes(phy, node=sample((n+1):(2*n-2))[1])
	phy.perm$tip.label <- phy$tip.label
	dist <- treedist(phy, phy.perm)
	sim$dist[i] <- dist[1]/2

	z <- summary(phylolm(Y ~ x, phy=phy.perm, model = "lambda"))
	sim$est[i] <- z$coef[2,1]
	sim$P[i] <- z$coef[2,4]
}
par(mfrow=c(1,2), mai=c(1,1,.2,.2))
plot(est ~ log(1 + dist), data=sim, xlab="log(1 + RF distance)", ylab="Estimate of the regression coefficient")
points(0,summary(phylolm(Y ~ x, phy=phy, model = "lambda"))$coef[2,1], col="red", pch=19)
plot(P ~ log(1 + dist), data=sim, xlab="log(1 + RF distance)", ylab="P-value for the regression coefficient")
points(0,summary(phylolm(Y ~ x, phy=phy, model = "lambda"))$coef[2,4], col="red", pch=19)

# This is the code for the power curve that you should skip unless you have a desparate need to do something like this. The code generates power curves for the test of H0:b1=0 with differing levels of tree uncertainty. Phylogenetic signal in x is set by lam.x. The data.frame "reject" keeps track of the RF distance between the original phylogeny and the permuation phylogeny, and also the depth of the node that is permuted.
n <- 100
b0 <- 0
b1range <- c(0, .1, .2, .3, .4, .5)
RFrange <- c(0, 1, 2, 4, 7, 15, 100)

lam.x <- 0
lam.e <- 1

nsims <- 2000
reject <- data.frame(b1=array(0,dim=nsims*length(b1range)*length(RFrange)), tree=NA, RF.cat=0, P=NA, depth=NA, RF=NA)

counter <- 0
i.tree <- 0
for(b1 in b1range) for(i in 1:nsims){
	i.tree <- i.tree + 1
	phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
	phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
	phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree
	
	x <- rTraitCont(phy.x, model = "BM", sigma = 1)
	e <- rTraitCont(phy.e, model = "BM", sigma = 1)
	x <- x[match(names(e), names(x))]
	Y <- b0 + b1 * x + e
	Y <- array(Y)
	rownames(Y) <- phy$tip.label
		
	for(i.RF in 1:length(RFrange)){
		counter <- counter + 1
		
		reject$tree[counter] <- i.tree
		reject$RF.cat[counter] <- RFrange[i.RF]
		
		if(RFrange[i.RF] == 0) {
			phy.perm <- phy
			reject$depth[counter] <- 0

			reject$RF[counter] <- 0
			reject$path.dif[counter] <- 0
		}
		if(RFrange[i.RF] > 0 & i.RF < length(RFrange)){
			flag <- F
			i.stop <- 1
			while(!flag & i.stop <= 100){
				i.stop <- i.stop + 1
				node <- sample(size=1, x=(n+1):(2*n-2))
				phy.perm <- rotateNodes(phy, node=node)
				phy.perm$tip.label <- phy$tip.label
				tdist <- treedist(phy, phy.perm)
				flag <- (RFrange[i.RF] <= tdist[1]/2) & (tdist[1]/2 < RFrange[i.RF+1])
			}
			if(i.stop < 100){
				nd <- node.depth(phy.perm, method=2)
				reject$depth[counter] <- (nd[node] - 1)/max(nd - 1)
	
				reject$RF[counter] <- tdist[1]/2
			}
		}

		z.pgls <- phylolm(Y ~ x, phy=phy.perm, model = "lambda")
		reject$b1[counter] <- b1
		reject$P[counter] <- summary(z.pgls)$coef[2,4]

		if(i.RF == length(RFrange)) {
			phy.perm <- phy
			reject$depth[counter] <- 0
			reject$RF[counter] <- 0

			z.lm <- lm(Y ~ x)
			reject$b1[counter] <- b1
			reject$P[counter] <- summary(z.lm)$coef[2,4]
		}

		show(reject[counter,])
	}
}
			
# This code took a long time to run, so I just saved the output and then reloaded it later. I've commented out these lines though.
# write.table(reject, file="tree uncertainty n=100 lam.x=0.csv", sep=",", row.names=F)
# reject <- read.csv(file="tree uncertainty n=100 lam.x=0.csv")

# Fig. 3.10
# Separate panels are produced by changing the number of species (n) and the phylogenetic signal in x (lam.x).
par(mfrow=c(1,1))
collist <- c(1:6,8)

reject <- reject[!is.na(reject$P),]
reject$fract <- reject$P < 0.05
power <- aggregate(reject$fract, by = list(reject$b1, reject$RF.cat), FUN = mean)
names(power) <- c("b1", "RF.cat", "rejected")

plot(rejected ~ b1, data=power[power$RF.cat == RFrange[1],], typ="l", xlim=c(0,1), ylim=c(0, 1), xlab="b1", ylab="Fraction rejected")

for(i in 2:length(RFrange)) {
	lines(rejected ~ b1, data=power[power$RF.cat == RFrange[i],], col=collist[i])
}
lines(c(0,10), c(.05,.05), lty=2)
legend(.6,.7,c("RF dist","   0", "   1", "   2-3", "   4-6", "   7-14", "   15+", "   Star"), col=c(0,collist), lty=1)

max(reject$RF, na.rm=T)

############################################################
# 3.8 Phylogenetic regression for binary data
############################################################

# This is the inverse logit function used to simulate binary data.
inv.logit <- function(x){
	1/(1 + exp(-x))
}

# The code below gives an example of a single simulated dataset that is then fit using phylogenetic GLM and phylogenetic logistic regression. It is probably worth just looking at this code, rather than proceeding with the full code for generating power curves (Fig. 3.11).
n <- 100
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

b1 <- 1
lam.x <- 1
lam.e <- 1
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree

x <- rTraitCont(phy.x, model = "BM", sigma = 1)
e <- rTraitCont(phy.e, model = "BM", sigma = 1)
x <- x[match(names(e), names(x))]
theta <- b1 * x + e
p <- inv.logit(theta)
Y <- rbinom(n=n, size=1, prob=p)
d <- data.frame(Y=Y, x=x, row.names=phy$tip.label)

# Fit the data using a binary PGLMM and phylogenetic logistic regression. Note that phylglm sometimes has convergence problems when the estimate of alpha is 54 (the maximum allowable value corresponding to no phylogenetic signal in the residuals). When the estimate of alpha > 54, you should just use the regression coefficient estimates from logistf().
binaryPGLMM(Y ~ x, phy=phy, data=d)
summary(phyloglm(Y ~ x, phy=phy, method = "logistic_MPLE", data=d))

# Fit the same data using methods that don't account for phylogeny.
summary(glm(Y ~ x, family="binomial", data=d))
summary(logistf(Y ~ x, data=d))

######
# 3.8.3 Type I errors and power for binary data

# Note that this subsection code is long and complicated, so you might want to skip it. Conceptually, there is nothing really new here, and I've only included it for completeness.

# Get critical values of the LLR for the parametric bootstrap around H0. These are placed in LLR.crit.
n <- 50
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
lam.e <- 1
phy.e <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.e))$tree

b0 <- 0
b1 <- 0

lam.x.list <- c(0,1)

nboot <- 2000
b1.crit <- data.frame(lam.x=lam.x.list, pglmm.max=NA, pglmm.min=NA, plog.max=NA, plog.min=NA)
lam.x.i <- 0
for(lam.x in lam.x.list){
	lam.x.i <- lam.x.i + 1
	phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree
	reject <- data.frame(pglmm=rep(NA, nboot), plog=NA)
	for(i in 1:nboot){
		x <- rTraitCont(phy.x, model = "BM", sigma = 1)
		e <- rTraitCont(phy.e, model = "BM", sigma = 1)
		x <- x[match(names(e), names(x))]
		q <- b1 * x + e
		p <- inv.logit(q)
		Y <- rbinom(n=n, size=1, prob=p)
		d <- data.frame(Y=Y, x=x, row.names=phy$tip.label)
	
		mod.pglmm <- binaryPGLMM(Y ~ x, phy=phy, data=d, cpp = TRUE)
		mod.plog <- phyloglm(Y ~ x, phy=phy, method = "logistic_MPLE", data=d)
		
		reject$pglmm[i] <- mod.pglmm$B[2]
		reject$plog[i] <- summary(mod.plog)$coef[2,1]
	}
	b1.crit$lam.x[lam.x.i] <- lam.x
	b1.crit$pglmm.max[lam.x.i] <- sort(reject$pglmm)[floor(.975*nboot)]
	b1.crit$pglmm.min[lam.x.i] <- sort(reject$pglmm)[ceiling(.025*nboot)]
	b1.crit$plog.max[lam.x.i] <- sort(reject$plog)[floor(.975*nboot)]
	b1.crit$plog.min[lam.x.i] <- sort(reject$plog)[ceiling(.025*nboot)]
}
round(b1.crit, digits=2)

# Compute power curves
b1range <- c(0, .2, .4, .6, .8, 1)

lam.x <- 0
phy.x <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.x))$tree

nsims <- 2000
reject <- data.frame(b1.true=array(NA,dim=nsims*length(b1range)), P.pglmm=NA, P.plog=NA, P.pglmm.boot0=NA, P.plog.boot0=NA, P.pgls=NA)

counter <- 0
for(b1 in b1range) for(i in 1:nsims){
	counter <- counter + 1
	
	x <- rTraitCont(phy.x, model = "BM", sigma = 1)
	e <- rTraitCont(phy.e, model = "BM", sigma = 1)
	x <- x[match(names(e), names(x))]
	q <- b1 * x + e
	p <- inv.logit(q)
	Y <- rbinom(n=n, size=1, prob=p)
	d <- data.frame(Y=Y, x=x, row.names=phy$tip.label)

	mod.pglmm <- binaryPGLMM(Y ~ x, phy=phy, data=d, cpp = TRUE)
	mod.plog <- phyloglm(Y ~ x, phy=phy, method = "logistic_MPLE", data=d)
	mod.pgls <- phylolm(Y ~ x, phy=phy, method = "lambda", data=d)
	
	reject$b1.true[counter] <- b1
	reject$P.pglmm[counter] <- mod.pglmm$B.pvalue[2]
	reject$P.plog[counter] <- summary(mod.plog)$coef[2,4]
	reject$P.pgls[counter] <- summary(mod.pgls)$coef[2,4]
	
	crit.min <- b1.crit$pglmm.min[b1.crit$lam.x == lam.x]
	crit.max <- b1.crit$pglmm.max[b1.crit$lam.x == lam.x]
	reject$P.pglmm.boot0[counter] <- (mod.pglmm$B[2] < crit.min | mod.pglmm$B[2] > crit.max)
	
	crit.min <- b1.crit$plog.min[b1.crit$lam.x == lam.x]
	crit.max <- b1.crit$plog.max[b1.crit$lam.x == lam.x]
	reject$P.plog.boot0[counter] <- (summary(mod.plog)$coef[2,1] < crit.min | summary(mod.plog)$coef[2,1] > crit.max)
}
# This code took a long time to run, so I just saved the output and then reloaded it later. I've commented out these lines though.
# write.table(reject, file="binary b1 power curves lam.x=0 n=50.csv", sep=",", row.names=F)
# reject <- read.csv(file="binary b1 power curves lam.x=0 n=50.csv")

# Plot power curves
reject$pglmm <- reject$P.pglmm < 0.05
power.pglmm <- aggregate(reject$pglmm, by = list(reject$b1.true), FUN = mean)
names(power.pglmm) <- c("b1", "rejected")

reject$plog <- reject$P.plog < 0.05
power.plog <- aggregate(reject$plog, by = list(reject$b1.true), FUN = mean)
names(power.plog) <- c("b1", "rejected")

power.pglmm.boot0 <- aggregate(reject$P.pglmm.boot0, by = list(reject$b1.true), FUN = mean)
names(power.pglmm.boot0) <- c("b1", "rejected")

power.plog.boot0 <- aggregate(reject$P.plog.boot0, by = list(reject$b1.true), FUN = mean)
names(power.plog.boot0) <- c("b1", "rejected")

reject$pgls <- reject$P.pgls < 0.05
power.pgls <- aggregate(reject$pgls, by = list(reject$b1.true), FUN = mean)
names(power.pgls) <- c("b1", "rejected")

# Fig. 3.11
par(mfrow=c(1,2))
plot(rejected ~ b1, data=power.pglmm, typ="l", main=expression(paste("No phylogenetic signal in ", italic(x))), ylab="Fraction rejected", ylim=c(0,1))
lines(rejected ~ b1, data=power.plog, col=2)
lines(rejected ~ b1, data=power.pglmm.boot0, col=3)
lines(rejected ~ b1, data=power.plog.boot0, col=4)
lines(rejected ~ b1, data=power.pgls, col=5)
lines(c(0,10), c(.05,.05), lty=2)
legend(0,1,c("PGLMM","PLOG","PGLMM.boot0","PLOG.boot0", "PGLS"), col=1:5, lty=1)

