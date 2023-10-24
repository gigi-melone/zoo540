# Zoology/Entomology 540: Problem Set 5
#Gigi Melone, Cassie Ceballos, and Arielle Link

# Due 24 October (midnight)

	# For the homework, ONLY TURN IN THE R CODE THAT YOU USED. Start with the code from correlated_data which is pasted below. You can then add your code, remove any of the code below that you don't need for the homework questions, and save it to a separate file.

# 1. What questions do you have about the material in Ch 3? What needs more explanation? I'm serious about asking this question, because I want to improve the book. (NOTE: This requires NO NEW R CODE.)
#how do we do 3? 
#do we always need to check our (phylogenetic) model with and without the phylogenetic signal? 
#what kind of traits could you measure with phylogenetic models? Does it have to be genetic data or can it be trait data like (seedset across species)? Is it delineated genetically or spatially? Can this be used to test how far apart populations are? 

#Can you use these models to see if spatial separation contributes to phylogenetic relatedness? 

# 2. Construct the right panel of figure 3.8 giving power curves for the detection of phylogenetic signal in simulated data using an OU transform. You can modify the code I provided for the left panel which simulates data using Pagel's λ transform which is below labeled "Fig. 3.8 left panel".

#~~~~~~~~~~~~~~~~~~~~ 
n <- 30
lam<-0
alpha<-50
nboot <- 200

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
nsims <- 200

reject <- data.frame(alpha=rep(NA, nsims*length(alpha.list)), lam.LRT=NA, OU.LRT=NA, lam.boot0=NA, OU.boot0=NA, K.perm=NA)
counter <- 0
for(alpha in alpha.list){
  phy.ou <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha=alpha))$tree
  for(i in 1:nsims){
    counter <- counter + 1
    reject$alpha[counter] <- alpha
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
    z.phylosig.K <- phylosig(Y.sim, tree=phy, method="K", test=TRUE, nsim=200)
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


par(mfrow=c(1,2))
plot(lam.LRT ~ alpha, data=power, typ="l", xlim=c(0,50), ylim=c(0, 1), xlab=expression(paste("Phylogenetic signal (", alpha, ")")), ylab="Fraction rejected")
lines(OU.LRT ~ alpha, data=power, col=2)
lines(lam.boot0 ~ alpha, data=power, col=3)
lines(OU.boot0 ~ alpha, data=power, col=4)
lines(K.perm ~ alpha, data=power, col=5)


lines(c(0,10), c(.05,.05), lty=2)
legend(.0,1,c("lam.LRT","OU.LRT", "lam.boot(H0)", "OU.boot(H0)", "K.perm"), col=1:5, lty=1)

#~~~~~~~~~~~~~~~~~~~~~~

# 3. The parametric bootstrap of H0 for phylogenetic signal (subsection 3.5.3) uses the log likelihood ratio (LLR) as the test statistic. It is also possible to use the phylogenetic signal parameter (λ or α) as the test statistic. This would involve the following: (i) Fit the model to the data and calculate the value of the phylogenetic signal parameter (λ or α). (ii) Refit the model with no phylogenetic signal and use the resulting parameter values to simulate a large number of datasets. (iii) Refit the model including the phylogenetic signal parameter (λ or α) for each dataset and collect these values. (iv) The resulting distribution of λ or α estimated from the simulated datasets approximates the distribution of the estimator of λ or α, allowing P-values to be calculated. Perform this bootstrap and compare the results to the bootstrap of H0 using the LLR.

#~~~~~~~~~~~~~~~~~~~~~~
n <- 30
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

# Simulate data with Pagel's lambda transform.
lam <- .5
phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
Y <- rTraitCont(phy.lam, model = "BM", sigma = 1)
#no signal
lam0 <- 0
phy.lam0 <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam0))$tree
Y1 <- rTraitCont(phy.lam0, model = "BM", sigma = 1)

# Estimate phylogenetic signal with Pagel's lambda and the OU transform.
z.lam <- phylolm(Y ~ 1, phy=phy, model = "lambda")
z.OU <- phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot")
summary(z.lam)
summary(z.OU)

# Perform the bootstrap using the parameters estimated under H0. The function simulate() is applied to the model z.0. The simulated datasets are fit with both no phylogeetic signals
nboot <- 2000
Y.sim <- Y1
LLR.lam <- 2*(z.lam$logLik - logLik(z.0)[1])
LLR.OU <- 2*(z.OU$logLik - logLik(z.0)[1])
boot <- data.frame(lam.lam=rep(NA,nboot), alpha.OU=NA)
for(i in 1:nboot){
  Y.sim <- simulate(z.0)[[1]]
  names(Y.sim) <- names(Y1)
  z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy.lam0, model = "lambda")
  z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy.lam0, model = "OUfixedRoot")
  z.0.sim <- lm(Y.sim ~ 1)
  boot$LLR.lam[i] <- 2*(z.lam.sim$logLik - logLik(z.0.sim))
  boot$LLR.OU[i] <- 2*(z.OU.sim$logLik - logLik(z.0.sim))
}

P.boot.LLR.lam <- mean(boot$LLR.lam > LLR.lam)
P.boot.LLR.OU <- mean(boot$LLR.OU > LLR.OU)
c(P.boot.LLR.lam, P.boot.LLR.OU)

#bootstrap around the lambda and alpha from the signal
nboot <- 2000
Y.sim <- Y
LLR.lam <- 2*(z.lam$logLik - logLik(z.0)[1])
LLR.OU <- 2*(z.OU$logLik - logLik(z.0)[1])
boot <- data.frame(lam.lam=rep(NA,nboot), alpha.OU=NA)
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



# Use the built-in bootstrap capability of phylolm() to bootstrap around the estimate.
nboot <- 20
z.lam.boot <- phylolm(Y ~ 1, phy=phy, model = "lambda", boot=nboot, full.matrix = TRUE)
summary(z.lam.boot) #lambda = 0.703
z.OU.boot <- phylolm(Y ~ 1, phy=phy, model = "OUfixedRoot", boot=nboot, full.matrix = TRUE)
summary(z.OU.boot) #alpha = 3.52

P.boot.lam <- mean(z.lam.boot$bootstrap[,3] < 0.01)
P.boot.OU <- mean(z.OU.boot$bootstrap[,3] > 49)

######
# 3.5.3 Parametric bootstrap of H0

# Use lm() to fit the null model (species are independent).
z.0 <- lm(Y ~ 1)

# Perform the bootstrap using the parameters estimated under H0. The function simulate() is applied to the model z.0. The simulated datasets are fit with both Pagel's lambda and the OU transform.
nboot <- 20
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
mean(boot$LLR.lam)
mean(boot$LLR.OU)


#~~~~~~~~~~~~~~~~~~~~~~

# 4. When investigating mixed models (chapters 1 and 2), we focused on testing hypotheses concerning 
#the fixed effects (the slope). Knowing what you know about testing for phylogenetic signal
#(sections 3.4 and 3.5), what methods could you use to test hypotheses regarding the random effects in
#a mixed model? You don't need to do this in R: I'm only asking for ideas. But you will impress me 
#if you can take the models from chapter 1 and perform a test of a null hypothesis concerning a random
#effect. 
#~~~~~~~~~~~~~~~~~~~~
#here are our ideas: 
#1. use REML 
#2. you will definitely have to do something with the variances of the random effects 
#3. you could use the lme4a package, and use the lme to test for random effects 
#4. something with a likelihood ratio test :) 
#~~~~~~~~~~~~~~~~~~~~
