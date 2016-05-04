library(cit)
library(reshape2)
library(ggplot2)
library(plyr)

#' 
#' Calculate fitted values directly
#'
#' @param y Response variable
#' @param  x Predictor
#' @export
#' @return Array
getFittedVals <- function(y, x)
{
	n <- length(x)
	bhat <- cov(y, x) / var(x)
	ahat <- mean(y) - bhat * mean(x)
	fitted <- ahat + x * bhat
	return(fitted)
}

#' Get residual values from OLS
#'
#'
#' @param y Response variable
#' @param x Predictor
#' @export
#' @return Array
getResiduals <- function(y, x)
{
	fitted <- getFittedVals(y, x)
	return(y - fitted)
}

# 		zA_adj <- zA
# 		zB_adj <- getResiduals(zB, zA)
# 	}
# 	return(list(zA = zA_adj, zB = zB_adj))
# }



#' Get pval directly
#' @param y Response variable
#' @param x Predictor
#' @export
#' @return p value
getPval <- function(y, x)
{
	n <- length(x)
	fitted <- getFittedVals(y, x)
	residuals <- y - fitted
	SSR <- sum((residuals - mean(residuals))^2)
	SSF <- sum((fitted - mean(fitted))^2)
	Fval <- (SSF) / (SSR/(n-2))
	pval <- pf(Fval, 1, n-2, lowe=F)
	return(pval)
}


#' Construct independent instruments
#'
#' Check that zA and zB are correlated
#' Calculate variance explained by zA and zB
#' Use the one that explains most variance as response
#' Use residuals as intrument
#' This will remove pleiotropy from instruments
#' It will also improve power
#'
#' @param A Trait A
#' @param B Trait B
#' @param zA Instrument for trait A
#' @param zB Instrument for trait B
#' @export
#' @return List of zA and zB
makeInstrumentsIndependent <- function(A, B, zA, zB)
{
	zB_adj <- getResiduals(zB, zA)
	zA_adj <- getResiduals(zA, zB)
	return(list(zA = zA_adj, zB = zB_adj))
}



#' Function to perform bi-direction 2SLS
#'
#' @param A Trait A
#' @param B Trait B
#' @param zA Instrument for trait A
#' @param zB Instrument for trait B
#' @param correction TRUE/FALSE for orthogonalising instruments
#' @export
#' @return Dataframe of 2sls
twoStageLS <- function(A, B, zA, zB, correction)
{
	if(correction)
	{
		instruments <- makeInstrumentsIndependent(A, B, zA, zB)
		zA <- instruments$zA
		zB <- instruments$zB
	}

	Ahat_zA <- getFittedVals(A, zA)
	Ahat_zB <- getFittedVals(A, zB)
	Bhat_zA <- getFittedVals(B, zA)
	Bhat_zB <- getFittedVals(B, zB)

	pvals <- array(0,4)

	pvals[1] <- getPval(A, Bhat_zA)
	pvals[2] <- getPval(A, Bhat_zB)
	pvals[3] <- getPval(B, Ahat_zA)
	pvals[4] <- getPval(B, Ahat_zB)

	nom <- c("ABA", "ABB", "BAA", "BAB")
	names(pvals) <- nom
	return(data.frame(as.list(pvals)))
}


#' Function to infer causal direction
#'
#' @param dat Output from twoStageLS
#' @param threshold Threshold
#' @export
#' @return String indicating causal direction
inferCausality <- function(dat, threshold = 0.01)
{
	dat <- subset(dat, select=c(ABA, ABB, BAA, BAB))
	sig <- unlist(dat) < threshold
	sig[is.na(sig)] <- FALSE
	AtoB <- c(TRUE, FALSE, TRUE, TRUE)
	BtoA <- c(TRUE, TRUE, FALSE, TRUE)
	conf <- c(TRUE, FALSE, FALSE, TRUE)

	if(all(sig == AtoB))
	{
		return("AtoB")
	} else if(all(sig == BtoA)) {
		return("BtoA")
	} else if(all(sig == conf)) {
		return("Confounder")
	} else {
		return("Unknown")
	}
}


makeProxy <- function(x, noise, bias)
{
	y <- x*bias + rnorm(length(x), 0, sd=noise)
	return(y)
}

testProxy <- function(x, Noise=seq(0.1, 1, by=0.1), Bias=seq(0.1, 1, by=0.1))
{
	require(lattice)
	dat <- expand.grid(noise=Noise, bias=Bias, rsq = NA)
	for(i in 1:nrow(dat))
	{
		dat$rsq[i] <- cor(x, makeProxy(x, dat$noise[i], dat$bias[i]))^2
	}
	print(wireframe(rsq ~ noise*bias, dat, drape=TRUE))
	return(dat)	
}

runSim <- function(N, NoiseT=seq(0.25, 1, by=0.25), BiasT=seq(0.25, 1, by=0.25), NoiseG=seq(0.25, 1, by=0.25), BiasG=seq(0.25, 1, by=0.25))
{
	require(plyr)
	d <- list()
	for(j in 1:length(N))
	{
		n <- N[j]
		L1 <- rbinom(n, 2, 0.5)
		L2 <- rbinom(n, 2, 0.5)
		G <- L1 + rnorm(n)
		T <- G + rnorm(n) + L2
		dat <- expand.grid(noiseT=NoiseT, biasT=BiasT, noiseG=NoiseG, biasG=BiasG, GT = NA, TG = NA, rsqG = NA, rsqT = NA, bdmr = NA)
		for(i in 1:nrow(dat))
		{
			message(i)
			G1 <- makeProxy(G, dat$noiseG[i], dat$biasG[i])
			T1 <- makeProxy(T, dat$noiseT[i], dat$biasT[i])
			dat$GT[i] <- cit.cp(L1, G1, T1)[1]
			dat$TG[i] <- cit.cp(L1, T1, G1)[1]
			dat$bdmr[i] <- inferCausality(twoStageLS(G1, T1, L1, L2, FALSE))
			dat$rsqG[i] <- cor(G1, G)^2
			dat$rsqT[i] <- cor(T1, T)^2
		}
		dat$n <- n
		d[[j]] <- dat
	}
	dat <- rbind.fill(d)
	return(dat)
}


dat <- runSim(c(100, 500, 1000, 5000, 10000))
save(dat, file="~/repo/MethylationIV/cit/results/20160504.RData")
