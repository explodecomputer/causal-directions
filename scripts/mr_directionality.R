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


make_geno <- function(n, p)
{
	stopifnot(all(p >= 0 & p <= 1))
	geno <- drop(sapply(p, function(x) rbinom(n, 2, x)))
	return(geno)
}


#' Create a phenotype using arbitrary number of known causal inputs
#'
#' For a set of \code{x} variables and effect sizes for each variable (\code{b}) 
#' \code{y} is constructed such that 
#' \code{y = Xb + e}
#' Given that the variance explained in \code{y} by \code{X}, is
#' \code{r^2 = sum(b * var(x) / (sqrt(x)*sqrt(y)))}
#' we can model \code{e ~ N(0, 1 - r^2)}
#'
#' @param effs Array of beta values for each input. Leave the vx and vy values to default (=1) to allow effs to be equal to the correlation between y and each x
#' @param indep Matrix of independent variables corresponding to effs
#' @param vy=1 The output variance of y
#' @param vx=rep(1, length(effs)) The desired scaled variance of x
#'
#' @export
#' @return Numeric array, simulated phenotype
#'
#' @examples \dontrun{
#' g1 <- make_geno(1000, 0.5)
#' g2 <- make_geno(1000, 0.3)
#' x1 <- rnorm(1000)
#' x2 <- rnorm(1000)
#' y <- make_phen(effs=c(0.2, 0.1, 0.15, 0.4), cbind(g1, g2, x1, x2))
#' 
#'}
make_phen <- function(effs, indep, vy=1, vx=rep(1, length(effs)))
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy))
	return(y)
}


make_system <- function(n, p, r_ab, r_za, noisea, noiseb)
{
	Z <- make_geno(n, p)
	A <- make_phen(r_za, Z)
	B <- make_phen(r_ab, A)
	Ap <- makeProxy(A, noisea, 1)
	Bp <- makeProxy(B, noiseb, 1)
	return(data.frame(Z, A, B, Ap, Bp))
}

# Test A causes B
TSLS <- function(A, B, Z)
{
	ahat <- getFittedVals(A, Z)
	summary(lm(B ~ ahat))
}


parameters <- expand.grid(
	n = c(100, 1000, 10000),
	p = 0.5,
	r_ab = sqrt(seq(0, 1, by=0.2)),
	r_za = c(sqrt(0.01), sqrt(0.05), sqrt(0.1)),
	noisea = sqrt(seq(0, 1, by=0.2)),
	noiseb = sqrt(seq(0, 1, by=0.2)),
	nsim = 1:100
)

for(i in 1:nrow(parameters))
{
	message(i)
	dat <- with(parameters[i,], make_system(n, p, r_ab, r_za, noisea, noiseb))
	parameters$cit_AB <- cit.cp(dat$Z, dat$Ap, dat$Bp)
	parameters$cit_BA <- cit.cp(dat$Z, dat$Bp, dat$Ap)
	parameters$cor_aap[i] <- cor(dat$A, dat$Ap)
	parameters$cor_bbp[i] <- cor(dat$B, dat$Bp)
	parameters$cor_abp[i] <- cor(dat$Ap, dat$Bp)
	parameters$cor_zap[i] <- cor(dat$Ap, dat$Z)
	parameters$cor_zbp[i] <- cor(dat$Bp, dat$Z)
	parameters$p_az[i] <- getPval(dat$Ap, dat$Z)
	parameters$p_bz[i] <- getPval(dat$Bp, dat$Z)
}

parameters$p_az <- -log10(parameters$p_az)
parameters$p_bz <- -log10(parameters$p_bz)
parameters$sig_mr <- parameters$p_bz > -log10(0.05)
parameters$correct_direction <- parameters$p_az > parameters$p_bz
parameters$bigna <- parameters$noisea > parameters$noiseb

save(parameters, file="~/repo/cit_measurement_error/results/p_comp_20160208.RData")
