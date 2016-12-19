library(TwoSampleMR)
library(ggplot2)
library(tidyr)

sim_y_from_x <- function(x, r)
{
	x <- scale(x)
	stopifnot(ncol(x) == length(r))
	stopifnot(sum(r^2) <= 1)
	e <- scale(rnorm(nrow(x))) * sqrt(1-sum(r^2))
	g <- x %*% r
	y <- g + e
	return(y)
}

make_vector_effects <- function(x, r2, pleiotropy_direction=NULL, pleiotropy_cor=NULL)
{
	if(!is.null(pleiotropy_cor))
	{
		stopifnot(pleiotropy_cor <= 1 & pleiotropy_cor >= -1)
		nsnp <- length(x)
		eff <- sim_y_from_x(x, pleiotropy_cor)
	} else {
		nsnp <- x
		eff <- rnorm(nsnp)
	}
	eff2 <- sign(eff) * eff^2
	eff2 <- eff2 / sum(abs(eff2)) * r2
	eff <- sign(eff2) * sqrt(abs(eff2))
	if(!is.null(pleiotropy_direction))
	{
		eff <- eff + pleiotropy_direction * sign(x)
	}
	return(eff)
}

make_proxy <- function(x, noise, bias)
{
	y <- x*bias + rnorm(length(x), 0, sd=noise)
	return(y)
}


make_system <- function(nid, nsnp, cor_A_B, cor_Z_A, cor_Z_B, pleiotropy_direction, cor_pleiotropy, noisea, noiseb)
{

	g <- matrix(rbinom(nid * nsnp, 2, 0.5), nid, nsnp)

	a <- make_vector_effects(nsnp, cor_Z_A^2)
	b <- make_vector_effects(a, cor_Z_B, pleiotropy_direction, cor_pleiotropy)
	cor(a,b)

	A <- sim_y_from_x(g, a)
	B <- sim_y_from_x(cbind(A, g), c(cor_A_B, b))

	Ap <- make_proxy(A, noisea, 1)
	Bp <- make_proxy(B, noiseb, 1)

	return(data.frame(
		Z = g,
		A = A,
		B = B,
		Ap = Ap,
		Bp = Bp
	))
}


get_regression <- function(x, y)
{
	n <- length(x)
	bhat <- cov(y, x) / var(x)
	ahat <- mean(y) - bhat * mean(x)
	fitted <- ahat + x * bhat
	residuals <- y - fitted
	SSR <- sum((residuals - mean(residuals))^2)
	SSF <- sum((fitted - mean(fitted))^2)
	Fval <- (SSF) / (SSR/(n-2))
	pval <- pf(Fval, 1, n-2, lowe=F)
	se <- sqrt(sum(residuals^2) / ((n - 2) * sum((x-mean(x))^2)))
	return(c(ahat, bhat, se, pval))
}

get_regressions <- function(X, y)
{
	nsnp <- ncol(X)
	mat <- matrix(0, nsnp, 4)
	for(i in 1:nsnp)
	{
		mat[i,] <- get_regression(X[,i], y)
	}
	mat <- as.data.frame(mat)
	names(mat) <- c("intercept", "b", "se", "pval")
	return(mat)
}

get_effects_from_system <- function(dat, reverse=FALSE)
{
	inst_names <- grep("Z.", names(dat))
	l <- list(
		A = get_regressions(dat[,inst_names], dat$A),
		B = get_regressions(dat[,inst_names], dat$B),
		Ap = get_regressions(dat[,inst_names], dat$Ap),
		Bp = get_regressions(dat[,inst_names], dat$Bp)
	)
	if(reverse)
	{
		names(l) <- gsub("A", "C", names(l))
		names(l) <- gsub("B", "A", names(l))
		names(l) <- gsub("C", "B", names(l))
	}
	return(l)
}

make_adjusted_plot <- function(out)
{	
	temp1 <- mr_ivw(out$A$b, out$B$b, out$A$se, out$B$se)
	temp2 <- mr_egger_regression(out$A$b, out$B$b, out$A$se, out$B$se)
	temp3 <- mr_weighted_median(out$A$b, out$B$b, out$A$se, out$B$se, default_parameters())

	d <- temp2$dat

	d1 <- data.frame(
		b_exp = d$b_exp,
		b_out = d$b_out,
		b_out_adj = d$b_exp * temp1$b,
		method = "IVW",
		id = 1:nrow(d)
	)	
	d2 <- data.frame(
		b_exp = d$b_exp,
		b_out = d$b_out,
		b_out_adj = d$b_exp * temp2$b + temp2$b_i,
		method = "MR Egger",
		id = 1:nrow(d)
	)	
	d3 <- data.frame(
		b_exp = d$b_exp,
		b_out = d$b_out,
		b_out_adj = d$b_exp * temp3$b,
		method = "Weighted median",
		id = 1:nrow(d)
	)	
	d <- rbind(d1, d2, d3)
	d <- gather(d, key, value, b_out, b_out_adj)

	s <- data.frame(intercept = c(0, temp2$b_i, 0),
		slope = c(temp1$b, temp2$b, temp3$b),
		pval = c(temp1$pval, temp2$pval, temp3$pval),
		method=c("IVW", "MR Egger", "Weighted median")
	)


	p <- ggplot(d, aes(x=b_exp, y=value)) +
	geom_line(aes(group=id)) +
	facet_grid(. ~ method) +
	geom_point(data=subset(d, key=="b_out")) +
	geom_abline(data=s, aes(intercept=intercept, slope=slope), linetype="dotted") +
	labs(x="Effect on exposure", y="Effect on outcome")
	return(list(s=s, p=p))
}


dat <- make_system(20000, 40, sqrt(0.4), sqrt(0.2), sqrt(0.001), -0.05, NULL, sqrt(0.5), sqrt(0.5))
out1 <- get_effects_from_system(dat, FALSE)
out2 <- get_effects_from_system(dat, TRUE)
make_adjusted_plot(out1)
make_adjusted_plot(out2)
dev.new()




head(dat)


r = cov(x,y) / sd(x)sd(y)
b = cov(x,y) / var(x)

pval + n -> r

b = sd(x)sd(y) / r / var(x)



