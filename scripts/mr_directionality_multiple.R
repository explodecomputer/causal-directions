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
	return(c(ahat, bhat, se, pval, n))
}

get_regressions <- function(X, y)
{
	nsnp <- ncol(X)
	mat <- matrix(0, nsnp, 5)
	for(i in 1:nsnp)
	{
		mat[i,] <- get_regression(X[,i], y)
	}
	mat <- as.data.frame(mat)
	names(mat) <- c("intercept", "b", "se", "pval", "n")
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


get_predicted_effects <- function(out)
{
	temp1 <- mr_ivw(out$A$b, out$B$b, out$A$se, out$B$se)
	temp2 <- mr_egger_regression(out$A$b, out$B$b, out$A$se, out$B$se)
	temp3 <- mr_weighted_median(out$A$b, out$B$b, out$A$se, out$B$se, default_parameters())

	d <- temp2$dat

	pred_ivw <- d$b_exp * temp1$b
	pred_egg <- d$b_exp * temp2$b + temp2$b_i
	pred_wme <- d$b_exp * temp3$b

	r_a <- get_r_from_pn(out$A$pval, out$A$n)
	r <- get_r_from_pn(out$B$pval, out$B$n)
	coef <- abs(r / out$B$b)
	r_ivw <- pred_ivw * coef
	r_egg <- pred_egg * coef
	r_wme <- pred_wme * coef

	return(data.frame(r_a^2, r^2, r_ivw^2, r_egg^2, r_wme^2))
}


get_r_from_pn <- function(p, n)
{
	get_p_from_r2n <- function(r2, n)
	{
		fval <- r2 * (n-2) / (1 - r2)
		pval <- pf(fval, 1, n-1, low=FALSE)
		return(pval)
	}
	optim.get_p_from_rn <- function(x, sample_size, pvalue)
	{
		abs(-log10(get_p_from_r2n(x, sample_size)) - -log10(pvalue))
	}

	if(length(p) > 1 & length(n) == 1)
	{
		message("Assuming n the same for all p values")
		n <- rep(n, length(p))
	}

	Fval <- qf(p, 1, n-1, low=FALSE)
	R2 <- Fval / (n - 2 + Fval)
	index <- !is.finite(Fval)
	if(any(index))
	{
		index <- which(index)
		for(i in 1:length(index))
		{
			R2[index[i]] <- optim(0.001, optim.get_p_from_rn, sample_size=n[index[i]], pvalue=p[index[i]])$par
		}
	}
	return(sqrt(R2))
}



dat <- make_system(20000, 40, sqrt(0.4), sqrt(0.2), sqrt(0.001), -0.05, NULL, sqrt(0.5), sqrt(0.5))
out1 <- get_effects_from_system(dat, FALSE)
out2 <- get_effects_from_system(dat, TRUE)
(res1 <- make_adjusted_plot(out1))
dev.new()
(res2 <- make_adjusted_plot(out2))



r1 <- get_predicted_effects(out1)
colSums(r1)
r2 <- get_predicted_effects(out2)
colSums(r2)

# This doesn't work...