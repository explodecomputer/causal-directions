## ---- load_up ----

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))


get_p_from_r2n <- function(r2, n)
{
	fval <- r2 * (n-2) / (1 - r2)
	pval <- pf(fval, 1, n-1, low=FALSE)
	return(pval)
}

get_r_from_pn <- function(p, n)
{
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


steiger_sensitivity <- function(rgx_o, rgy_o, ...)
{
	if(rgy_o > rgx_o)
	{
		a <- rgy_o
		b <- rgx_o
	} else {
		a <- rgx_o
		b <- rgy_o
	}

	d <- expand.grid(rxx_o=seq(rgx_o,1,length.out=50), ryy_o=seq(rgy_o,1,length.out=50), type=c("A","B"))
	d$rgy <- rgy_o / d$ryy_o
	d$rgx <- rgx_o / d$rxx_o
	d$z <- d$rgy - d$rgx
	d$z[d$type=="A"] <- 0
	temp <- wireframe(
		z ~ rxx_o * ryy_o, 
		groups=type, 
		data=d, 
		scales=list(arrows=FALSE), 
		col.groups = colorRampPalette(c("red", "blue"))(2), 
		drape=FALSE, 
		ylab=expression(rho[xx[o]]), 
		xlab=expression(rho[yy[o]]),
		zlab=expression(rho[gy]-rho[gx]),
		par.settings = list(axis.line=list(col="transparent")),
		...
	)

	vz <- a * log(a) - b * log(b) + a*b*(log(b)-log(a))
	vz0 <- -2*b - b * log(a) - a*b*log(a) + 2*a*b

	vz1 <- abs(vz - vz0)

	sensitivity <- vz0 / (2 * vz0 + abs(vz))
	sensitivity_ratio <- vz1 / vz0

	return(list(
		vz = vz,
		vz0 = vz0,
		vz1 = vz1,
		sensitivity = sensitivity,
		sensitivity_ratio = sensitivity_ratio,
		pl = temp
	))
}


mr_steiger <- function(p_exp, p_out, n_exp, n_out, r_xxo = 1, r_yyo=1, ...) 
{
	require(psych)
	index <- any(is.na(p_exp)) | any(is.na(p_out)) | any(is.na(n_exp)) | any(is.na(n_out))
	p_exp <- p_exp[!index]
	p_out <- p_out[!index]
	n_exp <- n_exp[!index]
	n_out <- n_out[!index]
	r_exp <- get_r_from_pn(p_exp, n_exp)
	r_out <- get_r_from_pn(p_out, n_out)

	r_exp_adj <- sqrt(r_exp^2 / r_xxo^2)
	r_out_adj <- sqrt(r_out^2 / r_yyo^2)

	sensitivity <- steiger_sensitivity(r_exp, r_out, ...)

	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	rtest_adj <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp_adj, r34 = r_out_adj)
	l <- list(
		r2_exp = r_exp^2, 
		r2_out = r_out^2, 
		r2_exp_adj = r_exp_adj^2, 
		r2_out_adj = r_out_adj^2, 
		correct_causal_direction = r_exp > r_out, 
		steiger_test = rtest$p,
		correct_causal_direction_adj = r_exp_adj > r_out_adj, 
		steiger_test_adj = rtest_adj$p,
		vz = sensitivity$vz,
		vz0 = sensitivity$vz0,
		vz1 = sensitivity$vz1,
		sensitivity = sensitivity$sensitivity,
		sensitivity_ratio = sensitivity$sensitivity_ratio,
		sensitivity_plot = sensitivity$pl
	)
	return(l)
}

mr_wald_ratio <- function(b_exp, b_out, se_exp, se_out, parameters) 
{
	if (length(b_exp) > 1) {
		return(list(b = NA, se = NA, pval = NA, nsnp = NA))
	}
	b <- b_out/b_exp
	se <- se_out/abs(b_exp)
	pval <- pnorm(abs(b)/se, lower.tail = F) * 2
	return(list(b = b, se = se, pval = pval, nsnp = 1))
}


## ---- read_data ----

load("~/repo/cit_measurement_error/results/mr_directionality_20160211.RData")

parameters <- subset(parameters, r_ab != 1)

parameters$p_az <- -log10(parameters$p_az)
parameters$p_bz <- -log10(parameters$p_bz)
parameters$r_ab <- parameters$r_ab^2
parameters$r_ab <- round(parameters$r_ab, 1)
parameters$r_za <- parameters$r_za^2
parameters$noisea <- parameters$noisea^2
parameters$noiseb <- parameters$noiseb^2
parameters$rhs <- abs(parameters$cor_ab * parameters$cor_bbp)
parameters$lhs <- abs(parameters$cor_aap)
parameters$rhs_lhs_diff <- parameters$lhs - parameters$rhs

# Calculate correlation test
a <- r.test(n=parameters$n, r12=parameters$cor_zap, r13=parameters$cor_zbp, r23=parameters$cor_abp)
parameters$cortest_t <- a$t
parameters$cortest_p <- -log10(a$p)
parameters$cortest_correct_direction <- sign(parameters$cortest_t) == 1
parameters$cortest_correct_direction[is.na(parameters$cortest_correct_direction)] <- FALSE


# Get the p-value for the MR analysis
# Use the cortest_t to determine which direction MR association to use
parameters$p_test <- parameters$p_bz
parameters$p_test[!parameters$cortest_correct_direction] <- parameters$p_az[!parameters$cortest_correct_direction]


# Get the CIT direction

get_cit_direction <- function(p1, p2, thresh)
{
	res <- array(0, length(p1))
	res[p1 <= thresh & p2 > thresh] <- 1
	res[p1 > thresh & p2 <= thresh] <- 2
	res[p1 <= thresh & p2 <= thresh] <- 3
	res[p1 > thresh & p2 > thresh] <- 4
	return(res)
}

parameters$cit_res <- get_cit_direction(parameters$cit_AB, parameters$cit_BA, 0.05)
parameters$cit_correct_direction <- parameters$cit_res == 1

# Get the p-value for the CIT direction
# ab < 0.05 & ba > 0.05 = correct
# ab > 0.05 & ba < 0.05 = incorrect
# ab < 0.05 & ba < 0.05 = no call
# ab > 0.05 & ba > 0.05 = incorrect


parameters$cit_p <- parameters$cit_AB
parameters$cit_p[parameters$cit_res == 2] <- parameters$cit_BA[parameters$cit_res == 2]
parameters$cit_p[parameters$cit_res == 3] <- 0.5
parameters$cit_p[parameters$cit_res == 4] <- 0.01

# Do we classify 4 as being significant wrong causal direction, or no causal inference?


parameters$cit_p <- -log10(parameters$cit_p)
parameters$cit_p[is.na(parameters$cit_AB) | is.na(parameters$cit_BA)] <- NA

# Gather
pl <- gather(parameters, test, direction_p_value, cit_p, cortest_p)
pl$correct_direction <- pl$cortest_correct_direction
pl$correct_direction[pl$test == "cit_p"] <- pl$cit_correct_direction[pl$test == "cit_p"]
pl$test_p_value <- pl$p_test
pl$test_p_value[pl$test == "cit_p"] <- pl$direction_p_value[pl$test == "cit_p"]
pl <- subset(pl, select=-c(cortest_correct_direction, cit_correct_direction, p_test))
pl$test[pl$test == "cit_p"] <- "CIT"
pl$test[pl$test == "cortest_p"] <- "MR"

# Fix NAs
pl$test_p_value[is.infinite(pl$test_p_value)] <- max(is.finite(pl$test_p_value))
pl$direction_p_value[is.infinite(pl$direction_p_value)] <- max(is.finite(pl$direction_p_value))

# Reduce
psum1 <- pl %>%
	group_by(n, p, r_ab, r_za, noisea, noiseb, test) %>%
	dplyr::summarise(
		prop_correct_direction=sum(correct_direction, na.rm=TRUE)/n(),
		prop_correct_direction_sig=sum(correct_direction & direction_p_value > -log10(0.05) & test_p_value > -log10(0.05), na.rm=TRUE)/n(),
		cor_aap=mean(cor_aap, na.rm=TRUE),
		cor_bbp=mean(cor_bbp, na.rm=TRUE),
		cor_abp=mean(cor_abp, na.rm=TRUE),
		cor_ab=mean(cor_ab, na.rm=TRUE),
		cor_zap=mean(cor_zap, na.rm=TRUE),
		cor_zbp=mean(cor_zbp, na.rm=TRUE),
		nsim=n(),
		test_sig = sum(test_p_value > -log10(0.05), na.rm=TRUE)/n(),
		direction_sig = sum(direction_p_value > -log10(0.05), na.rm=TRUE)/n(),
		rhs=mean(rhs),
		lhs=mean(cor_aap),
		rhs_lhs_diff=lhs-rhs,
		prop_sig_incorrect=sum(!correct_direction & direction_p_value > -log10(0.05) & test_p_value > -log10(0.05), na.rm=T)/n(),
		prop_sig_correct=sum(correct_direction & direction_p_value > -log10(0.05) & test_p_value > -log10(0.05), na.rm=T)/n(),
		prop_nonsig=sum(direction_p_value <= -log10(0.05) | test_p_value <= -log10(0.05), na.rm=T)/n()
	)
psum1$rhs_lhs_diff_bin <- cut(psum1$rhs_lhs_diff, breaks=seq(-0.2, 1.1, by=0.1))


## ---- cit_measurement_error_figure ----

load("~/repo/cit_measurement_error/results/20141114.RData")
dat <- gather(dat, eval, pval, GT, TG, factor_key=TRUE)
levels(dat$eval) <- c("Correct causal model", "Incorrect causal model")

ggplot(subset(dat, n %in% c(100, 1000, 10000)), aes(x=sqrt(rsq),y=-log10(pval))) +
geom_point(aes(colour=eval)) +
labs(y=expression(-log[10]*p), x=expression(cor(x, x[O])), colour="") +
facet_grid(n ~ ., scale="free_y")


## ---- d_relationship_figure ----

# There is something that i don't understand about this
# Replacing with steiger_sensitivity chunk below

ineq <- expand.grid(
	cora = seq(0, 1, by=0.02),
	corb = seq(0, 1, by=0.2),
	ab = c(0.1, 0.5, 0.9)
)

ineq$lhs <- ineq$cora
ineq$rhs <- ineq$corb * ineq$ab
ineq$d <- ineq$lhs - ineq$rhs
ineq$ablab <- paste0("cor(x,y) = ", ineq$ab)


area <- group_by(ineq, ablab) %>%
	do({

		temp1 <- subset(., corb == max(corb) & d <= 0)

		y1 <- min(temp1$d)
		x1 <- min(temp1$cora)

		y2 <- min(abs(temp1$d), 0)
		x2 <- max(temp1$cora)

		temp2 <- subset(., corb == min(corb))
		y3 <- min(abs(temp2$d), 0)
		x3 <- min(temp2$cora)

		data.frame(cora=c(x1, x2, x3), d=c(y1, y2, y3), f=1)
	})


p1 <- ggplot(ineq, aes(x=cora, y=d)) +
geom_polygon(data=area, fill="grey", aes(x=cora, y=d)) +
geom_line(aes(colour=corb, group=corb)) +
geom_hline(yintercept=0, linetype=2) +
facet_wrap(~ ablab) +
labs(x=TeX("$cor(x, x_o)$"), y="d", colour=TeX("$cor(y, y_o)$"))

p2 <- mr_steiger(
	p_exp = get_p_from_r2n(0.1, 1000),
	p_out = get_p_from_r2n(0.1 * 0.9, 1000),
	n_exp = 1000,
	n_out = 1000
)$sensitivity_plot
grid.arrange(p1, p2, ncol=1)


## ---- cit_mr_comparison_figure ----

psum2 <- gather(psum1, eval, value, prop_sig_correct, prop_sig_incorrect, prop_nonsig, factor_key=TRUE) %>%
	dplyr::group_by(n, rhs_lhs_diff_bin, r_za, eval, test) %>%
	dplyr::summarise(value=mean(value, na.rm=TRUE), nsim=n())

temp2 <- do.call(rbind, strsplit(as.character(psum2$rhs_lhs_diff_bin), split=","))
psum2$rhs_lhs_diff_bin_numeric <- as.factor(gsub("\\(", "", temp2[,1]))
psum2$rhs_lhs_diff_bin_lab <- as.factor(paste0("d = ", psum2$rhs_lhs_diff_bin_numeric))
psum2$rhs_lhs_diff_bin_lab <- factor(psum2$rhs_lhs_diff_bin_lab, levels=levels(psum2$rhs_lhs_diff_bin_lab)[order(as.numeric(as.character(levels(psum2$rhs_lhs_diff_bin_numeric))))])

levels(psum2$eval) <- c("Evidence for causality\nunder correct model", "Evidence for causality\nunder incorrect model", "No evidence for causality")

ggplot(subset(psum2, round(r_za,2)==0.01), aes(x=test, y=value)) +
geom_bar(stat="identity", aes(fill=eval)) +
facet_grid(n ~ rhs_lhs_diff_bin) +
scale_fill_brewer(type="qual") +
labs(x="Method", y="Proportion of simulations", fill="") +
theme(legend.key=element_rect(size=3), legend.key.size = unit(2, "lines"))



## ---- plot1 ----

temp <- subset(psum1, n==1000 & r_ab != 0)
temp$r_za <- paste0("cor(A,Z)=",temp$r_za^2)
temp$r_ab <- paste0("cor(A,B)=",temp$r_ab^2)
temp$noisea <- temp$noisea^2
temp$noiseb <- temp$noiseb^2
wireframe(correct_direction_sig ~ noisea*noiseb | as.factor(r_ab) + as.factor(r_za), data=temp, drape=TRUE)


## ---- plot2 ----

ineq <- expand.grid(
	noisea = seq(0, 1, by=0.02),
	noiseb = seq(0, 1, by=0.02)
)

ineq$ab <- ineq$noisea / ineq$noiseb
ineq$ab[ineq$ab > 1] <- NA

wireframe(ab ~ noisea * noiseb, ineq, drape=TRUE, screen = list(x=-30), scales = list(arrows = FALSE))


## ---- plot4 ----

# how does test significance change with rhs_lhs_diff

ggplot(subset(psum1, r_za==0.1 & r_ab != 0), aes(x=rhs_lhs_diff, y=test_sig)) +
geom_point(aes(colour=as.factor(test)), size=0.2) +
geom_smooth(se=FALSE, aes(colour=as.factor(test))) +
facet_grid(n ~ r_ab)



## ---- plot5 ----

# how does direction significance change with rhs_lhs_diff

ggplot(subset(psum1, r_za==0.1 & r_ab != 0), aes(x=rhs_lhs_diff, y=direction_sig)) +
geom_point(aes(colour=as.factor(test)), size=0.2) +
geom_smooth(se=FALSE, aes(colour=as.factor(test))) +
facet_grid(n ~ r_ab)




ggplot(subset(psum1, r_za==0.1 & r_ab != 0), aes(x=correct_direction, y=direction_sig)) +
geom_point(aes(colour=as.factor(test)), size=0.2) +
geom_smooth(se=FALSE, aes(colour=as.factor(test))) +
facet_grid(n ~ rhs_lhs_diff_bin)






## ---- plot7 ----


ggplot(subset(parameters, r_za==0.1 & r_ab==0.6 & noiseb==0), aes(x=cor_aap, y=cit_p)) +
geom_point(aes(colour=cit_correct_direction)) +
facet_grid(n ~ ., scale="free_y")


temp <- subset(parameters, r_za==0.1 & r_ab==0.6 & noiseb==0)










## ---- misc ----



ggplot(subset(psum1, r_za==0.1 & n == 1000), aes(y=prop_cortest_sign, x=prop_cortest_sig)) +
geom_point(aes(colour=as.factor(r_ab))) +
facet_grid(noisea ~ noiseb)



ggplot(subset(parameters, sig_mr), aes(x=rhs_lhs_diff, y=cortest_p)) +
geom_point(aes(colour=cortest_correct_direction)) +
facet_grid(n ~ r_za, scale="free_y")



ggplot(subset(parameters, sig_mr), aes(x=rhs_lhs_diff, y=-log10(cit_p))) +
geom_point(aes(colour=cit_correct_direction)) +
facet_grid(n ~ r_za, scale="free_y")


a <- subset(parameters, select=c(n, rhs_lhs_diff, cit_correct_direction, cit_p))
b <- subset(parameters, select=c(n, rhs_lhs_diff, cortest_correct_direction, cortest_p))
a$test <- "CIT"
b$test <- "MR"
names(a) <- names(b) <- c("n", "rhs_lhs_diff", "direction", "pval", "test")
ab <- rbind(a,b)

ggplot(ab, aes(x=rhs_lhs_diff, y=pval)) +
geom_point(aes(colour=direction)) +
facet_grid(test ~ n, scale="free_y")



psum2 <- parameters %>%
	filter(sig_mr) %>%
	group_by(n, p, r_ab, r_za, noisea, noiseb) %>%
	summarise(
		correct_direction=sum(correct_direction)/n(),
		sig_mr=sum(sig_mr)/n(),
		cor_aap=mean(cor_aap),
		cor_bbp=mean(cor_bbp),
		cor_abp=mean(cor_abp),
		sig_mr=sum(p_test > -log10(0.05))/n()
	)


test_cor_diff_1sample <- function(A, B, Z)
{
	r12 <- cor(A, Z)
	r13 <- cor(B, Z)
	r23 <- cor(A, B)
	
}

test_cor_diff_2sample <- function(A, B, Za, Zb)
{

}


n <- 10000
a <- rnorm(n)
b <- rnorm(n) + a

get_cor_from_pval(summary(lm(a ~ b))$coefficients[2,4], n)
cor(a,b)^2


qf(0.209, 1, 999, low=FALSE)


# compare cor_zap with cor_zbp
# see how often you get a significant result for the wrong direction?




ggplot(subset(parameters, r_za==0.1), aes(y=p_test, x=cor_aap-cor_bbp)) +
geom_point(aes(colour=correct_direction)) +
facet_grid(n ~ r_ab)



ggplot(subset(psum1, r_za == 0.1), aes(y=correct_direction, x=cor_abp)) +
geom_point(aes(colour=as.factor(n))) +
facet_grid(noisea ~ noiseb)


ggplot(subset(psum2, n==10000), aes(x=cor_abp, y=correct_direction)) +
geom_point(aes(colour=as.factor(r_za))) +
geom_line(aes(colour=as.factor(r_za))) +
facet_grid(noisea ~ noiseb, scale="free_y")



ggplot(parameters, aes(x=cor_abp, y=p_az - p_bz)) +
geom_point(aes(colour=as.factor(r_za))) +
geom_smooth(aes(colour=as.factor(r_za))) +
facet_grid(n ~ bigna, scale="free_y") +
scale_colour_brewer()






psum3 <- parameters %>%
	group_by(n, p, r_ab, r_za, noisea, noiseb) %>%
	summarise(
		correct_direction=sum(correct_direction)/n(),
		sig_mr=sum(sig_mr)/n(),
		cor_aap=mean(cor_aap),
		cor_bbp=mean(cor_bbp),
		cor_abp=mean(cor_abp),
		sig_mr_correct=sum(p_test > -log10(0.05) & correct_direction)/n(),
		sig_mr_incorrect=sum(p_test > -log10(0.05) & !correct_direction)/n(),
		non_sig_mr=sum(p_test <= -log10(0.05))/n()
	) %>%
	gather(outcome, value, sig_mr_correct, sig_mr_incorrect, non_sig_mr)

ggplot(subset(psum3, n == 10000 & r_ab > 0.23 & r_za==0.1), aes(y=value, x=cor_abp)) +
geom_point(aes(colour=outcome)) +
geom_line(aes(colour=outcome)) +
scale_fill_brewer(type="qual") +
facet_grid(noisea ~ noiseb)



## ---- causality_exists_tpr ----

pl$causality_exists <- TRUE
pl$causality_exists[pl$test == "MR"] <- pl$p_bz[pl$test == "MR"] > -log10(0.05)
pl$causality_exists[pl$test == "CIT"] <- pl$cit_res[pl$test == "CIT"] <= 3
pl$rhs_lhs_diff_bin <- cut(pl$rhs_lhs_diff, breaks=10)
temp2 <- do.call(rbind, strsplit(as.character(pl$rhs_lhs_diff_bin), split=","))
pl$rhs_lhs_diff_bin_numeric <- as.factor(gsub("\\(", "", temp2[,1]))
pl$rhs_lhs_diff_bin_lab <- as.factor(paste0("d = ", pl$rhs_lhs_diff_bin_numeric))
pl$rhs_lhs_diff_bin_lab <- factor(pl$rhs_lhs_diff_bin_lab, levels=levels(pl$rhs_lhs_diff_bin_lab)[order(as.numeric(as.character(levels(pl$rhs_lhs_diff_bin_numeric))))])

temp <- dplyr::group_by(pl, test, r_ab, r_za, n, noisea, noiseb) %>%
	dplyr::summarise(prop_causality_exists = sum(causality_exists) / n())
temp$r_za <- paste0("cor(x,g) = ", temp$r_za)
temp$n <- paste0("n = ", temp$n)
temp$noiseb <- paste0("cor(Y,Yo) = ", temp$noiseb)

ggplot(subset(temp, r_ab == 0.6 & r_za == "cor(x,g) = 0.01" & noiseb %in% c("cor(Y,Yo) = 0", "cor(Y,Yo) = 0.6")), aes(x=noisea, y = prop_causality_exists)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(noiseb ~ n) +
labs(y="True positive rate", x = expression(cor(X, X[O])), fill="Test")


## ----causality_exists_other ----

ggplot(subset(temp, r_ab == 0 & noiseb==0), aes(x=noisea, y = prop_causality_exists)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za)

ggplot(subset(temp, r_ab == 0 & noisea==0), aes(x=noiseb, y = prop_causality_exists)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za)

ggplot(subset(temp, r_ab == 0.2 & noisea==0), aes(x=noiseb, y = prop_causality_exists)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za)

ggplot(subset(temp, r_ab == 0.2 & noisea > 0.2), aes(x=noiseb, y = prop_causality_exists)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za)

ggplot(subset(temp, r_ab == 0.8), aes(x=rhs_lhs_diff_bin_numeric, y = prop_causality_exists)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za)



pl$causality_exists_correct <- TRUE
pl$causality_exists_correct[pl$test == "MR"] <- pl$p_bz[pl$test == "MR"] > -log10(0.05) & pl$correct_direction[pl$test == "MR"] & pl$direction_p_value[pl$test == "MR"] > -log10(0.05)
pl$causality_exists_correct[pl$test == "CIT"] <- pl$cit_res[pl$test == "CIT"] == 1


pl$cd <- TRUE
pl$cd[pl$test == "MR"] <- sign(pl$cortest_t[pl$test == "MR"]) == 1
pl$cd[pl$test == "CIT"] <- pl$cit_AB[pl$test == "CIT"] < pl$cit_BA[pl$test == "CIT"]

pl <- pl[order(pl$test, pl$r_ab, pl$r_za, pl$n, pl$noisea, pl$noiseb), ]
plmr <- subset(pl, test=="MR")
plcit <- subset(pl, test=="CIT")

all(plmr$r_ab == plcit$r_ab)
all(plmr$r_za == plcit$r_za)
all(plmr$n == plcit$n)
all(plmr$noisea == plcit$noisea)
all(plmr$noiseb == plcit$noiseb)

plmr$method_direction_comparison <- ""
plmr$method_direction_comparison[plmr$cd & plcit$cd] <- "Both"
plmr$method_direction_comparison[!plmr$cd & plcit$cd] <- "CIT"
plmr$method_direction_comparison[plmr$cd & !plcit$cd] <- "MR"
plmr$method_direction_comparison[!plmr$cd & !plcit$cd] <- "Neither"
plmr$method_direction_comparison_num <- plmr$cd - plcit$cd

table(plmr$method_direction_comparison)
with(plmr, tapply(rhs_lhs_diff, method_direction_comparison, range))
table(plmr$method_direction_comparison_num)

temp <- dplyr::group_by(plmr, r_ab, r_za, n, noisea, noiseb) %>%
	dplyr::summarise(prop_causality_exists_correct = mean(method_direction_comparison_num))
# temp$r_za <- paste0("cor(x,g) = ", temp$r_za)
# temp$n <- paste0("n = ", temp$n)


ggplot(subset(temp, round(r_za,2) == 0.01), aes(x = noisea, y = prop_causality_exists_correct)) +
geom_point(aes(colour = as.factor(noiseb))) +
facet_grid(n ~ r_ab) +
scale_colour_brewer(type="seq")



ggplot(subset(temp, r_ab == 0.4 & round(noiseb,2) == 0), aes(x=noisea, y = prop_causality_exists_correct)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za) +
labs(y="True positive rate", x = expression(cor(X, X[O])), fill="Test")

ggplot(subset(temp, r_ab == 0.4 & round(noisea,2) == 0), aes(x=noiseb, y = prop_causality_exists_correct)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(n ~ r_za) +
labs(y="True positive rate", x = expression(cor(X, X[O])), fill="Test")


ggplot(subset(temp, r_ab == 0.6 & n == "n = 1000"), aes(x=r_za, y = prop_causality_exists_correct)) +
geom_bar(stat="identity", position="dodge", aes(fill=as.factor(test))) +
facet_grid(noisea ~ noiseb) +
labs(y="True positive rate", x = expression(cor(X, X[O])), fill="Test")


## ---- shakhbazov ----


qtl_orig <- read.csv("../data/12864_2016_2498_MOESM16_ESM.csv")
qtl_orig <- qtl_orig[order(qtl_orig$PVALUE_expr), ]
qtl_orig <- subset(qtl_orig, !duplicated(paste0(meth_ind, exp_ind)))
cors <- read.csv("../data/12864_2016_2498_MOESM8_ESM.csv")

me <- expand.grid(rxx_o=c(0.5,0.75,1), ryy_o=c(0.5,0.75,1))
qtl <- group_by(me, rxx_o, ryy_o) %>%
	do({
		x <- .
		y <- qtl_orig
		y$rxx_o <- x$rxx_o[1]
		y$ryy_o <- x$ryy_o[1]
		return(y)
	})

qtl$r_exp <- NA
qtl$r_meth <- NA
qtl$dir <- NA
qtl$dir_p <- NA
qtl$mr_eff <- NA
qtl$mr_se <- NA
qtl$mr_p <- NA
for(i in 1:nrow(qtl))
{
	l <- mr_steiger(
		qtl$PVALUE_meth[i], 
		qtl$PVALUE_expr[i], 
		610, 862,
		qtl$rxx_o[i],
		qtl$ryy_o[i],
		screen = list(z = 70, x = -60, y = 3)
	)
	l$sensitivity_plot
	qtl$r_exp[i] <- l$r2_out_adj
	qtl$r_meth[i] <- l$r2_exp_adj
	qtl$dir[i] <- l$correct_causal_direction_adj
	qtl$dir_p[i] <- l$steiger_test_adj
	qtl$sensitivity[i] <- l$sensitivity

	if(qtl$AL1_meth[i] != qtl$AL1_expr[i])
	{
		qtl$EFFECT_meth[i] <- qtl$EFFECT_meth[i] * -1
	}

	if(qtl$dir[i])
	{
		a <- mr_wald_ratio(qtl$EFFECT_meth[i], qtl$EFFECT_expr[i], qtl$SE_meth[i], qtl$SE_expr[i])
		qtl$mr_eff[i] <- a$b
		qtl$mr_se[i] <- a$se
		qtl$mr_p[i] <- a$pval
	} else {
		a <- mr_wald_ratio(qtl$EFFECT_expr[i], qtl$EFFECT_meth[i], qtl$SE_expr[i], qtl$SE_meth[i])
		qtl$mr_eff[i] <- a$b
		qtl$mr_se[i] <- a$se
		qtl$mr_p[i] <- a$pval
	}
}


# Does expression or methylation most likely have the causal effect?

shakhbazov <- merge(qtl, subset(cors, select=c(meth_ind, exp_ind, pearson)), by=c("meth_ind", "exp_ind"), all.x=TRUE)
shakhbazov$mr_r <- sign(shakhbazov$mr_eff) * sqrt(abs(qnorm(shakhbazov$mr_p, low=FALSE))^2 / (abs(qnorm(shakhbazov$mr_p, low=FALSE))^2 + 400))
shakhbazov$dir <- as.character(shakhbazov$dir)
shakhbazov$dir[shakhbazov$dir=="TRUE"] <- "Methylation causes Expression"
shakhbazov$dir[shakhbazov$dir=="FALSE"] <- "Expression causes Methylation"
real_index <- shakhbazov$rxx_o == 1 & shakhbazov$ryy_o == 1


shakhsummary <- dplyr::group_by(shakhbazov, rxx_o, ryy_o, dir) %>%
	dplyr::summarise(
		n = n(),
		nsig = sum(dir_p < 0.05, na.rm=TRUE),
		pos = sum(mr_eff > 0) / n(),
		corr = cor(pearson^2, mr_r^2),
		sens = mean(sensitivity, na.rm=TRUE)
	) %>% as.data.frame()


## ---- shakhbazov_tests ----

# Is methylation's effect on expression more likely to be negative?

shakhtest1 <- fisher.test(
	table(sign(shakhbazov$mr_eff[real_index]), shakhbazov$dir[real_index])
)


# Is methylation's effect on expression more likely to be negative amongst dir < 0.05?

shakhtest2 <- fisher.test(
	table(
		sign(shakhbazov$mr_eff[real_index & shakhbazov$dir_p < 0.05]), 
		shakhbazov$dir[real_index & shakhbazov$dir_p < 0.05]
	)
)


# Is methylation more likely to cause expression amongst dor < 0.05

shakhtest3 <- binom.test(
	sum(shakhbazov$dir == "Methylation causes Expression" & shakhbazov$dir_p < 0.05 & real_index), 
	sum(shakhbazov$dir_p < 0.05 & real_index), 
	0.5
)


# Is methylation more likely to cause expression?

shakhtest4 <- binom.test(
	sum(shakhbazov$dir == "Methylation causes Expression" & real_index), 
	sum(real_index),
	0.5
)



## ---- shakhplot ----

# temp <- subset(shakhbazov, dir_p < 0.05 & real_index) %>%
# 	group_by(dir) %>%
# 	mutate(mr_eff = as.numeric(scale(mr_eff)), pos = ifelse(mr_eff >=0, "Positive effect", "Negative effect"))


# ggplot(temp, aes(x=mr_eff)) + 
# geom_density(aes(fill=dir), alpha=0.5, bw=0.04) + 
# # facet_grid(. ~ pos, scale="free_x") +
# labs(x="Causal effect (unit/unit change)", fill="Causal direction")

# temp2 <- table(temp$pos, temp$dir) %>% as.data.frame()
# ggplot(temp2, aes(y=Freq, x=Var2)) +
# geom_bar(aes(fill=Var1), position="stack", stat="identity")



temp <- gather(shakhsummary, key, value, n, nsig)
temp$rxx_o_lab <- paste0("cor(E,Eo) = ", temp$rxx_o)
temp$ryy_o_lab <- paste0("cor(M,Mo) = ", temp$ryy_o)

# ggplot(temp, aes(x=key, y=value)) +
# geom_bar(stat="Identity", aes(fill=dir), position="stack") +
# facet_grid(rxx_o_lab ~ ryy_o_lab) 


ggplot(subset(temp, key=="n"), aes(x=as.factor(rxx_o), y=value)) +
geom_bar(stat="Identity", aes(fill=dir), position="stack") +
geom_hline(yintercept = nrow(qtl_orig)/2, linetype="dotted") +
facet_grid(. ~ ryy_o_lab) +
labs(x="cor(E,Eo)", fill=NULL)



## ---- shakhtab ----

shakhtab <- shakhsummary
names(shakhtab) <- c("$\\rho_{x,x_o}$", "$\\rho_{y,y_o}$", "Causal direction", "Count", "$p_{Steiger} < 0.05$", "P(+ve effect)", "$cor(\\rho_{MR}, \\rho_{P})$", "Mean sensitivity")
kable(shakhtab)



## ---- shakhbazov_mr_vs_cor ----

cor(shakhbazov$mr_r, shakhbazov$pearson)
summary(lm(abs(mr_r) ~ abs(pearson), shakhbazov))
plot(abs(mr_r) ~ abs(pearson), shakhbazov)
abline(lm(abs(mr_r) ~ abs(pearson), shakhbazov))

ggplot(shakhbazov, aes(x=abs(pearson), y=abs(mr_r))) +
# geom_errorbar(aes(ymin = abs(mr_r) - mr_se * 1.96, ymax = abs(mr_r) + mr_se * 1.96), colour="grey") +
geom_point() +
stat_smooth(method="lm") +
facet_grid(. ~ dir) +
labs(x=TeX("$|\\rho_{P}|$"), y=TeX("$|\\rho_{MR}|$"))





## ---- steiger_sensitivity ----

sensitivity_parameters <- expand.grid(
	r_xy = seq(-1,1,0.05),
	r_gx = sqrt(c(0.01,0.05,0.1)),
	n=1000,
	sensitivity = NA
)

sensitivity_parameters$r_gy <- sensitivity_parameters$r_gx^2 * sensitivity_parameters$r_xy^2

sensitivity_parameters$r_gy <- sensitivity_parameters$r_xy * sensitivity_parameters$r_gx
sensitivity_parameters$p_gy <- get_p_from_r2n(sensitivity_parameters$r_gy^2, sensitivity_parameters$n)
sensitivity_parameters$p_gx <- get_p_from_r2n(sensitivity_parameters$r_gx^2, sensitivity_parameters$n)

for(i in 1:nrow(sensitivity_parameters))
{
	sensitivity_parameters$sensitivity[i] <- mr_steiger(sensitivity_parameters$p_gx[i], sensitivity_parameters$p_gy[i], sensitivity_parameters$n[i], sensitivity_parameters$n[i])$sensitivity_ratio
}


## ---- steiger_sensitivity_plot ----

example_sensitivity <- mr_steiger(
	get_p_from_r2n(sqrt(0.01), 10001),
	get_p_from_r2n(sqrt(0.01 * 0.1), 10001),
	10001, 10001,
	screen=list(z = 55, x = -60, y = 3), bty="n"
)

example_sensitivity

p1 <- example_sensitivity$sensitivity_plot

p2 <- ggplot(sensitivity_parameters, aes(y=sensitivity, x=r_xy, group=factor(r_gx^2))) +
	geom_point(aes(colour=factor(r_gx^2))) +
	geom_line(aes(colour=factor(r_gx^2))) +
	labs(x=expression(rho[xy]), y=expression(V[z<0]/V[z>=0]), colour=expression(rho[gx]^2)) +
	scale_colour_brewer(type="qual")
p2 + scale_y_log10()
grid.arrange(
	textGrob("a)", x=unit(0.1, "npc")), textGrob("b)", x=unit(0.1, "npc")),
	p1, p2, 
	ncol=2, heights=c(1,10)
)

