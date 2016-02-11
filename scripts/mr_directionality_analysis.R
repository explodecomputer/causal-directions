library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)

get_cor_from_pval <- function(p, n)
{
	Fval <- qf(p, 1, n-1, low=FALSE)
	R2 <- Fval / (n - 2 + Fval)
	return(R2)
}

load("~/repo/cit_measurement_error/results/p_comp_20160210.RData")

parameters$p_az <- -log10(parameters$p_az)
parameters$p_bz <- -log10(parameters$p_bz)
parameters$sig_mr <- parameters$p_bz > -log10(0.05)
parameters$correct_direction <- parameters$p_az > parameters$p_bz
parameters$bigna <- parameters$noisea > parameters$noiseb

parameters$p_test <- parameters$p_bz
parameters$p_test[!parameters$correct_direction] <- parameters$p_az[!parameters$correct_direction]
parameters$r_ab <- parameters$r_ab^2
parameters$r_za <- parameters$r_za^2
parameters$noisea <- parameters$noisea^2
parameters$noiseb <- parameters$noiseb^2

parameters$rhs <- sqrt(parameters$r_ab) * parameters$cor_bbp
parameters$lhs <- parameters$cor_aap
parameters$rhs_lhs_diff <- parameters$lhs - parameters$rhs

a <- r.test(n=parameters$n, r12=parameters$cor_zap, r13=parameters$cor_zbp, r23=parameters$cor_abp)
parameters$cortest_t <- a$t
parameters$cortest_p <- -log10(a$p)
parameters$cortest_sign_correct <- sign(parameters$cortest_t) == 1

table(parameters$rhs > parameters$cor_aap)

parameters <- subset(parameters, ! r_ab %in% c(0, 1))

parameters$cit_direction <- parameters$cit_AB < parameters$cit_BA
table(parameters$cit_direction, sign(parameters$cortest_t))

parameters$cit_p <- parameters$cit_AB
parameters$cit_p[!parameters$cit_direction] <- parameters$cit_BA[!parameters$cit_direction]
parameters$cit_p <- -log10(parameters$cit_p)


ps <- parameters %>%
	gather()



psum1 <- parameters %>%
	group_by(n, p, r_ab, r_za, noisea, noiseb) %>%
	summarise(
		correct_direction=sum(correct_direction)/n(),
		correct_direction_sig=sum(correct_direction & p_test > -log10(0.05))/n(),
		cor_aap=mean(cor_aap),
		cor_bbp=mean(cor_bbp),
		cor_abp=mean(cor_abp),
		cor_zap=mean(cor_zap),
		cor_zbp=mean(cor_zbp),
		nsim=n(),
		prop_cortest_sign = sum(sign(cortest_t) == 1, na.rm=T)/n(),
		prop_cortest_sig = sum(cortest_p > -log10(0.05),na.rm=T)/n(),
		cortest_t = mean(cortest_t, na.rm=T),
		sig_mr=sum(p_test > -log10(0.05))/n(),
		bigna=first(bigna),
		rhs=mean(rhs),
		lhs=mean(cor_aap),
		rhs_lhs_diff=lhs-rhs
	)


ggplot(subset(psum1, r_za==0.1 & n == 1000), aes(y=prop_cortest_sign, x=prop_cortest_sig)) +
geom_point(aes(colour=as.factor(r_ab))) +
facet_grid(noisea ~ noiseb)



ggplot(subset(parameters, sig_mr), aes(x=rhs_lhs_diff, y=cortest_p)) +
geom_point(aes(colour=cortest_sign_correct)) +
facet_grid(n ~ r_za, scale="free_y")



ggplot(subset(parameters, sig_mr), aes(x=rhs_lhs_diff, y=-log10(cit_p))) +
geom_point(aes(colour=cit_direction)) +
facet_grid(n ~ r_za, scale="free_y")


a <- subset(parameters, select=c(rhs_lhs_diff, cit_direction, cit_p))
b <- subset(parameters, select=c(rhs_lhs_diff, cortest_sign_correct, cortest_p))
a$test <- "CIT"
b$test <- "MR"
names(a) <- names(b) <- c("rhs_lhs_diff", "direction", "pval", "test")
ab <- rbind(a,b)

ggplot(ab, aes(x=rhs_lhs_diff, y=pval)) +
geom_point(aes(colour=direction)) +
facet_grid(test ~ ., scale="free_y")



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


compare cor_zap with cor_zbp
see how often you get a significant result for the wrong direction?




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


library(lattice)
temp <- subset(psum1, n==1000 & r_ab != 0)
temp$r_za <- paste0("cor(A,Z)=",temp$r_za^2)
temp$r_ab <- paste0("cor(A,B)=",temp$r_ab^2)
temp$noisea <- temp$noisea^2
temp$noiseb <- temp$noiseb^2
wireframe(correct_direction_sig ~ noisea*noiseb | as.factor(r_ab) + as.factor(r_za), data=temp, drape=TRUE)





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




ineq <- expand.grid(
	noisea = seq(0, 1, by=0.02),
	noiseb = seq(0, 1, by=0.02)
)

ineq$ab <- ineq$noisea / ineq$noiseb
ineq$ab[ineq$ab > 1] <- NA

summary(ineq$ab)

library(lattice)
wireframe(ab ~ noisea * noiseb, ineq, drape=TRUE, screen = list(x=-60,y=180), scales = list(arrows = FALSE))


