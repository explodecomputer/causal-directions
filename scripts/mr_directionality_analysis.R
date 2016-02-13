## ---- load_up ----

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(lattice))

get_cor_from_pval <- function(p, n)
{
	Fval <- qf(p, 1, n-1, low=FALSE)
	R2 <- Fval / (n - 2 + Fval)
	return(R2)
}


## ---- read_data ----

load("~/repo/cit_measurement_error/results/mr_directionality_20160211.RData")

parameters$p_az <- -log10(parameters$p_az)
parameters$p_bz <- -log10(parameters$p_bz)
parameters$r_ab <- parameters$r_ab^2
parameters$r_ab <- round(parameters$r_ab, 1)
parameters$r_za <- parameters$r_za^2
parameters$noisea <- parameters$noisea^2
parameters$noiseb <- parameters$noiseb^2
parameters$rhs <- parameters$cor_ab * parameters$cor_bbp
parameters$lhs <- parameters$cor_aap
parameters$rhs_lhs_diff <- parameters$lhs - parameters$rhs

parameters$correct_direction <- parameters$p_az > parameters$p_bz

parameters$p_test <- parameters$p_bz
parameters$p_test[!parameters$correct_direction] <- parameters$p_az[!parameters$correct_direction]

# Calculate correlation test
a <- r.test(n=parameters$n, r12=parameters$cor_zap, r13=parameters$cor_zbp, r23=parameters$cor_abp)
parameters$cortest_t <- a$t
parameters$cortest_p <- -log10(a$p)
parameters$cortest_correct_direction <- sign(parameters$cortest_t) == 1


parameters$cit_correct_direction <- parameters$cit_AB < parameters$cit_BA
parameters$cit_correct_direction[is.na(parameters$cit_correct_direction)] <- FALSE

parameters$cit_p <- parameters$cit_AB
parameters$cit_p[!parameters$cit_correct_direction] <- parameters$cit_BA[!parameters$cit_correct_direction]
parameters$cit_p <- -log10(parameters$cit_p)
parameters$cit_p[is.na(parameters$cit_AB) | is.na(parameters$cit_BA)] <- NA

pl <- gather(parameters, test, direction_p_value, cit_p, cortest_p)
pl$correct_direction <- pl$cortest_correct_direction
pl$correct_direction[pl$test == "cit_p"] <- pl$cit_correct_direction[pl$test == "cit_p"]
pl$test_p_value <- pl$p_test
pl$test_p_value[pl$test == "cit_p"] <- pl$direction_p_value[pl$test == "cit_p"]
pl <- subset(pl, select=-c(cortest_correct_direction, cit_correct_direction, p_test))
pl$test[pl$test == "cit_p"] <- "CIT"
pl$test[pl$test == "cortest_p"] <- "MR"

pl$test_p_value[is.infinite(pl$test_p_value)] <- max(is.finite(pl$test_p_value))
pl$direction_p_value[is.infinite(pl$direction_p_value)] <- max(is.finite(pl$direction_p_value))


psum1 <- pl %>%
	group_by(n, p, r_ab, r_za, noisea, noiseb, test) %>%
	summarise(
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
		prop_sig_incorrect=sum(!correct_direction & direction_p_value > -log10(0.05), na.rm=T)/n(),
		prop_sig_correct=sum(correct_direction & direction_p_value > -log10(0.05), na.rm=T)/n(),
		prop_nonsig=sum(direction_p_value <= -log10(0.05), na.rm=T)/n()
	)
psum1$rhs_lhs_diff_bin <- cut(psum1$rhs_lhs_diff, breaks=10)


## ---- cit_measurement_error_figure ----

load("~/repo/cit_measurement_error/results/20141114.RData")
dat <- gather(dat, eval, pval, GT, TG, factor_key=TRUE)
levels(dat$eval) <- c("Correct inference", "Incorrect inference")

ggplot(dat, aes(x=sqrt(rsq),y=-log10(pval))) +
geom_point(aes(colour=eval)) +
labs(y=expression(-log[10]*p), x=expression(cor(x, x[O])), colour="") +
facet_grid(n ~ ., scale="free_y")



## ---- cit_mr_comparison_figure ----

psum2 <- gather(psum1, eval, value, prop_sig_correct, prop_sig_incorrect, prop_nonsig, factor_key=TRUE) %>%
	group_by(n, rhs_lhs_diff_bin, r_za, eval, test) %>%
	summarise(value=mean(value, na.rm=TRUE))	

temp <- do.call(rbind, strsplit(as.character(psum2$rhs_lhs_diff_bin), split=","))
psum2$rhs_lhs_diff_bin_numeric <- as.factor(gsub("\\(", "", temp[,1]))
psum2$rhs_lhs_diff_bin_lab <- as.factor(paste0("d = ", psum2$rhs_lhs_diff_bin_numeric))
psum2$rhs_lhs_diff_bin_lab <- factor(psum2$rhs_lhs_diff_bin_lab, levels=levels(psum2$rhs_lhs_diff_bin_lab)[order(as.numeric(as.character(levels(psum2$rhs_lhs_diff_bin_numeric))))])

levels(psum2$eval) <- c("Evidence for causality with correct direction", "Evidence for causality with incorrect direction", "No evidence for causality")

ggplot(subset(psum2, r_za==0.1), aes(x=test, y=value)) +
geom_bar(stat="identity", aes(fill=eval)) +
facet_grid(n ~ rhs_lhs_diff_bin_lab) +
scale_fill_brewer(type="qual") +
labs(x="Method", y="Proportion of simulations")



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


## ---- plot3 ----

ineq <- expand.grid(
	noisea = seq(0, 1, by=0.02),
	noiseb = seq(0, 1, by=0.02),
	ab = seq(0, 1, by=0.2)
)

ineq$lhs <- ineq$noisea
ineq$rhs <- ineq$noiseb * ineq$ab
ineq$d <- ineq$lhs - ineq$rhs

ggplot(ineq, aes(x=noisea, y=noiseb)) +
geom_tile(aes(fill=d)) +
facet_wrap(~ ab)


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




