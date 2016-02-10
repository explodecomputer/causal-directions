library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)


load("~/repo/cit_measurement_error/results/p_comp_20160208.RData")

parameters$p_test <- parameters$p_bz
parameters$p_test[!parameters$correct_direction] <- parameters$p_az[!parameters$correct_direction]
parameters$r_ab <- parameters$r_ab^2
parameters$r_za <- parameters$r_za^2
parameters$noisea <- parameters$noisea^2
parameters$noiseb <- parameters$noiseb^2

parameters$rhs <- sqrt(parameters$r_ab) * parameters$cor_bbp

table(parameters$rhs > parameters$cor_aap)
table(parameters$correct_direction & parameters$sig_mr, parameters$rhs < parameters$cor_aap & parameters$sig_mr)


psum1 <- parameters %>%
	group_by(n, p, r_ab, r_za, noisea, noiseb) %>%
	summarise(
		correct_direction=sum(correct_direction)/n(),
		correct_direction_sig=sum(correct_direction & p_test > -log10(0.05))/n(),
		sig_mr=sum(sig_mr)/n(),
		cor_aap=mean(cor_aap),
		cor_bbp=mean(cor_bbp),
		cor_abp=mean(cor_abp),
		sig_mr=sum(p_test > -log10(0.05))/n(),
		bigna=first(bigna)
	)

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



