library(cit)
library(reshape2)
library(ggplot2)
library(plyr)

load("~/repo/MethylationIV/cit/results/20141114.RData")


datmrp <- ddply(datmr, .(noise, bias, n), function(x){
	x <- mutate(x)
	a <- sum(x$bdmr == "AtoB") / nrow(x)
	ap <- sum(x$bdmr_p == "AtoB") / nrow(x)
	b <- sum(x$bdmr == "BtoA") / nrow(x)
	bp <- sum(x$bdmr_p == "BtoA") / nrow(x)
	con <- sum(x$bdmr == "Confounder") / nrow(x)
	conp <- sum(x$bdmr_p == "Confounder") / nrow(x)
	u <- sum(x$bdmr == "Unknown") / nrow(x)
	up <- sum(x$bdmr_p == "Unknown") / nrow(x)
	rsq <- mean(x$rsq)
	y <- expand.grid(cause = c("Correct inference", "Incorrect inference", "Unknown"), variable = c("True", "Proxy"))
	y$rsq <- rsq
	y$prob <- c(a,b,u,ap,bp,up)
	return(y)
})


dat2 <- ddply(dat, .(n), function(x)
{
	reshape(mutate(x), varying=list(c("GT", "TG")), direction="long", times=c("Correct inference", "Incorrect inference"))
})

dat2$tim <- dat2$time
dat3 <- dat2
dat3$GT[dat3$GT < 1e-100] <- 1e-100

ggplot(subset(dat3, sqrt(rsq) > 0.75 & tim=="Incorrect inference"), aes(x=sqrt(rsq),y=-log10(GT))) +
stat_smooth(aes(colour=factor(n)), se=FALSE) +
facet_grid(tim ~ .) +
labs(y=expression(-log[10]*p), x=expression(cor(X, X[p])), colour="")
ggsave(file="~/repo/MethylationIV/cit/images/samplesize_zoomed.pdf")


ggplot(subset(dat2), aes(x=sqrt(rsq),y=-log10(GT))) +
geom_point(aes(colour=time)) +
stat_smooth(aes(colour=time), se=FALSE) +
labs(y=expression(-log[10]*p), x=expression(cor(X, X[p])), colour="") +
facet_grid(. ~ n)
ggsave("~/repo/MethylationIV/cit/images/proxy_cit.pdf")

ggplot(subset(dat2), aes(x=sqrt(rsq),y=-log10(GT))) +
geom_point(aes(colour=time)) +
labs(y=expression(-log[10]*p), x=expression(cor(X, X[p])), colour="") +
facet_grid(n ~ ., scale="free_y")
ggsave("~/repo/MethylationIV/cit/images/proxy_cit_n.pdf", width=7, height=10)


levels(datmrp$cause) <- c("Correct inference", "Incorrect inference", "Inconclusive")
ggplot(subset(datmrp, variable=="Proxy"), aes(y=prob, x=sqrt(rsq))) +
geom_line(aes(colour=cause)) +
geom_point(aes(colour=cause)) +
scale_colour_brewer(type="qual") +
facet_grid(. ~ n) +
labs(y="Frequency", x=expression(cor(X, X[p])), colour="")
ggsave("~/repo/MethylationIV/cit/images/proxy_bdmr.pdf")


