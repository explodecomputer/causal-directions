library(cit)
library(tidyverse)

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}


# n <- 10000
# u <- rnorm(n)
# g <- rbinom(n, 2, 0.5)
# x <- g + u + rnorm(n)
# y <- x + u + rnorm(n)

# i <- 18
# u <- rnorm(param$n[i])
# g <- rbinom(param$n[i], 2, 0.5)
# x <- scale(g) * sqrt(param$vg[i]) + u * sqrt(param$vu[i]) + rnorm(param$n[i], sd = sqrt(1 - param$vg[i] - param$vu[i]))
# y <- x * sqrt(param$vx[i]) + u * sqrt(param$vu[i]) + rnorm(param$n[i], sd = sqrt(1 - param$vx[i] - param$vu[i]))

# cor(g,x)^2
# cor(x,y)^2
# cor(x,u)^2
# cor(u,y)^2


# cit.cp(G = x, L = g, T = y)


# n <- 10000
# u <- rnorm(n)
# g <- rbinom(n, 2, 0.5)
# x <- g + u + rnorm(n)
# y <- x + u*0 + rnorm(n)


# cit.cp(G = x, L = g, T = y)

param <- expand.grid(
	n=c(1000,3000,6000),
	vu=runif(300, 0, 0.5),
	vg=c(0.2),
	vx=c(0.2)
)

for(i in 1:nrow(param))
{
	message(i, " of ", nrow(param))
	u <- rnorm(param$n[i])
	g <- rbinom(param$n[i], 2, 0.5)
	x <- scale(g) * sqrt(param$vg[i]) + u * sqrt(param$vu[i]) + rnorm(param$n[i], sd = sqrt(1 - param$vg[i] - param$vu[i]))
	y <- x * sqrt(param$vx[i]) + u * sqrt(param$vu[i]) + rnorm(param$n[i], sd = sqrt(1 - param$vx[i] - param$vu[i]))
	param$cit[i] <- cit.cp(G = x, L = g, T = y)[1]
	param$mr[i] <- fastAssoc(y, g)$pval
}

save(param, file="../results/cit_collider.rdata")

###

library(tidyverse)
load("../results/cit_collider.rdata")

res <- gather(param, key=test, value=pval, cit, mr)
res$test <- as.factor(res$test)
levels(res$test) <- c("CIT", "y ~ g")

ggplot(subset(res), aes(x=vu, y=-log10(pval))) +
geom_point(aes(colour=test)) +
geom_smooth(aes(colour=test)) +
scale_colour_brewer(type="qual", palette=2) +
facet_grid(n ~ .) +
labs(x="Confounder effect", colour="Test")

ggsave("../images/cit_confounder.pdf")



#####

library(cit)
library(dplyr)
library(ggplot2)

nid <- 5000

eff <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5)
eff <- expand.grid(conf=c(eff, -eff, 0), geff=c(eff, -eff, 0))
dim(eff)
l <- list()
for(i in 1:nrow(eff))
{
	message(i)
	g <- rnorm(nid)
	u <- rnorm(nid)
	a <- g * eff$geff[i] + u * eff$conf[i] + rnorm(nid)
	b <- g * eff$geff[i] + u * eff$conf[i] + rnorm(nid)
	x <- as.data.frame(rbind(
		cit.cp(L=g, G=a, T=b),
		cit.cp(L=g, G=b, T=a)
	))
	x$sim <- i
	x$dir <- c("AB", "BA")
	l[[i]] <- x
}



l1 <- bind_rows(l)
l1$conf <- rep(eff$conf, each=2)
l1$geff <- rep(eff$geff, each=2)

l2 <- group_by(l1, sim) %>%
	summarise(nsig = sum(p_cit < 0.05), conf=conf[1], geff=geff[1])

save(l1, l2, file="../results/cit_confounder.rdata")


load("../results/cit_confounder.rdata")
ggplot(l1, aes(x=conf, y=geff)) +
geom_point(aes(colour=p_cit < 0.05)) +
facet_grid(. ~ dir)


l2$nsig <- as.factor(l2$nsig)
levels(l2$nsig) <- c("None", "One direction", "Both directions")

ggplot(l2, aes(x=conf, y=geff)) +
geom_point(aes(colour=as.factor(nsig)), size=4) +
scale_colour_brewer(type="qual") +
labs(x="Confounder effect", y="Genetic effect", colour="CIT significance")
ggsave("../images/cit_confounder2.pdf")
