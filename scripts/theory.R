nsim <- 100

c_gy <- array(0, nsim)
c_gyhat <- array(0, nsim)
c_gy_exp <- array(0, nsim)
c_gyhat_exp <- array(0, nsim)

c_gy0 <- array(0, nsim)
c_gy0hat <- array(0, nsim)
c_gy0_exp <- array(0, nsim)
c_gy0hat_exp <- array(0, nsim)


n <- 100000

for(i in 1:nsim)
{
	message(i)
	g <- rbinom(n, 2, 0.5)

	a_g <- runif(1)
	b_g <- runif(1)
	e_g <- runif(1)

	a_m <- runif(1)
	b_m <- runif(1)
	e_m <- runif(1)

	a_x <- runif(1)
	b_x <- runif(1)
	e_x <- runif(1)

	a_n <- runif(1)
	b_n <- runif(1)
	e_n <- runif(1)

	x <- a_g + b_g * g + rnorm(n, sd=sqrt(e_g))
	x0 <- a_m + b_m * x + rnorm(n, sd=sqrt(e_m))
	y <- a_x + b_x * x + rnorm(n, sd=sqrt(e_x))
	y0 <- a_n + b_n * y + rnorm(n, sd=sqrt(e_n))
	yhat <- fitted.values(lm(y ~ x0))

	c_gy[i] <- cov(g, y)
	c_gy_exp[i] <- b_x * b_g * var(g)

	c_gyhat[i] <- cov(g, yhat)
	c_gyhat_exp[i] <- b_x * b_g * var(g) * b_m^2 * var(x) / (b_m^2 * var(x) + e_m)


	y0hat <- fitted.values(lm(y0 ~ x0))

	c_gy0[i] <- cov(g, y0)
	c_gy0_exp[i] <- b_n * b_x * b_g * var(g)

	c_gy0hat[i] <- cov(g, y0hat)
	c_gy0hat_exp[i] <- b_n * b_x * b_g * var(g) * b_m^2 * var(x) / (b_m^2 * var(x) + e_m)

}


dat <- rbind(
	data.frame(obs=c_gy, exp=c_gy_exp, what="cov_gy"),
	data.frame(obs=c_gyhat, exp=c_gyhat_exp, what="cov_gyhat")
)

dat0 <- rbind(
	data.frame(obs=c_gy0, exp=c_gy0_exp, what="cov_gy"),
	data.frame(obs=c_gy0hat, exp=c_gy0hat_exp, what="cov_gyhat")
)

library(ggplot2)

ggplot(dat, aes(x=exp, y=obs)) +
geom_point() +
facet_grid(. ~ what)

ggplot(dat0, aes(x=exp, y=obs)) +
geom_point() +
facet_grid(. ~ what)

