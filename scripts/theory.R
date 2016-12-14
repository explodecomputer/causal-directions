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
	e_m <- runif(1)*10

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

lm(obs ~ exp, dat)
lm(obs ~ exp, dat0)


##

n <- 10000
x <- rnorm(n)

bm <- 0.5
am <- 5
em <- rnorm(n) / 5
em <- rep(0, n)

bx <- 2
ax <- 3
ex <- rnorm(n)

bn <- 50
an <- 10
en <- rnorm(n)/10

x0 <- am + bm * x + em
y <- ax + bx * x + ex
y0 <- an + bn * y + en

py <- fitted.values(lm(y ~ x))
py0 <- fitted.values(lm(y0 ~ x))
py00 <- fitted.values(lm(y0 ~ x0))

plot(py ~ py0)
plot(py00 ~ py0)

summary(lm(py0 ~ py00))
summary(lm(py0 ~ py))

cov(py0, py00)
bx^2 * bn^2 * bm * var(x)

lm(y ~ x)
lm(y0 ~ x)

var(py0)
var(py00)

pred <- 10 + 0.1 * (3 + 2 * x)

plot(pred ~ py0)

x <- rnorm(n)
e <- rnorm(n)/2
b <- 10
x0 <- 3 + b * x + e

cov(x, x0)
b * var(x)



n <- 1000
u <- rnorm(n)
a <- u + rnorm(n)
b <- u + rnorm(n)
summary(lm(a ~ b + u))



# two variables, instrumented...

# if x causes y then use steiger 

we have r_iv^2 and steiger test statistic

steiger requires rgx, rgy, ngx, ngy


n <- 30000
g <- rbinom(n, 2, 0.5)
x <- rnorm(n) + g
y <- 2 * x + rnorm(n)

yo <- y + rnorm(n) * 5
xo <- x + rnorm(n)

bgx <- cov(g,x) / var(g)
bgy <- cov(g,y) / var(g)
(biv <- bgy / bgx)

bgxo <- cov(g,xo) / var(g)
bgyo <- cov(g,yo) / var(g)
(bivo <- bgyo / bgxo)

cor(g,x)^2
cor(g,xo)^2
cor(x,xo)^2 * cor(g,x)^2


get_r_from_pn <- function(p, n)
{
	Fval <- qf(p, 1, n-1, low=FALSE)
	R2 <- Fval / (n - 2 + Fval)
	return(R2)
}


mr_steiger <- function(p_exp, p_out, n_exp, n_out) 
{
	require(psych)
	index <- any(is.na(p_exp)) | any(is.na(p_out)) | any(is.na(n_exp)) | any(is.na(n_out))
	p_exp <- p_exp[!index]
	p_out <- p_out[!index]
	n_exp <- n_exp[!index]
	n_out <- n_out[!index]
	r_exp <- get_r_from_pn(p_exp, n_exp)
	r_out <- get_r_from_pn(p_out, n_out)
	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	l <- list(r2_exp = r_exp^2, r2_out = r_out^2, correct_causal_direction = r_exp > r_out, steiger_test = rtest$p)
	return(l)
}

get_p_from_rn(0.46,890)


optim()


get_p_from_rn <- function(r, n)
{
	fval <- r * (n-2) / (1 - r)
	pval <- pf(fval, 1, n-1, low=FALSE)
	return(pval)
}


get_p_from_rn(get_r_from_pn(0.0001, 100), 100)





rgx_o <- rgx * rxx_o
rgx_o <- rgx * rxx_o


real rgx can range from rgx to 1

rgx_o <- 0.3

rxx_o <- seq(rgx_o,1,length.out=100)
rgx <- rgx_o / rxx_o
# rxx_o[rxx_o > 1] <- NA

plot(rgx ~ rxx_o)



b * ((a-1) * log(b*x/a))


rgy_o <- 0.04

ryy_o <- seq(rgy_o,1,length.out=100)
rgy <- rgy_o / ryy_o
# ryy_o[ryy_o > 1] <- NA

plot(rgy ~ ryy_o)
points(rgx ~ rxx_o)


x = measurement error
y = true r2
c = measured r2

c = xy


library(lattice)


rgx_o <- 0.7
rgy_o <- 0.2

d <- expand.grid(rxx_o=seq(rgx_o,1,length.out=70), ryy_o=seq(rgy_o,1,length.out=70), type=c("A","B"))
d$rgy <- rgy_o / d$ryy_o
d$rgx <- rgx_o / d$rxx_o
d$z <- d$rgy - d$rgx
d$z[d$type=="A"] <- 0
d$col <- d$z > 0
temp <- wireframe(z ~ ryy_o * rxx_o, groups=type, data=d, scales=list(arrows=FALSE), col.groups = colorRampPalette(c("red", "blue"))(2), drape=FALSE, pretty=TRUE)

0.6*0.8

temp <- subset(d, round(z, 2) == 0 & type=="B")
head(temp)
dim(temp)
dim(d)
plot(ryy_o ~ rxx_o, data=temp)

nsim <- 100000
ryy_o <- runif(nsim, min=rgy_o, max=1)
rxx_o <- runif(nsim, min=rgx_o, max=1)

rgy <- rgy_o / ryy_o
rgx <- rgx_o / rxx_o
table(rgy > rgx)
d <- data.frame(z=rgy - rgx, x=rxx_o, y=ryy_o)
d$col <- d$z > 0
cloud(z ~ y * x, data=d, scales=list(arrows=FALSE, col="black"), groups=col, col.groups=c("red", "blue"))



rgy_o * (2 * rgx_o - log(rgy_o) - 2) + rgy_o * log(rgy_o / rgx_o)

rgy_o * (2 * rgx_o - log())

rgx_o * rgy_o * log(rgy_o/rgx_o) - rgy_o * log(rgy_o)



rgy_o * log(rgy_o) + rgx_o*(rgy_o*(log(rgx_o) - log(rgy_o)) - log(rgx_o))
rgx_o * log(rgx_o) + rgy_o*(rgx_o*(log(rgy_o) - log(rgx_o)) - log(rgy_o))


cloud(Sepal.Length ~ Petal.Length * Petal.Width, 
                 data = iris, cex = .8, 
                 groups = Species, 
                 main = "Stereo",
                 screen = list(z = 20, x = -70, y = 3),
                 # par.settings = par.set,
                 scales = list(col = "black"))



d <- expand.grid(x=seq(-1,))

a <- 0.2
x <- seq(a,1,length=100)
b <- 0.01
y <- seq(b,1,length=100)

act <- sqrt(b^2 / y^4 + a^2 / x^4 + 1)
bin <- a/x^2 + 0.5*x^2 / a * (b^2/y^4 + 1)

plot(act ~ bin)



get_r_from_pn <- function(p, n)
{
	Fval <- qf(p, 1, n-1, low=FALSE)
	# print(Fval)
	if(!is.finite(Fval))
	{
		get_p_from_rn <- function(r, n)
		{
			fval <- r * (n-2) / (1 - r)
			pval <- pf(fval, 1, n-1, low=FALSE)
			return(pval)
		}
		optim.get_p_from_rn <- function(x, sample_size, pvalue)
		{
			abs(-log10(get_p_from_rn(x, sample_size)) - -log10(pvalue))
		}
		R2 <- optim(0.1, optim.get_p_from_rn, sample_size=n, pvalue=p)$par
		# print(R2)
		return(R2)
	}
	R2 <- Fval / (n - 2 + Fval)
	return(R2)
}

steiger_sensitivity <- function(rgx_o, rgy_o)
{
	if(rgy_o > rgx_o)
	{
		a <- rgy_o
		b <- rgx_o
	} else {
		a <- rgx_o
		b <- rgy_o
	}

	d <- expand.grid(rxx_o=seq(rgx_o,1,length.out=70), ryy_o=seq(rgy_o,1,length.out=70), type=c("A","B"))
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
		xlab=expression(rho[xx[o]]), 
		ylab=expression(rho[yy[o]]),
		zlab=expression(rho[gy]-rho[gx])
	)


	vz <- a * log(a) - b * log(b) + a*b*(log(b)-log(a))
	vz0 <- -2*b - b * log(a) - a*b*log(a) + 2*a*b

	sensitivity <- vz0 / (2 * vz0 + abs(vz))

	return(list(
		vz = vz,
		vz0 = vz0,
		sensitivity = sensitivity,
		pl = temp
	))
}


mr_steiger <- function(p_exp, p_out, n_exp, n_out) 
{
	require(psych)
	index <- any(is.na(p_exp)) | any(is.na(p_out)) | any(is.na(n_exp)) | any(is.na(n_out))
	p_exp <- p_exp[!index]
	p_out <- p_out[!index]
	n_exp <- n_exp[!index]
	n_out <- n_out[!index]
	r_exp <- get_r_from_pn(p_exp, n_exp)
	r_out <- get_r_from_pn(p_out, n_out)

	sensitivity <- steiger_sensitivity(r_exp, r_out)

	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	l <- list(
		r2_exp = r_exp^2, 
		r2_out = r_out^2, 
		correct_causal_direction = r_exp > r_out, 
		steiger_test = rtest$p,
		vz = sensitivity$vz,
		vz0 = sensitivity$vz0,
		sensitivity = sensitivity$sensitivity,
		sensitivity_plot = sensitivity$pl
	)
	return(l)
}


mr_steiger(1e-120, 1e-30, 900, 300)

temp <- steiger_sensitivity(0.3, 0.05)
pdf("")