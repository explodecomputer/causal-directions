n <- 100000
vu <- 1
vg <- 1
bxy <- 1
bux <- 1
buy <- -5
bgx <- 05
vey <- 1
vex <- 1
u <- rnorm(n, sd=sqrt(vu))
g <- rnorm(n, sd=sqrt(vg))
x <- bgx * g + bux * u + rnorm(n, sqrt(vex))
y <- bxy * x + buy * u + rnorm(n, sqrt(vey))


cor(g,x)^2
bgx^2 / var(x)


cor(u,x)^2
bux^2/var(x)

# cor(x,y)^2

cor(g,y)^2


cor(g,x)^2 * cor(y,x)^2
cor(g,x)^2 * bxy^2 / var(y)

cov(g, y)

cov(g, bxy * (bgx * g + bux * u) + buy * u)
cov(g, bxy * bgx * g) / var(g)
cov(g,y) / var(g)

(cov(x,y)/var(x))^2 * bgx^2 * var(g) + (cov(x,y)/var(x))^2 * bux^2 * var(u) + buy^2 * var(u) + 1
var(y)

(bxy^2 * bgx^2 * var(g)^2) / ((cov(x,y)/var(x))^2 * bgx^2 * var(g) + (cov(x,y)/var(x))^2 * bux^2 * var(u) + buy^2 * var(u) + 1)

cor(g,y)^2
cor(g,x)^2

cor(x,y)^2

lm(y ~ x)

var(x)
var(y)

lm(y ~ x)$coefficients[2]^2 * var(x)/var(y)



get_bxy <- function(bxy, bgx, bux, buy, vg, vu, vex, vey)

get_calcs <- function(bxy, bgx, bux, buy, vg, vu, vex, vey)
{
	bxyo <- ((bgx^2*bxy*vg + bux^2*bxy*vu + bux*buy*vu + bxy) / (vg*bgx^2 + vu*bux^2 + vex))

	vx <- bgx^2 * vg + bux^2 * vu + vex
	vy <- bxyo^2 * vx + buy^2 * vu + vey

	rsqxyo <- bxyo^2 * vx / vy
	rsqxyos <- rsqxyo * sign(bxyo)
	rsqxy <- bxy^2 * vx / vy
	rsqxys <- rsqxy * sign(bxy)
	rsqgx <- bgx^2 / (vg * bgx^2 + vu * bux^2 + vex)
	rsqgy <- (bxy^2 * bgx^2 * vg^2) / (bxyo^2 * bgx^2 * vg + bxyo^2 * bux^2 * vu + buy^2 * vu + vey)

	rsqux <- bux^2 / (vg*bgx^2 + vu*bux^2 + vex)
	rsquy <- buy^2 / (bxyo^2 * bgx^2 * vg + bxyo^2 * bux^2 * vu + buy^2 * vu + vey)
	rsquxs <- rsqux * sign(bux)
	rsquys <- rsquy * sign(buy)
	return(data.frame(bxyo=bxyo, rsqgx=rsqgx, rsqgy=rsqgy, rsqux=rsqux, rsquy=rsquy, rsqxy=rsqxy, rsqxyo=rsqxyo, rsqxyos=rsqxyos, rsquxs=rsquxs, rsquys=rsquys, rsqxys=rsqxys))
}

library(tidyverse)

dat <- expand.grid(
	bxy=seq(-2, 2, length.out=11),
	bux=seq(-2, 2, length.out=110),
	buy=seq(-2, 2, length.out=110),
	bgx=sqrt(0.1),
	vu=1,
	vg=1,
	vex=1,
	vey=1
)

res <- get_calcs(dat$bxy, dat$bgx, dat$bux, dat$buy, dat$vg, dat$vu, dat$vex, dat$vey)
dat <- cbind(dat, res)

dat$rsqxysog <- cut(dat$rsqxyos, breaks=10)
dat$rsqxysg <- cut(dat$rsqxys, breaks=10)

ggplot(dat, aes(x=rsquxs, y=rsquys)) +
# geom_tile(aes(fill=as.factor((rsqgx-rsqgy)>0))) +
geom_point(aes(colour=as.factor((rsqgx-rsqgy)>0))) +
facet_grid(. ~ rsqxysog) +
labs(x=expression(Signed~R[ux]^2), y=expression(Signed~R[uy]^2), colour="Correct Steiger\ntest result")
ggsave("../images/unmeasured_conf_steiger.pdf", width=10, height=6)


ggplot(dat, aes(x=bux, y=buy)) +
geom_point(aes(colour=as.factor(citres))) +
facet_grid(bxy ~ bgx)

dat$citres <- NA
for(i in 1:nrow(dat))
{
	message(i)
	dat$citres[i] <- do_cit(dat$bxy[i], dat$bgx[i], dat$bux[i], dat$buy[i], dat$vg[i], dat$vu[i], dat$vex[i], dat$vey[i])
}

library(cit)
do_cit <- function(bxy, bgx, bux, buy, vg, vu, vex, vey, n=500)
{
	u <- rnorm(n, sd=sqrt(vu))
	g <- rnorm(n, sd=sqrt(vg))
	x <- bgx * g + bux * u + rnorm(n, sqrt(vex))
	y <- bxy * x + buy * u + rnorm(n, sqrt(vey))

	a <- cit.cp(g, x, y)
	b <- cit.cp(g, y, x)

	return(a[1] > b[1])
}




cov(x,y)/var(x)
cor(g,x)^2
cor(g,y)^2

get_rgy <- function()



cov(g,y)^2



cov(g, y)
var(g) * bgx * cov(x,y) / var(x)
var(g) * bgx * bxy


cov(x,y)
bgx^2*bxy*var(g) + bux^2*bxy*var(u) + bux*buy*var(u) + bxy

var(g)*bgx^2 + var(u)*bux^2 + 1

(bgx^2*bxy*var(g) + bux^2*bxy*var(u) + bux*buy*var(u) + bxy) / (var(g)*bgx^2 + var(u)*bux^2 + 1)


cov(x,y) / var(x)



bxy + (bux / bxy * buy/bxy)




(var(g)*bxy*bgx^2 + bux^2*bxy*var(u) + bux*buy*var(u)) / var(x)

(bxy*bgx^2*var(g) + (bux^2 * bxy + bux * buy) * var(u)) / var(x)



### mediator and collider

n <- 100000
bgx <- 5
bxy <- -2
bux <- -4
buy <- 0
vg <- 1
vu <- 1
vey <- 1
vex <- 1

u <- rnorm(n, sd=sqrt(vu))
g <- rnorm(n, sd=sqrt(vg))

x <- g * bgx + u * bux + rnorm(n, sd=sqrt(vex))
y <- x * bxy + u * buy + rnorm(n, sd=sqrt(vey))

# summary(lm(g ~ u))
# summary(lm(g ~ u + x))


uhat <- u - fitted.values(lm(x ~ g + u))
# summary(lm(g ~ uhat))


xhat <- fitted.values(lm(x ~ g))
yres1 <- residuals(lm(y ~ xhat))


yres2 <- residuals(lm(y ~ x))

# summary(lm(y ~ g))

# summary(lm(y ~ g + x))
# summary(lm(yres1 ~ g))
summary(lm(yres2 ~ g))


bxyo <- ((bgx^2*bxy*vg + bux^2*bxy*vu + bux*buy*vu + bxy) / (vg*bgx^2 + vu*bux^2 + vex))

summary(lm(yres2 ~ g))


(bxy - bxyo) * bgx

