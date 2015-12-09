n <- 10000
gx <- rbinom(n, 2, 0.3)
gz <- rbinom(n, 2, 0.3)
u1 <- rnorm(n)
u2 <- rnorm(n)
u3 <- rnorm(n)
ex <- rnorm(n)
ez <- rnorm(n)
ey <- rnorm(n)
x <- 0.3 * gx + u1 + u2 + ex
z <- 0.5 * gz + x + u1 + u3 + ez
y <- x + z + u2 + u3 + ey




standardise <- function(x) 
{
	a <- mean(x)
	return((x - mean(x)) / sd(x) + a)
}


n <- 100000
g <- rbinom(n, 2, 0.3)
e1 <- rnorm(n)
e2 <- rnorm(n)
e3 <- rnorm(n)
x1 <- standardise(0.4 * g + e1)
x2 <- standardise(1 * x1 + e2)
y <- 0.7 * x1 + e3

cov(g, y)

0.4 * 0.7 * var(g)


cov(g, x1)
0.4 * var(g)


cov(g, x2)

0.4 * 0.3 * var(g)



yhat1 <- residuals(lm(y ~ x1))
yhat2 <- residuals(lm(y ~ x2))

cov(yhat1, yhat2)


cov(g, yhat1)
cov(g, yhat2)


cov(g, x1)

var(g) * 0.4

cov(g, x2)
var(g) * 0.4

cov(g, 0.4 * g + e1 + e2)
cov(g, 0.4 * g + e1)

cor(x2, I(0.4 * g + e1 + e2))



lm(y ~ x1)
lm(y ~ x2)

cov(g, x1) / var(g)

cov(g, x2)

cov(g, y)

cov(x1, y)
cov(x2, y)



cov(x1, x2)
cor(x1, x2)


0.3*var(g) + 0.3*var(e1)


library(fGarch)
x <- rsnorm(n, xi=5)
x <- x - min(x) + 1
hist(x)

x1 <- log(x)

cor(x1, x)

library(GenABEL)


a <- rbeta(n, 0.4, 3)
hist(a)

b <- rntransform(a)

plot(a ~ b)
cor(a, b)

