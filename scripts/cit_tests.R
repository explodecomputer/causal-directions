## ---- cit_tests ----

suppressWarnings(suppressPackageStartupMessages(library(network)))
suppressWarnings(suppressPackageStartupMessages(library(sna)))
suppressWarnings(suppressPackageStartupMessages(library(ggnetwork)))
suppressWarnings(suppressPackageStartupMessages(library(gridExtra)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

nid <- 10000
l <- m <- list()

# 1

mat <- matrix(c(
	0,0,0,
	1,0,0,
	0,1,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[1]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 1")

g <- rnorm(nid)
a <- g + rnorm(nid)
b <- a + rnorm(nid)
m[[1]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[1]] <- as.data.frame(m[[1]])
m[[1]]$model <- 1
m[[1]]$direction <- c("Testing A to B", "Testing B to A")


# 2

mat <- matrix(c(
	0,0,0,
	0,0,1,
	1,0,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[2]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 2")

g <- rnorm(nid)
b <- g + rnorm(nid)
a <- b + rnorm(nid)
m[[2]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[2]] <- as.data.frame(m[[2]])
m[[2]]$model <- 2
m[[2]]$direction <- c("Testing A to B", "Testing B to A")



# 3

mat <- matrix(c(
	0,0,0,
	1,0,0,
	1,0,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[3]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 3")

g <- rnorm(nid)
a <- g + rnorm(nid)
b <- g + rnorm(nid)
m[[3]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[3]] <- as.data.frame(m[[3]])
m[[3]]$model <- 3
m[[3]]$direction <- c("Testing A to B", "Testing B to A")



# 4

mat <- matrix(c(
	0,0,0,
	1,0,1,
	0,0,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[4]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 4")

g <- rnorm(nid)
b <- rnorm(nid)
a <- g + b + rnorm(nid)
m[[4]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[4]] <- as.data.frame(m[[4]])
m[[4]]$model <- 4
m[[4]]$direction <- c("Testing A to B", "Testing B to A")


# 5

mat <- matrix(c(
	0,0,0,
	0,0,0,
	1,1,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[5]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 5")

g <- rnorm(nid)
a <- rnorm(nid)
b <- g + a + rnorm(nid)
m[[5]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[5]] <- as.data.frame(m[[5]])
m[[5]]$model <- 5
m[[5]]$direction <- c("Testing A to B", "Testing B to A")


# 6

mat <- matrix(c(
	0,0,0,
	1,0,0,
	1,1,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[6]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 6")

g <- rnorm(nid)
a <- g + rnorm(nid)
b <- a + g + rnorm(nid)
m[[6]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[6]] <- as.data.frame(m[[6]])
m[[6]]$model <- 6
m[[6]]$direction <- c("Testing A to B", "Testing B to A")


# 7

mat <- matrix(c(
	0,0,0,
	1,0,1,
	1,0,0
), 3, 3)

n <- network(mat, directed=TRUE)
n %v% "what" <- c("G", "A", "B")
l[[7]] <- ggplot(ggnetwork(n, layout="circle", arrow.gap=0.05), aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_nodes(size = 11, aes(color=what)) +
  geom_edges(arrow = arrow(type = "closed")) +
  geom_nodelabel(aes(label = what), fontface = "bold") +
  theme_blank() +
  theme(legend.position="none") +
  labs(title="Model 7")

g <- rnorm(nid)
b <- g + rnorm(nid)
a <- b + g + rnorm(nid)
m[[7]] <- rbind(
	cit.cp(L=g, G=a, T=b),
	cit.cp(L=g, G=b, T=a)
)
m[[7]] <- as.data.frame(m[[7]])
m[[7]]$model <- 7
m[[7]]$direction <- c("Testing A to B", "Testing B to A")



#####

res <- bind_rows(m)
res <- gather(res, key=test, value=pval, p_cit, p_TassocL, p_TassocGgvnL, p_GassocLgvnT, p_LindTgvnG)
res$sig <- res$pval < 0.0001
res$test <- as.factor(res$test)
levels(res$test) <- c("Omnibus CIT", "1. G ~ L|T", "4. L !~ T|G", "2. T ~ G|L", "3. T ~ L")
res$test <- as.character(res$test)

dots <- ggplot(res, aes(x=test, y=as.factor(model))) +
geom_point(aes(colour=sig), size=5) +
facet_grid(. ~ direction) +
theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) +
labs(x="CIT test components", y="Causal model", colour="P-value < 0.05")


models <- grid.arrange(grobs=l, ncol=3, newpage=FALSE)

save(models, dots, file="../results/cit_tests.RData")


####


load("../results/cit_tests.RData")
grid.arrange(models, dots, ncol=2)


