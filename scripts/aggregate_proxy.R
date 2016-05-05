n <- c(100, 500, 1000, 5000, 10000)
l <- list()
for(i in 1:length(n))
{
	load(paste0("~/repo/cit_measurement_error/results/20160504_", n[i], ".RData"))
	l[[i]] <- dat
}

dat <- do.call(rbind, l)
save(dat, file="~/repo/cit_measurement_error/results/20160504.RData")
