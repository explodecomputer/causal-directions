library(dplyr)

splits <- 100
savefile <- "~/repo/cit_measurement_error/results/mr_directionality_20170126.RData"

l <- list()
for(i in 1:splits)
{
	cat(i, "\n")
	filename <- paste0("~/repo/cit_measurement_error/scratch/resultswc_", i, ".RData")
	if(file.exists(filename))
	{
		load(filename)
		l[[i]] <- parameters
	} else {
		message("Missing ", filename)
	}
}

parameters <- rbind_all(l)
save(parameters, file = savefile)