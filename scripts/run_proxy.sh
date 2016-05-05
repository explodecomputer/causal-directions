#!/bin/bash

module add libraries/gnu_builds/gsl-1.16

Rscript proxy.R 100 &
Rscript proxy.R 500 &
Rscript proxy.R 1000 &
Rscript proxy.R 5000 &
Rscript proxy.R 10000 &
wait

Rscript aggregate_proxy.R

