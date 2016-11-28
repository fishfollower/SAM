context("nsher example")
setwd("../nsher/")

stopifnot(require("testthat"),
		require("stockassessment"))
source("../compare.fits.R")
source("script.R")
load("fit.expected.Rdata")#reads in fit.exp
compare.fits(fit, fit.exp)

