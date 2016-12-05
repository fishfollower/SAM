run.compare <<- function(name){
	#name is a string
	context(paste0(name, " example"))
	setwd(paste0("../", name, "/"))
	stopifnot(require("testthat"),
		require("stockassessment"))
	source("../compare.fits.R")
	source("script.R")
	load("fit.expected.Rdata")#reads in fit.exp
	compare.fits(fit, fit.exp)
}