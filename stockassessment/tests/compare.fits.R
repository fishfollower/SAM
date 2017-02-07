compare.fits <<- function(fit.new, fit.exp)
{
	#fit.new and fit.exp are sam objects
	#tolerances were copied from adcomp/tmb_examples/tools/unittest.R
	test_that("same objective", {expect_equal(fit.new$opt$objective, fit.exp$opt$objective, tol=1e-4)})
	test_that("same parameters", {expect_equal(fit.new$opt$par, fit.exp$opt$par, tol=1e-5)})
	test_that("same cov.fixed", {expect_equal(fit.new$sdrep$cov.fixed, fit.exp$sdrep$cov.fixed, tol=1e-4)})
	test_that("same diag.cov", {expect_equal(fit.new$sdrep$diag.cov.random, fit.exp$sdrep$diag.cov.random, tol=1e-7)})
	test_that("same gradient", {expect_equal(fit.new$sdrep$gradient.fixed, fit.exp$sdrep$gradient.fixed, tol=1e-4)})
	test_that("same par.fixed", {expect_equal(fit.new$sdrep$par.fixed, fit.exp$sdrep$par.fixed, tol=1e-5)})
	test_that("same par.random", {expect_equal(fit.new$sdrep$par.random, fit.exp$sdrep$par.random, tol=1e-6)})
	test_that("same sdrep$value", {expect_equal(fit.new$sdrep$value, fit.exp$sdrep$value, tol=1e-6)})
}