compare.fits <<- function(new.fit, exp.fit)
{
	#new.fit and exp.fit are sam objects
	#tolerances were copied from adcomp/tmb_examples/tools/unittest.R
	test_that("same objective", {expect_equal(new.fit$opt$objective, exp.fit$opt$objective, tol=1e-4)})
	test_that("same parameters", {expect_equal(new.fit$opt$par, exp.fit$opt$par, tol=1e-5)})
	test_that("same cov.fixed", {expect_equal(new.fit$sdrep$cov.fixed, exp.fit$sdrep$cov.fixed, tol=1e-4)})
	test_that("same diag.cov", {expect_equal(new.fit$sdrep$diag.cov.random, exp.fit$sdrep$diag.cov.random, tol=1e-7)})
	test_that("same gradient", {expect_equal(new.fit$sdrep$gradient.fixed, exp.fit$sdrep$gradient.fixed, tol=1e-4)})
	test_that("same par.fixed", {expect_equal(new.fit$sdrep$par.fixed, exp.fit$sdrep$par.fixed, tol=1e-5)})
	test_that("same par.random", {expect_equal(new.fit$sdrep$par.random, exp.fit$sdrep$par.random, tol=1e-6)})
	test_that("same sdrep$value", {expect_equal(new.fit$sdrep$value, exp.fit$sdrep$value, tol=1e-6)})
}