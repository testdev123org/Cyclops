library("testthat")

#
# These tests are believed to be broken; they need confirmation and fixes
#

test_that("Find covariate by name and number", {
	counts <- c(18,17,15,20,10,20,25,13,12)
	outcome <- gl(3,1,9)
	treatment <- gl(3,3)
	tolerance <- 1E-4
	
	glmFit <- glm(counts ~ outcome + treatment, family = poisson()) # gold standard	
	
	dataPtr <- createCcdDataFrame(counts ~ outcome + treatment, 
																modelType = "pr")
	
	ccdFit <- fitCcdModel(dataPtr,
												prior = prior("laplace",																			
												exclude = c("(Intercept)", "outcome2", "outcome3")),
												control = control(noiseLevel = "silent"))

	# Shrinkage on treatment-effects
	expect_less_than(coef(ccdFit)[4], coef(glmFit)[4])
	expect_less_than(coef(ccdFit)[5], coef(glmFit)[5])	
	
	dataPtr2 <- createCcdDataFrame(counts ~ outcome + treatment, 
																modelType = "pr")
	
	
	ccdFit2 <- fitCcdModel(dataPtr2,
												prior = prior("laplace", 
																			exclude = c(1:3)),
												control = control(noiseLevel = "silent"))	
	# Check c(i:j) notation
	expect_equal(coef(ccdFit), coef(ccdFit2))
})

test_that("Simple reductions", {
	counts <- c(18,17,15,20,10,20,25,13,12)
	outcome <- gl(3,1,9)
	treatment <- gl(3,3)
	tolerance <- 1E-4
	
	dataPtr <- createCcdDataFrame(counts ~ outcome + treatment, 
																modelType = "pr")
	
	expect_error(reduce(dataPtr, 0))
	expect_error(reduce(dataPtr, "BAD"))
	expect_equal(reduce(dataPtr, c(1,2)), c(9,3))
	
})

test_that("Error when covariate not found", {
	counts <- c(18,17,15,20,10,20,25,13,12)
	outcome <- gl(3,1,9)
	treatment <- gl(3,3)
	tolerance <- 1E-4
	
	dataPtr <- createCcdDataFrame(counts ~ outcome + treatment, 
																modelType = "pr")
	
	expect_error(
	fitCcdModel(dataPtr,
												prior = prior("laplace", 
																			exclude = c("BAD", "outcome2", "outcome3")),
												control = control(noiseLevel = "silent")))
	
	dataPtr2 <- createCcdDataFrame(counts ~ outcome + treatment, 
																 modelType = "pr")
	
	expect_error(
	fitCcdModel(dataPtr2,
												 prior = prior("laplace", 
												 							exclude = c(10,1:3)),
												 control = control(noiseLevel = "silent")))

})

test_that("Preclude profiling regularized coefficients", {
	counts <- c(18,17,15,20,10,20,25,13,12)
	outcome <- gl(3,1,9)
	treatment <- gl(3,3)
	tolerance <- 1E-4
		
	dataPtr <- createCcdDataFrame(counts ~ outcome + treatment, 
																 modelType = "pr")
	
	ccdFit <- fitCcdModel(dataPtr,
												prior = prior("laplace", exclude = "(Intercept)"),
												control = control(noiseLevel = "silent"))
	
	expect_true(
		!is.null(confint(ccdFit, "(Intercept)")) # not regularized
	)
	expect_error(
		confint(ccdFit, "outcome2") # regularized
	)
	
})

# test_that("Check profile conditional posterior vs likelihood", {
# })

# test_that("Check default regularization variance", {
# })

# test_that("Check starting regularization with cross validation", {
# })

# test_that("Standardize covariates", {
# })

# test_that("Check correct dimensions in matrices in createCcdDataFrame", {
# })

# test_that("Fail to convergence", {
# })

# test_that("Use shared_ptr to handle most data", {
# })

# test_that("Return data summary statistics", {
# })

# test_that("Return segmented reductions", {
# })