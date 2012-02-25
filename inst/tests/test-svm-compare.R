context("SVM comparisons")

data("promoters", package="MLplay")

test_that("custom spectrum kernel == kernlab,stringkernel", {
  sfeatures <- SpectrumFeatures(promoters, length=4)
  y <- values(promoters)$class

  m1 <- ksvm(sfeatures, y, kernel="vanilla", C=5, type="C-svc")
  m2 <- ksvm(promoters, y, kernel="spectrum", length=4, type="C-svc", C=5)
})


test_that("kernlab spectrum kernel == shogun", {

})

test_that("custom spectrum kernel == shogun", {

})
