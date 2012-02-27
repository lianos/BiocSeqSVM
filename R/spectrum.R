spectrumFeatures <- function(strings, degree=4) {
  features <- lapply(degree, function(deg) {
    xx <- oligonucleotideFrequency(strings, 4)
    xx[, colSums(xx) > 0, drop=FALSE]
  })
  if (length(features) > 1) {
    ans <- do.call(cbind, features)
  } else {
    ans <- features[[1]]
  }
  ans
}

## Manual spectrum kernel
if (FALSE) {

}
