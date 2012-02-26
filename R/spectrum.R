spectrumFeatures <- function(strings, degree=4) {
  features <- lapply(degree, function(deg) {
    xx <- oligonucleotideFrequency(strings, 4)
    xx[colSums(xx) > 0]
  })
  do.call(cbind, features)
}

## Manual spectrum kernel
if (FALSE) {

}
