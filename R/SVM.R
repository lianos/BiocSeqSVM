
################################################################################
## Tutorial

##' Utility function to plot data in the dead-center of the plot
##'
##' @param X 2D data matrix
##' @param y Numeric labels vector, as long as their are rows as X
simplePlot <- function(X, y, ...) {
  stopifnot(ncol(X) == 2)
  stopifnot(nrow(X) == length(y))

  xlim <- c(min(X[,1]) -1, max(X[,1]) + 1)
  ylim <- c(min(X[,2]) -1, max(X[,2]) + 1)

  is.pos <- y == 1
  matplot(X[is.pos, 1], X[is.pos, 2], pch="+", col="blue", cex=1.5, xlim=xlim,
          ylim=ylim, xlab="", ylab="")
  matplot(X[!is.pos, 1], X[!is.pos, 2], pch="-", col="red", cex=1.5, add=TRUE)
}


##' Creates a grid of values to cover the surface of a plot
##'
##' Taken from shogun-toolbox's svm_classification.R function found in
##' the r_static interface
##'
##' @param a vector of points to tile along the rows
##' @param a vector of points to tile along the columns
meshgrid <- function(a,b) {
  list(x=outer(b*0,a,FUN="+"), y=outer(b,a*0,FUN="+"))
}

##' Plots the decision surface of an SVM on 2D data.
##'
##' Adapted from shogun-toolbox's svm_classification.R function found in
##' the r_static interface
##'
##' @param X A 2D data matrix. Rows are observations, columns are features
##' @param y Label vector: has as many entries as there are rows in X
##' @param C The cost parameter for the SVM
##' @param ... Arguments to pass to the SVM contstructor
##' @param .pacakge Either 'shikken' or 'kernlab' -- determines which
##' SVM package to use to build the classifier.
showSVM <- function(X, y, kernel='linear', C=1, ...,
                    .package=c("shikken", "kernlab")) {
  .package <- match.arg(.package)
  args <- list(...)
  if (.package == "shikken") {
    svmf <- SVM
    args$type <- NULL
    args <- c(list(X, y, kernel=kernel, C=C), args)
  } else if (.package == "kernlab") {
    svmf <- ksvm
    if (kernel == 'linear') {
      kernel <- 'vanilla'
    } else if (kernel == 'gaussian') {
      kernel <- "rbf"
      sigma <- args$width
      if (is.null(sigma)) {
        sigma <- 1
      }
      args$kpar <- list(sigma=1/sigma)
    }
    args <- c(list(X, y, kernel=kernel, C=C, type="C-svc"), args)
  }
  ## model <- SVM(X, y, kernel=kernel, C=C, ...)
  model <- do.call(svmf, args)

  ## Calculate grid of points to draw decision values over that will cover
  ## the range of our data
  xlim <- c(min(X[,1] - 1), max(X[,1] + 1))
  ylim <- c(min(X[,2] - 1), max(X[,2] + 1))
  x1 <- seq(xlim[1] - 1, xlim[2] + 1, length.out=100)
  x2 <- seq(ylim[1] - 1, ylim[2] + 1, length.out=100)
  Xtest <- meshgrid(x1, x2)
  Xtest <- matrix(c(Xtest$x, Xtest$y), ncol=2)

  ## Get the decision values from the SVM for the the test data.
  out <- predict(model, Xtest, "decision")

  ## Now it's time to make pretty things
  z <- t(matrix(out, 100, 100))
  cols <- terrain.colors(1000)

  ## image(x1, x2, z, col=cols)
  ## image.plot in fields package gives us a handy color bar/legend
  image.plot(x1, x2, z, col=terrain.colors(50))
  contour(x1, x2, z, add=TRUE)

  ## Get indices to support vectors
  svs <- SVindex(model)

  posSVs <- X[y ==  1 & 1:nrow(X) %in% svs,, drop=FALSE]
  negSVs <- X[y == -1 & 1:nrow(X) %in% svs,, drop=FALSE]

  pos <- X[y ==  1 & !1:nrow(X) %in% svs, ]
  neg <- X[y == -1 & !1:nrow(X) %in% svs, ]

  matplot(posSVs[,1], posSVs[,2], pch="+", col="red", add=TRUE, cex=1.5)
  matplot(negSVs[,1], negSVs[,2], pch="-", col="red", add=TRUE, cex=1.5)

  matplot(pos[,1], pos[,2], pch="+", col="black",add=TRUE, cex=1)
  matplot(neg[,1], neg[,2], pch="-", col="black",add=TRUE, cex=1)
  title(paste("Decision surface for", kernel, "SVM"))

  invisible(model)
}