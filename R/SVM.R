
################################################################################
## Tutorial
simplePlot <- function(X, y, ...) {
  xlim <- c(min(X[,1]) -1, max(X[,1]) + 1)
  ylim <- c(min(X[,2]) -1, max(X[,2]) + 1)

  is.pos <- y == 1
  matplot(X[is.pos, 1], X[is.pos, 2], pch="+", col="blue", cex=1.5, xlim=xlim,
          ylim=ylim, xlab="", ylab="")
  matplot(X[!is.pos, 1], X[!is.pos, 2], pch="-", col="red", cex=1.5, add=TRUE)
}

meshgrid <- function(a,b) {
  list(x=outer(b*0,a,FUN="+"), y=outer(b,a*0,FUN="+"))
}

showSVM <- function(X, y, kernel='linear', C=1, ...) {
  model <- SVM(X, y, kernel=kernel, C=C, ...)

  ## Get a grid of points to draw decision values over
  ## x1 <- x2 <- seq(-5, 5, length.out=100)

  xlim <- c(min(X[,1] - 1), max(X[,1] + 1))
  ylim <- c(min(X[,2] - 1), max(X[,2] + 1))

  xrange <- range(X)
  x1 <- seq(xlim[1] -1, xlim[2] +1, length.out=100)
  x2 <- seq(ylim[1] -1, ylim[2] +1, length.out=100)

  Xtest <- meshgrid(x1, x2)
  Xtest <- matrix(c(Xtest$x, Xtest$y), ncol=2)

  ## Xtest <- meshgrid(x1, x2)
  ## Xtest <- matrix(c(Xtest$x, Xtest$y), ncol=2)


  out <- predict(model, Xtest, "decision")
  z <- t(matrix(out, 100, 100))
  # cols <- topo.colors(1000)
  cols <- terrain.colors(1000)
  ## browser()
  image(x1, x2, z, col=cols)
  contour(x1, x2, z, add=TRUE)

  ## Get indices to support vectors
  svs <- SVindex(model)
  b <- bias(model)

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
## Manual
if (FALSE) {

}
