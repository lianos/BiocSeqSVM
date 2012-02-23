## Play here after you're done with prep.R
load('../data/pasillaData.rda')

plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black",
                       sep=0.5, ...) {
  if (inherits(x, "GenomicRanges")) {
    stopifnot(length(unique(seqnames(x))) == 1 & length(unique(strand(x)) == 1))
    x <- ranges(x)
  }
  height <- 1
  if (is(xlim, "Ranges")) {
    xlim <- c(min(start(xlim)), max(end(xlim)))
  }
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins) * (height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main) +	axis(1)
}


## Let's make our own string kernel
if (FALSE) {

}
