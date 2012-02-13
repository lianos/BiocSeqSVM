if (FALSE) {
  library(rtracklayer)

  ## gff <- system.file("extdata", "Drosophila_melanogaster.BDGP5.25.62.gtf.gz",
  ##                    package="MLplay")
  ## gff <- system.file("extdata", 'Dmel.BDGP5.25.62.gtf.rds', package="MLplay")
  gff <- file.path("/Users/stavros/cBio/bioc/MLplay/inst/extdata",
                   'Dmel.BDGP5.25.62.gtf.rds')
  annos <- readRDS(gff)
  summary.gff <- system.file("extdata", "Dmel.BDGP5.25.62.DEXSeq.chr.gff",
                             package="pasilla")
}

## From Marc
## But if you wanted to do such a thing, I would recommend parsing the GFF,
## and then calling the generalized makeTranscriptDb() function, and then
## after that looking at the generalized makeTxDbPackage() function.  From
## there, I think you can see how it should be pretty easy to extend things
## to make something like a pair of makeTranscriptDbFromGFF() and a
## makeTxDbPackageFromGFF() functions.

.parseGroupLineList <- function(x, max.n, var.names, sep=";") {
  x <- strsplit(x, "\\s", perl=TRUE)
  x <- lapply(x, function(xx) xx[nchar(xx) > 0])

  if (!missing(max.n)) {
    if (length(x) > max.n) {
      stop("Illegal number of columns")
    }
  }

  vals <- lapply(x, function(xx) gsub("'", "", gsub('"', "", xx[[2L]])))
  names(vals) <- sapply(x, '[[', 1L)

  if (!missing(var.names)) {
    if (!all(vnames %in% var.names)) {
      stop("Unexpected variable name")
    }
    vals <- vals[var.names]
  }

  vals
}

.parseGroup <- function(group, sep=";", N=50) {
  ##############################################################################
  ## Figure out number of columns
  n.cols <- sapply(gregexpr(sep, head(group, N)), length)
  max.n <- max(n.cols)
  info <- strsplit(group, sep, fixed=TRUE)

  max.parse <- .parseGroupLineList(info[[max.n]])
  ans <- lapply(max.parse, function(x) character(length(info)))
  var.names <- names(max.parse)
  names(ans) <- var.names

  for (i in 1:length(info)) {
    this <- .parseGroupLineList(info[[i]])
    for (name in var.names) {
      ans[[name]][i] <- this[[name]]
    }
  }

  as.data.frame(ans)
}

makeTranscriptDbFromGFF <- function(gff, version=c("1", "2", "3"),
                                    genome="hg19",
                                    circ_seqs=DEFAULT_CIRC_SEQS) {
  if (!file.exists(gff)) {
    stop("Unable to read file: ", gff)
  }
  version <- match.arg(version)
  annos <- import.gff(gff, version, genome)

  for (chr in names(annos)) {
    a <- as(annos[chr], "GRanges")
    meta.list <- strsplit(values(a)$group, ";", fixed=TRUE)

    #############################################################
    ## Figure out number of columns

    parsed <- lapply(meta.list, function(x) {
      x <- strsplit(x, "\\s", perl=TRUE)
      x <- lapply(x, function(xx) xx[nchar(xx) > 0])
      vals <- lapply(x, function(xx) gsub("'", "", gsub('"', "", xx[[2L]])))
      names(vals) <- sapply(x, '[[', 1L)
      as.data.frame(vals)
    })
  }
}
