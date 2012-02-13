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

##' Parses the 'group' column in a gff/gtf file into a data.frame
##'
##' This function assumes that the group entry is like so:
##'
##'   attrib1 "value"; attrib2 "value"; attrib3 "value"; ...
##'
##' @param group The "group" column from a GFF/GTF file
##' @param group The separator used two separate attributes
##' @param N The number of lines to read from the head of the GFF to guess the
##' max number of columns the output data.frame will have.
.parseGroup <- function(group, sep=";", N=100) {
  n <- sapply(gregexpr(sep, head(group, N)), length)
  N.cols <- max(n)

  group <- strsplit(group, sep, fixed=TRUE)
  example <- .parseGroupLineList(group[[which(n == N.cols)[1L]]])
  var.names <- names(example)

  ans <- lapply(example, function(x) character(length(info)))

  for (i in 1:length(group)) {
    this <- .parseGroupLineList(group[[i]])
    use <- names(this)[names(this) %in% var.names]
    for (name in use) {
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

  ## import.gff column renaming (version 1?):
  ##   feature -> type
  ##   frame -> phase
  annos <- import.gff(gff, version, genome)

  for (chr in names(annos)) {
    a <- as(annos[chr], "GRanges")
    group <- .parseGroup(values(a)$group)
  }
}
