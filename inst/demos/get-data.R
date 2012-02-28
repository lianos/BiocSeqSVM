library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ggplot2)

libary(devtools)
install_github('shikken', 'lianos', 'staticwrap')

get.me <- c(
  "mm9.anno.rds",      ## has compressed mm9 genome annotation
  "NOVA.mm9.rda",      ## post-processed NOVA hits-clip file
  "dm3.anno.rds",      ## has compressed dm3 genome annotation
  "pasillaData.rda")   ## dm.exons

base.url <- 'http://cbio.mskcc.org/~lianos/files/bioc2012'

system('mkdir svm')

for (get in get.me) {
  download.file(paste(base.url, get, sep="/"),
                paste('svm', get, sep="/"))
}

ag.mm9 <- readRDS('svm/mm9.anno.rds')
ag.dm3 <- readRDS('svm/dm3.anno.rds')
load('svm/NOVA.mm9.rda')
load('svm/pasillaData.rda')

