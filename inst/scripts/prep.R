################################################################################
## Prepare some data
library(rtracklayer)
library(pasilla)

data(pasillaGenes) ## count data for all genes
data(pasillaExons) ## count data for subset of exons

## Parses the pasilla 'Dmel.BDGP5.25.62.DEXSeq.chr.gff' file and saves it
## into inst/extdata/Dmel.DEXSeq.exons.rda
if (FALSE) {
  dm.exons <- parseFlattenedGFF(as.data.table=FALSE)
  save(dm.exons, file='../data/Dmel.DEXSeq.exons.rda')
}

## Load the pasilla hit data for later saving
if (FALSE) {
  count.files <- c('untreated1fb.txt', 'untreated2fb.txt', 'untreated3fb.txt',
                   'untreated4fb.txt',
                   'treated1fb.txt', 'treated2fb.txt', 'treated3fb.txt')
  count.files <- system.file('extdata', count.files, package="pasilla")

  ## The count files have comments at the end about number of ambiguous,
  ## or unaliged, etc. reads there are. These lines start with '_'
  counts <- lapply(count.files, read.table, header=FALSE, comment.char="_",
                   stringsAsFactors=FALSE)

  ## Make sure the IDs match
  ids <- counts[[1]][[1]]
  ids.match <- lapply(counts[-1], function(x) ids == x[[1]])
  ids.match <- do.call(cbind, ids.match)
  all(ids.match)

  ids <- Reduce(intersect, lapply(counts, '[[', 1L))
  nrow(counts[[1]]) == length(ids)

  ## Make the table
  ctable <- do.call(cbind, lapply(counts, '[[', 2L))
  colnames(ctable) <- gsub('.txt', '', basename(count.files))

  gene.info <- strsplit(as.character(counts[[1]][[1]]), ':', fixed=TRUE)
  gene.id <- sapply(gene.info, '[[', 1L)
  exon.id <- sapply(gene.info, function(x) as.integer(x[2]))

  p.counts <- cbind(data.frame(gene_id=gene.id, bin=exon.id),
                    as.data.frame(ctable))
  save(p.counts, file="../inst/extdata/pasilla.counts.rda")
}

## Put the dm.exons and p.counts together in one data file, save them
if (FALSE) {
  save(p.counts, dm.exons, file="../data/pasillaData.rda")
}

## Parses the dataset_S2_Nova_CLIP_cluster_mm9.bed dataset from
## Zhang, et al. Argonaute HITS-CLIP decodes microRNA–mRNA interaction maps
hits.dir <- '/Users/stavros/cBio/bioc/data/HITS-CLIP'
if (FALSE) {
  library(SeqTools)
  dl('gc', 'seqs', 'tags')
  mm.peaks <- file.path(hits.dir, 'dataset_S2_Nova_CLIP_cluster_mm9.bed')
  mm.peaks <- import.bed(mm.peaks)
  mm.peaks <- as(mm.peaks, "GRanges")
  mm.anno <- readRDS("/Users/stavros/cBio/projects/TagSeq/inst/extdata/annotated.genome.mm9.rds")
  pa <- annotateReads(mm.peaks, mm.anno)
  nova.peaks <- pa
  save(nova.peaks, file='../NOVA.mm9.rda')
}

## Parse wigFix file for dm3 and mm9 conservation
if (FALSE) {
  ##############################################################################
  ## mm9
  ##
  ## All wigFix files were downloaded from:
  ## http://hgdownload.cse.ucsc.edu/goldenPath/mm9/phyloP30way/vertebrate/
  mm9.cons <- dir('/Users/stavros/cBio/data/conservation/mm9/phyloP27/vertebrate',
                  full.names=TRUE)
  for (file in dir(mm9.dir)) {
    info <- readLines(gzfile(file))
    if (!length(grep("step=1", info[1]))) {
      warning(basename(file), " isn't step=1")
    } else {
      outname <- gsub(".wigFix.gz", ".rds", file)
      xx <- as.double(info[-1])
    }
  }


  ##############################################################################
  ## dm3
  ##
  ## Files were d/led from:
  ## http://hgdownload.cse.ucsc.edu/goldenPath/dm3/phastCons15way/
}

## Data for toy string kernels
##
## We'll use the promoter-gene and splice-junction gene sequences?
##   ftp://ftp.ics.uci.edu//pub/machine-learning-databases/molecular-biology
if (FALSE) {
  pr.data <- read.table('../inst/extdata/promoters.data', sep=',',
                        stringsAsFactors=FALSE)
  pr.data[[3]] <- toupper(gsub("\\s+", "", pr.data[[3]]))
  pr.data[[1]] <- ifelse(pr.data[[1]] == "+", 1L, -1L)
  names(pr.data) <- c('class', 'name', 'sequence')
  promoters <- DNAStringSet(pr.data[[3]])
  values(promoters) <- DataFrame(pr.data[,-3])
  save(promoters, file='../data/promoters.rda')
}
