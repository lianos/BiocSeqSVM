library(rtracklayer)
library(data.table)

##' Parses a DEXSeq-flattened GFF into data.frame
parseFlattenedGFF <- function(efile, as.data.table=FALSE) {
  if (missing(efile)) {
    efile <- system.file('extdata', 'Dmel.BDGP5.25.62.DEXSeq.chr.gff',
                         package="pasilla")
  }
  all.exons <- as(import.gff(efile), "GRanges")
  all.group <- .parseGroup(values(all.exons)$group)

  exons <- cbind(as(all.exons, 'data.frame'), all.group)
  exons <- subset(exons, type == "exonic_part")
  exons$bin <- as.integer(as.character(exons$exonic_part_number))

  ## To reconstruct exon__part_number you do:
  ## formatC(e$bin, 2, flag='0')
  col.order <- c('gene_id', 'bin', 'seqnames', 'start', 'end', 'width',
                 'strand', 'transcripts')
  exons <- exons[, col.order]
  exons <- data.table(exons, key=c("seqnames", "strand", "start"))

  ## The current exon is split into two bins.
  ## The exon on the following row can be a continuiation of this one.
  run.on.fwd <- with(exons, {
    ifelse(head(strand, -1L) == '+' & head(end, -1L) + 1L == tail(start, -1L),
           TRUE, FALSE)
  })
  run.on.fwd <- c(run.on.fwd, FALSE)

  run.on.rev <- with(exons, {
    ifelse(tail(strand, -1L) == '-' & tail(start, -1L) - 1L == head(end, -1L),
           TRUE, FALSE)
  })
  run.on.rev <- c(FALSE, run.on.rev)
  exons$run.on <- run.on.fwd | run.on.rev

  setkeyv(exons, c("seqnames", "strand", "gene_id", "bin"))

  if (!as.data.table) {
    exons <- as.data.frame(exons)
  }

  exons
}

if (FALSE) {
  dm.exons <- parseFlattenedGFF(as.data.table=FALSE)
  save(dm.exons, file='../inst/extdata/Dmel.DEXSeq.exons.rda')


  ## Parse the pasilla flattened annotations
  afile <-
  groups <- values(dm.exons)$group
  dm.group <-
  values(dm.exons)$source <- NULL
  values(dm.exons)$group <- NULL
  values(dm.exons) <- cbind(values(dm.exons), DataFrame(dm.group))

  exons <- subset(as.data.frame(dm.exons), type == "exonic_part")


  gb <- exons[, list(seqnames=seqnames[1L], start=min(start),
                     end=max(end), strand=strand[1L]),
              by=gene_id]
  gb <- transform(gb, width=end - start + 1L)
  gbounds <- subset(as.data.frame(dm.exons), type == "aggregate_gene")
  gbounds <- data.table(gbounds, key="gene_id")
}
library(pasilla)

## pasillaGenes has count data for all genes
data(pasillaGenes)

## pasillaExons has count data for
data(pasillaExons)

gcounts <- counts(pasillaGenes)
ecounts <- counts(pasillaExons)


