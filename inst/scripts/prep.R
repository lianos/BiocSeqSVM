################################################################################
## Prepare some data
library(pasilla)

data(pasillaGenes) ## count data for all genes
data(pasillaExons) ## count data for subset of exons

## Parses the pasilla 'Dmel.BDGP5.25.62.DEXSeq.chr.gff' file and saves it
## into inst/extdata/Dmel.DEXSeq.exons.rda
if (FALSE) {
  dm.exons <- parseFlattenedGFF(as.data.table=FALSE)
  save(dm.exons, file='../inst/extdata/Dmel.DEXSeq.exons.rda')
}
