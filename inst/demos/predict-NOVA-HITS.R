library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(shikken)

################################################################################
## Local setup
if (FALSE) {
  ag.mm9 <- readRDS('~/cBio/projects/TagSeq/inst/extdata/annotated.genome.mm9.rds')
  ag.dm3 <- readRDS('/Users/stavros/cBio/projects/GenomicCache/GenomicCache.dm3.ensGene/cache/annotated.chromosomes/annotated.collapse-cover.up-500.down-5000.cds-cover-max.rds')

  load('/Users/stavros/cBio/bioc/BiocSeqSVM/data/Dmel.DEXSeq.exons.rda')
  load('/Users/stavros/cBio/bioc/BiocSeqSVM/data/NOVA.mm9.rda')
  load('/Users/stavros/cBio/bioc/BiocSeqSVM/data/pasillaData.rda')

  dex <- read.table('/Users/stavros/cBio/bioc/BiocSeqSVM/data/ps.diff-exons.txt',
                    stringsAsFactors=FALSE, header=TRUE, sep="\t")
  library(GenomicCache)
  library(SeqTools)
  gcm <- GenomicCache("/Users/stavros/cBio/projects/GenomicCache/GenomicCache.mm9.knownGene")
  gcd <- GenomicCache('/Users/stavros/cBio/projects/GenomicCache/GenomicCache.dm3.ensGene')
  n1 <- GFGene("Nova1", gcm)
  n2 <- GFGene("Nova2", gcm)

  dl('seqstore')

  bams <- c(wt1='/Users/stavros/cBio/bioc/data/pasilla/untreated-1/tophat_out/accepted_hits.bam',
            wt2='/Users/stavros/cBio/bioc/data/pasilla/untreated-4/tophat_out/accepted_hits.bam',
            ps1='/Users/stavros/cBio/bioc/data/pasilla/ps-si-1/tophat_out/accepted_hits.bam')
  bams <- as.list(bams)
  bams <- lapply(bams, BamFile)

  ## Fly
  ## http://www-huber.embl.de/pub/DEXSeq/psfb/testForDEU.html
  tenm <- GFGene("Ten-m", gcd) ## exon 7: ~ 22362702
  tranges <- ag.mm9[values(ag.mm9)$symbol == "Ten-m"]
  br <- GFGene("br", gcd) ## exon 12 1537033 chrX [1531317, 1531367]
}

################################################################################
## AMI setup
if (FALSE) {
  ag.mm9 <- readRDS('/home/steve/ml-data/annotated.genome.mm9.rds')
  ag.dm3 <- readRDS('/home/steve/ml-data/annotated.genome.dm3.rds')
  base <- 'http://cbio.mskcc.org/~lianos/files/bioc2012'
  fetch <- c('Dmel.DEXSeq.exons.rda', 'NOVA.mm9.rda', 'pasillaData.rda')
  for (get in fetch) {
    download.file(paste(base, get, sep="/"), get)
    load(get)
  }

  dex <- system.file("data", "ps.diff-exons.txt", package="BiocSeqSVM")
}

################################################################################
## Use spectrum kernel to learn preferred binding landscape of NOVA in mousr
head(nova.peaks)
dt.mm9 <- as(ag.mm9, 'data.table')
dt.nova <- as(nova.peaks, 'data.table')

## Where does NOVA like to bind on the genome?
nova.summary <- dt.nova[, list(score=sum(score)), by='exon.anno']

gg.angle <- opts(axis.text.x=theme_text(angle=-45, hjust=0, vjust=1))
g <- ggplot(nova.summary, aes(x=exon.anno, y=score, fill=exon.anno)) +
  geom_bar(stat="identity") + theme_bw() + gg.angle +
  opts(title="NOVA binding site regions in mm9")
print(g)

## Let's remove intergnic hits
dt.nova <- subset(dt.nova, exon.anno != 'intergenic')
dt.mm9 <- subset(dt.mm9, entrez.id %in% dt.nova$entrez.id)

cbound <- subset(dt.nova, exon.anno == "cds")
cbound <- cbound[order(cbound$score, decreasing=TRUE),]

ibound <- subset(dt.nova, exon.anno == "intron")
ibound <- ibound[order(ibound$score, decreasing=TRUE),]
summary(ibound$width[1:500])
summary(ibound$score[1:500])

ubound <- subset(dt.nova, exon.anno == "utr3")
ubound <- ubound[order(ubound$score, decreasing=TRUE),]
summary(ubound$score[1:500])

## I know I should be using ggplot2, but ...
## Does binding in intron vs 3'utr look different?
plot.densities(trim.data(cbound$width, 0.05),
               trim.data(ibound$width, 0.05),
               trim.data(ubound$width, 0.05),
               legend=c('cds', 'intron', 'utr3'),
               main="Binding widths")

## How about their scores?
plot.densities(trim.data(cbound$score, 0.05),
               trim.data(ibound$score, 0.05),
               trim.data(ubound$score, 0.05),
               legend=c('cds', 'intron', 'utr3'),
               main="Score distro")

ipos <- subset(ibound, width > 20 & width < 100 & score > 10)
cpos <- subset(cbound, width > 20 & width < 100 & score > 10)
pos <- as(rbind(ipos, cpos), "GRanges")
posdna <- getSeq(Mmusculus, pos)

gexpr <- dt.nova$entrez.id
ineg <-
################################################################################
## dm3
dt.dm3 <- as(ag.dm3, 'data.table')
flybase <- mget(ifelse(is.na(dt.dm3$entrez.id), 'XXX', dt.dm3$entrez.id),
                org.Dm.egFLYBASE, ifnotfound=NA)
dt.dm3$fb <- sapply(flybase, '[[', 1L)




