library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(shikken)

################################################################################
## Use spectrum kernel to learn preferred binding landscape of NOVA in mouse
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

## Let's remove intergnic hits from the HITS data
dt.nova <- subset(dt.nova, exon.anno != 'intergenic')

## Subset our annotated genome to only include genes that are bound
## by NOVA -- the idea is that we only want to work with genes that are
## expressed, but we don't have any expression data, so this is the
## best we can do for now.
dt.mm9 <- subset(dt.mm9, entrez.id %in% dt.nova$entrez.id & !is.na(entrez.id))

## filter annotation to "expressed" genic regions only
ag.mm9 <- as(dt.mm9, "GRanges")

## Lets break down genic binding sites into different categories.
## Ones that bind in:
##   - introns (maybe the affect up/downstream splicing)
##   - cds  (maybe the affect their own splicing
##   - utr3 (you can take a guess at why it might bind here)
cbound <- subset(dt.nova, exon.anno == "cds")
cbound <- cbound[order(cbound$score, decreasing=TRUE),]

ibound <- subset(dt.nova, exon.anno == "intron")
ibound <- ibound[order(ibound$score, decreasing=TRUE),]
summary(ibound$width[1:500])
summary(ibound$score[1:500])

ubound <- subset(dt.nova, exon.anno == "utr3")
ubound <- ubound[order(ubound$score, decreasing=TRUE),]
summary(ubound$score[1:500])

## What do binding sites in these different regions look like?
## Are the peaks different widths?
## Are their "scores" different?

## I know I should be using ggplot2, but ...
## Does binding in intron vs 3'utr look different?
plot.densities(trim.data(cbound$width, 0.05),
               trim.data(ibound$width, 0.05),
               trim.data(ubound$width, 0.05),
               legend=c('cds', 'intron', 'utr3'),
               main="Binding widths")

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## This is why you usually want to start w/ the raw reads your self
## and process it in a way that makes sense -- or at least in a way
## that makes sense to you.

## What does their score distribution look like?
plot.densities(trim.data(cbound$score, 0.05),
               trim.data(ibound$score, 0.05),
               trim.data(ubound$score, 0.05),
               legend=c('cds', 'intron', 'utr3'),
               main="Score distro")
## Wonderful -- what's up with intronic binding sites?
## We may never know

################################################################################
## OK, let's take a closer look at what intronic binding events look like

## Lets get some "good" peaks that land in introns
ipos <- subset(ibound, width > 20 & width < 100 & score > 10)

## where are NOVA binding sites located in the intron? Are they biased
## towards the start? end? random?
igr <- as(ipos, "GRanges")
ihost <- subsetByOverlaps(ag.mm9[values(ag.mm9)$exon.anno == 'intron'], igr)
ixref <- match(igr, ihost)
iinfo <- cbind(ipos, as(ihost, 'data.table')[ixref,])
head(iinfo)

## Sanity check
all(iinfo$symbol == iinfo$symbol.1)


## How far along the intron is the binding site?
## We assume the protein binds at center of peak.
ifar.along <- with(iinfo, {
  pos <- ifelse(strand == '+', end, start)
  ifelse(strand == '+',
         (pos - start.1) / width.1,
         (end.1 - pos) / width.1)
})
## peaks that straddle intron/exon boundaries are problematic
## Look at IGV to convince yourself numbers > 1 are due to this
head(subset(iinfo, ifar.along > 1))
ifar.along[ifar.along < 0] <- 0
ifar.along[ifar.along > 1] <- 1
iinfo$far.along <- ifar.along

plot(density(subset(iinfo, width.1 > 1000)$far.along))
## Geez -- that was as surprise.


## I'm going to assert that NOVA only affect splicing when it binds
## close to an exon
icansplice <- subset(iinfo, far.along < .1 | far.along > .9)
## is there a bias now?
plot(density(icansplice$far.along))
icansplice.dna <- getSeq(Mmusculus, as(icansplice, "GRanges"))

## NOVA seems to prefer binding down stream of exons -- ~66 percent
## of the time.
pr.ds <- sum(icansplice$far.along < .4) / nrow(icansplice)

## -----------------------------------------------------------------------------
## Get a negative set
##
## How about we get introns that are expressed but are not bound?
inonova <- ag.mm9[values(ag.mm9)$exon.anno == 'intron' &
                  values(ag.mm9)$symbol %in% iinfo$symbol]
inonova <- inonova[countOverlaps(inonova, nova.peaks) == 0]

## make sure we got what we asked for
all(values(inonova)$symbol %in% ipos$symbol)
length(subsetByOverlaps(inonova, igr)) == 0

## randomize just for fun
inonova <- inonova[sample(length(inonova))]

## sample upstream/downstream regions of exons according to same distro
## as positive set
take.ds <- rbinom(length(inonova), 1, pr.ds)
resize.fix <- ifelse(take.ds == 1, 'start', 'end')
neg.gr <- resize(inonova, width=70, fix=resize.fix)

## get sequences upstream and downstream of exons
neg.dna <- getSeq(Mmusculus, neg.gr)

X <- c(icansplice.dna, neg.dna)
y <- c(rep(1, length(icansplice.dna)), rep(-1, length(neg.dna)))

## C.neg and C are swapped in shikken now (bug) -- wil fix!
## This training take a long time
m <- SVM(X, y, kernel="spectrum", degree=4, C.neg=100000, C=10, threads=2)

## Some sequences have N -- this is illegal, The index of these
## sequences have been returned to m

X <- X[-m]
y <- y[-m]

## TODO: train SVM w/o imbalance C's

m <- SVM(X, y, kernel="spectrum", degree=4, C.neg=100000, C=1, threads=2)
preds <- predict(m, X, type="decision")
sum(sign(preds) == y) / length(y)               ## 0.856
sum(sign(preds) == 1 & y == 1) / sum(y == 1)    ## 1!
sum(sign(preds) == -1 & y == -1) / sum(y == -1) ## 0.848

## predict on all peaks inside introns now
nintrons <- nova.peaks[values(nova.peaks)$exon.anno == "intron"]
iall.seqs <- getSeq(Mmusculus, nintrons)
iall.preds <- predict(m, iall.seqs, type="decision")
ineg.preds <- preds[y == -1]

plot(iall.preds, values(nintrons)$score, pch=19, col="#3f3f3fff")
title("Decision value for intronic peaks")

## Note that we only trained on 1506 positive examples and 26,432
##
## We are trying to predict on 151,775 positive examples!
plot.densities(iall.preds, ineg.preds, legend=c('peak', 'nopeak'),
               main="Decision value for SVM")

################################################################################
## dm3
## Can we transfer SVM to fly?
library(BSgenome.Dmelanogaster.UCSC.dm3)
dex <- read.table('http://cbio.mskcc.org/~lianos/files/bioc2012/ps.diff-exons-some.txt',
                  stringsAsFactors=FALSE, header=TRUE)
dex <- as.data.frame(dex)
dex <- transform(dex, key=paste(gene_id, formatC(bin, 2, flag='0'), sep=":"))


## intersect results with annotation
dm.exons <- as.data.frame(dm.exons)
dm.exons <- transform(dm.exons,
                      key=paste(gene_id, formatC(bin, 2, flag='0'), sep=":"))

dm.expr <- subset(dm.exons, gene_id %in% dex$gene_id)
dm.expr$is.dex <- dm.expr$key %in% dex$key

dm.pos <- as(subset(dm.expr, is.dex), "GRanges")
dm.pos.up <- flank(dm.pos, 60, start=TRUE) + 5
dm.pos.up.seq <- getSeq(Dmelanogaster, dm.pos.up)

preds.pos.up <- predict(m, dm.pos.up.seq, type="decision")


dm.pos.dn <- flank(dm.pos, 60, start=FALSE) + 5
dm.pos.dn.seq <- getSeq(Dmelanogaster, dm.pos.dn)

preds.pos.down <- predict(m, dm.pos.dn.seq, type="decision")

## get max decision value
preds.pos <- pmax(preds.pos.up, preds.pos.down)


dm.neg <- as(subset(dm.expr, !is.dex), "GRanges")
dm.neg.up <- flank(dm.neg, 70, start=TRUE)
dm.neg.up.seq <- getSeq(Dmelanogaster, dm.neg.up, as.character=FALSE)
preds.neg.up <- predict(m, dm.neg.up.seq, type="decision")

dm.neg.dn <- flank(dm.neg, 70, start=FALSE)
dm.neg.dn.seq <- getSeq(Dmelanogaster, dm.neg.dn, as.character=FALSE)
preds.neg.dn <- predict(m, dm.neg.dn.seq, type="decision")

preds.neg <- pmax(preds.neg.up, preds.neg.dn)

boxplot(list(pos=preds.pos, neg=preds.neg))

## call expression of genes
## expt.names <- names(p.counts)[-(1:2)]
## p.counts <- data.table(p.counts, key=c('gene_id', 'bin'))
## dm.expr <- p.counts[, {
##   list(wt=sum(untreated1fb, untreated2fb, untreated3fb, untreated4fb),
##        ps=sum(treated1fb, treated2fb, treated3fb))
## }, by="gene_id"]

## extract genes from dm.exons that have *a* differentially spliced exon
dm.expr <- subset()



## The DEXSeq results are indexed by flybase gene id, we need to xref
## those to entrz.id in my genome annotation table
dt.dm3 <- as(ag.dm3, 'data.table')

flybase <- mget(ifelse(is.na(dt.dm3$entrez.id), 'XXX', dt.dm3$entrez.id),
                org.Dm.egFLYBASE, ifnotfound=NA)
dt.dm3$geneID <- sapply(flybase, '[[', 1L)

## Find the genes that have differentially spliced exons
dex.genes <- subset(dt.dm3, !is.na(geneID) & geneID %in% dex$geneID)

## Identify which cds exons are not differentially spliced
no.dex <- as(subset(dex.genes, exon.anno == 'cds'), "GRanges")
no.dex <- no.dex[countOverlaps(no.dex, dex.gr) == 0]



##################################################################################
## Scratch
################################################################################
## Local setup
# if (FALSE) {
#   ag.mm9 <- readRDS('~/cBio/projects/TagSeq/inst/extdata/annotated.genome.mm9.rds')
#   ag.dm3 <- readRDS('/Users/stavros/cBio/projects/GenomicCache/GenomicCache.dm3.ensGene/cache/annotated.chromosomes/annotated.collapse-cover.up-500.down-5000.cds-cover-max.rds')
#
#   ## dm3 is duplicated!
#   if (FALSE) {
#     dt <- as(ag.dm3, 'data.table')
#     setkeyv(dt, c('seqnames', 'strand', 'start'))
#     u <- unique(dt)
#     ag.dm3 <- as(u, 'GRanges')
#     saveRDS(ag.dm3, '/Users/stavros/cBio/projects/GenomicCache/GenomicCache.dm3.ensGene/cache/annotated.chromosomes/annotated.collapse-cover.up-500.down-5000.cds-cover-max.rds')
#   }
#
#   ## dm.exons
#   ## load('/Users/stavros/cBio/bioc/BiocSeqSVM/data/Dmel.DEXSeq.exons.rda')
#
#   ## nova.peaks
#   load('/Users/stavros/cBio/bioc/BiocSeqSVM/data/NOVA.mm9.rda')
#
#   ## p.counts, dm.exons
#   load('/Users/stavros/cBio/bioc/BiocSeqSVM/data/pasillaData.rda')
#
#   dex <- read.table('/Users/stavros/cBio/bioc/BiocSeqSVM/data/ps.diff-exons-some.txt',
#                     stringsAsFactors=FALSE, header=TRUE, sep="\t")
#
#
#   library(GenomicCache)
#   library(SeqTools)
#   gcm <- GenomicCache("/Users/stavros/cBio/projects/GenomicCache/GenomicCache.mm9.knownGene")
#   gcd <- GenomicCache('/Users/stavros/cBio/projects/GenomicCache/GenomicCache.dm3.ensGene')
#   n1 <- GFGene("Nova1", gcm)
#   n2 <- GFGene("Nova2", gcm)
#
#   dl('seqstore')
#
#   bams <- c(wt1='/Users/stavros/cBio/bioc/data/pasilla/untreated-1/tophat_out/accepted_hits.bam',
#             wt2='/Users/stavros/cBio/bioc/data/pasilla/untreated-4/tophat_out/accepted_hits.bam',
#             ps1='/Users/stavros/cBio/bioc/data/pasilla/ps-si-1/tophat_out/accepted_hits.bam')
#   bams <- as.list(bams)
#   bams <- lapply(bams, BamFile)
#
#   ## Fly
#   ## http://www-huber.embl.de/pub/DEXSeq/psfb/testForDEU.html
#   tenm <- GFGene("Ten-m", gcd) ## exon 7: ~ 22362702
#   tranges <- ag.mm9[values(ag.mm9)$symbol == "Ten-m"]
#   br <- GFGene("br", gcd) ## exon 12 1537033 chrX [1531317, 1531367]
# }
#
# ################################################################################
# ## AMI setup
# if (FALSE) {
#   # ag.mm9 <- readRDS('/home/steve/ml-data/annotated.genome.mm9.rds')
#   # ## ag.dm3 <- readRDS('/home/steve/ml-data/annotated.genome.dm3.rds')
#   # base <- 'http://cbio.mskcc.org/~lianos/files/bioc2012'
#   #
#   # fetch <- c('Dmel.DEXSeq.exons.rda',
#   #            'NOVA.mm9.rda',  ## nova.peaks
#   #            'pasillaData.rda', ## dm.exons
#   #            'dm3.anno.rds') ## dm3 genome annotation
#   # for (get in fetch) {
#   #   download.file(paste(base, get, sep="/"), get)
#   #   if (length(grep("rds$", get)) == 0) {
#   #     load(get)
#   #   }
#   # }
#   # ag.dm3 <- readRDS('dm3.anno.rds')
#   #
#   # dex <- system.file("data", "ps.diff-exons.txt", package="BiocSeqSVM")
# }
#
