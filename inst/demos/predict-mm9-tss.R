library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicRanges)
library(shikken)

ag.mm9 <- readRDS('/home/steve/ml-data/annotated.genome.mm9.rds')
## ag.mm9 <- readRDS('~/cBio/projects/TagSeq/inst/extdata/annotated.genome.mm9.rds')

chr1.prom <- ag.mm9[seqnames(ag.mm9) == 'chr1' &
                    values(ag.mm9)$exon.anno == 'utr5*']
summary(width(chr1.prom))
## remove outliers
pos <- chr1.prom[width(chr1.prom) > 500]
pos <- resize(pos, width=1, fix='end') + 100

neg <- ag.mm9[seqnames(ag.mm9) == 'chr1' &
              values(ag.mm9)$exon.anno == 'utr3*']

neg <- neg[width(neg) > 500]
neg <- resize(neg, width=1, fix='start') + 100

N <- 500
posdna <- getSeq(Mmusculus, pos, as.character=FALSE)[1:N]
negdna <- getSeq(Mmusculus, neg, as.character=FALSE)[1:N]

X <- c(posdna, negdna)
y <- c(rep(1, length(posdna)), rep(-1, length(negdna)))

m <- SVM(X, y, kernel="spectrum", degree=4, C=100)
## X <- X[-m]
## y <- y[-m]

preds <- predict(m, X)
sum(preds == y) / length(y)

## Qustions:
##   - can you find better values for C & degree?
##   - if we take a better random set for the background, does
##     it make the problem harder?

## Fetch promoters from different chromosome and predict!
pos.test <- ag.mm9[seqnames(ag.mm9) == 'chr17' & values(ag.mm9)$exon.anno == 'utr5*']
summary(width(chr17.pos))
pos.test <- pos.test[width(pos.test) > 500]
pos.test <- resize(pos.test, width=1, fix='end') + 100
pt.dna <- getSeq(Mmusculus, pos.test, as.character=FALSE)
p.preds <- predict(m, pt.dna)
sum(p.preds == 1) / length(pt.dna)

nt.test <- ag.mm9[seqnames(ag.mm9) == 'chr17' & values(ag.mm9)$exon.anno == 'utr3*']
summary(width(nt.test))
nt.test <- resize(nt.test, width=1, fix="start") + 100
nt.dna <- getSeq(Mmusculus, nt.test, as.character=FALSE)
n.preds <- predict(m, nt.dna)
sum(n.preds == -1) / length(nt.dna)

## Question:
##   You see that the accuracy isn't as high as we observed from when we
##   trained our model.
##     - Why?
##     - What should we do about it?
