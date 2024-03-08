#Creating Manhattan Plot
#Reference: https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotManhattan/PlotManhattan.html


library(regioneR)
set.seed(123456)

createDataset <- function(num.snps=20000, max.peaks=5) {
  hg19.genome <- filterChromosomes(getGenome("hg19"))
  snps <- sort(createRandomRegions(nregions=num.snps, length.mean=1, 
                                   length.sd=0, genome=filterChromosomes(getGenome("hg19"))))
  names(snps) <- paste0("rs", seq_len(num.snps))
  snps$pval <- rnorm(n = num.snps, mean = 0.5, sd = 1)
  snps$pval[snps$pval<0] <- -1*snps$pval[snps$pval<0]
  #define the "significant peaks"
  peaks <- createRandomRegions(runif(1, 1, max.peaks), 8e6, 4e6)
  peaks
  for(npeak in seq_along(peaks)) {
    snps.in.peak <- which(overlapsAny(snps, peaks[npeak]))
    snps$pval[snps.in.peak] <- runif(n = length(snps.in.peak), 
                                     min=0.1, max=runif(1,6,8))
  }
  snps$pval <- 10^(-1*snps$pval)
  return(list(peaks=peaks, snps=snps))
}
ds <- createDataset()
ds$snps
 
a


library(karyoploteR)

kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = ds$peaks, points.cex = 0.8)



chr1.snps <- ds$snps[seqnames(ds$snps)=="chr1"]
chr1.top.snp <- which.min(chr1.snps$pval)
chr1.top.snp

kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = chr1.top.snp)


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set1", r0=autotrack(1,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "2blues", r0=autotrack(2,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "greengray", r0=autotrack(3,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "rainbow", r0=autotrack(4,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = c("orchid", "gold", "orange"), r0=autotrack(5,5))

transf.pval <- -log10(ds$snps$pval)
points.col <- colByValue(transf.pval, colors=c("#BBBBBB00", "orange"))
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = points.col)


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, 
                      highlight = ds$peaks, highlight.col = "orchid",
                      suggestive.col="orange", suggestive.lwd = 3,
                      genomewide.col = "red", genomewide.lwd = 6)


transf.pval <- -log10(ds$snps$pval)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, pval = transf.pval, logp = FALSE )

kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)
kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax)


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)
ymax <- kp$latest.plot$computed.values$ymax
ticks <- c(0, seq_len(floor(ymax)))
kpAxis(kp, ymin=0, ymax=ymax, tick.pos = ticks)

kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)


snps <- kp$latest.plot$computed.values$data
suggestive.thr <- kp$latest.plot$computed.values$suggestiveline
#Get the names of the top SNP per chr
top.snps <- tapply(seq_along(snps), seqnames(snps), function(x) {
  in.chr <- snps[x]
  top.snp <- in.chr[which.max(in.chr$y)]
  return(names(top.snp))
})
#Filter by suggestive line
top.snps <- top.snps[snps[top.snps]$y>suggestive.thr]
#And select all snp information based on the names
top.snps <- snps[top.snps]

top.snps



kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)

kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3)


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)

kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=4, cex=1.6, col="red")
kpPoints(kp, data = top.snps, pch=1, cex=1.6, col="red", lwd=2, ymax=10)




kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)

kpPlotMarkers(kp, data=top.snps, labels=names(top.snps), srt=45, y=0.9, 
              ymax=10, r0=0.8, line.color="red")
kpSegments(kp, data=top.snps, y0=top.snps$y, y1=8, ymax=10, col="red")


reg <- extendRegions(ds$peaks, 15e6, 15e6)
kp <- plotKaryotype(plot.type=4)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, r0=autotrack(1,4))
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = ds$peaks, r0=autotrack(1,4))
kpAddLabels(kp, labels = "Trait 2", srt=90, pos=3, r0=autotrack(2,4))
kp <- kpPlotManhattan(kp, data=createDataset()$snps,  r0=autotrack(2,4))
kpAddLabels(kp, labels = "Trait 3", srt=90, pos=3, r0=autotrack(3,4))
kp <- kpPlotManhattan(kp, data=createDataset()$snps,  r0=autotrack(3,4))
kpAddLabels(kp, labels = "Trait 4", srt=90, pos=3, r0=autotrack(4,4))
kp <- kpPlotManhattan(kp, data=createDataset()$snps,  r0=autotrack(4,4))

kpRect(kp, data=reg, y0=0, y1=1, col=NA, border="red", lwd=3)



kp <- plotKaryotype(plot.type=4)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10, r0=0.5)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = ds$peaks, r0=0.5, r1=1, ymax=10)
kpAddLabels(kp, labels = "Trait 2", srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10, r0=0.5, r1=0, tick.pos = c(5,10))
kp <- kpPlotManhattan(kp, data=createDataset()$snps, r0=0.5, r1=0, ymax=10, points.col = "2blues")

genes <- createRandomRegions(nregions = 5, genome = hg19.genome)
kp <- plotKaryotype(plot.type=4)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set3", ymax=10)

kpPlotMarkers(kp, data = genes, labels = paste0("Gene", seq_along(genes)), cex=1.8, y=0.85)


kp <- plotKaryotype(plot.type=4)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10, r1=0.8)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set3", ymax=10, r1=0.8)

kpPlotDensity(kp, data=snps, col="lightskyblue1", r0=0.82, window.size = 10e6)
kpAddLabels(kp, labels = "Density", srt=90, pos=3, cex=1.2, label.margin = 0.025, r0=0.82)


kp <- plotKaryotype(plot.type=1)
kpAxis(kp, ymin=0, ymax=10)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set3", ymax=10)
kpPlotMarkers(kp, data = genes, labels = paste0("Gene", seq_along(genes)), 
              text.orientation = "horizontal", cex=1.2, y=0.85, pos=4)


## Zooming

kp <- plotKaryotype(plot.type=4, zoom="chr1:50e6-90e6")
kpAddBaseNumbers(kp, add.units = TRUE, cex=1)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10, r1=0.8)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set3", ymax=10, r1=0.8)

kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3, cex=1.6, col="red", r1=0.8)
kpPoints(kp, data = top.snps, pch=1, cex=1.6, col="red", lwd=2, ymax=10, r1=0.8)

kpPlotDensity(kp, data=snps, col="lightskyblue1", r0=0.82, window.size = 10e6)
kpAddLabels(kp, labels = "Density", srt=90, pos=3, cex=1.2, label.margin = 0.025, r0=0.82)

 
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")


library(TxDb.Hsapiens.UCSC.hg19.knownGene)

kp <- plotKaryotype(plot.type=4, zoom="chr1:68e6-73e6")
kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e6)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10, r0=0.2)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set3", ymax=10, r0=0.2)

kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3, cex=1.6, col="red", r0=0.2)
kpPoints(kp, data = top.snps, pch=1, cex=1.6, col="red", lwd=2, ymax=10, r0=0.2)

genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8, gene.name.position = "left")

