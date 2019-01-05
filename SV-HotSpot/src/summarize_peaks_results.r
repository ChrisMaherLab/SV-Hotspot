#!/gapp/x64linux/opt/R3.1.2/bin/Rscript

# Summarize results for all peaks 
# Created by: Abdallah Eteleeb <eteleeb@gmial.com>

library(plyr)

args = commandArgs(T)
out.dir = args[1]
reg.of.int = args[2]

is.roi.avail = FALSE
if(reg.of.int !=0 ) { is.roi.avail = TRUE }

####################################################################################################################
##### function to compute percentage of sv types 
all.bp <- read.table(paste0(out.dir,'/processed_data/peaks_overlap_bp.tsv'), header =F, stringsAsFactors = F)
colnames(all.bp) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample', 'sv.type')
all.bp <- unique(all.bp)
sv.type.counts <- count(all.bp, c('p.name', 'sv.type'))

pct.sv.fun <- function (pk, w){
   p.sv <- sv.type.counts[sv.type.counts$p.name==pk, ]
   dom.svtype <- p.sv[p.sv$freq == max(p.sv$freq),"sv.type"]
   dom.svtype <- paste(dom.svtype, collapse = "|")
   num.sam <- length(unique(all.bp[all.bp$p.name==pk, 'sample']))
   density <- num.sam/w
   num.sv.types<- sum(p.sv$freq)
   p.sv$pct.sv <- paste0(p.sv$sv.type, '(', signif((p.sv$freq/num.sam)*100, digits = 4), '%)')
   dd <- data.frame(p.name=pk, num.samples = num.sam, density, pct.sv.types = paste(p.sv$pct.sv, collapse = "|"), dom.svtype)
    
   return (dd)
}
####################################################################################################################

##### read all peaks 
all.peaks <- read.table(paste0(out.dir,'/peaks/all_peaks.bed'), header =F, sep="\t", stringsAsFactors=F)
colnames(all.peaks) <- c('p.chr', 'p.start', 'p.stop','p.name','p.id', 'num.samples', 'pct.samples', 'sample')

##### read peaks overlap/nearby genes  
if (file.exists(paste0(out.dir,'/peaks_with_overlap_nearby_genes.tsv'))) {
   p.with.genes <- read.table(paste0(out.dir,'/peaks_with_overlap_nearby_genes.tsv'), header =F, sep="\t", stringsAsFactors=F)
   colnames(p.with.genes) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample',
                              'g.chr', 'g.start', 'g.stop', 'gene', 'g.score', 'g.strand', 'dist','g.pos')
} else {
  stop ('peaks_with_overlap_nearby_genes.tsv file was not found!. Makr sure all previous steps were run correctly.')
}

##### read peaks overlap/nearby region of interest
if (file.exists(paste0(out.dir,'/peaks_overlap_with_region_of_interest.tsv'))) {
  #p.with.roi <- read.table(paste0(out.dir,'/peaks_with_overlap_nearby_region_of_interest.tsv'), header =F, sep="\t", stringsAsFactors=F)
  p.with.roi <- read.table(paste0(out.dir,'/peaks_overlap_with_region_of_interest.tsv'), header =F, sep="\t", stringsAsFactors=F)
  colnames(p.with.roi) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample',
                           'roi.chr', 'roi.start', 'roi.stop', 'roi.name', 'dist')
}

annot.peaks <- NULL
for (i in 1:nrow(all.peaks)) {
  
  pk <- all.peaks$p.name[i]
  p.res <- all.peaks[all.peaks$p.name==pk, ]
  p.res$p.id <- NULL
  p.width <- abs(p.res$p.stop- p.res$p.start)
  
  #### compute percentage of sv types 
  pct.sv <- pct.sv.fun(pk, p.width)  
  p.res <- merge(p.res, pct.sv)
  p.res <- p.res[, c("p.name", "p.chr", "p.start", "p.stop", "num.samples", "pct.samples", "sample", "pct.sv.types", "density", "dom.svtype")]

  ### extract overlapped/nearby genes 
  ov.genes <- unique(p.with.genes[p.with.genes$p.name==pk & p.with.genes$g.pos=="overlap", 'gene'])
  nearby.genes <- unique(p.with.genes[p.with.genes$p.name==pk & p.with.genes$g.pos=="nearby", 'gene'])
  
  ### extract overlapped region of interes 
  if (is.roi.avail) {
    #ov.roi <- unique(p.with.roi[p.with.roi$p.name==pk & p.with.roi$roi.pos=="overlap", 'roi.name'])
    #nearby.roi <- unique(p.with.roi[p.with.roi$p.name==pk & p.with.roi$roi.pos=="nearby", 'roi.name']) 
    ov.roi <- unique(p.with.roi[p.with.roi$p.name==pk, 'roi.name'])

  #  d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"),
  #                overlap.roi=paste(ov.roi, collapse = "|"), nearby.roi=paste(nearby.roi, collapse = "|"))
    d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"), overlap.roi=paste(ov.roi, collapse = "|"))

  } else {
    d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"))
  }

  annot.peaks <- rbind(annot.peaks, d)
  
} ### end of chromosomes 

write.table(annot.peaks, file=paste0(out.dir,'/processed_data/annotated_peaks_summary.tsv'), sep="\t", row.names=F, quote = F)

