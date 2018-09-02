#!/gapp/x64linux/opt/R3.1.2/bin/Rscript

# Summarize results for all peaks 
# Created by: Abdallah Eteleeb <eteleeb@gmial.com>

library(plyr)

args = commandArgs(T)
out.dir = args[1]
reg.of.int = args[2]

if(reg.of.int !=0 ) { is.roi.avail = TRUE }

####################################################################################################################
##### function to compute percentage of sv types 
all.bp <- read.table(paste0(out.dir,'/temp/peaks_overlap_bp.tsv'), header =F, stringsAsFactors = F)
colnames(all.bp) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'pct.samples', 'sample', 'sv.type')
all.bp <- unique(all.bp)
sv.type.counts <- count(all.bp, c('p.name', 'sv.type'))
pct.sv.fun <- function (pk){
   p.sv <- sv.type.counts[sv.type.counts$p.name==pk, ]
   num.sam <- length(unique(all.bp[all.bp$p.name==pk, 'sample']))
   #num.sv.types<- sum(p.sv$freq)
   p.sv$pct.sv <- paste0(p.sv$sv.type, '(', signif((p.sv$freq/num.sam)*100, digits = 4), '%)')
   dd <- data.frame(p.name=pk, num.samples = num.sam, pct.sv.types = paste(p.sv$pct.sv, collapse = "|"))
   return (dd)
}
####################################################################################################################
##### read all peaks 
all.peaks <- read.table(paste0(out.dir,'/peaks/all_peaks.bed'), header =F, sep="\t", stringsAsFactors=F)
colnames(all.peaks) <- c('p.chr', 'p.start', 'p.stop','p.name','p.id', 'pct.samples')

##### read peaks overlap/nearby genes  
if (file.exists(paste0(out.dir,'/peaks_with_overlap_nearby_genes.tsv'))) {
   p.with.genes <- read.table(paste0(out.dir,'/peaks_with_overlap_nearby_genes.tsv'), header =F, sep="\t", stringsAsFactors=F)
   colnames(p.with.genes) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'pct.samples','g.chr', 'g.start', 'g.stop', 'gene', 'g.score', 'g.strand', 'dist','g.pos')
} else {
  stop ('peaks_with_overlap_nearby_genes.tsv was not found!. Makr usre you run previous steps correctly.')
}

if (file.exists(paste0(out.dir,'/peaks_overlap_nearby_region_of_interest.tsv'))) {
  p.with.roi <- read.table(paste0(out.dir,'/peaks_overlap_nearby_region_of_interest.tsv'), header =F, sep="\t", stringsAsFactors=F)
  colnames(p.with.roi) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'pct.samples','roi.chr', 'roi.start', 'roi.stop', 'roi.name', 'dist','roi.pos')
}


annot.peaks <- NULL
for (i in 1:nrow(all.peaks)) {
  p <- all.peaks$p.name[i]
  p.res <- all.peaks[ all.peaks$p.name==p, ]
  p.res$p.id <- NULL
  
  #### compute percentage of sv types 
  pct.sv <- pct.sv.fun(p)  
  p.res <- merge(p.res, pct.sv)
  
  ov.genes <- unique(p.with.genes[p.with.genes$p.name==p & p.with.genes$g.pos=="overlap", 'gene'])
  nearby.genes <- unique(p.with.genes[p.with.genes$p.name==p & p.with.genes$g.pos=="nearby", 'gene'])
  
  if (is.roi.avail) {
   ov.roi <- unique(p.with.roi[p.with.roi$p.name==p & p.with.roi$roi.pos=="overlap", 'roi.name'])
   nearby.roi <- unique(p.with.roi[p.with.roi$p.name==p & p.with.roi$roi.pos=="nearby", 'roi.name']) 

  d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"),
                  overlap.roi=paste(ov.roi, collapse = "|"), nearby.roi=paste(nearby.roi, collapse = "|"))
  } else {

  d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"))
  }

  annot.peaks <- rbind(annot.peaks, d)
  
} ### end of chromosomes 

write.table(annot.peaks, file=paste0(out.dir,'/temp/annotated_peaks_summary.tsv'), sep="\t", row.names=F, quote = F)

