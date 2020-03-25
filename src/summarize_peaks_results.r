#!/usr/bin/Rscript

# Summarize results for all peaks 
# Created by: Abdallah Eteleeb <eteleeb@gmial.com>

library(plyr)

args = commandArgs(T)

out.dir = args[1]
reg.of.int = args[2]
roi.lbl = args[3]
exp.file = args[4]

### set the region of interest parameters
is.roi.avail = FALSE
if(reg.of.int !=0 ) { 
   is.roi.avail = TRUE 
   num.of.roi.files =  length(unlist(strsplit(reg.of.int,  ",")))
   #if (roi.lbl == 0) {
   #    roi.lbl  = gsub("\\..*", "", basename(reg.of.int))
   #}
}

####################################################################################################################
##### function to compute percentage of sv types 
all.bp <- read.table(paste0(out.dir,'/processed_data/peaks_overlap_bp.tsv'), header =F, stringsAsFactors = F)
colnames(all.bp) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample', 'sv.type')
all.bp <- unique(all.bp)
sv.type.counts <- count(all.bp, c('p.name', 'sv.type'))

pct.sv.fun <- function (pk,w,num.sam){
   p.sv <- sv.type.counts[sv.type.counts$p.name==pk, ]
   dom.svtype <- p.sv[p.sv$freq == max(p.sv$freq),"sv.type"]
   dom.svtype <- paste(dom.svtype, collapse = "|")
   #num.sam <- length(unique(all.bp[all.bp$p.name==pk, 'sample']))
   #density <- num.sam/w
   num.sv.types<- sum(p.sv$freq)
   p.sv$pct.sv <- paste0(p.sv$sv.type, '(', signif((p.sv$freq/num.sam)*100, digits = 4), '%)')
   dd <- data.frame(p.name=pk, pct.sv.types = paste(p.sv$pct.sv, collapse = "|"), dom.svtype)
    
   return (dd)
}
####################################################################################################################

##### read all peaks 
all.peaks <- read.table(paste0(out.dir,'/peaks/all_peaks.bed'), header =F, sep="\t", stringsAsFactors=F)
colnames(all.peaks) <- c('p.chr', 'p.start', 'p.stop','p.name','p.id', 'num.samples', 'pct.samples', 'sample')

##### read peaks overlap/nearby genes  
if (file.exists(paste0(out.dir,'/processed_data/peaks_with_overlap_nearby_genes.tsv'))) {
   p.with.genes <- read.table(paste0(out.dir,'/processed_data/peaks_with_overlap_nearby_genes.tsv'), header =F, sep="\t", stringsAsFactors=F)
   colnames(p.with.genes) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample',
                              'g.chr', 'g.start', 'g.stop', 'gene', 'g.score', 'g.strand', 'dist','g.pos')
} else {
  stop ('peaks_with_overlap_nearby_genes.tsv file was not found!. Make sure all previous steps were run correctly.')
}

##### read peaks overlap/nearby region of interest
if (file.exists(paste0(out.dir,'/processed_data/peaks_overlap_with_region_of_interest.tsv'))) {
  #p.with.roi <- read.table(paste0(out.dir,'/processed_data/peaks_with_overlap_nearby_region_of_interest.tsv'), header =F, sep="\t", stringsAsFactors=F)
  p.with.roi <- read.table(paste0(out.dir,'/processed_data/peaks_overlap_with_region_of_interest.tsv'), header =T, sep="\t", stringsAsFactors=F)
  #colnames(p.with.roi) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample',
  #                         'roi.chr', 'roi.start', 'roi.stop', 'roi.name', 'dist')
  ### extract region of interest name(s)
  roi.name.cols = colnames(p.with.roi)[grepl("name", colnames(p.with.roi)) & colnames(p.with.roi) !="p.name"]  
}
 

annot.peaks <- NULL
for (i in 1:nrow(all.peaks)) {
  
  pk <- all.peaks$p.name[i]
  p.res <- all.peaks[all.peaks$p.name==pk, ]
  p.res$p.id <- NULL
  p.width <- abs(p.res$p.stop- p.res$p.start)
  
  #### compute percentage of sv types 
  pct.sv <- pct.sv.fun(pk, p.width, p.res$num.samples) 
  p.res <- merge(p.res, pct.sv)
  p.res <- p.res[, c("p.name", "p.chr", "p.start", "p.stop", "num.samples", "pct.samples", "sample", "pct.sv.types", "dom.svtype")]

  ### extract overlapped/nearby genes 
  ov.genes <- unique(p.with.genes[p.with.genes$p.name==pk & p.with.genes$g.pos=="overlap", 'gene'])
  nearby.genes <- unique(p.with.genes[p.with.genes$p.name==pk & p.with.genes$g.pos=="nearby", 'gene'])
  
  ### extract overlapped region of interest
  if (is.roi.avail) {
     ov.roi = NULL
     for (k in 1:length(roi.name.cols)) { 
          c.roi = unique(p.with.roi[p.with.roi$p.name==pk, roi.name.cols[k]])
          c.roi = c.roi[!is.na(c.roi)] 
          c.roi = paste(c.roi, collapse="|")
          
          if (k == 1 ) { 
             ov.roi = as.data.frame(c.roi)
          } else {
             ov.roi = as.data.frame(cbind(ov.roi , c.roi))
          }
     }
     colnames(ov.roi) <- gsub(".name", "", roi.name.cols)

    ### combine all 
    d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"), ov.roi)

  } else {
    d <- data.frame(p.res, overlap.genes=paste(ov.genes, collapse = "|"), nearby.genes=paste(nearby.genes, collapse = "|"))
  }

  annot.peaks <- rbind(annot.peaks, d)
  
} ### end of peaks

### set column names 
if (is.roi.avail) {
  colnames(annot.peaks) <- c('Peak.name','Chr','Start','End','Number.SV.samples','Percentage.SV.samples','SV.sample','Percentage.SV.types',
                             'Dominant.svtype','Overlapped.genes','Nearby.genes',paste0('Overlapped.',gsub(".name", "", roi.name.cols)))
} else {
  colnames(annot.peaks) <- c('Peak.name','Chr','Start','End','Number.SV.samples','Percentage.SV.samples','SV.sample','Percentage.SV.types',
                             'Dominant.svtype','Overlapped.genes','Nearby.genes')
}

### write final results 
annot.peaks[annot.peaks==""] <- NA
annot.peaks <- annot.peaks[order(annot.peaks$Percentage.SV.samples, decreasing = T), ] 
if (exp.file != 0) {
  write.table(annot.peaks, file=paste0(out.dir,'/processed_data/annotated_peaks_summary.tsv'), sep="\t", row.names=F, quote = F)
} else {
  write.table(annot.peaks, file=paste0(out.dir,'/annotated_peaks_summary.tsv'), sep="\t", row.names=F, quote = F)
}

cat('Done.\n')

