#!/usr/bin/Rscript

# Summarize results for all peaks 
# Created by: Abdallah Eteleeb <eteleeb@gmial.com>

library(plyr)

args = commandArgs(T)
out.dir = args[1]

roi.files = list.files(paste0(out.dir,"/processed_data/ROIs"))

all.roi.files = NULL
for (i in 1:length(roi.files)) {
  f = read.table(paste0(out.dir,"/processed_data/ROIs/",roi.files[i]), header =F, sep="\t", stringsAsFactor = F)
  roi.name = gsub(".*peaks_overlap_with_", "",roi.files[i])
  roi.name =  gsub("\\..*", "", roi.name)
  colnames(f) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample',
                  paste0(roi.name,'.chr'), paste0(roi.name,'.start'), paste0(roi.name,'.stop'), paste0(roi.name,'.name'), paste0(roi.name,'.dist'))
  if (i == 1 ) { 
    all.roi.files = f
  } else {
    all.roi.files = merge(all.roi.files , f, all = T)   
  }

}

### write the results 
write.table(all.roi.files, file=paste0(out.dir,'/processed_data/peaks_overlap_with_region_of_interest.tsv'), sep="\t", row.names=F, quote = F)

### remove ROIs folder 
unlink(paste0(out.dir,"/processed_data/ROIs"), recursive = TRUE)

