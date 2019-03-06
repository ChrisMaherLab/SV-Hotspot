#!#!/usr/bin/Rscript

# Created by: Ha X. Dang <haxdang@gmail.com>
# Modified by: Abdallah Eteleeb <eteleeb@gmail.com>

args = commandArgs(T)

sv.file = args[1]
out.dir = args[2]

### read structural variants file 
if (file.exists((sv.file))) {
  sv <- read.table(sv.file, header =T, sep="\t", stringsAsFactors = F, check.names=F)
  sv$sample = sub('/.*$', '', sv$name)
  sv$svtype = sub('^.*/', '', sv$name)
  sv = sv[sv$svtype %in% c('DUP', 'DEL'),]
  sv$pos1 = floor((sv$start1 + sv$end1)/2)
  sv$pos2 = ceiling((sv$start2 + sv$end2)/2)
  ### extract the final file 
  sv = sv[,c('chrom1', 'pos1', 'pos2', 'name', 'score', 'strand1')]
} else {
  stop('structural variants file was not found!')
}

### remove wide peaks (< = 3M)
#sv$sv.width = abs(sv$pos1 - sv$pos2)
#sv = sv[sv$sv.width <= 10000000, ]
#sv$sv.width = NULL

### write the results to be overlapped with segments 
#sv <- sv[,c('chrom1', 'pos1', 'pos2', 'sample', 'svtype')]
write.table(sv, file=paste0(out.dir, "/processed_data/del_dup_sv.bed"), sep ="\t", quote = F, row.names = F, col.names = F)

