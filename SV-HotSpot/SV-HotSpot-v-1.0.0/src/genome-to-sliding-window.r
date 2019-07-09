#!/usr/local/bin/Rscript

# Read chr sizes and generate sliding windows in bed format
# Created by: Ha X. Dang <haxdang attt gmail dottt com>

args = commandArgs(T)

chr.size.file = args[1]
win.size = as.integer(args[2])
step.size = as.integer(args[3])
out.dir = args[4]
chr.of.int = args[5]

### read all breakpoints
bps = read.table(paste0(out.dir,"/processed_data/all_bp.bed"), header =F, stringsAsFactors=F, sep='\t')
colnames(bps) = c('chr','start','stop','name', 'score', 'strand')

### read chromosomes size file 
chr = read.table(chr.size.file, header=T, stringsAsFactors=F, sep='\t')

### extract chromsome names
if (chr.of.int == "ALL") {
    chr.names = chr$chrom
    chr = chr[chr$chrom %in% unique(bps$chr) & chr$chr %in% chr.names,]
    chr.names.found = chr$chr
    cat('Chromosomes to segment:\n')
    print (chr.names.found)
} else {
    chr.names.found = unlist(strsplit(chr.of.int,  ",")) 
    cat('Chromosomes to segment:\n')
    print (chr.names.found)
}


bb = NULL
for (c in chr.names.found){
    cat(c, '\n')
    n = chr$size[chr$chr == c]
    b = data.frame(chr=c, start=seq(0,n,step.size))
    b$stop = b$start + win.size
    b$stop[nrow(b)] = n
    if (is.null(bb)){bb = b}else{bb = rbind(bb, b)}
}

bb$start = as.integer(bb$start)
bb$stop = as.integer(bb$stop)

write.table(bb, file=paste0(out.dir,'/processed_data/genome.segments.bed'), row.names=F, col.names=F, quote=F, sep='\t')

