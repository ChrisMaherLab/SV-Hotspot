#!/gapp/x64linux/opt/R3.1.2/bin/Rscript

# Read chr sizes and generate sliding windows in bed format
# Created by: Ha X. Dang <haxdang attt gmail dottt com>

args = commandArgs(T)
#if (length(args) != 4){
#    stop('<thisfile> <chrom-size.tsv> <window-size-in-bases> <step-in-bases> <out.bed.file>\n')
#}
# D - sliding window size; d - sliding distance

chr.size.file = args[1]
D = as.integer(args[2])
d = as.integer(args[3])
out.dir = args[4]

### read all breakpoints
bps = read.table(paste0(out.dir,"/temp/all_bp.bed"), header =F, stringsAsFactors=F, sep='\t')
colnames(bps) = c('chr','start','stop','name', 'score', 'strand')

chr = read.table(chr.size.file, header=F, stringsAsFactors=F, sep='\t')
colnames(chr) = c('chr', 'size')
chr.names = paste0('chr', c(1:22, 'X', 'Y'))
chr = chr[chr$chr %in% unique(bps$chr) & chr$chr %in% chr.names,]
chr.names.found = chr$chr
cat('Chromosomes to segment:\n')
print (chr.names.found)

bb = NULL
for (c in chr.names.found){
    cat(c, '\n')
    n = chr$size[chr$chr == c]
    b = data.frame(chr=c, start=seq(0,n,d))
    b$stop = b$start + D
    b$stop[nrow(b)] = n
    if (is.null(bb)){bb = b}else{bb = rbind(bb, b)}
}

bb$start = as.integer(bb$start)
bb$stop = as.integer(bb$stop)

write.table(bb, file=paste0(out.dir,'/temp/genome.segments.bed'), row.names=F, col.names=F, quote=F, sep='\t')

