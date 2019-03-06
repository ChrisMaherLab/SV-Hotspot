#!/usr/bin/Rscript

# Summarize sample count per sliding window from bed overlap
# Created by: Ha X. Dang <haxdang@gmail.com>
# Modified by: Abdallah Eteleeb <eteleeb@gmail.com>

args = commandArgs(T)

chr.file = args[1]
out.dir = args[2]

chr.name = gsub(".bed","", chr.file)
cat (paste("\nSummarizing sample counts for chromosme",chr.name, "\n"))

#read window vs break point overlap 
cat('Reading window data...\n')
x = read.table(paste0(out.dir,'/processed_data/segments_with_bps_per_chr/',chr.file), header=F, sep='\t', stringsAsFactors=F)
x = x[, c(1,2,3,7)]
colnames(x) = c('chr', 'start', 'stop', 'sample')
x$svtype = gsub('^.*/', '', x$sample)
x$sample = gsub('/.*', '', x$sample)
svtypes = setdiff(unique(x$svtype), '.')

### extract the total number of samples 
bp = read.table(paste0(out.dir,'/processed_data/all_bp.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(bp) = c('chr', 'start', 'stop', 'name', 'score', 'strand')
bp$sample = gsub('/.*$', '', bp$name)
bp$svtype = gsub('^.*/', '', bp$name)
### extract total number of samples 
total.samples <- length(unique(bp$sample))

#### collect all windows (should have everything including windows w/o break points)
w = x[, c('chr', 'start', 'stop')]
w = w[!duplicated(w),]
rownames(w) = paste(w$chr, w$start, w$stop, sep='_')
cat('Found total', nrow(w), 'unique windows\n')

# count number of samples having at least one break point in each window
countSamples <- function(x, svtype='ALL'){
    cat('Counting samples for', svtype, ' structural variants ...')
    if (svtype != 'ALL'){
        z = x[x$svtype == svtype,]
    }else{
        z = x
    }
    z = z[, 1:4]
    
    # make sure each sample is reported only once for each window
    z = z[!duplicated(z),]

    z1 = aggregate(sample ~ chr+start+stop, data=z, FUN=paste, collapse=',')
    z2 = aggregate(sample ~ chr+start+stop, data=z, FUN=length)
    colnames(z2)[4] = 'num.samples'
    zz = merge(z2, z1)
    zz$num.samples[zz$sample == '.'] = 0
    if (svtype != 'ALL'){
        missing.windows = w[!(rownames(w) %in% paste(zz$chr, zz$start, zz$stop, sep='_')),]
        missing.windows$num.samples = as.integer(0)
        missing.windows$sample = '.'
        zz = rbind(zz, missing.windows)
    }

    zz$pct.samples = zz$num.samples/total.samples*100
    zz$pos = (zz$start + zz$stop)/2
    zz$svtype = svtype
    cat(nrow(zz), 'windows reported (including windows w/o count).\n')
    if(nrow(zz) != nrow(w)){stop('ERR: something went wrong, number of reported windows != total windows.\n')}
    return(zz)
}

# summarize by types separately, then combine. This is slow!
cat('Summarizing...\n')
cts = countSamples(x)
for (svtype in svtypes){
     ### run the counting function for current sv type
     cts = rbind(cts, countSamples(x, svtype))
 }

### write final counts 
saveRDS(cts, file=paste0(out.dir, '/processed_data/counts/', chr.name, '.counts.rds'), compress=FALSE)
#saveRDS(all.cts, file=file.path(out.dir, 'counts.rds'), compress=FALSE)


