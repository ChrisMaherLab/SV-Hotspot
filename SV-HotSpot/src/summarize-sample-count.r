#!/gapp/x64linux/opt/R3.1.2/bin/Rscript

# Summarize sample count per sliding window from bed overlap
# Created by: Ha X. Dang <haxdang@gmail.com>
# Modified by: Abdallah Eteleeb <eteleeb@gmail.com>

args = commandArgs(T)

out.dir = args[1]
plotit = as.logical(args[2])

if (is.na(plotit)) {plotit = F}

library(grid)
#library(gridBase)
library(gridExtra)
library(ggplot2)
library(reshape2)

#read window vs break point overlap 
cat('Reading window data...\n')
x = read.table(paste0(out.dir,'/temp/genome.segments.with.bps.bed'), header=F, sep='\t', stringsAsFactors=F)
x = x[, c(1,2,3,7)]
colnames(x) = c('chr', 'start', 'stop', 'sample')
x$svtype = gsub('^.*/', '', x$sample)
x$sample = gsub('/.*', '', x$sample)
#svtypes = c('BND', 'DEL', 'DUP', 'NS', 'INV')
svtypes = setdiff(unique(x$svtype), '.')
#x = x[x$svtype %in% c(svtypes, '.'),]

samples = setdiff(unique(x$sample), '.')
total.samples = length(samples)

#chrs = c('chr21')
#x = x[x$chr %in% chrs,]

# collect all windows (should have everything including windows w/o break points)
w = x[, c('chr', 'start', 'stop')]
w = w[!duplicated(w),]
rownames(w) = paste(w$chr, w$start, w$stop, sep='_')
cat('Found total', nrow(w), 'unique windows\n')

###
#system(paste0("cut -f1-4 | sort | uniq | groupBy -g 1,2,3 -c 4 -o count > ", out.dir, "/counts.rds")) 
###

# count number of samples having at least one break point in each window
# return data.frame of sample count and details for all windows
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
    if(nrow(zz) != nrow(w)){stop('ERR: something wrong, number of reported windows != total windows.\n')}
    return(zz)
}

# summarize by types separately, then combine. This is slow!
cat('Summarizing...\n')
aaa = countSamples(x)
for (svtype in svtypes){
    aaa = rbind(aaa, countSamples(x, svtype))
}
zz = aaa

### write final counts 
saveRDS(zz, file=file.path(out.dir, 'count.rds'))


