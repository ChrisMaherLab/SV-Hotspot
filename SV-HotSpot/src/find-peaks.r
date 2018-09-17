#!/usr/bin/Rscript


args = commandArgs(T)
#args = c('output/count.rds', 'ALL', 'chr14,chr21,chrX', '100', '5', 'output/peaks')
#args = c('output/count.rds', 'ALL', 'ALL', '100', '5', 'output/peaks')

library(peakPick)
library(ggplot2)
library(gridExtra)

#count.file = args[1]
svtype = args[1]
chrs = args[2]
peakPick.win = as.numeric(args[3])
peakPick.minsd = as.numeric(args[4])
pct.samples.cutoff = as.numeric(args[5])
out.dir = args[6]
dist = as.numeric(args[7])
sv.path=args[8]
genes.of.int <- args[9]

if (file.exists(paste0(out.dir,'/counts.rds')) ){
  cat('Reading sliding window sample count...\n')
  x0 = readRDS(paste0(out.dir,'/counts.rds'))
} else {
  stop ("File \"count.rds\" was not found, make sure this file exists.\n")
}

#chrs = c('chr21', 'chrX', 'chr14')
if (chrs[1] == "ALL"){
    chrs = paste0('chr', c(1:22, 'X', 'Y'))
}
chrs = intersect(chrs, unique(x0$chr))
cat('Chromosomes to analyze:\n')
print(chrs)

# gene info
g = read.table(genes.of.int, header=F, stringsAsFactors=F)
colnames(g) = c('chr', 'start', 'stop', 'gene')
g$pos = (g$start + g$stop)/2
g = g[, c('chr', 'pos', 'gene')]
g = unique(g)

# group peaks nearby (distance < d)
####################### START OF FUNCTION (Group peaks) ######################
groupPeaks <- function(x, d, m){
    x = x[order(x$pos),]
    peaks =  which(x$peak)
    last.peak = peaks[1]
    x$group = NA
    group = 1
    for (pk in peaks){
        if (x$pos[pk] > x$pos[last.peak] + d){
            group = group + 1
        }
        x$group[pk] = group
        last.peak = pk
    }

 num.peak.groups = length(unique(x$group[!is.na(x$group)]))
 cat('Number of peak groups: ', num.peak.groups, '\n')
    
    # produce bed file
    gr = NULL
    if (num.peak.groups > 0){
        z = x[order(x$pos),]
        start = z[!is.na(z$group) & !duplicated(z$group), c('chr', 'start', 'group')]
        z = z[order(z$pos, decreasing=T),]
        stop = z[!is.na(z$group) & !duplicated(z$group), c('chr', 'stop', 'group')]
        gr = merge(start, stop)
        gr = gr[order(gr$group),]
        z = z[order(z$pct.samples, decreasing=T),]
        z = z[!is.na(z$group) & !duplicated(z$group),]
        gr = merge(gr, z[, c('group', 'pct.samples')])
        gr = gr[, c('chr', 'start', 'stop', 'group', 'pct.samples')]
    }

    return(list(x=x, group=gr))
}
####################### END OF FUNCTION (Group peaks) ######################

# call peaks for each chr for chosen svtype and plot results
out.dir <- paste0(out.dir,'/peaks')
dir.create(out.dir)

for (chr in chrs){
    cat(chr, '\n')
    x = x0[x0$chr == chr & x0$svtype == svtype,]
    x = x[order(x$pos),]

    #print ("Identifying peaks based on peak calling method ...")
    mat = as.matrix(data.frame(a=1, b=x$pct.samples))
    # call peakPick
    # peakPick.win = 100; peakPick.minsd = 5
    h <- detect.spikes(mat, roi=c(peakPick.win+1, nrow(x)-peakPick.win-1), winlen=peakPick.win, spike.min.sd=peakPick.minsd)
      
    # post processing peakPick output 
    x$peak = F
    x$peak[which(h[,2])] = T
    #filter peaks with low % of samples
    x$peak[x$pct.samples < pct.samples.cutoff] = F
    
    ## group peaks 
    #print(table(x$peak))
    z = groupPeaks(x, dist)

    write.table(z$group, file=paste0(out.dir, '/', chr, '.peak.group.bed'), sep='\t', quote=F, row.names=F)
    write.table(z$x, file=paste0(out.dir, '/', chr, '.peak.bed'), sep='\t', quote=F, row.names=F)

    g1 = g[g$chr == chr,]
    fake = data.frame(chr=chr, pos=-1000, gene='', stringsAsFactors=F)
    if (nrow(g1) == 0){g1 = fake}
    g1$y = -2

    png(paste0(out.dir, '/', chr, '.png'), units='in', width=10, height=5, res=200)
    p1 = (ggplot(x[x$pct.samples > 0,], aes(x=pos, y=pct.samples))
          + geom_bar(stat='identity', fill='black', size=1)
          + scale_x_continuous(limits=c(0,max(x$pos)))
          + theme_bw()
          + ggtitle('Windows')
    )
    if (is.null(z$group)){#no peak or peak groups
        #p2 = ggplot(NULL) + geom_point()
        p2 = (ggplot(g1)
                + geom_text(data=g1, aes(x=pos, y=y, label=gene), angle=90, hjust=1,color='blue', size=1.5)
                + scale_x_continuous(limits=c(0,max(x$pos)))
                + scale_y_continuous(limits=c(-10, max(x$pct.samples)+10), breaks=seq(0,100,10))
                + theme_bw()
                + ggtitle('Peaks (grouped)')
                + ylab('pct.samples')
        )
    }else{
        p2 = (ggplot(x[x$peak,], aes(x=pos, y=pct.samples))
                + geom_bar(stat='identity', size=0, color='red', fill='red')
                + geom_point(shape=1, size=0.5, alpha=0.5)
                + theme_bw()
                + geom_text(data=g1, aes(x=pos, y=y, label=gene), angle=90, hjust=1,color='blue', size=1.5)
                + geom_segment(data=z$group, aes(x=start, y=(pct.samples+5), xend=stop, yend=(pct.samples+5)), color='blue')
                + geom_text(data=z$group, aes(x=(start+stop)/2, y=(pct.samples+7), label=group), size=2, color='blue', angle=90, hjust=0)
                + scale_x_continuous(limits=c(0,max(x$pos)))
                + scale_y_continuous(limits=c(-10, max(x$pct.samples)+10), breaks=seq(0,100,10))
                + ggtitle('Peaks (grouped)')
        )
    }

    grid.arrange(p1,p2,nrow=2)
    dev.off()

}


