#!/usr/bin/Rscript

library(RCircos)
library(ggplot2)

usage = '

Plot whole genome circos and individual chromsome plots of sample counts
of SVs (all and individual)

Rscript plot_whole_genome.r
    <SV-Hotspot_result_dir>
    <annotation_dir>
    <genome_assembly_version, eg. hg38>
    <plot_circos, eg. TRUE, FALSE>
    <chromosomes_to_plot, eg. "ALL", "chr1,chrX"> 
    <genes_to_show, eg. "ERG,PTEN,ETV1">
    <color_genes_by_association_direction_with_SV, eg. "TRUE", "FALSE">
    <output_figure_format, eg. "png", "pdf">
    <output_figure_dir>

    eg. Rscript ../plot_whole_genome.r \
          ../output/100-30-50-prank/sv-hotspot-output \
          ../annotations hg38 TRUE chr1,chrX,chr21 ERG,AR,PTEN TRUE png out

'


args = commandArgs(T)

if (length(args) != 9){stop(usage)}

svh.res.dir = args[1]
annot.dir = args[2]
genome = args[3]
plotCircos = args[4] == 'TRUE'
chroms.to.plot = unique(unlist(strsplit(args[5], ',')))
if (chroms.to.plot[1] == 'ALL'){chroms.to.plot = paste0('chr', c(1:22,'X','Y'))}
genes.to.show = unique(unlist(strsplit(args[6], ',')))
colorGeneByDE = args[7] == 'TRUE'
out.format = args[8]
out.dir = args[9]

# reprocess window count data for more efficient plotting/figure
# by creating bigger windows, and grouping smaller windows within
# a bigger window and use the max value of the smaller windows
# for the bigger window.
rebinWindowCount <- function(x, num.win.per.chr=100, win=NULL,
                             chr.col='chr', start.col='start',
                             stop.col='stop', val.col=NULL){
  #if(length(unique(x[[chr.col]])) != 1){return(x)}
  chrs = unique(x[[chr.col]])
  if (is.null(val.col)){val.col = colnames(x)[4]}
  x$value = x[[val.col]]
  x$chr = x[[chr.col]]
  y = NULL
  for (chrom in chrs){
    z = x[x[[chr.col]] == chrom,]
    z = z[order(z[[start.col]]),]

    # create bigger bins/windows
    xmin = min(z[[start.col]])
    xmax = max(z[[stop.col]])
    if(is.null(win)){
        win = floor((xmax - xmin)/num.win.per.chr)
    }
    ww = seq(xmin, xmax, win)
    if (ww[length(ww)] < xmax){ww = c(ww,xmax)} else {ww[length(ww)] = xmax}

    # group smaller original windows completely contained within
    # a bigger windows, and take the max value of those smaller windows
    wstarts = ww[-length(ww)]
    wstops = ww[-1]-1
    z$wstart = z[[start.col]]
    z$wstop = z[[stop.col]]
    for (i in 1:length(wstarts)){
        sel = z[[start.col]] >= wstarts[i] & z[[stop.col]] <= wstops[i]
        z$wstart[sel] = wstarts[i]
        z$wstop[sel] = wstops[i]
    }
    # keep max value
    chp = aggregate(value ~ chr + wstart + wstop, data=z, FUN=max)
    colnames(chp) = c('chr', 'start', 'stop', val.col)
    chp$pos = (chp$start + chp$stop)/2
    if (is.null(y)){y = chp}else{y=rbind(y,chp)}
  }
  return(y)
}

cat('Preparing data...\n')

annot = read.table(file.path(annot.dir, 'hg38/genes.bed'),
                   header=T, stringsAsFactors=F, sep='\t')
colnames(annot) = c('chr', 'start', 'stop', 'gene', 'score', 'strand')

# add AR enhancer to genome annotation
if (genome == 'hg38'){
    ehc = data.frame(chr='chrX', start=66896703, stop=66927905,
                 gene='AR_enhcr', strand='.', score=1000, stringsAsFactors=F)
    annot = rbind(annot, ehc)
}


built.in.genomes = c('hg18', 'hg19', 'hg38', 'mm9', 'mm10', 'dm3', 'dm6', 'rn4', 'rn5','rn6')
if (!(genome %in% built.in.genomes)){
    stop(paste0('ERROR: genome', genome, ' not supported. Quit!\n'))
}


#Load the genome cytoband data
if (genome == "hg38") {
  data(UCSC.HG38.Human.CytoBandIdeogram)
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
} else if (genome == "hg19") {
  data(UCSC.HG19.Human.CytoBandIdeogram)
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
} else if (genome == "mm10") {
  data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
  cyto.info <- UCSC.Mouse.GRCm38.CytoBandIdeogram
} else if (genome == "rn4") {
  data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram)
  cyto.info <- UCSC.Baylor.3.4.Rat.cytoBandIdeogram
} else {
  file.path= file.path(annot.dir, "cytoband_info",
                       paste0(genome, ".cytoBand.txt.gz"))
  if (file.exists(file.path)) {
    cyto.info <- read.table(file.path,
                            colClasses = c("character", "numeric", "numeric",
                                           "character", "character"),
                            sep = "\t", stringsAsFactors = FALSE)  
    colnames(cyto.info) = c('Chromosome','chromStart','chromEnd',
                            'Name','Stain')
  } else {
    warning(paste0("Circos plot failed.
                   Cytoband information not available for ", genome,".\n"))
  }
  
}

# read all counts for SV types
cnt = readRDS(paste0(svh.res.dir, '/processed_data/counts.rds'))

# make sure count windows are within genome coordinate
chrlen = aggregate(chromEnd ~ Chromosome, data=cyto.info, FUN=max)
colnames(chrlen) = c('chr', 'chr.len')
cnt = merge(cnt, chrlen)
long = cnt$stop > cnt$chr.len
cnt$stop[long] = cnt$chr.len[long]

# read peaks & gene associations
pkall = read.table(file.path(svh.res.dir,'annotated_peaks_summary.filtered.tsv'),
                   header=T, stringsAsFactors=F)
pk = pkall
gg = read.table(paste0(svh.res.dir, '/genes.associated.with.SVs.filtered.tsv'),
                header=T, stringsAsFactors=F, sep='\t')
gg = gg[gg$Gene %in% genes.to.show, c('Gene', 'Rank', 'Peak.name',
                                      'Dominant.DE.Direction.WithSV')]
colnames(gg)[1] = 'gene'
gg = merge(annot[annot$gene %in% genes.to.show,], gg, all.x=T)
gg$PlotColor = 'black'
if (colorGeneByDE){
    gg$PlotColor[which(gg$Dominant.DE.Direction.WithSV == 'Upregulated')] = 'red'
    gg$PlotColor[which(gg$Dominant.DE.Direction.WithSV == 'Downregulated')] = 'blue'
}
#gg$PlotColor[gg$gene %in% c('TMPRSS2')] = 'black'
#gg$PlotColor[gg$gene %in% c('ERG')] = 'red'
#gg$grp = NA
#gg$grp[gg$PlotColor == 'red'] = 1
#gg$grp[gg$PlotColor == 'blue'] = 2
#gg$grp[gg$gene %in% c('ERG', 'TMPRSS2')] = 3

gg = gg[!is.na(gg$Peak.name) | gg$gene == 'AR_enhcr',]
cols = c('chr', 'start', 'stop', 'gene')
gg = gg[, c(cols, setdiff(colnames(gg), cols))]


pk = pk[pk$Peak.name %in% gg$Peak.name,]
pk = pk[, 1:4]
colnames(pk) = c('peak', 'chr', 'start', 'stop')

# mark Windows overlapping with peaks for coloring when plot
# currently not used because these look a bit ugly due to
# windows are overlapping
markWinOverlapPeaks <- function(x, pk){
  x$overlap = F
  for (i in 1:nrow(pk)){
    pki = pk[i,]
    sel = x$chr == pki$chr & (x$start >= pki$start & x$start <= pki$stop |
                 x$stop >= pki$start & x$stop <= pki$stop |
                 x$start <= pki$start & x$stop >= pki$stop)
    if(any(sel)){ x$overlap[sel] = T }
  }
  return(x)
}
#z = markWinOverlapPeaks(svall, pk)

# Prepare some colors for SV types
colrs=c('All SVs'='deeppink', 'DUP'='#b53f4d', 'DEL'='#2c7fb8', 'BND'='#2ca25f','INS'='#fec44f', 'INV'='#c994c7')

### extract data
svall = cnt[cnt$svtype == 'ALL', c('chr', 'start', 'stop', 'pct.samples')]
svall = rebinWindowCount(svall)
svall = markWinOverlapPeaks(svall, pk)
#svall = markWinOverlapPeaks(svall, svall[svall$overlap,])
svall$PlotColor <- colrs['All SVs']
#svall$PlotColor[svall$overlap] = 'red'
bnd = cnt[cnt$svtype == 'BND', c('chr', 'start', 'stop', 'pct.samples')]
bnd = rebinWindowCount(bnd)
bnd$PlotColor = colrs['BND']
dup = cnt[cnt$svtype == 'DUP', c('chr', 'start', 'stop', 'pct.samples')]
dup = rebinWindowCount(dup)
dup$PlotColor = colrs['DUP']
del = cnt[cnt$svtype == 'DEL', c('chr', 'start', 'stop', 'pct.samples')]
del = rebinWindowCount(del)
del$PlotColor = colrs['DEL']
inv = cnt[cnt$svtype == 'INV', c('chr', 'start', 'stop', 'pct.samples')]
inv = rebinWindowCount(inv)
inv$PlotColor = colrs['INV']
ins = cnt[cnt$svtype == 'INS', c('chr', 'start', 'stop', 'pct.samples')]
ins = rebinWindowCount(ins)
ins$PlotColor = colrs['INS']

max.val = max(svall$pct.samples)

if (plotCircos){

    cat('Plotting genome circos plot...\n')
      
    RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,
                                tracks.inside=5, tracks.outside=4)
    ##### modify plot parameters 
    para <- RCircos.Get.Plot.Parameters()
    para$text.size = 0.7
    para$track.height <- 0.175
    para$track.background = 'white'
    para$chrom.width = 0.05
    para$grid.line.color = 'white'
    #para$plot.radius = 3
    #para$track.out.start = 1.1
    #para$chr.name.po = 0.8
    #para$chr.ideo.pos = 0.9
    RCircos.Reset.Plot.Parameters(para)

    dir.create(out.dir, recursive=T)
    if (out.format == 'pdf'){
        pdf(file.path(out.dir, "circos.pdf"),  height=7, width=7)
    }else{
        png(file.path(out.dir, "circos.png"),  height=7, width=7,
            unit='in', res=600)
    }
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()

    # draw tracks of sample counts
    RCircos.Histogram.Plot(svall, data.col=4, track.num=1, "out", min.value=0,
                           max.value=max.val, inside.pos=1.9,outside.pos=2.3)
    RCircos.Histogram.Plot(dup, data.col=4, track.num=1, "in",
                           min.value=0, max.value=max.val)
    RCircos.Histogram.Plot(del, data.col=4, track.num=2, "in",
                           min.value=0, max.value=max.val)
    RCircos.Histogram.Plot(bnd, data.col=4, track.num=3, "in",
                           min.value=0, max.value=max.val)
    RCircos.Histogram.Plot(ins, data.col=4, track.num=4, "in",
                           min.value=0, max.value=max.val)
    RCircos.Histogram.Plot(inv, data.col=4, track.num=5, "in",
                           min.value=0, max.value=max.val)

    # draw genes of interest names, RCircos gene coloring is unstable if
    # multiple colors are mixed, hence, draw genes in batches of diff. colors
    # note: close genes with different colors won't be drawn nicely
    for (colr in unique(gg$PlotColor)){
        RCircos.Gene.Connector.Plot(gg[gg$PlotColor == colr,], track.num=2,
                                   side="out", inside.pos=2.3, outside.pos=2.4)
        RCircos.Gene.Name.Plot(gg[gg$PlotColor == colr,], name.col=4,
                     track.num=3, side="out", inside.pos=2.41, outside.pos=2.7)
    }
    legend(1.9,2.4, legend=names(colrs), fill=colrs,  horiz = F, bty="n",
             border="white", cex=0.8, x.intersp=0.5)
    mtext(paste0('max = ',round(max.val),"%"), at=0, cex=0.7, padj=40)
    dev.off()
}

# plot sample counts for all and individual SV types per chromosome
# plotSVSampleCount(chrom='chrX', countData=cnt, peakData=pkall,
#    gene=gg, out.dir=out.dir)
plotSVSampleCount <- function(chroms, countData, peakData, genes=NULL,
                              out.format='png', png.dpi=300, out.dir){
    for (chrom in chroms){
        cat('Plotting ', chrom, '...\n')
        # subsetting count and peak data for chrom to plot
        # count for individual SV types
        #z = countData[countData$chr == chrom & countData$svtype != 'ALL',]
        #z$svtype = factor(z$svtype, levels=rev(c('INV', 'INS', 'BND', 'DEL', 'DUP')))
        # count for all SV types
        #za = cnt[cnt$svtype == 'ALL' & cnt$chr == chrom,]
        #ymax = max(za$pct.samples)
        # peaks
        pk = peakData[peakData$Chr == chrom,]
        pk$pos = (pk$Start + pk$End)/2
        pk$pct.samples = pk$Percentage.SV.samples
        # genes to show
        g = NULL
        if (!is.null(genes)){
            g = genes[genes$chr == chrom,]
            g$pos = (g$start + g$stop)/2
        }

        za = countData[countData$chr == chrom & countData$pct.samples > 0,]
        za$pct.samples[za$svtype != 'ALL'] = -za$pct.samples[za$svtype != 'ALL']
        za$svtype[za$svtype == 'ALL'] = 'All SVs'
        za$svtype = factor(za$svtype,
                      levels=c('All SVs', 'DUP', 'DEL', 'BND', 'INS', 'INV'))
        pchr = (ggplot(za) + geom_bar(aes(x=pos, y=pct.samples, fill=svtype),
                                      stat='identity', size=1)
          + geom_segment(data=pk, aes(x=Start, xend=End, y=Percentage.SV.samples,
                     yend=Percentage.SV.samples), color='black', size=0.5)
          + geom_point(data=pk, aes(x=pos,y=Percentage.SV.samples), color='black',
                                    shape=4, size=0.25)
          + theme_bw()
          + theme(panel.grid=element_blank(),
                  legend.key.size = unit(2.5,'mm'),
                  legend.title=element_blank(),
                  legend.text=element_text(size=9),
                  #axis.text.x=element_blank(),
                  #legend.position='none'
          )
          #+ ggtitle(paste('Percentage of samples across', chrom, 'segments (windows)'))
          + scale_fill_manual(name="SV Type:", values=colrs)
          + scale_color_manual(values=c('blue', 'red'))
          + scale_x_continuous(labels = scales::comma)
          + scale_y_continuous(labels=abs)
          + xlab(paste0(chrom, ' position'))
          + ylab('Pct. samples')
        )

        out.file = file.path(out.dir, paste0(chrom, '.', out.format))
        if (out.format == 'png'){
            ggsave(pchr, file=out.file, width=7, height=2, dpi=png.dpi)
        }else{
            ggsave(pchr, file=out.file, width=7, height=2, useDingbats=F, title='')
        }
    }
}

cat('Plotting individual chromosomes...\n')
plotSVSampleCount(chrom=chroms.to.plot, countData=cnt, peakData=pkall,
                  genes=gg, out.dir=out.dir)

