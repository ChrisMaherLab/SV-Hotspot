#!/usr/bin/Rscript

# Find regions/peaks whose SVs altered expression of nearby genes
# Written by Abdallah Eteleeb & Ha Dang



usage = '

Plot peak regions and expression of associated genes

usage : Rscipt plot_peak_region.r
    <csv-of-peaks, eg. pX.55.1,p1.2.3>
    <SV-HotSpot_result_dir>
    <bedpe_file_of_SVs>
    <output_figure_dir>
    <tsv_file_of_exp_matrix>
    <tsv_file_of_cna_segments>
    <bed_file_of_regulatory_or_other_nongene_annotations>
    <bedGraph_file_of_ChIP-Seq_coverage>
    <cna_gain_cutoff>
    <cna_loss_cutoff>
    <ChIP-Seq_label>
    <peak_plot_left_extend_in_bases>
    <peak_plot_right_extend_in_bases>
    <csv_of_genes_to_plot_and_show, if not given or ="-", show all genes in region>
    <format, eg. "wide" = exp. at bottom, if not given, exp on right>

    eg. Rscipt plot_peak_region.r pX.55.1 ./sv-hotspot-output data/sv.bedpe
        outfigs exp.tsv cna.seg.tsv GeneHancers.bed H3K27ac.bg 2.99 1.35
        H3K27ac 500000 500000 - wide


'


library(ggplot2)
library(reshape2)
library(grid)
#library(gridBase)
library(gridExtra)
library(gtable)
library(ggsignif)
library(stringr)

args = commandArgs(T)

#args = unlist(strsplit('pX.55.1 output/100-30-50-prank/sv-hotspot-output data/sv_filtered_2020.bedpe output/100-30-50-prank/figs data/su2c_exp_unique.tsv data/cna.seg.tsv data/GeneHancers.bed data/H3K27ac.chrX.cut.bg 2.99 1.35 H3K27ac 500000 500000 AR,AR_enhancer wide', ' '))

if (length(args) < 13){stop(usage)}

pks.to.plot = args[1]
res.dir = args[2]
sv.file = args[3]
out.dir = args[4]
exp.file = args[5]
cn.file = args[6]
roi.file = args[7]
chip.seq = args[8]
t.amp = as.numeric(args[9])
t.del = as.numeric(args[10])
chip.cov.lbl= args[11]
#roi.lbl = args[12]
left.ext = as.numeric(args[12])
right.ext = as.numeric(args[13])
genes.to.show = unlist(strsplit(args[14], ','))
if (genes.to.show=='none' || genes.to.show == '-'){genes.to.show = NULL}
layout = args[15]
if(is.na(layout)){layout = 'narrow'}


colrs=c('DUP'='#b53f4d', 'DEL'='#2c7fb8', 'BND'='#2ca25f','INS'='#fec44f', 'INV'='#c994c7')
# how much the region plot will expand to right to encompass legends
pctRightExpand = 0.15
if(layout == 'wide'){pctRightExpand = pctRightExpand/2}

#### function to align figures 
AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  plots.grobs <- lapply(list(...), ggplotGrob)
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],unit(abs(diff(c(LegendWidth(x),max.legends.width))), "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

# reprocess chip-seq for more efficient plotting/figure
# overlap with bigger window, and take max of overlapping small bin
# also convert cov to pct for convenient plot and align
binChip <- function(chip, win=NULL, convert.to.pct=F){
    if(length(unique(chip$chrom)) != 1){return(chip)}
    chip = chip[order(chip$pos),]
    xmin = min(chip$start)
    xmax = max(chip$end)
    if(is.null(win)){
        win = floor((xmax - xmin)/500)
    }
    ww = seq(xmin, xmax, win)
    if (ww[length(ww)] < xmax){ww = c(ww,xmax)}else{ww[length(ww)] = xmax}
    wstarts = ww[-length(ww)]
    wstops = ww[-1]-1
    chip$wstart = chip$start; chip$wstop = chip$end
    for (i in 1:length(wstarts)){
        sel = with(chip, start >= wstarts[i] & end <= wstops[i])
        chip$wstart[sel] = wstarts[i]
        chip$wstop[sel] = wstops[i]
    }
    chp = aggregate(cov ~ chrom + wstart + wstop, data=chip, FUN=max)
    colnames(chp) = c('chrom', 'start', 'end', 'cov')
    chp$pos = (chp$start + chp$end)/2
    max.chip.cov = max(chp$cov)
    if (convert.to.pct){
        chp$cov = chp$cov/max.chip.cov*100
    }
    return(list(chip=chp,max.cov=max.chip.cov))
}


# Function to pile up DUP and DEL calls in region (chrom:left-right)
pileUp <- function(x, chrom, left, right){
    # select events that overlap with region to plot (left, right)
    x = x[x$chrom1 == chrom,]
    x = x[x$pos1 <= right & x$pos1 >= left | 
          x$pos2 <= right & x$pos2 >= left |
          x$pos1 <= left & x$pos2 >= right, ]
    # order/pile up
    x$sign = ifelse(x$svtype == 'DUP', 1, 0)
    x = x[order(x$svtype, x$sign*x$pos1),]
    x$samp = paste0(x$svtype, '/', x$sample, '/', x$pos1)
    x$samp = factor(x$samp, levels=unique(x$samp))
    return(x)
}

# Function to pile up DUP, DEL, INS, INV, BND  calls in region (chrom:left-right)
pileUp2 <- function(x, chrom, left, right, svtypes=c('DUP','DEL')){
    # select events that overlap with region to plot (left, right)
    x = x[x$chrom1 == chrom | x$chrom2 == chrom,]
    #print(c(chrom, left, right))
    # DUP, DEL, INS, INV (same chromsome breakpoints)
    sel = x$chrom1 == chrom & x$chrom2 == chrom &
           (x$pos1 <= right & x$pos1 >= left | 
           x$pos2 <= right & x$pos2 >= left |
           x$pos1 <= left & x$pos2 >= right)
    # translocation
    sel.bnd.1 = x$svtype == 'BND' & x$chrom1 == chrom & x$pos1 <= right & x$pos1 >= left
    sel.bnd.2 = x$svtype == 'BND' & x$chrom2 == chrom & x$pos2 <= right & x$pos2 >= left
    x$pos1[x$svtype == 'BND' & !sel.bnd.1] = NA
    x$pos2[x$svtype == 'BND' & !sel.bnd.2] = NA

    x = x[sel | sel.bnd.1 | sel.bnd.2,]
    #print(table(x$svtype))
    #print(table(x$chrom1, x$chrom2))

    # order/pile up
    x$sign = ifelse(x$svtype == 'DUP', 1, 0)
    x = x[order(x$svtype, x$sign*x$pos1),]
    x$samp = paste0(x$svtype, '/', x$sample, '/', x$pos1)
    x$samp = factor(x$samp, levels=unique(x$samp))
    return(x)
}



################################# FUNCTION TO PLOT SVs EXPRESSION #########################################
plot.exp <- function (g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, gene, pk) {
  
  ##################### plot the expression for aLL SV samples (SVs types vs non-SVs) ###############################
  if (length(g.exp[g.exp$sample %in% nonSV.pats, 'gene.exp']) > 0 ) {
    SVs.exp = data.frame(grp = "SVs", exp = g.exp[g.exp$sample.status=="SVs", 'gene.exp'])
    nonSVs.exp = data.frame(grp = "non-SVs", exp = g.exp[g.exp$sample %in% nonSV.pats, 'gene.exp'])
    
    if (length(BND.pats) > 0 & length(g.exp[g.exp$sample %in% BND.pats, 'gene.exp']) > 0) {BND.exp.n = data.frame(grp = "BND", exp = g.exp[g.exp$sample %in% BND.pats, 'gene.exp'])} else {BND.exp.n=NULL}
    if (length(DUP.pats) > 0 & length(g.exp[g.exp$sample %in% DUP.pats, 'gene.exp']) > 0) {DUP.exp.n = data.frame(grp = "DUP", exp = g.exp[g.exp$sample %in% DUP.pats, 'gene.exp'])} else {DUP.exp.n=NULL}
    if (length(INS.pats) > 0 & length(g.exp[g.exp$sample %in% INS.pats, 'gene.exp']) > 0) {INS.exp.n = data.frame(grp = "INS", exp = g.exp[g.exp$sample %in% INS.pats, 'gene.exp'])} else {INS.exp.n=NULL}
    if (length(DEL.pats) > 0 & length(g.exp[g.exp$sample %in% DEL.pats, 'gene.exp']) > 0) {DEL.exp.n = data.frame(grp = "DEL", exp = g.exp[g.exp$sample %in% DEL.pats, 'gene.exp'])} else {DEL.exp.n=NULL}
    if (length(INV.pats) > 0 & length(g.exp[g.exp$sample %in% INV.pats, 'gene.exp']) > 0) {INV.exp.n = data.frame(grp = "INV", exp = g.exp[g.exp$sample %in% INV.pats, 'gene.exp'])} else {INV.exp.n=NULL}
    
    svtype.exp = do.call('rbind', list(nonSVs.exp, SVs.exp, BND.exp.n, DUP.exp.n, INS.exp.n, DEL.exp.n, INV.exp.n))
    svtype.exp$grp <- factor(svtype.exp$grp, levels=c('non-SVs', 'SVs', 'BND','DUP', 'INS','DEL', 'INV'))
  } else {
    svtype.exp <- NULL
  }
  
  sv.grps.n = as.data.frame(combn(as.character(unique(svtype.exp$grp)),2), stringsAsFactors = F)
  idxs <- as.data.frame(which(sv.grps.n =="non-SVs", arr.ind=TRUE))
  sv.grps.n <- sv.grps.n[,idxs$col]
  svCMPlist = as.list(sv.grps.n[,1:ncol(sv.grps.n)])
  
  title.size1 = length(levels(factor(svtype.exp$grp))) * 3.5
  if (title.size1 <= 7) { title.size = 12 }
  
  lbls.n = c(paste0("non-SVs\n(n=",nrow(nonSVs.exp),")"), paste0("SVs\n(n=",nrow(SVs.exp),")"), 
             paste0("BND\n(n=",nrow(BND.exp.n),")"), 
             paste0("DUP\n(n=",nrow(DUP.exp.n),")"), paste0("INS\n(n=",nrow(INS.exp.n),")"),
             paste0("DEL\n(n=",nrow(DEL.exp.n),")"), paste0("INV\n(n=",nrow(INV.exp.n),")"))
  lbls.n = grep(paste(as.character(unique(svtype.exp$grp)), collapse="|"), lbls.n, value=TRUE)
  
  max.exp = max(log2(svtype.exp$exp+1))
  e1 <- ggplot(svtype.exp, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp), outlier.shape=NA, outlier.size=0.5) 
  e1 <- e1 + geom_jitter(position=position_jitter(0.3), shape=1, size=0.75) + theme_bw() 
  e1 <- e1 + labs(x='', y=paste(g, 'expression')) + ggtitle('All samples')
  e1 <- e1 + theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
                   axis.text.y=element_text(size=12, color="black"), 
                   axis.title.y=element_text(size=12),
                   #panel.background=element_rect(color="black"),
                   plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
                   legend.position="none",
                   panel.border = element_rect(linetype='solid', color='black'),
                   plot.margin=unit(c(1,1,1,5), 'mm')
  )
  e1 = e1 + scale_fill_manual(name="", values =c("non-SVs"="gray", "SVs"="orange2", "BND"="#2ca25f", "DUP"="#b53f4d", "INS"="#fec44f","DEL"="#2c7fb8","INV"="#c994c7"))
  e1 = e1 + scale_x_discrete(labels=lbls.n)
  e1 = e1 + scale_y_continuous(expand=c(0.05,0.05*max.exp))
  e1 = e1 + geom_signif(comparisons= svCMPlist, step_increase=0.1, textsize = 4, map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  
  # pval = (wilcox.test(g.exp[g.exp$sample.status=="SVs", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp']))$p.value
  # #pval= sprintf(pval, fmt="%#.5f")
  # pval.txt = sprintf("p = %.1g", pval)
  # pval = signif(pval, digits=3)
  # 
  # g.exp$sample.status <- factor(g.exp$sample.status, levels=c('non-SVs', 'SVs'))
  # max.exp = max(log2(g.exp$gene.exp+1))
  # e1 <- ggplot(g.exp, aes(x=sample.status, y=log2(gene.exp+1))) 
  # e1 <- e1 + geom_boxplot(aes(fill=sample.status), outlier.shape=NA, outlier.size=0.5) + geom_jitter(position=position_jitter(0.3), shape=1, size=0.75) + theme_bw() 
  # e1 <- e1 + labs(x='', y=paste(gene, 'expression'))# + ggtitle(paste0('\n',gene, ' expression in SV\nand non-SV samples'))
  # e1 <- e1 + ggtitle('All samples')
  # e1 <- e1 + theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
  #                  axis.text.y=element_text(size=12, color="black"), 
  #                  axis.title.y=element_text(size=12),
  #                  #panel.background=element_blank(),
  #                  plot.title = element_text(size = 12, hjust=0.5, color="black", face="plain"),
  #                  legend.position="none",
  #                  panel.border = element_rect(linetype='solid', color='black'),
  #                  plot.margin=unit(c(1,1,1,5), 'mm')
  # )
  # e1 = e1 + scale_fill_manual(name="", values =c("non-SVs"="gray", "SVs"="orange2"))
  # e1 = e1 + scale_x_discrete(labels=c(paste0("None\n(n=",nrow(g.exp[g.exp$sample.status=="non-SVs",]),")"), paste0("SVs\n(n=",nrow( g.exp[g.exp$sample.status=="SVs",]),")")))
  # #e1 = e1 + geom_signif(comparisons=list(c('non-SVs','SVs')))
  # e1 = e1 + annotate("text", x = 1.5, y = max.exp*1.1, label = pval.txt, cex=4, hjust=0.5)
  # e1 = e1 + scale_y_continuous(expand=c(0.05,0.05*max.exp))
  
  ###################### plot the expression for SV samples (SV types vs non-SVs) ########################
  if (length(g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="neut", 'gene.exp']) > 0 ) {
    SVs.exp = data.frame(grp = "SVs", exp = g.exp[g.exp$sample.status=="SVs" & g.exp$gene.cn.status=="neut", 'gene.exp'])
    nonSVs.exp = data.frame(grp = "non-SVs", exp = g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])
    
    if (length(BND.pats) > 0 & length(g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {BND.exp.n = data.frame(grp = "BND", exp = g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="neut" , 'gene.exp'])} else {BND.exp.n=NULL}
    if (length(DUP.pats) > 0 & length(g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {DUP.exp.n = data.frame(grp = "DUP", exp = g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {DUP.exp.n=NULL}
    if (length(INS.pats) > 0 & length(g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {INS.exp.n = data.frame(grp = "INS", exp = g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {INS.exp.n=NULL}
    if (length(DEL.pats) > 0 & length(g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {DEL.exp.n = data.frame(grp = "DEL", exp = g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {DEL.exp.n=NULL}
    if (length(INV.pats) > 0 & length(g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {INV.exp.n = data.frame(grp = "INV", exp = g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {INV.exp.n=NULL}
    
    svtype.exp = do.call('rbind', list(nonSVs.exp, SVs.exp, BND.exp.n, DUP.exp.n, INS.exp.n, DEL.exp.n, INV.exp.n))
    svtype.exp$grp <- factor(svtype.exp$grp, levels=c('non-SVs', 'SVs', 'BND','DUP', 'INS','DEL', 'INV'))
  } else {
    svtype.exp <- NULL
  }
  
  ## plot 
  if (!is.null(svtype.exp)) {
    ### prepare comarison list
    sv.grps.n = as.data.frame(combn(as.character(unique(svtype.exp$grp)),2), stringsAsFactors = F)
    idxs <- as.data.frame(which(sv.grps.n =="non-SVs", arr.ind=TRUE))
    sv.grps.n <- sv.grps.n[,idxs$col]
    svCMPlist = as.list(sv.grps.n[,1:ncol(sv.grps.n)])
    
    title.size1 = length(levels(factor(svtype.exp$grp))) * 3.5
    if (title.size1 <= 7) { title.size = 12 }
    
    lbls.n = c(paste0("non-SVs\n(n=",nrow(nonSVs.exp),")"), paste0("SVs\n(n=",nrow(SVs.exp),")"), 
               paste0("BND\n(n=",nrow(BND.exp.n),")"), 
               paste0("DUP\n(n=",nrow(DUP.exp.n),")"), paste0("INS\n(n=",nrow(INS.exp.n),")"),
               paste0("DEL\n(n=",nrow(DEL.exp.n),")"), paste0("INV\n(n=",nrow(INV.exp.n),")"))
    lbls.n = grep(paste(as.character(unique(svtype.exp$grp)), collapse="|"), lbls.n, value=TRUE)
    
    #neut.title = paste0(g, ' expression by SV type in peak ', pk,'\n(CN neutral samples only)')
    max.exp = max(log2(svtype.exp$exp+1))
    e2 <- ggplot(svtype.exp, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp), outlier.shape=NA, outlier.size=0.5) 
    e2 <- e2 + geom_jitter(position=position_jitter(0.3), shape=1, size=0.75) + theme_bw() 
    e2 <- e2 + labs(x='', y=paste(g, 'expression')) + ggtitle('Gene is CN-neutral')
    e2 <- e2 + theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
                     axis.text.y=element_text(size=12, color="black"), 
                     axis.title.y=element_text(size=12),
                     #panel.background=element_rect(color="black"),
                     plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
                     legend.position="none",
                     panel.border = element_rect(linetype='solid', color='black'),
                     plot.margin=unit(c(1,1,1,5), 'mm')
    )
    e2 = e2 + scale_fill_manual(name="", values =c("non-SVs"="gray", "SVs"="orange2", "BND"="#2ca25f", "DUP"="#b53f4d", "INS"="#fec44f","DEL"="#2c7fb8","INV"="#c994c7"))
    e2 = e2 + scale_x_discrete(labels=lbls.n)
    e2 = e2 + scale_y_continuous(expand=c(0.05,0.05*max.exp))
    e2 = e2 + geom_signif(comparisons= svCMPlist, step_increase=0.1, textsize = 4, map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  } else {
    e2 <- NULL
  }
  # if (length(g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="neut", 'gene.exp']) > 0 ) {
  #   nonSVs.exp = data.frame(grp = "non-SVs", exp = g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])
  #   if (length(BND.pats) > 0 & length(g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {BND.exp.n = data.frame(grp = "BND", exp = g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="neut" , 'gene.exp'])} else {BND.exp.n=NULL}
  #   if (length(DUP.pats) > 0 & length(g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {DUP.exp.n = data.frame(grp = "DUP", exp = g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {DUP.exp.n=NULL}
  #   if (length(INS.pats) > 0 & length(g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {INS.exp.n = data.frame(grp = "INS", exp = g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {INS.exp.n=NULL}
  #   if (length(DEL.pats) > 0 & length(g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {DEL.exp.n = data.frame(grp = "DEL", exp = g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {DEL.exp.n=NULL}
  #   if (length(INV.pats) > 0 & length(g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {INV.exp.n = data.frame(grp = "INV", exp = g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {INV.exp.n=NULL}
  #   
  #   svtype.exp = do.call('rbind', list(nonSVs.exp, BND.exp.n, DUP.exp.n, INS.exp.n, DEL.exp.n, INV.exp.n))
  #   svtype.exp$grp <- factor(svtype.exp$grp, levels=c('non-SVs', 'BND','DUP', 'INS','DEL', 'INV'))
  # } else {
  #   svtype.exp <- NULL
  # }
  # 
  # 
  # ######### plot SVs types for neutral samples only ######
  # if (!is.null(svtype.exp)) {
  #   ### prepare comarison list
  #   sv.grps.n = as.data.frame(combn(as.character(unique(svtype.exp$grp)),2), stringsAsFactors = F)
  #   idxs <- as.data.frame(which(sv.grps.n =="non-SVs", arr.ind=TRUE))
  #   sv.grps.n <- sv.grps.n[,idxs$col]
  #   svCMPlist = as.list(sv.grps.n[,1:ncol(sv.grps.n)])
  #   
  #   title.size1 = length(levels(factor(svtype.exp$grp))) * 3.5
  #   if (title.size1 <= 7) { title.size = 12 }
  #   
  #   lbls.n = c(paste0("non-SVs\n(n=",nrow(nonSVs.exp),")"), paste0("BND\n(n=",nrow(BND.exp.n),")"), 
  #              paste0("DUP\n(n=",nrow(DUP.exp.n),")"), paste0("INS\n(n=",nrow(INS.exp.n),")"),
  #              paste0("DEL\n(n=",nrow(DEL.exp.n),")"), paste0("INV\n(n=",nrow(INV.exp.n),")"))
  #   lbls.n = grep(paste(as.character(unique(svtype.exp$grp)), collapse="|"), lbls.n, value=TRUE)
  #   neut.title = paste0(gene, ' expression by SV type in peak ', pk,'\n(CN neutral samples only)')
  #   neut.title = 'Gene is CN-neutral'
  #   max.exp = max(log2(svtype.exp$exp+1))
  #   e2 <- ggplot(svtype.exp, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp), outlier.shape=NA, outlier.size=0.5) 
  #   e2 <- e2 + geom_jitter(position=position_jitter(0.3), shape=1, size=0.75) + theme_bw() 
  #   e2 <- e2 + labs(x='', y=paste(gene, 'expression')) + ggtitle(neut.title)
  #   e2 <- e2 + theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
  #                    axis.text.y=element_text(size=12, color="black"), 
  #                    axis.title.y=element_text(size=12),
  #                    #panel.background=element_rect(color="black"),
  #                    plot.title = element_text(size = 12, hjust=0.5, color="black", face="plain"),
  #                    legend.position="none",
  #                    panel.border = element_rect(linetype='solid', color='black'),
  #                    plot.margin=unit(c(1,1,1,5), 'mm')
  #              )
  #   e2 = e2 + scale_fill_manual(name="", values =c("non-SVs"="gray", "BND"="#2ca25f", "DUP"="#b53f4d", "INS"="#fec44f","DEL"="#2c7fb8","INV"="#c994c7"))
  #   e2 = e2 + scale_x_discrete(labels=sub('non-SVs', 'None', lbls.n))
  #   e2 = e2 + scale_y_continuous(expand=c(0.05,0.05*max.exp))
  #   e2 = e2 + geom_signif(comparisons= svCMPlist, step_increase=0.1, textsize = 4, map_signif_level=function(p)sprintf("p = %.1g", p), show.legend=F)
  # } else {
  #   e2 <- NULL
  # }
   
  ####################### plot gene expression incorporating copy number #########################################
  #if (is.cn.avail) {
    g.exp$group = "NC"
    g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="neut", "group"] = "GNPN"
    g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="neut", "group"] = "GAPN"
    g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="amp", "group"] = "GNPA"
    g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="amp", "group"] = "GAPA"
    g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="neut", "group"] = "GDPN"
    g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="del", "group"] = "GNPD"
    g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="del", "group"] = "GDPD"
    
    g.exp$group <- factor(g.exp$group, levels=c('GNPN', 'GNPA', 'GAPN','GAPA','GDPN','GNPD','GDPD','NC'))
    
    lbls = c("GNPN"=paste0("GeneNeut\npeakNeut\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="neut", "group"]),")"), 
             "GAPN"=paste0("GeneAmp\nPeakNeut\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="neut", "group"]),")"),
             "GNPA"=paste0("GeneNeut\nPeakAmp\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="amp", "group"]),")"),
             "GAPA"=paste0("GeneAmp\nPeakAmp\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="amp", "group"]),")"), 
             "GDPN"=paste0("GeneDel\nPeakNeut\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="neut", "group"]),")"),
             "GNPD"=paste0("GeneNeut\nPeakDel\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="del", "group"]),")"),
             "GDPD"=paste0("GeneDel\nPeakDel\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="del", "group"]),")"), 
             "NC"=paste0("Other\n(n=",nrow(g.exp[g.exp$group=="NC",])))
  
#     lbls = c("GNPN"=paste0("GN/PN\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="neut", "group"]),")"), 
#              "GAPN"=paste0("GA/PN\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="neut", "group"]),")"),
#              "GNPA"=paste0("GN/PA\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="amp", "group"]),")"),
#              "GAPA"=paste0("GA/PA\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="amp", "group"]),")"), 
#              "GDPN"=paste0("GD/PN\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="neut", "group"]),")"),
#              "GNPD"=paste0("GN/PD\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="del", "group"]),")"),
#              "GDPD"=paste0("GD/PD\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="del", "group"]),")"), 
#              "NC"=paste0("NC\n(n=",nrow(g.exp[g.exp$group=="NC",])))

    ### construct the list for all possible values 
    amp.grp = as.character(unique(g.exp$group[g.exp$group %in% c("GNPN","GAPN", "GNPA","GAPA")])) 
    del.grp = as.character(unique(g.exp$group[g.exp$group %in% c("GNPN", "GDPN", "GNPD", "GDPD")])) 
    
    if (length(amp.grp) > 2) { 
      amp.cmps = data.frame(combn(amp.grp,2), stringsAsFactors = F)
      ampCMPlist = as.list(amp.cmps[,1:ncol(amp.cmps)]) 
    } else { ampCMPlist = NULL }
    
    if (length(del.grp) > 2) { 
      del.cmps = data.frame(combn(del.grp,2), stringsAsFactors = F)
      delCMPlist = as.list(del.cmps[,1:ncol(del.cmps)]) 
    } else { delCMPlist = NULL}
    
    myCMPlist = c(ampCMPlist, delCMPlist)
    
    if (length(amp.grp) == 2) { myCMPlist[[length(myCMPlist)+1]] = amp.grp }
    if (length(del.grp) ==2) { myCMPlist[[length(myCMPlist)+1]] = del.grp }
    
    title.size2 = length(levels(factor(g.exp$group))) * 3.5  
    if (title.size2 <= 7) { title.size2 = 12 }
    
    tt = paste0(gene,' expression with the presence and\nabsence of CN at ', gene, '/peak regions')
    tt = 'By peak/gene CN status'

    max.exp = max(log2(g.exp$gene.exp+1))
    e3 <- ggplot(g.exp, aes(x=group, y=log2(gene.exp+1))) + geom_boxplot(aes(fill=group), outlier.shape=NA, outlier.size=0.5)
    e3 <- e3 + geom_jitter(position=position_jitter(0.3), shape=1, size=0.75) + theme_bw() 
    e3 <- e3 + labs(x='', y=paste(gene, 'expression')) + ggtitle(tt)
    e3 <- e3 + theme(axis.text.x=element_text(size=11, vjust=0.5, color="black", angle=60),
                     axis.text.y=element_text(size=12, color="black"),
                     axis.title=element_text(size=12), 
                     #panel.background=element_blank(),
                     plot.title = element_text(size = 12, hjust=0.5, color="black", face="plain"),
                     legend.position="none",
                     panel.border = element_rect(linetype='solid', color='black'),
                     plot.margin=unit(c(1,1,1,5), 'mm')
               )
    e3 = e3 + scale_fill_manual(name="", values = c("GNPN"="gray", "GAPN"="#f03b20", "GNPA"="#b53f4d", "GAPA"="salmon", 
                                                    "GDPN"="#a6bddb", "GNPD"="#2c7fb8", "GDPD"="skyblue2")) 
    e3 = e3 + scale_x_discrete(labels=lbls)
    e3 = e3 + scale_y_continuous(expand=c(0.05,0.05*max.exp))
    e3 = e3 + geom_signif(comparisons=myCMPlist , step_increase=0.15, textsize = 4, map_signif_level=function(p)sprintf("p = %.1g", p), show.legend=F)
  #}    
  
  ### preapre and return results 
  #width1 = length(levels(svtype.exp$grp))*0.5
  #width2 = length(levels(factor(g.exp$group)))*0.5
  #if (width1==1) { width1 =1.5}
  #if (width2==1) { width2 =1.5}
  col.width = max(length(levels(svtype.exp$grp)), length(levels(factor(g.exp$group))))*0.5
  if (col.width <= 1) { col.width = 2 }
  e1 = e1 + theme(legend.position='none')
  e2 = e2 + theme(legend.position='none')
  e3 = e3 + theme(legend.position='none')
  exp.plots <- list(e1,e2,e3, col.width)
  names(exp.plots) <- c( "e1", "e2","e3", "ColWidth")
  return (exp.plots)
  
}
##################################################################################################################


############################## FUNCTION TO PLOT SVs EXPRESSION FOR AMP/DEL #######################################
plot.exp.amp.del <- function (g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, out.dir, gene, pk) {
  
  if (length(g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="amp", 'gene.exp']) > 0 ) {
    nonSVs.exp.a = data.frame(grp = "non-SVs", exp = g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="amp", 'gene.exp'])
    if (length(BND.pats) > 0 & length(g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="amp", 'gene.exp']) > 0) {BND.exp.a = data.frame(grp = "BND", exp = g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="amp", 'gene.exp'])} else {BND.exp.a=NULL}
    if (length(DUP.pats) > 0 & length(g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="amp", 'gene.exp']) > 0) {DUP.exp.a = data.frame(grp = "DUP", exp = g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="amp", 'gene.exp'])} else {DUP.exp.a=NULL}
    if (length(INS.pats) > 0 & length(g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="amp", 'gene.exp']) > 0) {INS.exp.a = data.frame(grp = "INS", exp = g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="amp", 'gene.exp'])} else {INS.exp.a=NULL}
    if (length(DEL.pats) > 0 & length(g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="amp", 'gene.exp']) > 0) {DEL.exp.a = data.frame(grp = "DEL", exp = g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="amp", 'gene.exp'])} else {DEL.exp.a=NULL}
    if (length(INV.pats) > 0 & length(g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="amp", 'gene.exp']) > 0) {INV.exp.a = data.frame(grp = "INV", exp = g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="amp", 'gene.exp'])} else {INV.exp.a=NULL}
    svtype.exp.a = do.call('rbind', list(nonSVs.exp.a,BND.exp.a, DUP.exp.a, INS.exp.a, DEL.exp.a, INV.exp.a))
    svtype.exp.a$grp <- factor(svtype.exp.a$grp, levels=c('non-SVs', 'BND','DUP', 'INS','DEL', 'INV'))
  } else {
    svtype.exp.a <- NULL
  }
   
  ######### plot SVs types for amplified samples only ######
  if (!is.null(svtype.exp.a) & length(unique(svtype.exp.a$grp)) >=2 ) {
   
    ### prepare comarison list
    sv.grps.a = as.data.frame(combn(as.character(unique(svtype.exp.a$grp)),2), stringsAsFactors = F)
    idxs <- as.data.frame(which(sv.grps.a =="non-SVs", arr.ind=TRUE))
    sv.grps.a <- sv.grps.a[,idxs$col]
    svCMPlist = as.list(sv.grps.a[,1:ncol(sv.grps.a)])
    
    title.size1 = length(levels(factor(svtype.exp.a$grp))) * 3.5
    if (title.size1 <= 7 | title.size1 > 12) { title.size1 = 12 }
    
    lbls.a = c(paste0("non-SVs\n(n=",nrow(nonSVs.exp.a),")"), paste0("BND\n(n=",nrow(BND.exp.a),")"), 
               paste0("DUP\n(n=",nrow(DUP.exp.a),")"), paste0("INS\n(n=",nrow(INS.exp.a),")"),
               paste0("DEL\n(n=",nrow(DEL.exp.a),")"), paste0("INV\n(n=",nrow(INV.exp.a),")"))
    lbls.a = grep(paste(as.character(unique(svtype.exp.a$grp)), collapse="|"), lbls.a, value=TRUE)
    amp.title = paste0(gene, ' expression by SV type in peak ', pk,'\n(amplified samples only)')
    
    a <- ggplot(svtype.exp.a, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp))
    a <- a + geom_jitter(position=position_jitter(0.3)) + theme_bw() 
    a <- a + labs(x='', y=paste(gene, 'expression')) + ggtitle(amp.title)
    a <- a + theme(axis.text.x=element_text(size=12, vjust=0.5, color="black"),
                       axis.text.y=element_text(size=12, color="black"), 
                       axis.title.y=element_text(size=14), panel.background=element_blank(),
                       plot.title = element_text(size = title.size1, hjust=0.5, color="black", face="bold"),
                       legend.position="none")
    a = a + scale_fill_manual(name="", values =c("non-SVs"="gray", "BND"="#2ca25f", "DUP"="#b53f4d", "INS"="#fec44f","DEL"="#2c7fb8","INV"="#c994c7"))
    a = a + scale_x_discrete(labels=lbls.a)
    a = a + geom_signif(comparisons= svCMPlist, step_increase=0.1)
    
    png(paste0(out.dir, "_amp_svtypes.png"))
    print(a)
    dev.off()
  }
  

  ######### plot SVs types for deleted samples only ######
  if (length(g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="del", 'gene.exp']) > 0 ) {
    nonSVs.exp.d = data.frame(grp = "non-SVs", cn.call="del", exp = g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="del", 'gene.exp'])
    if (length(BND.pats) > 0 & length(g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="del", 'gene.exp'])) {BND.exp.d = data.frame(grp = "BND", cn.call="del", exp = g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="del", 'gene.exp'])} else {BND.exp.d=NULL}
    if (length(DUP.pats) > 0 & length(g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="del", 'gene.exp'])) {DUP.exp.d = data.frame(grp = "DUP", cn.call="del", exp = g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="del", 'gene.exp'])} else {DUP.exp.d=NULL}
    if (length(INS.pats) > 0 & length(g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="del", 'gene.exp'])) {INS.exp.d = data.frame(grp = "INS", cn.call="del", exp = g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="del", 'gene.exp'])} else {INS.exp.d=NULL}
    if (length(DEL.pats) > 0 & length(g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="del", 'gene.exp'])) {DEL.exp.d = data.frame(grp = "DEL", cn.call="del", exp = g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="del", 'gene.exp'])} else {DEL.exp.d=NULL}
    if (length(INV.pats) > 0 & length(g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="del", 'gene.exp'])) {INV.exp.d = data.frame(grp = "INV", cn.call="del", exp = g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="del", 'gene.exp'])} else {INV.exp.d=NULL}
    svtype.exp.d = do.call('rbind', list(nonSVs.exp.d,BND.exp.d, DUP.exp.d, INS.exp.d, DEL.exp.d, INV.exp.d))
    svtype.exp.d$grp <- factor(svtype.exp.d$grp, levels=c('non-SVs', 'BND','DUP', 'INS','DEL', 'INV'))
  } else {
    svtype.exp.d <- NULL 
  }
  
  if (!is.null(svtype.exp.d) & length(unique(svtype.exp.d$grp)) >=2) {
    
    ### prepare comarison list
    sv.grps.d = as.data.frame(combn(as.character(unique(svtype.exp.d$grp)),2), stringsAsFactors = F)
    idxs <- as.data.frame(which(sv.grps.d =="non-SVs", arr.ind=TRUE))
    sv.grps.d <- sv.grps.d[,idxs$col]
    svCMPlist = as.list(sv.grps.d[,1:ncol(sv.grps.d)])
    
    title.size2 = length(levels(factor(svtype.exp.d$grp))) * 3.5
    if (title.size2 <= 7 | title.size2 > 12) { title.size2 = 12 }
    
    lbls.d = c(paste0("non-SVs\n(n=",nrow(nonSVs.exp.d),")"), paste0("BND\n(n=",nrow(BND.exp.d),")"), 
               paste0("DUP\n(n=",nrow(DUP.exp.d),")"), paste0("INS\n(n=",nrow(INS.exp.d),")"),
               paste0("DEL\n(n=",nrow(DEL.exp.d),")"), paste0("INV\n(n=",nrow(INV.exp.d),")"))
    lbls.d = grep(paste(as.character(unique(svtype.exp.d$grp)), collapse="|"), lbls.d, value=TRUE)
    del.title = paste0(gene, ' expression by SV type in peak ', pk,'\n(deleted samples only)')
                       
    d <- ggplot(svtype.exp.d, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp))
    d <- d + geom_jitter(position=position_jitter(0.3)) + theme_bw() 
    d <- d + labs(x='', y=paste(gene, 'expression')) + ggtitle(del.title)
    d <- d + theme(axis.text.x=element_text(size=12, vjust=0.5, color="black"),
                       axis.text.y=element_text(size=12, color="black"), 
                       axis.title.y=element_text(size=14), panel.background=element_blank(),
                       plot.title = element_text(size = title.size2, hjust=0.5, color="black", face="bold"),
                       legend.position="none")
    d = d + scale_fill_manual(name="", values =c("non-SVs"="gray", "BND"="#2ca25f", "DUP"="#b53f4d", "INS"="#fec44f","DEL"="#2c7fb8","INV"="#c994c7"))
    d = d + scale_x_discrete(labels=lbls.d)
    d = d + geom_signif(comparisons= svCMPlist, step_increase=0.1)
    
    png(paste0(out.dir, "_del_svtypes.png"))
    print(d)
    dev.off()
  }
 
}
##################################################################################################################


################################### FUNCTION TO PLOT PEAKS REGIONS ################################################
plot.region <- function(pk, pk.corr, gene, genes.in.p, p.roi, D=NULL){

   #construct region coordinates 
   right = max(pk.corr$Start, pk.corr$End, genes.in.p$g.start, genes.in.p$g.stop)
   left =  min(pk.corr$Start, pk.corr$End, genes.in.p$g.start, genes.in.p$g.stop)
   width = abs(pk.corr$End - pk.corr$Start)/1000

   ### add left and right extensions if provided 
   left <- left - left.ext     
   right <- right + right.ext
   if(gene %in% c('ETV1', 'FOXA1')){
       ww = pk.corr$End - pk.corr$Start  + 1
       left = pk.corr$End + ww*3
       right = pk.corr$Start - ww*3
   }
   D = right - left
   #scale binwidth accordingly based on region width
   binwidth = D/75
   #genes within region
   g.corr = genes.in.p[genes.in.p$gene ==gene, ]
 
   ### extract SVs data
   x = cts[cts$chr == pk.corr$Chr & cts$pos > left & cts$pos < right,]  
   x = x[x$svtype !="ALL",]
   
   ### for DUP and DEL only 
   x2 = x 
   x2 = x2[x2$svtype  %in% c("DUP","DEL"),]
   x2[x2$svtype == "DEL", "num.samples"] <-  x2[x2$svtype == "DEL", "num.samples"] * -1

   ### make the title 
   #title = paste0('Associated Gene: ',gene,' (Peak locus: ',pk.corr$Chr, ':',  pk.corr$Start, '-',  pk.corr$End,')')

   ### compute the DUP and DEL pileup of SV
   dup_del = pileUp(sv, pk.corr$Chr, left, right)
   
   dup_del$pos1[dup_del$pos1 < left] = left
   dup_del$pos2[dup_del$pos2 > right] = right

   rwidth = right - left + 1
   #print(rwidth)
   rstep = max(round(rwidth/5/5/10^5)*5*10^5, 2*10^5)
   #print(rstep)
   brks = seq(0,10^9,rstep)
   brks = brks[brks >= left & brks <= right]
   xlabs = as.character(brks/(10^6))
   xlabs[length(xlabs)] = paste0(xlabs[length(xlabs)], 'Mb')

   #print(brks)
   #print(labs)
   

  ################################## plot region copy number ###############################################
  #if (is.cn.avail) {
     reg.width = (right - left)+1
     reg.cn = cn_data[cn_data$chrom == pk.corr$Chr & (cn_data$pos > left | cn_data$pos < right),]
     reg.cn = reg.cn[tolower(reg.cn$cn.call) %in% c("amp", "del"), ]
     
     if (nrow(reg.cn) !=0) { 
       
       s = round(reg.width/200)
       imin = min(reg.cn$start)
       imax = max(reg.cn$end)
       chr.name=unique(reg.cn$chrom)
       reg.win = data.frame(chr=chr.name, start=seq(imin,imax,s))
       reg.win$stop = reg.win$start + s - 1

       #TODO: randomize file name and use system tmp for tmp output such as reg_win.tsv, reg.cn.tsv
       write.table(reg.cn, file=paste0(res.dir,'/processed_data/reg.cn.tsv'), quote=F, row.names=F, sep="\t", col.names=F)
       write.table(reg.win, file=paste0(res.dir,'/processed_data/reg.win.tsv'), quote=F, row.names=F, sep="\t", col.names=F)
       
       ### overlap windows with copy number 
       #cat ("Overlapping windows with copy number ..")
       #TODO: randomize file name and use system tmp for tmp output such as win_cn.tsv
       system(paste0("intersectBed -wo -a ",res.dir,"/processed_data/reg.win.tsv -b ", res.dir,"/processed_data/reg.cn.tsv | cut -f2,3,7,9 | sort | uniq | sort -k1,1 -k2,2 | groupBy -full -g 1,2 -c 4 -o count > ", res.dir, "/processed_data/win_cn.tsv")) 
       
       win.data = read.table(paste0(res.dir,"/processed_data/win_cn.tsv"), header = F, sep ="\t")
       colnames(win.data) = c('start', 'stop','sample', 'cn.call', 'num.samples')
       win.data$pos = (win.data$start + win.data$stop)/2
       ### add dummy data for visualization purpuses 
       win.data=rbind(win.data, data.frame(start=-1,stop=-2,sample="dummp1",cn.call="amp",num.samples=0, pos=-1))
       win.data=rbind(win.data, data.frame(start=-1,stop=-2,sample="dummp1",cn.call="del",num.samples=0, pos=-1))
       
       # manually create breaks for y axis, and change 1st to '000' so plots will align
       nn = max(win.data$num.samples); ii = max(round(nn/40)*10,5)
       bb = seq(0,nn,ii); ll=bb; ll[1]='000'

       p0 = ggplot(win.data, aes(x=pos,y=num.samples, fill=cn.call)) + geom_bar(stat="identity")
       p0 = p0 + scale_fill_manual(name="", values=c("amp"="#b53f4d", "del"="#2c7fb8"), labels=c("amp"="Gain", "del"="Loss"))
       p0 = p0 + theme_bw() + xlab('') + ylab('CNA freq.')
       p0 = p0 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                       plot.title=element_text(size=12, hjust=0.5),
                       #axis.ticks = element_blank(),
                       #axis.text.x=element_blank(),
                       axis.text.x=element_text(size=10, color="black"),
                       axis.text.y=element_text(size=12, color="black"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_text(size=12, color="black"),
                       legend.key.size = unit(2.5,'mm'),
                       legend.title=element_blank(),
                       legend.text=element_text(size=9),
                       legend.position=c(1-pctRightExpand/2,0.5),
                       legend.background=element_blank(),
                       plot.margin = unit(c(1,1,1,1), 'mm'),
                       panel.border = element_rect(color='black', linetype='solid')
                 )
       #p0 = p0 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
       p0 = p0 + scale_x_continuous(breaks=brks, labels = xlabs , limits=c(left, right), expand=c(0.05,0.05), position='top') + scale_y_continuous(breaks=bb, labels=ll)

       p0 = p0 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
       p0 = p0 + ggtitle(paste0(pk, '(', p.corr$Chr, ':',  p.corr$Start, '-',  p.corr$End, ')'))
       p0 = p0 + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))
       
     } else {
        p0 <- NULL
      }    
    
  #}

  ################################## Plot DUP $ DEL Freq ####################################
    dup_del$svtype <- factor(dup_del$svtype, levels=c('DUP', 'DEL'))

    p1 = (ggplot(dup_del)
    + geom_segment(aes(x=pos1, xend=pos2, y=samp, yend=samp, color=svtype))
    + theme_bw()
    + ylab('Dup/del pileup')
    + scale_x_continuous(breaks=brks, labels=xlabs, expand=c(0.05,0.05))
    + scale_y_discrete(breaks=dup_del$samp[c(1, nrow(dup_del))], labels=c('000',''))
    + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.title=element_blank(), #axis.ticks.y = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=12, color='white'),
            axis.title.x=element_blank(), axis.title.y=element_text(size=12, color="black"),
            legend.key.size = unit(2.5,'mm'),
            legend.title=element_blank(),
            legend.text=element_text(size=9),
            legend.position=c(1-pctRightExpand/2,0.5),
            legend.background=element_blank(),
            plot.margin = unit(c(1,1,1,1), 'mm'),
            panel.border = element_rect(color='black', linetype='solid')
    )
    + scale_color_manual(name="", values=colrs)
    + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
    + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))
    )
  
  ################################## plot SVs (DUP and DEL only) ###############################################
  x2$svtype <- factor(x2$svtype, levels=c('DUP','BND','INS','INV','DEL'))
  nmax = max(x2$num.samples); nmin = min(x2$num.samples)
  # manually create breaks for y axis, and change 1st to '000' so plots will align
  nmax = max(x2$num.samples); nmin = min(x2$num.sample)
  nn = nmax - nmin + 1; ii = max(round(nn/40)*10,5)
  bbdup = seq(0,nmax,ii); lldup=bbdup; lldup[1]='000'
  bbdel = seq(0,nmin,-ii); bbdel = bbdel[-1]
  bb = c(rev(bbdel), bbdup); ll = c(bbdel, lldup); ll = gsub('-', '', ll)


  p2 = ggplot(x2, aes(x=pos, y=num.samples, fill=svtype)) + geom_bar(stat="identity")
  p2 = p2 + theme_bw() + xlab('') + ylab('Num. samples')
  p2 = p2 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  p2 = p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 plot.title=element_blank(), #axis.ticks = element_blank(),
                 axis.text.x=element_blank(), axis.text.y=element_text(size=12, color="black"),
                 axis.title.x=element_blank(),
                 axis.title.y=element_text(size=12, color="black"),
                 legend.key.size = unit(2.5,'mm'),
                 legend.title=element_blank(),
                 legend.text=element_text(size=9),
                 legend.position=c(1-pctRightExpand/2,0.5),
                 legend.background=element_blank(),
                 plot.margin = unit(c(1,1,1,1), 'mm'),
                 panel.border = element_rect(color='black', linetype='solid')
            )
   #p2 = p2 + scale_y_continuous(labels=abs)
   p2 = p2 + scale_y_continuous(breaks=bb, labels=ll)
   p2 = p2 + scale_x_continuous(breaks=brks, limits=c(left, right), expand=c(0.05,0.05))
   p2 = p2 + scale_fill_manual(values=colrs)
   p2 = p2 + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))

 ################################## plot SVs (all) ###############################################
  x$svtype <- factor(x$svtype, levels=c('DUP','BND','INS','INV','DEL'))

  # manually create breaks for y axis, and change 1st to '000' so plots will align
  zz = aggregate(num.samples ~ pos, data=x, FUN=sum)
  nn = max(zz$num.samples); ii = max(round(nn/40)*10,5)
  bb = seq(0,nn,ii); ll=bb; ll[1]='000'

 
  p22 = ggplot(x, aes(x=pos, y=num.samples, fill=svtype)) + geom_bar(stat="identity")
  p22 = p22 + theme_bw() + xlab('') + ylab('Num. samples')
  p22 = p22 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  p22 = p22 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 plot.title=element_text(size=16, hjust=0.5, face="bold"), 
                 panel.background=element_rect(color="black"),
                 #axis.text.x=element_text(size=12, color="black"),
                 axis.text.x=element_blank(),
                 axis.text.y=element_text(size=12, color="black"),
                 #axis.title.x=element_text(size=14, color="black"),
                 axis.title.x=element_blank(),
                 axis.title.y=element_text(size=12, color="black"),
                 legend.key.size = unit(2.5,'mm'),
                 legend.title=element_blank(),
                 legend.text=element_text(size=9),
                 legend.position=c(1-pctRightExpand/2,0.5),
                 legend.background=element_blank(),
                 plot.margin = unit(c(1,1,1,1), 'mm'),
                 panel.border = element_rect(color='black', linetype='solid')
              )
   #p22 = p22 + scale_fill_manual(name="SV type", values=c('BND'='#2ca25f','INS'='#fec44f', 'INV'='#c994c7', 'DUP'='#b53f4d', 'DEL'='#2c7fb8'), labels=c('BND'='BND','INS'='INS', 'INV'='INV', 'DUP'='DUP', 'DEL'='DEL'), guide = guide_legend(override.aes = list(size = 7)))
   p22 = p22 + scale_fill_manual(values=colrs)
   #p22 = p22 + scale_x_continuous(labels = scales::comma, limits=c(left, right), expand=c(0.05,0.05))
   p22 = p22 + scale_x_continuous(breaks=brks, labels = xlabs , limits=c(left, right), expand=c(0.05,0.05)) + scale_y_continuous(breaks=bb, labels=ll)
   p22 = p22 + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))

  ######################################## plot chip-seq data ###########################################
  if (is.chip.avail) {
     reg.chip = chip_seq[chip_seq$chrom == pk.corr$Chr & 
                         chip_seq$pos >= left & 
                         chip_seq$pos <= right,]


     reg.chip0 <<- reg.chip
     pk.corr0 <<- pk.corr

     # down sampling chip data for plotting
     zz = binChip(reg.chip, convert.to.pct=T)
     reg.chip = zz$chip
     max.chip.cov = zz$max.cov
     chip.max.cov = sprintf('%.2g', zz$max.cov)

     p3 = ggplot (reg.chip,aes(x=pos, y=cov))# + geom_bar(stat="identity")
     p3 = p3 + geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=cov),
                         color='black', size=0.1, fill='black')
     p3 = p3 + theme_bw() + xlab('') + ylab(paste0(chip.cov.lbl))
     p3 = p3 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End),
                          color='black', linetype='dashed')
     p3 = p3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_text(size=12, color="black"),
                     axis.title.x=element_blank(),
                     axis.title.y=element_text(size=12, color="black"),
                     plot.margin = unit(c(1,1,1,1), 'mm'),
                     panel.border = element_rect(color='black', linetype='solid')
               )
     p3 = p3 + scale_x_continuous(breaks=brks, limits=c(left, right), expand=c(0.05,0.05))
     p3 = p3 + annotate('text', x=right, y=50, angle=90, size=3, hjust=0.5,vjust=1,
                        label=paste0(' \nMax = ', chip.max.cov))
     p3 = p3 + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))
     #pdf('tmpfigs/tmp.pdf'); print(p3); dev.off()


  } else {
    p3 = NULL
  }

  ################################## plot region of interest annotation ###############################################
  if (!is.null(p.roi) && nrow(p.roi) !=0 ) {
      p.roii <<- p.roi
      pk.corrr <<- pk.corr
      ll <<- left
      rr <<- right
      p.roi = p.roi[p.roi$chrom == pk.corr$Chr & p.roi$start > left & p.roi$end < right,]
  }
  if (!is.null(p.roi) && nrow(p.roi) !=0 ) {
      rtypes = unique(p.roi$roi.type)
      p4 = (ggplot(p.roi)
            + geom_segment(aes(x=start, xend=end, y=roi.type, yend=roi.type,
                               color=roi.type), size=6)
            + theme_bw()
            + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                  axis.title.y=element_text(size=12, color="black"),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_text(size=12, color="white"),
                  #legend.key.size = unit(0.1,'line'),
                  legend.title=element_blank(),
                  legend.text=element_text(size=9),
                  legend.position=c(0.8,0.8),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  plot.margin = unit(c(1,1,1,1), 'mm'),
                  panel.border = element_rect(color='black', linetype='solid')
            )
            + ylab('Annot')
            + scale_color_manual(values=c('black', 'darkgreen', 'darkred'))
            + scale_x_continuous(breaks=brks, limits=c(left, right), expand=c(0.05,0.05))
            + scale_y_discrete(breaks=rtypes[1], labels='000')
            + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
            + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))
            + guides(color = guide_legend(override.aes = list(size = 2)))
      )

  } else {
     p4 = NULL
  }

 ################################## plot gene annotation ###############################################
  # add enhancer
  #ARE = 'hg38: chrX:66896703-66927905'
  #erg = data.frame(g.chr='chr21', g.start=41464305, g.stop=41508158, gene='TMPRSS2', g.strand='-', stringsAsFactors=F)
  #ehc = data.frame(g.chr='chrX', g.start=66896703, g.stop=66927905, gene='Enhancer', g.strand='.', stringsAsFactors=F)
  #genes.in.p = rbind(genes.in.p, erg, ehc)

  if (!is.null(genes.to.show)){genes.in.p = genes.in.p[genes.in.p$gene %in% genes.to.show,]}
   
  p5 = (ggplot(genes.in.p)
        #+ geom_segment(aes(x=g.start, xend=g.stop, y=2, yend=2),
        #       color=with(genes.in.p,ifelse(g.strand=="+",'darkgreen',
        #       ifelse(g.strand=='-', 'blue', 'red'))), size=5)
        + geom_rect(aes(xmin=g.start, xmax=g.stop, ymin=0, ymax=0.8, fill=g.strand))
        + theme_bw()
        #+ geom_text(data=genes.in.p, aes(x=(g.start+g.stop)/2, y=1.25, yend=1.25,label=gene),
        #              color=with(genes.in.p, ifelse(g.strand=="+",'darkgreen',
        #              ifelse(g.strand=='-', 'blue', 'red'))), size=4, hjust=0.5)
        + geom_text(aes(x=((g.start+g.stop)/2), y=-0.15, label=gene, color=g.strand),
                    vjust=1, hjust=0.5, size=5, fontface='italic')
        + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                  axis.ticks=element_blank(),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.y=element_text(size=12, color='white'),
                  axis.title.y=element_text(size=12, color="black"),
                  plot.margin = unit(c(1,1,1,1), 'mm'),
                  #panel.border = element_rect(color='black', linetype='solid'),
                  panel.border = element_blank(),
                  legend.position='none'
        )
        + ylab('')
        + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
        + scale_y_continuous(breaks=c(1), labels=c('000'), limits=c(-1, 1))
        + scale_size_identity()
        + scale_color_manual(values=c('+'='darkblue', '-'='darkgreen', '.'='black'))
        + scale_fill_manual(values=c('+'='darkblue', '-'='darkgreen', '.'='black'))
        + coord_cartesian(xlim=c(left, right + (right-left)*pctRightExpand))
  )
 
 
  ####### make plots list and return the result 
  plots <- list(p0,p1,p2,p22,p3,p4,p5)
  #names(plots) <- c( "p0", "p1","p2","p22","p3",p4,"p5")
  plots <- plots[!sapply(plots, is.null)]
  return (plots)
}
############################################# END OF PLOT REGION FUNCTION ###################################################################

### read structural variants file 
cat('Reading all breakpoints...')
if (file.exists((sv.file))) {
  sv <- read.table(sv.file, header =T, sep="\t", stringsAsFactors = F, check.names=F)
  sv$sample = sub('/.*$', '', sv$name)
  sv$svtype = sub('^.*/', '', sv$name)
  sv = sv[sv$svtype %in% c('DUP', 'DEL'),]
  sv$pos1 = (sv$start1 + sv$end1)/2
  sv$pos2 = (sv$start2 + sv$end2)/2
} else {
  stop('structural variants file was not found!')
}

### read all break points
bp = read.table(paste0(res.dir,'/processed_data/all_bp.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(bp) = c('chr', 'start', 'stop', 'name', 'score', 'strand')
bp$pos = (bp$start+bp$stop)/2
bp$sample = gsub('/.*$', '', bp$name)
bp$svtype = gsub('^.*/', '', bp$name)
### extract total number of samples 
samples.with.sv <- unique(bp$sample)
cat('done.\n')

### read windows counts file 
cat('Reading sliding window sample count...')
if (file.exists(paste0(res.dir,'/processed_data/counts.rds')) ){
  cts = readRDS(paste0(res.dir,'/processed_data/counts.rds'))
} else {
  stop (paste("File \"",res.dir, "/counts.rds\" was not found!.\n", sep=""))
}
cat('done.\n')

### read annotated peaks summary file 
res = read.table(paste0(res.dir, '/annotated_peaks_summary.tsv'), header=T, sep='\t', stringsAsFactors=F)
### filter peaks with effected genes only 
res = res[!is.na(res$Associated.genes), ]

### read copy number file 
cat('Reading copy number data for genes and peaks...')
if (file.exists((cn.file))) {
  cn_data <- read.table(cn.file, header =T, check.names=F, sep="\t", stringsAsFactors = F, comment.char="")
  cn_data$pos = (cn_data$start+cn_data$end)/2
  ### check if cn.call is exists 
  if (!"cn.call" %in% colnames(cn_data)) {
    cn_data$cn.call <- "neut"
    cn_data[cn_data$cn > t.amp, 'cn.call'] <- "amp"
    cn_data[cn_data$cn < t.del, 'cn.call'] <- "del"
  }
} else {
  stop(paste0('Copy number file \"', cn.file, '\" was not found!'))
}
#cat('done.\n')

### read peaks with copy number data
if (file.exists(paste0(res.dir,'/processed_data/peaks_with_cn.bed'))) {
  pks.cn <- read.table(paste0(res.dir,'/processed_data/peaks_with_cn.bed'), header =F, check.names=F, sep="\t", stringsAsFactors = F)
  colnames(pks.cn) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'samples',
                        'cn.chr', 'cn.start', 'cn.stop', 'sample', 'seg.cn', 'cn.call','dist', 'cn.value')
} else {
  stop (paste0("Peaks copy number file \"", res.dir, "/processed_data/peaks_with_cn.bed\" was not found!.\n"))
}
#cat('done.\n')

### read genes with copy number data
if (file.exists(paste0(res.dir,'/processed_data/genes_with_cn.bed'))) {
  genes.cn <- read.table(paste0(res.dir,'/processed_data/genes_with_cn.bed'), header =F, check.names=F, sep="\t", stringsAsFactors = F)
  colnames(genes.cn) <- c('g.chr', 'g.start', 'g.stop', 'gene', 'score', 'strand', 'cn.chr', 'cn.start', 'cn.stop', 
                          'sample', 'seg.cn', 'cn.call','dist', 'cn.value')
} else {
  stop (paste0("Genes copy number file \"", res.dir, "/processed_data/genes_with_cn.bed\" was not found!.\n"))
}
cat('done.\n')

### read chip-seq coverage file
#chip.seq = 0
if (chip.seq != 0) {
  is.chip.avail = TRUE 
  chiph = 3 # previous value = 2
  cat('Reading chip-seq coverage data...') 
  if (file.exists(chip.seq) ) {
    chip_seq <- read.table(chip.seq, header =T, sep="\t", stringsAsFactors = F, check.names=F)
    chip_seq$pos = (chip_seq$start+chip_seq$end)/2
  } else {
    stop (paste0("Avergaed chip coverage file \"", chip.seq, "\" was not found!.\n"))
  }
  cat('done.\n')
} else {
  is.chip.avail = FALSE
  chiph = NULL
}

### read annotation file  
annot <- read.table(paste0(res.dir,"/processed_data/genes.bed"), header =T, sep="\t", check.names=F, comment.char = "$")
colnames(annot) <- c('chr', 'start', 'stop', 'gene', 'score', 'strand') 
#feature.col <- colnames(annot)[4]

### read raw peak calls
cat('Reading peaks with breakpoints...')
pks = read.table(paste0(res.dir, '/processed_data/peaks_overlap_bp.tsv'), header=F, stringsAsFactors=F, sep='\t')
colnames(pks) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample', 'sv.type')
cat('done.\n')

### read gene overlap/nearby peaks 
genes.and.peaks <- read.table(paste0(res.dir,'/processed_data/peaks_with_overlap_nearby_genes.tsv'), header =F, sep="\t", stringsAsFactors=F, check.names=F)
colnames(genes.and.peaks) = c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples','sample', 
                              'g.chr', 'g.start', 'g.stop', 'gene', 'g.score', 'g.strand', 'dist','g.pos')

### read expression data
cat('Reading expression data...')
if (file.exists((exp.file))) {
  exp <- read.table(exp.file, header =T, check.names=F, sep="\t", stringsAsFactors = F)
  exp.data.cols <- colnames(exp)[-1]
} else {
  stop (paste0("Expression file \"", exp.file, "\" was not found!.\n"))
}
cat('done.\n')

### extract samples that have expression but do have SVs
samples.with.no.SVs <- unique(colnames(exp)[!exp.data.cols %in% samples.with.sv])

#### read region of interest file(s) 
if (roi.file !=0) {
  cat('Reading region(s) of interest data...')
  is.roi.avail = TRUE
  roi.files = unlist(strsplit(roi.file, ","))
  roi.col.names <- paste0('Overlapped.', gsub("\\..*", "", basename(roi.files)))
  for (j in 1:length(roi.files)){
    roi.ff <- roi.files[j]
    roi.name <- gsub("\\..*", "", basename(roi.ff))
    
    if (file.exists(roi.ff)) {
      roi = read.table(roi.ff, header =T, sep="\t", stringsAsFactors=F)
      roi = roi[,1:4]
    } else {
      stop (paste0("Region of interest file \"", roi.ff, "\" was not found!.\n"))
    }
    assign(paste(roi.name,"roi", sep = "."), roi)
  }
  cat('done.\n')
} else {
  is.roi.avail = FALSE
}

######################################### plot peaks #############################################
#### create directory for plots 
dir.create(paste0(out.dir, '/peaks-plots'), recursive=T, showWarnings = FALSE) 

### extract peaks
pks.to.plot = unlist(strsplit(pks.to.plot,  ","))

### check if the availability of peaks  
if (!any(pks.to.plot %in% res$Peak.name)) {
  stop('None of the peaks were found the in the results file.\n') 
}
cat('Peaks to plot:', pks.to.plot[pks.to.plot %in% res$Peak.name],'\n')
if (length(pks.to.plot[!pks.to.plot %in% res$Peak.name]) !=0) {
  cat('Peaks not found in the results file:', pks.to.plot[!pks.to.plot %in% res$Peak.name], '\n')
}

for (k in 1:length(pks.to.plot)) {
  pk = pks.to.plot[k]
  cat('\n','Plotting peak', pk, '\n')
  ### extract peak coordinates
  p.corr <- res[res$Peak.name==pk, c('Chr','Start','End', 'Percentage.SV.samples', 'Number.SV.samples')]
  
  ### extract effected genes 
  assoc.genes <- unlist(strsplit(res[res$Peak.name==pk, 'Associated.genes'],  "\\|"))
  ### keep genes that have expression data
  assoc.genes <-  assoc.genes[assoc.genes %in% exp[,1]]
  
  if (length(assoc.genes) ==0 ) { stop('No genes found to be associated with this peak.') }
  
  ### extract locus information for effected genes 
  genes.in.peak <- unique(genes.and.peaks[genes.and.peaks$p.name==pk,c('g.chr','g.start','g.stop','gene','g.strand')])

  ### extract peak data
  pp = pks[pks$p.name==pk & pks$sample %in% samples.with.sv, ]
  pp.cn = pks.cn[pks.cn$p.name==pk & pks.cn$sample %in% samples.with.sv, c('p.name','sample', 'cn.value', 'cn.call')]

  #### include samples with no SVs as neutral samples 
  if (length(samples.with.no.SVs) > 0) {
    pp.cn = rbind(pp.cn, data.frame(p.name=pk, sample=samples.with.no.SVs, cn.value=0, cn.call="neut"))
  }
  
  ### extract peak copy number samples
  pk.neut.samples <- unique(pp.cn[pp.cn$cn.call=="neut",'sample'])
  pk.amp.samples <-  unique(pp.cn[pp.cn$cn.call=="amp",'sample'])
  pk.del.samples <-  unique(pp.cn[pp.cn$cn.call=="del",'sample'])
  
  ### extract sv/other samples 
  sv.pats <- unique(pp$sample)
  DUP.pats <- unique(pp[pp$sv.type=="DUP", 'sample'])    
  BND.pats <- unique(pp[pp$sv.type=="BND", 'sample'])    
  INS.pats <- unique(pp[pp$sv.type=="INS", 'sample'])    
  DEL.pats <- unique(pp[pp$sv.type=="DEL", 'sample'])
  INV.pats <- unique(pp[pp$sv.type=="INV", 'sample'])    
  nonSV.pats <- c(unique(samples.with.sv[!samples.with.sv %in% sv.pats]), samples.with.no.SVs) 
  
  ### extract region of interest results for current peak 
  if (is.roi.avail) {
    #p.roi.res <- res[res$Peak.name==pk, 'overlap.roi']
    p.roi.res <- as.data.frame(res[res$Peak.name==pk, roi.col.names])
    colnames(p.roi.res) = roi.col.names
    roi.annot <- NULL
    for (k in roi.col.names) {
      roi.f <- eval(parse(text = paste0(gsub("Overlapped.", "", k), ".roi")))
      roi.r <- roi.f[roi.f$name %in% unlist(strsplit(as.character(p.roi.res[,k]), "\\|")),]
      #disable ROI overlapping peaks
      roi.r <- roi.f
      if (nrow(roi.r) !=0) {
      	  roi.r$roi.type <- gsub("Overlapped.", "", k)
      	  roi.annot <- rbind(roi.annot, roi.r)
      }
    }

    #p.roi.res <- roi[roi$name %in% unlist(strsplit(p.roi.res, "\\|")),]
    if (nrow(roi.annot) == 0 || is.null(roi.annot)) { roih = NULL } else { roih = length(unique(roi.annot$roi.type))+0.5 }
  } else {
    roi.annot <- NULL
    roih = NULL
  }
 
  ### loop through genes in the peak 
  p.genes.res <- NULL
  for (j in 1:length(assoc.genes)) {
    g = assoc.genes[j]

    # check if genes to show were provided
    if (!is.null(genes.to.show)) {
      if (!is.null(genes.to.show) && !(g %in% genes.to.show)) {next}
    } 
    
    ### extract gene copy number samples  
    g.neut.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="neut",'sample'])
    g.amp.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="amp",'sample'])
    g.del.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="del",'sample'])
    
    ### extract expression 
    g.exp <- as.data.frame(t(exp[exp[,1]==g, exp.data.cols]))
    g.exp$sample <- rownames(g.exp)
    rownames(g.exp) <- NULL
    colnames(g.exp) <- c('gene.exp', 'sample')
    g.exp <- g.exp[,c('sample', 'gene.exp')]
    
    #### a dd sample status column
    g.exp$sample.status <- 'non-SVs'
    g.exp[g.exp$sample %in% sv.pats, 'sample.status'] <- 'SVs'
    
    #### add gene copy number status 
    g.exp$gene.cn.status <- "neut"
    g.exp[g.exp$sample %in% g.amp.samples, "gene.cn.status"] <- "amp"
    g.exp[g.exp$sample %in% g.del.samples, "gene.cn.status"] <- "del"
    
    #### add peak copy number status 
    g.exp$pk.cn.status <- "neut"
    g.exp[g.exp$sample %in% pk.amp.samples, "pk.cn.status"] <- "amp"
    g.exp[g.exp$sample %in% pk.del.samples, "pk.cn.status"] <- "del"
    
    ############### run the function to plot expression  #################
    plot.ex <- plot.exp(g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, g, pk)
    #plot.exp.amp.del(g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, paste0(out.dir,'/peaks-plots/',g,'_',pk), g, pk)
    
    ############### run the function to plot the peak regin #################
    p.reg = suppressWarnings( plot.region(pk, p.corr, g, genes.in.peak, roi.annot) )
    n = length(p.reg) 
    
    p.reg[[length(p.reg)+1]] = plot.ex[["e1"]]
    p.reg[[length(p.reg)+1]] = plot.ex[["e2"]]
    p.reg[[length(p.reg)+1]] = plot.ex[["e3"]]
    
    ### align all plots 
    all.plots <- suppressWarnings( do.call(AlignPlots, p.reg) )

    ### set the layout matrix (version 1 with three columns)
    #mat = matrix(ncol=3, nrow=n+1)
    #mat[, 1] = 1:(n+1)
    #mat[, 2] = c(1:n, n+2)
    #mat[, 3] = c(1:n, n+3)
    ### set the layout matrix (version 1 with two columns)
    mat = matrix(ncol=2, nrow=n)
    mat[, 1] = 1:(n)
    mat[, 2] = c(rep(n+1,2), rep(n+2,2), rep(n+3, n-4))


    ### set the height 
    cnh =3 #4
    ddh =3 #4
    svh1 =3 #4
    svh2 =3 #4
    geneh = 1 
    #boxp = 5
    myheights = c(cnh, ddh, svh1,svh2, chiph, roih, geneh)
    mywidths = c(3, plot.ex[["ColWidth"]])
    #mywidths = c(1.5, plot.ex[["w1"]], plot.ex[["w2"]])
    #mywidths[2] = mywidths[2]*1.1
    plot.height = sum(myheights)*0.60
    plot.width = 8 
    plot.ncols = 2

    # layout for expression plot at bottom (wide format)
    if (!is.na(layout) && layout == 'wide'){
      mat = matrix(ncol=3, nrow=(n+1))
      mat[1:n,] = matrix(rep(1:n,3), ncol=3)
      mat[n+1,] = seq(n+1,n+3)
      #myheights = c(myheights, 8)
      myheights = c(4,4,4,4,2, roih, geneh, 8)
      mywidths = c(2,2,2)
      p.reg[[n+2]] = p.reg[[n+2]] + ylab(NULL)
      p.reg[[n+3]] = p.reg[[n+3]] + ylab(NULL)
      # fit a letter size paper
      plot.width = 12 #8
      plot.height = 11
      plot.ncols = 3
    }


    #### plot all 
    g.title = paste0('Peak locus: ',p.corr$Chr, ':',  p.corr$Start, '-',  p.corr$End,'\n')
    mytitle=textGrob(g.title, gp=gpar(fontsize=20,fontface="bold"))

    pdf(paste0(out.dir,'/peaks-plots/',g,'_',pk,'.pdf'),
        width=plot.width, height=plot.height, title='', useDingbats=F, onefile=FALSE)
    #grid.arrange(grobs=all.plots, nrow=n,ncol=2, layout_matrix=mat, heights = myheights, widths=mywidths)#, top = mytitle)
    grid.arrange(arrangeGrob(grobs=p.reg, ncol=plot.ncols,
                  layout_matrix=mat, heights = myheights, widths=mywidths))
    dev.off()
    
  }  ### end of genes in the peak 
  
}    ### end of peaks 

### remove temporary files 
unlink(paste0(res.dir, "/processed_data/win_cn.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.win.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.cn.tsv"))   


