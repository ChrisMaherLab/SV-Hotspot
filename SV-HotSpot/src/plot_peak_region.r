#!/usr/bin/Rscript

# Find regions/peaks whose SVs altered expression of nearby genes
# Written by Abdallah Eteleeb & Ha Dang

library(ggplot2)
library(reshape2)
library(grid)
#library(gridBase)
library(gridExtra)
library(gtable)
library(ggsignif)

args = commandArgs(T)

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

################################# FUNCTION TO PLOT SVs EXPRESSION #########################################
plot.exp <- function (g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, gene, pk) {
  
  ##################### plot the expression for SV samples (SVs vs non-SVs) ###############################
  pval = (wilcox.test(g.exp[g.exp$sample.status=="SVs", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp']))$p.value
  pval= sprintf(pval, fmt="%#.5f")
  g.exp$sample.status <- factor(g.exp$sample.status, levels=c('non-SVs', 'SVs'))
  e1 <- ggplot(g.exp, aes(x=sample.status, y=log2(gene.exp+1))) + geom_boxplot(aes(fill=sample.status)) + theme_bw() 
  e1 <- e1 + labs(x='', y='Log2(expression)') + ggtitle(paste0('\n',gene, ' expression in SV\nand non-SV samples'))
  e1 <- e1 + theme(axis.text.x=element_text(size=12, vjust=0.5, color="black"),
                   axis.text.y=element_text(size=14, color="black"), 
                   axis.title.y=element_text(size=16), panel.background=element_blank(),
                   plot.title = element_text(size = 16, hjust=0.5, color="black", face="plain"),
                   legend.position="none")
  e1 = e1 + scale_fill_manual(name="", values =c("non-SVs"="gray", "SVs"="orange2"))
  e1 = e1 + scale_x_discrete(labels=c(paste0("non-SVs\n(n=",nrow(g.exp[g.exp$sample.status=="non-SVs",]),")"), paste0("SVs\n(n=",nrow( g.exp[g.exp$sample.status=="SVs",]),")")))
  #e1 = e1 + geom_signif(comparisons=list(c('non-SVs','SVs')))
  e1 = e1 + annotate("text", x = 0.8, y = max(log2(g.exp$gene.exp+1)), label = paste0('p=', pval), cex=4, fontface="bold")
  
  ###################### plot the expression for SV samples (SV types vs non-SVs) ########################
  if (length(g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="neut", 'gene.exp']) > 0 ) {
    nonSVs.exp = data.frame(grp = "non-SVs", exp = g.exp[g.exp$sample %in% nonSV.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])
    if (length(BND.pats) > 0 & length(g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {BND.exp.n = data.frame(grp = "BND", exp = g.exp[g.exp$sample %in% BND.pats & g.exp$gene.cn.status=="neut" , 'gene.exp'])} else {BND.exp.n=NULL}
    if (length(DUP.pats) > 0 & length(g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {DUP.exp.n = data.frame(grp = "DUP", exp = g.exp[g.exp$sample %in% DUP.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {DUP.exp.n=NULL}
    if (length(INS.pats) > 0 & length(g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {INS.exp.n = data.frame(grp = "INS", exp = g.exp[g.exp$sample %in% INS.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {INS.exp.n=NULL}
    if (length(DEL.pats) > 0 & length(g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {DEL.exp.n = data.frame(grp = "DEL", exp = g.exp[g.exp$sample %in% DEL.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {DEL.exp.n=NULL}
    if (length(INV.pats) > 0 & length(g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="neut" , 'gene.exp']) > 0) {INV.exp.n = data.frame(grp = "INV", exp = g.exp[g.exp$sample %in% INV.pats & g.exp$gene.cn.status=="neut", 'gene.exp'])} else {INV.exp.n=NULL}
    
    svtype.exp = do.call('rbind', list(nonSVs.exp, BND.exp.n, DUP.exp.n, INS.exp.n, DEL.exp.n, INV.exp.n))
    svtype.exp$grp <- factor(svtype.exp$grp, levels=c('non-SVs', 'BND','DUP', 'INS','DEL', 'INV'))
  } else {
    svtype.exp <- NULL
  }
 
  
  ######### plot SVs types for neutral samples only ######
  if (!is.null(svtype.exp)) {
    ### prepare comarison list
    sv.grps.n = as.data.frame(combn(as.character(unique(svtype.exp$grp)),2), stringsAsFactors = F)
    idxs <- as.data.frame(which(sv.grps.n =="non-SVs", arr.ind=TRUE))
    sv.grps.n <- sv.grps.n[,idxs$col]
    svCMPlist = as.list(sv.grps.n[,1:ncol(sv.grps.n)])
    
    title.size1 = length(levels(factor(svtype.exp$grp))) * 3.5
    if (title.size1 <= 7) { title.size = 12 }
    
    lbls.n = c(paste0("non-SVs\n(n=",nrow(nonSVs.exp),")"), paste0("BND\n(n=",nrow(BND.exp.n),")"), 
               paste0("DUP\n(n=",nrow(DUP.exp.n),")"), paste0("INS\n(n=",nrow(INS.exp.n),")"),
               paste0("DEL\n(n=",nrow(DEL.exp.n),")"), paste0("INV\n(n=",nrow(INV.exp.n),")"))
    lbls.n = grep(paste(as.character(unique(svtype.exp$grp)), collapse="|"), lbls.n, value=TRUE)
    neut.title = paste0(gene, ' expression by SV type in peak ', pk,'\n(CN neutral samples only)')
    
    e2 <- ggplot(svtype.exp, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp)) + theme_bw() 
    e2 <- e2 + labs(x='', y='Log2(expression)') + ggtitle(neut.title)
    e2 <- e2 + theme(axis.text.x=element_text(size=12, vjust=0.5, color="black"),
                     axis.text.y=element_text(size=14, color="black"), 
                     axis.title.y=element_text(size=16), panel.background=element_blank(),
                     plot.title = element_text(size = 16, hjust=0.5, color="black", face="plain"),
                     legend.position="none")
    e2 = e2 + scale_fill_manual(name="", values =c("non-SVs"="gray", "BND"="#2ca25f", "DUP"="#b53f4d", "INS"="#fec44f","DEL"="#2c7fb8","INV"="#c994c7"))
    e2 = e2 + scale_x_discrete(labels=lbls.n)
    e2 = e2 + geom_signif(comparisons= svCMPlist, step_increase=0.1)
  } else {
    e2 <- NULL
  }
  
  
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
    
 #   lbls = c("GNPN"=paste0("GeneNeut\npeakNeut\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="neut", "group"]),")"), 
 #            "GAPN"=paste0("GeneAmp\nPeakNeut\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="neut", "group"]),")"),
 #            "GNPA"=paste0("GeneNeut\nPeakAmp\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="amp", "group"]),")"),
 #            "GAPA"=paste0("GeneAmp\nPeakAmp\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="amp", "group"]),")"), 
 #            "GDPN"=paste0("GeneDel\nPeakNeut\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="neut", "group"]),")"),
 #            "GNPD"=paste0("GeneNeut\nPeakDel\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="del", "group"]),")"),
 #            "GDPD"=paste0("GeneDel\nPeakDel\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="del", "group"]),")"), 
 #            "NC"=paste0("Other\n(n=",nrow(g.exp[g.exp$group=="NC",])))
    
 
    lbls = c("GNPN"=paste0("GN/PN\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="neut", "group"]),")"), 
             "GAPN"=paste0("GA/PN\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="neut", "group"]),")"),
             "GNPA"=paste0("GN/PA\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="amp", "group"]),")"),
             "GAPA"=paste0("GA/PA\n(n=",length(g.exp[g.exp$gene.cn.status=="amp" & g.exp$pk.cn.status=="amp", "group"]),")"), 
             "GDPN"=paste0("GD/PN\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="neut", "group"]),")"),
             "GNPD"=paste0("GN/PD\n(n=",length(g.exp[g.exp$gene.cn.status=="neut" & g.exp$pk.cn.status=="del", "group"]),")"),
             "GDPD"=paste0("GD/PD\n(n=",length(g.exp[g.exp$gene.cn.status=="del" & g.exp$pk.cn.status=="del", "group"]),")"), 
             "NC"=paste0("NC\n(n=",nrow(g.exp[g.exp$group=="NC",])))

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
    
    tt = paste0(gene,' expression with the presence and\nabsence of CN at ', gene, '/peak regions\n(G=gene, P=Peak, N,A,D=CN events)')
    
    e3 <- ggplot(g.exp, aes(x=group, y=log2(gene.exp+1))) + geom_boxplot(aes(fill=group)) + theme_bw()
    #e3 <- e3 + stat_summary(fun.data = give.n, geom = "text", size=5) 
    e3 <- e3 + labs(x='', y='Log2(expression)') + ggtitle(tt)
    e3 <- e3 + theme(axis.text.x=element_text(size=12, vjust=0.5, color="black"),
                     axis.text.y=element_text(size=14, color="black"),
                     axis.title=element_text(size=16), 
                     panel.background=element_blank(),
                     plot.title = element_text(size = 16, hjust=0.5, color="black", face="plain"),
                   legend.position="none")
    e3 = e3 + scale_fill_manual(name="", values = c("GNPN"="gray", "GAPN"="#f03b20", "GNPA"="#b53f4d", "GAPA"="salmon", 
                                                    "GDPN"="#a6bddb", "GNPD"="#2c7fb8", "GDPD"="skyblue2")) 
    e3 = e3 + scale_x_discrete(labels=lbls)
    e3 = e3 + geom_signif(comparisons=myCMPlist , step_increase=0.1)
  #}    
  
  ### preapre and return results 
  width1 = length(levels(svtype.exp$grp))*0.5
  width2 = length(levels(factor(g.exp$group)))*0.5
  if (width1==1) { width1 =1.5}
  if (width2==1) { width2 =1.5}
  exp.plots <- list(e1,e2,e3, width1, width2)
  names(exp.plots) <- c( "e1", "e2","e3", "w1", "w2")
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
    
    a <- ggplot(svtype.exp.a, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp)) + theme_bw() 
    a <- a + labs(x='', y='Log2(expression)') + ggtitle(amp.title)
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
                       
    d <- ggplot(svtype.exp.d, aes(x=grp, y=log2(exp+1))) + geom_boxplot(aes(fill=grp)) + theme_bw() 
    d <- d + labs(x='', y='Log2(expression)') + ggtitle(del.title)
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
   left <- left - left.ext - 10000     
   right <- right + right.ext  +10000
   D = right - left
   #r.width = pk.corr$End - pk.corr$Start
   #scale binwidth accordingly based on region width
   binwidth = D/75
   #genes within region
   g.corr = genes.in.p[genes.in.p$gene ==gene, ]
 
   ### extract SVs data
   #x = bp[bp$chr == pk.corr$p.chr & bp$pos > left & bp$pos < right,]
   x = cts[cts$chr == pk.corr$Chr & cts$pos > left & cts$pos < right,]  ### using counts 
   x = x[x$svtype !="ALL",]
   
   ### for DUP and DEL only 
   x2 = x 
   x2 = x2[x2$svtype  %in% c("DUP","DEL"),]
   x2[x2$svtype == "DEL", "num.samples"] <-  x2[x2$svtype == "DEL", "num.samples"] * -1

   ### make the title 
   title = paste0('Associated Gene: ',gene,'\n (Peak locus: ',pk.corr$Chr, ':',  pk.corr$Start, '-',  pk.corr$End, '; Peak name=', pk,'; Peak width=', width,'kb)')

   ### compute the DUP and DEL pileup of SV
   dup_del = pileUp(sv, pk.corr$Chr, left, right)
   
   dup_del$pos1[dup_del$pos1 < left] = left
   dup_del$pos2[dup_del$pos2 > right] = right
 
  ################################## plot region copy number ###############################################
  #if (is.cn.avail) {
     reg.cn = cn_data[cn_data$chr == pk.corr$Chr & (cn_data$pos > left | cn_data$pos < right),]
     reg.cn = reg.cn[tolower(reg.cn$cn.call) %in% c("amp", "del"), ]
     #reg.width = pk.corr$p.stop - pk.corr$p.start + 1

     reg.width = (right - left)+1
     s = round(reg.width/500)
     imin = min(reg.cn$start)
     imax = max(reg.cn$end)
     chr.name=unique(reg.cn$chr)
     reg.win = data.frame(chr=chr.name, start=seq(imin,imax,s))
     reg.win$stop = reg.win$start + s - 1
     write.table(reg.cn, file=paste0(res.dir,'/processed_data/reg.cn.tsv'), quote=F, row.names=F, sep="\t", col.names=F)
     write.table(reg.win, file=paste0(res.dir,'/processed_data/reg.win.tsv'), quote=F, row.names=F, sep="\t", col.names=F)

     ### overlap windows with copy number 
     #cat ("Overlapping windows with copy number ..")
     system(paste0("intersectBed -wo -a ",res.dir,"/processed_data/reg.win.tsv -b ", res.dir,"/processed_data/reg.cn.tsv | cut -f2,3,7,9 | sort | uniq | sort -k1,1 -k2,2 | groupBy -full -g 1,2 -c 4 -o count | cut -f2-6 > ", res.dir, "/processed_data/win_cn.tsv")) 

     win.data = read.table(paste0(res.dir,"/processed_data/win_cn.tsv"), header = F, sep ="\t")
     colnames(win.data) = c('start', 'stop','sample', 'cn.call', 'num.samples')
     win.data$pos = (win.data$start + win.data$stop)/2
     ### add dummy data for visualization purpuses 
     win.data=rbind(win.data, data.frame(start=-1,stop=-2,sample="dummp1",cn.call="amp",num.samples=0, pos=-1))
     win.data=rbind(win.data, data.frame(start=-1,stop=-2,sample="dummp1",cn.call="del",num.samples=0, pos=-1))
     
     p0 = ggplot(win.data, aes(x=pos,y=num.samples, fill=cn.call)) + geom_bar(stat="identity")
     p0 = p0 + scale_fill_manual(name="", values=c("amp"="#b53f4d", "del"="#2c7fb8"), labels=c("amp"="Gain", "del"="Loss"), drop=FALSE,
                                 guide = guide_legend(override.aes = list(size = 7)))
     p0 = p0 + theme_bw() + xlab('') + ylab('Copy number\nalterations') + ggtitle(title)
     p0 = p0 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                   plot.title=element_text(size=16, hjust=0.5, face="bold"), axis.ticks = element_blank(),
                   axis.text.x=element_blank(), axis.text.y=element_text(size=12, color="black"),
                   axis.title.x=element_text(size=14, color="black"), axis.title.y=element_text(size=12, color="black"),
                   legend.key.size = unit(0.8,"cm"), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))
     #p0 = p0 + xlim(left, right )
     p0 = p0 + scale_x_continuous(limits=c(left, right), expand=c(0.02,0.02))
     p0 = p0 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
     #p0 = p0 + geom_vline(xintercept=67711690, color='black', linetype='dashed', size=0.3)
     #p1 = p1 + scale_y_continuous(breaks = seq(min(reg.cn$cn), max(reg.cn$cn),10))
  #}

  ################################## Plot DUP $ DEL Freq ####################################
  #dup_del = pileUp(sv, pk.corr$p.chr, left, right)  
  #prepare breaks and labs in Mb
   brks = seq(0,10^9,10^6)
   brks = brks[brks >= left & brks <= right]
   labs = brks/(10^6)
   
   dup_del$svtype <- factor(dup_del$svtype, levels=c('DUP', 'DEL'))

    p1 = (ggplot(dup_del)
    + geom_segment(aes(x=pos1, xend=pos2, y=samp, yend=samp, color=svtype))
    + theme_bw(base_size=8)
    + coord_cartesian(xlim=c(left, right))
    + ylab('Duplication &\ndeletion events')
    + scale_x_continuous(breaks=brks, expand=c(0.02,0.02))
    + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.title=element_text(size=16, hjust=0.5, face="bold"), axis.ticks = element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(), legend.key=element_rect(fill=NA),
            axis.title.x=element_blank(), axis.title.y=element_text(size=12, color="black"),
            legend.key.size = unit(0.8,"cm"), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))
    + scale_color_manual(name="", values=c('DUP'='#b53f4d', 'DEL'='#2c7fb8'), guide = guide_legend(override.aes = list(size = 7)))
    + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
    #+ geom_vline(xintercept=67711690, color='black', linetype='dashed', size=0.3)
    #+ geom_rect(xmin=g.corr$g.start, xmax=g.corr$g.stop, ymin=0, ymax=30, color="white", alpha=0.005)
    #+ xlim(left, right)
    )
  
  ################################## plot SVs (DUP and DEL only) ###############################################
  x2$svtype <- factor(x2$svtype, levels=c('DUP','BND','INS','INV','DEL'))
  p2 = ggplot(x2, aes(x=pos, y=num.samples, fill=svtype)) + geom_bar(stat="identity")
  p2 = p2 + theme_bw() + xlab('') + ylab('Number of\nsamples')
  #p2 = p2 + geom_text(label=ifelse(abs(x2$num.samples) >=9 & x2$svtype=="DEL", x2$sample,''), align=90)
  p2 = p2 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  #p2 = p2 + geom_vline(xintercept=67711690, color='black', linetype='dashed', size=0.3)
  p2 = p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 plot.title=element_text(size=16, hjust=0.5, face="bold"), axis.ticks = element_blank(),
                 axis.text.x=element_blank(), axis.text.y=element_text(size=12, color="black"),
                 axis.title.x=element_text(size=14, color="black"), axis.title.y=element_text(size=12, color="black"),
                 legend.key.size = unit(0.8,"cm"), legend.title=element_blank(), legend.text=element_text(size=12))
   #p2 = p2 + scale_y_continuous(breaks = sv.brks, labels = sv.lbls, limits=c(min(x2$num.samples), max(x2$num.samples)))
   p2 = p2 + scale_y_continuous(labels=abs)
   #p2 = p2 + xlim(left, right )
   p2 = p2 + scale_x_continuous(limits=c(left, right), expand=c(0.02,0.02))
   p2 = p2 + scale_fill_manual(name="SV type", values=c('BND'='#2ca25f','INS'='#fec44f', 'INV'='#c994c7', 'DUP'='#b53f4d', 'DEL'='#2c7fb8'), 
                                        labels=c('BND'='BND','INS'='INS', 'INV'='INV', 'DUP'='DUP', 'DEL'='DEL'),
                                        guide = guide_legend(override.aes = list(size = 7)))

 ################################## plot SVs (all) ###############################################
  x$svtype <- factor(x$svtype, levels=c('DUP','BND','INS','INV','DEL'))
 
  p22 = ggplot(x, aes(x=pos, y=num.samples, fill=svtype)) + geom_bar(stat="identity")
  p22 = p22 + theme_bw() + xlab('') + ylab('Number of\nsamples')
  p22 = p22 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  p22 = p22 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 plot.title=element_text(size=16, hjust=0.5, face="bold"), 
                 #axis.text.x=element_blank(), axis.text.y=element_text(size=12, color="black"),
                 axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"),
                 axis.title.x=element_text(size=14, color="black"), axis.title.y=element_text(size=12, color="black"),
                 legend.key.size = unit(0.8,"cm"), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))
   p22 = p22 + scale_fill_manual(name="SV type", values=c('BND'='#2ca25f','INS'='#fec44f', 'INV'='#c994c7', 'DUP'='#b53f4d', 'DEL'='#2c7fb8'), 
                                        labels=c('BND'='BND','INS'='INS', 'INV'='INV', 'DUP'='DUP', 'DEL'='DEL'),
                                        guide = guide_legend(override.aes = list(size = 7)))
   p22 = p22 + scale_x_continuous(labels = scales::comma, limits=c(left, right), expand=c(0.02,0.02))

  ######################################## plot chip-seq data ###########################################
  if (is.chip.avail) {
     reg.chip = chip_seq[chip_seq$chrom == pk.corr$Chr & chip_seq$pos > left | chip_seq$pos < right,]

     p3 = ggplot (reg.chip,aes(x=pos, y=cov)) + geom_bar(stat="identity", width = D/300)
     p3 = p3 + theme_bw() + xlab('') + ylab(chip.cov.lbl)
     #p3 = p3 + geom_rect(xmin=g.corr$g.start, xmax=g.corr$g.stop, ymin=0, ymax=30, color="white", alpha=0.005)
     p3 = p3 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
     #p3 = p3 + geom_vline(xintercept=67711690, color='black', linetype='dashed', size=0.3)
     p3 = p3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks = element_blank(),
                  #axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"),
                  axis.text.x=element_blank(), axis.text.y=element_text(size=12, color="black"),
                  axis.title.x=element_text(size=14, color="black"), axis.title.y=element_text(size=12, color="black"))
     #p3 = p3 + xlim(min(x$pos), max(x$pos) )
     #p3 = p3 + xlim(left, right )
     p3 = p3 + scale_x_continuous(limits=c(left, right), expand=c(0.02,0.02))
  } else {
    p3 = NULL
  }


  ################################## plot region of interest annotation ###############################################
  if (!is.null(p.roi) && nrow(p.roi) !=0 ) {
      p4 = ggplot(p.roi) + geom_segment(aes(x=start, xend=end, y=roi.type, yend=roi.type), color='black', size=6) + theme_bw(base_size = 12)
      p4 = p4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                      axis.ticks=element_blank(),
                      axis.text.y=element_text(hjust=0.5, size=12, angle=ifelse(length(unique(p.roi$roi.type)) >1 ,45,90)),
                      axis.title.y=element_text(size=12, color="black"),
                      axis.title.x=element_blank(),axis.text.x=element_blank(),
                      panel.background=element_blank(), panel.border=element_blank(),
                      axis.line = element_line(colour = "black")
                      ) 
  #   p4 = ggplot(p.roi) + geom_segment(aes(x=start, xend=end, y=0, yend=0), color='black', size=6) + theme_bw(base_size = 12)
  #   p4 = p4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
  #                   axis.ticks=element_blank(), 
  #                   axis.title.y=element_text(size=12, color="black"),
  #                   panel.background=element_blank(), panel.border=element_blank(),
  #                   axis.line = element_line(colour = "black"),
  #                   axis.text=element_blank(), axis.title.x=element_blank())
      p4 = p4 + ylab("Region(s) of\ninterest")
  #   #p4 = p4 + xlim(min(x$pos), max(x$pos) )
  #   #p4 = p4 + xlim(left, right )
      p4 = p4 + scale_x_continuous(limits=c(left, right), expand=c(0.02,0.02))
      p4 = p4 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  #   #p4 = p4 + scale_y_continuous(breaks=c(-0.5,0), limits=c(-0.5, 0))
  
  } else {
     p4 = NULL
  }

 ################################## plot gene annotation ###############################################
  p5 = ggplot(genes.in.p) + geom_segment(aes(x=g.start, xend=g.stop, y=2, yend=2), color=ifelse(genes.in.p$g.strand=="+", 'red', 'blue'), size=5) + theme_bw() 
  p5 = p5 + geom_text(data=genes.in.p, aes(x=(g.start+g.stop)/2, y=1.7, yend=1.7,label=paste0(gene, ' (',g.strand,')')), color=ifelse(genes.in.p$g.strand=="+",'red','blue'), size=4, hjust=0.5)
  p5 = p5 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                  axis.ticks=element_blank(), axis.title.y=element_text(size=12, color="black"),
                  panel.background=element_blank(), panel.border=element_blank(),
                  axis.text=element_blank(), axis.title.x=element_blank())
  p5 = p5 + ylab('')
  #p5 = p5 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  #p5 = p5 + xlim(left, right )
  p5 = p5 + scale_x_continuous(limits=c(left, right), expand=c(0.02,0.02))
  p5 = p5 + scale_y_continuous(breaks=c(1,2), limits=c(1, 2))
  p5 = p5 + scale_size_identity()
 
 
  ####### make plots list and return the result 
  plots <- list(p0,p1,p2,p22,p3,p4, p5)
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
res = read.table(paste0(res.dir, '/annotated_peaks_summary_final.tsv'), header=T, sep='\t', stringsAsFactors=F)
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

### read peaks with copy number data
if (file.exists(paste0(res.dir,'/processed_data/peaks_with_cn.bed'))) {
  pks.cn <- read.table(paste0(res.dir,'/processed_data/peaks_with_cn.bed'), header =F, check.names=F, sep="\t", stringsAsFactors = F)
  colnames(pks.cn) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'samples',
                        'cn.chr', 'cn.start', 'cn.stop', 'sample', 'seg.cn', 'cn.call','dist', 'cn.value')
} else {
  stop (paste0("Peaks copy number file \"", res.dir, "/processed_data/peaks_with_cn.bed\" was not found!.\n"))
}

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
if (chip.seq != 0) {
  is.chip.avail = TRUE 
  chiph = 2
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

### extract samples without SVs
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

# if (roi.file !=0) { 
#   is.roi.avail = TRUE
#   num.roi.files = unlist(strsplit(roi.file, ","))
#   if (file.exists(roi.file)) {
#     roi = read.table(roi.file, header =T, sep="\t", stringsAsFactors=F)
#     roi = roi[,1:4] 
#   } else {
#     stop (paste0("Region of interest file \"", roi.file, "\" was not found!.\n"))     
#   }
# } else {
#   is.roi.avail = FALSE
# }

######################################### plot peaks #############################################
#### create directory for plots 
dir.create(paste0(out.dir, '/peaks-plots'))

### extract peaks
pks.to.plot = unlist(strsplit(pks.to.plot,  ","))

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
  pp.cn = pks.cn[pks.cn$p.name==pk & pks.cn$sample %in% samples.with.sv, ]

  #### include samples with no SVs as neutral samples 
  if (length(samples.with.no.SVs) > 0) {
    pp.cn = rbind(pp.cn, data.frame(p.name=pk, sample=samples.with.no.SVs, cn.value=0, cn.call="neut"))
  }
  
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
      if (nrow(roi.r) !=0) {
      	  roi.r$roi.type <- gsub("Overlapped.", "", k)
      	  roi.annot <- rbind(roi.annot, roi.r)
      }
    }

    #p.roi.res <- roi[roi$name %in% unlist(strsplit(p.roi.res, "\\|")),]
    if (nrow(roi.annot) == 0 || is.null(roi.annot)) { roih = NULL } else { roih = length(unique(roi.annot$roi.type)) }
  } else {
    roi.annot <- NULL
    roih = NULL
  }
 
  ### loop through genes in the peak 
  p.genes.res <- NULL
  for (j in 1:length(assoc.genes)) {
    g = assoc.genes[j]

    ### extract gene copy number samples  
    g.neut.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="neut",'sample'])
    g.amp.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="amp",'sample'])
    g.del.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="del",'sample'])
    
    ### extract peak copy number samples
    pk.neut.samples <- unique(pp.cn[pp.cn$cn.call=="neut",'sample'])
    pk.amp.samples <-  unique(pp.cn[pp.cn$cn.call=="amp",'sample'])
    pk.del.samples <-  unique(pp.cn[pp.cn$cn.call=="del",'sample'])
    
    ### extract expression 
    g.exp <- as.data.frame(t(exp[exp[,1]==g, exp.data.cols]))
    g.exp$sample <- rownames(g.exp)
    rownames(g.exp) <- NULL
    colnames(g.exp) <- c('gene.exp', 'sample')
    g.exp <- g.exp[,c('sample', 'gene.exp')]
    
    #### a dd sample status column
    g.exp$sample.status <- 'non-SVs'
    g.exp[g.exp$sample %in% sv.pats, 'sample.status'] <- 'SVs'
    
    #### add tandem duplication status 
    #g.exp$td.status <- 'non-TD'
    #g.exp[g.exp$sample %in% DUP.pats, 'td.status'] <- 'TD'
    
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
    plot.exp.amp.del(g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, paste0(out.dir,'/peaks-plots/',g,'_',pk), g, pk)
    
    ############### run the function to plot the peak regin #################
    p.reg = plot.region(pk, p.corr, g, genes.in.peak, roi.annot)
    n = length(p.reg) 
    
    p.reg[[length(p.reg)+1]] = plot.ex[["e1"]]
    p.reg[[length(p.reg)+1]] = plot.ex[["e2"]]
    p.reg[[length(p.reg)+1]] = plot.ex[["e3"]]
    
    ### align all plots 
    all.plots <- do.call(AlignPlots, p.reg)

    ### set the layout matrix 
    mat = matrix(ncol=3, nrow=n+1)
    mat[, 1] = 1:(n+1)
    mat[, 2] = c(1:n, n+2)
    mat[, 3] = c(1:n, n+3)
    #mat[, 2] = c(rep(n+1, k), rep(n+2, n-k))
    ### set the height 
    cnh = 2
    ddh = 2
    svh1 = 2
    svh2 = 2
    geneh = 1 
    boxp = 5
    myheights = c(cnh, ddh, svh1,svh2, chiph, roih, geneh, boxp)
    mywidths = c(1.5, plot.ex[["w1"]], plot.ex[["w2"]])
    #mywidths = c(2, 3.5, 3)  
    
    #### plot all 
    pdf(paste0(out.dir,'/peaks-plots/',g,'_',pk,".pdf"), width=20, height=sum(myheights), title='', useDingbats=F, onefile=FALSE)
    grid.arrange(grobs=all.plots, nrow=n+1,ncol=3, layout_matrix=mat, heights = myheights, widths=mywidths)
    dev.off()
    
  }  ### end of genes in the peak 
  
}    ### end of peaks 

### remove temporary files 
unlink(paste0(res.dir, "/processed_data/win_cn.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.win.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.cn.tsv"))   


