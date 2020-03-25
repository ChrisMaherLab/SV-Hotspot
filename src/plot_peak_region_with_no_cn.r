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
#cn.file = args[6]
roi.file = args[6]
chip.seq = args[7]
#t.amp = as.numeric(args[9])
#t.del = as.numeric(args[10])
chip.cov.lbl= args[8]
#roi.lbl = args[12]
left.ext = as.numeric(args[9])
right.ext = as.numeric(args[10])

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
  
  pval = (wilcox.test(g.exp[g.exp$sample.status=="SVs", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp']))$p.value
  #pval= sprintf(pval, fmt="%#.5f")
  pval = signif(pval, digits=2)

  g.exp$sample.status <- factor(g.exp$sample.status, levels=c('non-SVs', 'SVs'))
  e1 <- ggplot(g.exp, aes(x=sample.status, y=log2(gene.exp+1))) 
  e1 <- e1 + geom_boxplot(aes(fill=sample.status)) + geom_jitter(position=position_jitter(0.3)) + theme_bw() 
  e1 <- e1 + labs(x='', y=paste(gene, 'expression')) + ggtitle(paste0('\n',gene, ' expression in SV\nand non-SV samples'))
  e1 <- e1 + theme(axis.text.x=element_text(size=14, vjust=0.5, color="black"),
                   axis.text.y=element_text(size=14, color="black"), 
                   axis.title.y=element_text(size=16), panel.background=element_blank(),
                   plot.title = element_text(size = 16, hjust=0.5, color="black", face="plain"),
                   legend.position="none")
  e1 = e1 + scale_fill_manual(name="", values =c("non-SVs"="gray", "SVs"="orange2"))
  e1 = e1 + scale_x_discrete(labels=c(paste0("non-SVs\n(n=",nrow(g.exp[g.exp$sample.status=="non-SVs",]),")"), paste0("SVs\n(n=",nrow( g.exp[g.exp$sample.status=="SVs",]),")")))
  #e1 = e1 + geom_signif(comparisons=list(c('non-SVs','SVs')))
  e1 = e1 + annotate("text", x = 0.8, y = max(log2(g.exp$gene.exp+1)), label = paste0('p=', pval), cex=5, fontface="bold")
  
  return (e1)
  
}
##################################################################################################################

################################### FUNCTION TO PLOT PEAKS REGIONS ################################################
plot.region <- function(pk, pk.corr, gene, genes.in.p, p.roi, D=NULL){

   #construct region coordinates 
   right = max(pk.corr$Start, pk.corr$End, genes.in.p$g.start, genes.in.p$g.stop)
   left =  min(pk.corr$Start, pk.corr$End, genes.in.p$g.start, genes.in.p$g.stop)
   width = abs(pk.corr$End - pk.corr$Start)/1000

   ### add left and right extensions if provided 
   left <- left - left.ext - 100000     
   right <- right + right.ext  +100000
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
 
  ################################## Plot DUP $ DEL Freq ####################################
   brks = seq(0,10^9,10^6)
   brks = brks[brks >= left & brks <= right]
   labs = brks/(10^6)
   
   dup_del$svtype <- factor(dup_del$svtype, levels=c('DUP', 'DEL'))

    p1 = (ggplot(dup_del)
    + geom_segment(aes(x=pos1, xend=pos2, y=samp, yend=samp, color=svtype))
    + theme_bw(base_size=8)
    + coord_cartesian(xlim=c(left, right))
    + ylab('Duplication &\ndeletion pileup')
    + scale_x_continuous(breaks=brks, expand=c(0.05,0.05))
    + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.title=element_blank(), axis.ticks = element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(), legend.key=element_rect(fill=NA),
            axis.title.x=element_blank(), axis.title.y=element_text(size=16, color="black"),
            legend.key.size = unit(1,"cm"), legend.text=element_text(size=14),
            panel.background=element_rect(color="black"))
    + scale_color_manual(name="", values=c('DUP'='#b53f4d', 'DEL'='#2c7fb8'), guide = guide_legend(override.aes = list(size = 7)))
    + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
    )
  
  ################################## plot SVs (DUP and DEL only) ###############################################
  x2$svtype <- factor(x2$svtype, levels=c('DUP','BND','INS','INV','DEL'))
  p2 = ggplot(x2, aes(x=pos, y=num.samples, fill=svtype)) + geom_bar(stat="identity")
  p2 = p2 + theme_bw() + xlab('') + ylab('Number of\nsamples')
  p2 = p2 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  p2 = p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                 plot.title=element_blank(), axis.ticks = element_blank(),
		 panel.background=element_rect(color="black"),
                 axis.text.x=element_blank(), axis.text.y=element_text(size=14, color="black"),
                 axis.title.x=element_text(size=16, color="black"), axis.title.y=element_text(size=16, color="black"),
                 legend.key.size = unit(1,"cm"), legend.title=element_blank(), legend.text=element_text(size=14))
   p2 = p2 + scale_y_continuous(labels=abs)
   p2 = p2 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
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
                 panel.background=element_rect(color="black"),
                 axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=14, color="black"),
                 axis.title.x=element_text(size=14, color="black"), axis.title.y=element_text(size=16, color="black"),
                 legend.key.size = unit(1,"cm"), legend.title=element_text(size=16, face="bold"), legend.text=element_text(size=14))
   p22 = p22 + scale_fill_manual(name="SV type", values=c('BND'='#2ca25f','INS'='#fec44f', 'INV'='#c994c7', 'DUP'='#b53f4d', 'DEL'='#2c7fb8'), 
                                        labels=c('BND'='BND','INS'='INS', 'INV'='INV', 'DUP'='DUP', 'DEL'='DEL'),
                                        guide = guide_legend(override.aes = list(size = 7)))
   p22 = p22 + scale_x_continuous(labels = scales::comma, limits=c(left, right), expand=c(0.05,0.05))

  ######################################## plot chip-seq data ###########################################
  if (is.chip.avail) {
     reg.chip = chip_seq[chip_seq$chrom == pk.corr$Chr & 
                         chip_seq$pos >= left & 
                         chip_seq$pos <= right,]

     p3 = ggplot (reg.chip,aes(x=pos, y=cov)) + geom_bar(stat="identity", width = D/300)
     p3 = p3 + theme_bw() + xlab('') + ylab(chip.cov.lbl)
     p3 = p3 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
     p3 = p3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks = element_blank(),
                     panel.background=element_rect(color="black"),
                     axis.text.x=element_blank(), axis.text.y=element_text(size=14, color="black"),
                     axis.title.x=element_blank(), axis.title.y=element_text(size=16, color="black"))
     p3 = p3 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
  } else {
    p3 = NULL
  }

  ################################## plot region of interest annotation ###############################################
  if (!is.null(p.roi) && nrow(p.roi) !=0 ) {
      p4 = ggplot(p.roi) + geom_segment(aes(x=start, xend=end, y=roi.type, yend=roi.type), color='black', size=6) + theme_bw(base_size = 12)
      p4 = p4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                      axis.ticks=element_blank(),
                      axis.text.y=element_text(hjust=0.5, size=12, angle=ifelse(length(unique(p.roi$roi.type)) >1 ,45,90),
                                               color= ifelse(length(unique(p.roi$roi.type)) >1 ,"black","NA")),
                      axis.title.y=element_text(size=16, color="black"),
                      axis.title.x=element_blank(),axis.text.x=element_blank(),
                      panel.background=element_rect(color="black")
                      #axis.line = element_line(colour = "black")
                      ) 
      p4 = p4 + ylab(ifelse(length(unique(p.roi$roi.type)) >1, "Regions of\ninterest", p.roi$roi.type))
      p4 = p4 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
      p4 = p4 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
  
  } else {
     p4 = NULL
  }

 ################################## plot gene annotation ###############################################
  p5 = ggplot(genes.in.p) + geom_segment(aes(x=g.start, xend=g.stop, y=2, yend=2), color=ifelse(genes.in.p$g.strand=="+", 'red', 'blue'), size=5) + theme_bw() 
  p5 = p5 + geom_text(data=genes.in.p, aes(x=(g.start+g.stop)/2, y=1.7, yend=1.7,label=paste0(gene, ' (',g.strand,')')), 
                      color=ifelse(genes.in.p$g.strand=="+",'red','blue'), size=4, hjust=0.5)
  p5 = p5 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                  axis.ticks=element_blank(), axis.title.y=element_text(size=12, color="black"),
                  panel.background=element_blank(), panel.border=element_blank(),
                  axis.text=element_blank(), axis.title.x=element_blank())
  p5 = p5 + ylab('')
  p5 = p5 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
  p5 = p5 + scale_y_continuous(breaks=c(1,2), limits=c(1, 2))
  p5 = p5 + scale_size_identity()
 
 
  ####### make plots list and return the result 
  plots <- list(p1,p2,p22,p3,p4, p5)
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
if (exp.file !=0) {
  is.exp.avail = TRUE
  cat('Reading expression data...')
  if (file.exists((exp.file))) {
    exp <- read.table(exp.file, header =T, check.names=F, sep="\t", stringsAsFactors = F)
    exp.data.cols <- colnames(exp)[-1]
    ### extract samples that have expression but do have SVs
    samples.with.no.SVs <- unique(colnames(exp)[!exp.data.cols %in% samples.with.sv])
  } else {
    stop (paste0("Expression file \"", exp.file, "\" was not found!.\n"))
  }
  cat('done.\n')
  
} else {
  is.chip.avail = FALSE
}


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
if(!dir.exists(paste0(out.dir, '/peaks-plots'))){ dir.create(paste0(out.dir, '/peaks-plots')) }

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

  #### include samples with no SVs as neutral samples 
  #if (length(samples.with.no.SVs) > 0) {
  #  pp.cn = rbind(pp.cn, data.frame(p.name=pk, sample=samples.with.no.SVs, cn.value=0, cn.call="neut"))
  #}
  
  ### extract peak copy number samples
  #pk.neut.samples <- unique(pp.cn[pp.cn$cn.call=="neut",'sample'])
  #pk.amp.samples <-  unique(pp.cn[pp.cn$cn.call=="amp",'sample'])
  #pk.del.samples <-  unique(pp.cn[pp.cn$cn.call=="del",'sample'])
  
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
    if (nrow(roi.annot) == 0 || is.null(roi.annot)) { roih = NULL } else { roih = length(unique(roi.annot$roi.type))+0.5 }
  } else {
    roi.annot <- NULL
    roih = NULL
  }
 
  ### loop through genes in the peak 
  p.genes.res <- NULL
  for (j in 1:length(assoc.genes)) {
    g = assoc.genes[j]

    ### extract expression 
    g.exp <- as.data.frame(t(exp[exp[,1]==g, exp.data.cols]))
    g.exp$sample <- rownames(g.exp)
    rownames(g.exp) <- NULL
    colnames(g.exp) <- c('gene.exp', 'sample')
    g.exp <- g.exp[,c('sample', 'gene.exp')]
    
    #### add sample status column
    g.exp$sample.status <- 'non-SVs'
    g.exp[g.exp$sample %in% sv.pats, 'sample.status'] <- 'SVs'
    
    ############### run the function to plot expression  #################
    if (is.exp.avail) {
      plot.ex <- plot.exp(g.exp, BND.pats,DUP.pats,INS.pats,DEL.pats,INV.pats, g, pk)
    }
    

    ############### run the function to plot the peak regin #################
    p.reg = suppressWarnings( plot.region(pk, p.corr, g, genes.in.peak, roi.annot) )
    n = length(p.reg) 
    
    if (is.exp.avail) {
      p.reg[[length(p.reg)+1]] = plot.ex
    }
    
    ### align all plots 
    all.plots <- suppressWarnings( do.call(AlignPlots, p.reg) )

    ### set the layout matrix 
    if (is.exp.avail) {
      mat = matrix(ncol=2, nrow=n)
      mat[, 1] = 1:(n)
      mat[, 2] = c( rep(n+1, 2), rep(n+2, n-2) )
    } else {
      mat = matrix(ncol=1, nrow=n)
      mat[, 1] = 1:(n)
    }
   
    ### set the height 
    ddh = 4
    svh1 = 4
    svh2 = 4
    geneh = 1 
    #boxp = 5
    myheights = c(ddh, svh1,svh2, chiph, roih, geneh)
    mywidths = c(3, 1.5)
   
    #### plot all 
    g.title = paste0('Peak locus: ',p.corr$Chr, ':',  p.corr$Start, '-',  p.corr$End,'\n')
    mytitle=textGrob(g.title, gp=gpar(fontsize=20,fontface="bold"))
    
    pdf(paste0(out.dir,'/peaks-plots/',g,'_',pk,'.pdf'), width=18, height=sum(myheights), title='', useDingbats=F, onefile=FALSE)
    grid.arrange(grobs=all.plots, nrow=n,ncol=2, layout_matrix=mat, heights = myheights, widths=mywidths, top = mytitle)
    dev.off()
    
  }  ### end of genes in the peak 
  
}    ### end of peaks 

### remove temporary files 
unlink(paste0(res.dir, "/processed_data/win_cn.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.win.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.cn.tsv"))   



