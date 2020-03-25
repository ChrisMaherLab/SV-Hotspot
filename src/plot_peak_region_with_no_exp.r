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
cn.file = args[5]
roi.file = args[6]
chip.seq = args[7]
t.amp = as.numeric(args[8])
t.del = as.numeric(args[9])
chip.cov.lbl= args[10]
left.ext = as.numeric(args[11])
right.ext = as.numeric(args[12])

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


################################### FUNCTION TO PLOT PEAKS REGIONS ################################################
plot.region <- function(pk, pk.corr, genes.in.p, p.roi, D=NULL){

   #construct region coordinates 
   right = max(pk.corr$Start, pk.corr$End, genes.in.p$g.start, genes.in.p$g.stop)
   left =  min(pk.corr$Start, pk.corr$End, genes.in.p$g.start, genes.in.p$g.stop)
   width = abs(pk.corr$End - pk.corr$Start)/1000

   ### add left and right extensions if provided 
   left <- left - left.ext - 100000     
   right <- right + right.ext  + 100000
   D = right - left
   #scale binwidth accordingly based on region width
   binwidth = D/75
   #genes within region
   #g.corr = genes.in.p[genes.in.p$gene ==gene, ]
   
   ### extract SVs data
   x = cts[cts$chr == pk.corr$Chr & cts$pos > left & cts$pos < right,]  
   x = x[x$svtype !="ALL",]
   
   ### for DUP and DEL only 
   x2 = x 
   x2 = x2[x2$svtype  %in% c("DUP","DEL"),]
   x2[x2$svtype == "DEL", "num.samples"] <-  x2[x2$svtype == "DEL", "num.samples"] * -1

   ### compute the DUP and DEL pileup of SV
   dup_del = pileUp(sv, pk.corr$Chr, left, right)
   
   dup_del$pos1[dup_del$pos1 < left] = left
   dup_del$pos2[dup_del$pos2 > right] = right
 
  ################################## plot region copy number ###############################################
  if (is.cn.avail) {
     reg.width = (right - left) + 1
     reg.cn = cn_data[cn_data$chrom == pk.corr$Chr & (cn_data$pos > left | cn_data$pos < right),]
     reg.cn = reg.cn[tolower(reg.cn$cn.call) %in% c("amp", "del"), ]
     
     if (nrow(reg.cn) !=0) {
       s = round(reg.width/500)
       imin = min(reg.cn$start)
       imax = max(reg.cn$end)
       
       chr.name=unique(reg.cn$chrom)
       reg.win = data.frame(chr=chr.name, start=seq(imin,imax,s))
       reg.win$stop = reg.win$start + s - 1
       write.table(reg.cn, file=paste0(res.dir,'/processed_data/reg.cn.tsv'), quote=F, row.names=F, sep="\t", col.names=F)
       write.table(reg.win, file=paste0(res.dir,'/processed_data/reg.win.tsv'), quote=F, row.names=F, sep="\t", col.names=F)
       
       ### overlap windows with copy number 
       #cat ("Overlapping windows with copy number ..")
       system(paste0("intersectBed -wo -a ",res.dir,"/processed_data/reg.win.tsv -b ", res.dir,"/processed_data/reg.cn.tsv | cut -f2,3,7,9 | sort | uniq | sort -k1,1 -k2,2 | groupBy -full -g 1,2 -c 4 -o count > ", res.dir, "/processed_data/win_cn.tsv")) 
       
       win.data = read.table(paste0(res.dir,"/processed_data/win_cn.tsv"), header = F, sep ="\t")
       colnames(win.data) = c('start', 'stop','sample', 'cn.call', 'num.samples')
       win.data$pos = (win.data$start + win.data$stop)/2
       ### add dummy data for visualization purpuses 
       win.data=rbind(win.data, data.frame(start=-1,stop=-2,sample="dummp1",cn.call="amp",num.samples=0, pos=-1))
       win.data=rbind(win.data, data.frame(start=-1,stop=-2,sample="dummp1",cn.call="del",num.samples=0, pos=-1))
       
       p0 = ggplot(win.data, aes(x=pos,y=num.samples, fill=cn.call)) + geom_bar(stat="identity")
       p0 = p0 + scale_fill_manual(name="", values=c("amp"="#b53f4d", "del"="#2c7fb8"), labels=c("amp"="Gain", "del"="Loss"), drop=FALSE,
                                   guide = guide_legend(override.aes = list(size = 7)))
       p0 = p0 + theme_bw() + xlab('') + ylab('Copy number\nalteration frequency')
       p0 = p0 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                       plot.title=element_text(size=16, hjust=0.5, face="bold"), axis.ticks = element_blank(),
                       axis.text.x=element_blank(), axis.text.y=element_text(size=14, color="black"),
                       axis.title.x=element_blank(), axis.title.y=element_text(size=16, color="black"),
                       legend.key.size = unit(1,"cm"), legend.text=element_text(size=14),
                       panel.background=element_rect(color="black"))
       p0 = p0 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
       p0 = p0 + geom_vline(xintercept=c(pk.corr$Start, pk.corr$End), color='black', linetype='dashed')
     } else {
       p0 <- NULL 
     }
  } else {
    p0 <- NULL
  }

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
     reg.chip = chip_seq[chip_seq$chrom == pk.corr$Chr & chip_seq$pos > left | chip_seq$pos < right,]

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
  p5 = p5 + geom_text(data=genes.in.p, aes(x=(g.start+g.stop)/2, y=1.5, yend=1.5,label=paste0(gene, ' (',g.strand,')')), 
                      color=ifelse(genes.in.p$g.strand=="+",'red','blue'), size=4, hjust=0.5, angle=90)
  p5 = p5 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                  axis.ticks=element_blank(), axis.title.y=element_text(size=12, color="black"),
                  panel.background=element_blank(), panel.border=element_blank(),
                  axis.text=element_blank(), axis.title.x=element_blank())
  p5 = p5 + ylab('')
  p5 = p5 + scale_x_continuous(limits=c(left, right), expand=c(0.05,0.05))
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
res = read.table(paste0(res.dir, '/annotated_peaks_summary.tsv'), header=T, sep='\t', stringsAsFactors=F)
### filter peaks with effected genes only 
#res = res[!is.na(res$Associated.genes), ]

### read copy number file 
if (cn.file !=0) {
  is.cn.avail <- TRUE
  cnh = 4
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
  cat('done.\n')
} else {
  is.cn.avail <- FALSE
  cnh = NULL
}

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

### loop through peaks to plot
for (k in 1:length(pks.to.plot)) {
  pk = pks.to.plot[k]
  cat('\n','Plotting peak', pk, '\n')
  ### extract peak coordinates
  p.corr <- res[res$Peak.name==pk, c('Chr','Start','End', 'Percentage.SV.samples', 'Number.SV.samples')]
  
  ### extract locus information for genes nearby the peak 
  genes.in.peak <- unique(genes.and.peaks[genes.and.peaks$p.name==pk,c('g.chr','g.start','g.stop','gene','g.strand')])

  ### extract peak data
  pp = pks[pks$p.name==pk & pks$sample %in% samples.with.sv, ]
  
  ### extract sv/other samples 
  sv.pats <- unique(pp$sample)
  DUP.pats <- unique(pp[pp$sv.type=="DUP", 'sample'])    
  BND.pats <- unique(pp[pp$sv.type=="BND", 'sample'])    
  INS.pats <- unique(pp[pp$sv.type=="INS", 'sample'])    
  DEL.pats <- unique(pp[pp$sv.type=="DEL", 'sample'])
  INV.pats <- unique(pp[pp$sv.type=="INV", 'sample'])    
  nonSV.pats <- unique(samples.with.sv[!samples.with.sv %in% sv.pats])
  
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
 
    ############### run the function to plot the peak regin #################
    p.reg = suppressWarnings ( plot.region(pk, p.corr, genes.in.peak, roi.annot) )
    n = length(p.reg) 
    ### align all plots 
    all.plots <- suppressWarnings ( do.call(AlignPlots, p.reg) )

    ### set the layout matrix (version 1 with three columns)
    mat = matrix(ncol=1, nrow=n)
    mat[, 1] = 1:(n)

    ### set the height 
    ddh = 4
    svh1 = 4
    svh2 = 4
    geneh = 1.5 
    #boxp = 5
    myheights = c(cnh, ddh, svh1,svh2, chiph, roih, geneh)
    
    #### plot all
    g.title = paste0('Peak locus: ',p.corr$Chr, ':',  p.corr$Start, '-',  p.corr$End,'\n')
    mytitle=textGrob(g.title, gp=gpar(fontsize=20,fontface="bold"))
    
    pdf(paste0(out.dir,'/peaks-plots/', pk,'.pdf'), width=18, height=sum(myheights), title='', useDingbats=F, onefile=FALSE)
    grid.arrange(grobs=all.plots, nrow=n,ncol=1, layout_matrix=mat, heights = myheights, top=mytitle)
    dev.off()
    
}    ### end of peaks 

### remove temporary files  
unlink(paste0(res.dir, "/processed_data/win_cn.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.win.tsv"))   
unlink(paste0(res.dir, "/processed_data/reg.cn.tsv"))   



