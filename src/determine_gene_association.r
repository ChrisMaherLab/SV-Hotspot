#!/usr/bin/Rscript

# Find regions/peaks whose SVs altered expression of the downstream genes
# Written by Abdallah Eteleeb <eteleeb@gmail.com> 

library(plyr)
library(RCircos)
require(data.table)

args = commandArgs(T)

out.dir = args[1]
exp.file = args[2]
cn.file = args[3]
t.amp = as.numeric(args[4])
t.del = as.numeric(args[5])
pval = as.numeric(args[6])
stat_test = args[7]
GsOI <- args[8]
genome = args[9]
tool.path = args[10]
#pre.fitler = args[8]

### determine t-test used 
if (!stat_test %in% c('wilcox.test', 't.test')) {
  stop(paste0("Statistical test ", stat_test, " is not available, please choose either \'wilcox.test'\ or \'t.test\'"))
}
s_test = eval(parse(text = stat_test))

############################ Fisher's method to combine p-values ##############################
# Fisher's combination
# Fisher p
fisher.p <- function(pvals, max.p=1){
  pvals = pvals[!is.na(pvals)]
  return(pchisq( -2*sum(log(pvals)), df=length(pvals), lower.tail=FALSE))
}

############################ Function to select the top peaks for each gene ##############################
pickTopPeaks <- function (peaks, genes, total_samples) {
  ### extract all genes
  assoc.genes <- as.character(unique(genes$Gene))
  
  ### function to extract the top peak
  selectTopPeak <- function (pks, by='pvalue') {
    
    if (by == 'pvalue'){
      gg = genes[genes$Peak.name %in% pks,]
      gg = gg[order(gg$Fisher.FDR),]
      #gg = gg[order(gg$Min.pval),]
      # if there are two peaks with same p-value (which is rare), only one is chosen
      # TODO: also sort by pct. samples after p-value
      return(gg$Peak.name[1])
    }
    
    # if not ranking by fdr, rank by percent samples
    top.peak <- peaks[peaks$Peak.name %in% pks, c('Peak.name', 'Start', 'End', 'Percentage.SV.samples')]
    top.peak2 <- top.peak[top.peak$Percentage.SV.samples == max(top.peak$Percentage.SV.samples), 'Peak.name']
    
    if (length(top.peak2) > 1) {
      top.peak$len <- top.peak$End - top.peak$Start
      top.peak2 <- top.peak[top.peak$len == min(top.peak$len), 'Peak.name']
      if (length(top.peak2) > 1) {
        top.peak2 <- sample(top.peak$Peak.name, size = 1)
      }
    }
    return (top.peak2)
  }
  
  filtered.res <- NULL
  for (i in 1:length(assoc.genes)) {
    g <- assoc.genes[i]
    ### extract peaks assoicated with the gene
    gene.peaks <- as.character(unique(genes[genes$Gene==g, 'Peak.name']))
    if (length(gene.peaks) == 1) {
      dd <- data.frame(Gene=g, Peak.name = gene.peaks, Peak.family = gene.peaks)
      filtered.res <- rbind(filtered.res, dd)
      next
    }
    ### loop through al peaks and compare
    all.pairs = as.data.frame(t(combn(as.character(gene.peaks),2)), stringsAsFactors = F)
    colnames(all.pairs) <- c('peak1', 'peak2')
    all.pairs$pval <- 0
    all.pairs$status <- NA
    for (j in 1:nrow(all.pairs)) {
      p1.sample <- unique(unlist(strsplit(peaks[peaks$Peak.name==all.pairs$peak1[j], 'SV.sample'], ",")))
      p2.sample <- unique(unlist(strsplit(peaks[peaks$Peak.name==all.pairs$peak2[j], 'SV.sample'], ",")))
      ov <- length(intersect(p1.sample, p2.sample))
      ss <- total_samples - (length(p1.sample) - ov) - (length(p2.sample) - ov) - ov
      
      test.mat <-matrix(c(ov, length(p1.sample) - ov,
                          length(p2.sample) - ov,  ss), nrow = 2,
                        dimnames = list(Peak1 = c("yes", "no"),
                                        Peak2 = c("yes", "no")))
      test.pval <- fisher.test(test.mat, alternative = "two.sided")$p.value
      
      all.pairs$pval[j] <- test.pval
      if (test.pval < 0.05) {
        all.pairs$status[j] <- 'D'
      } else {
        all.pairs$status[j] <- 'I'
      }
      
    }  ## end of all pairs
    
    ### determine the status of the peaks
    if (nrow(all.pairs[all.pairs$status=="I", ])==0) {
      topPeak <- selectTopPeak(gene.peaks)
      d <- data.frame(Gene=g, Peak.name = topPeak, Peak.family = paste(gene.peaks, collapse = ","))
      filtered.res <- rbind(filtered.res, d)
    } else {
      dep.peaks <- all.pairs[all.pairs$status=="D", c('peak1', 'peak2')]
      dep.peaks <- unique(c(dep.peaks$peak1, dep.peaks$peak2))
      topPeak <- selectTopPeak(gene.peaks)
      d <- data.frame(Gene=g, Peak.name = topPeak, Peak.family = paste(gene.peaks, collapse = ","))
      filtered.res <- rbind(filtered.res, d)
      
      indep.peaks <- all.pairs[all.pairs$status=="I", c('peak1', 'peak2')]
      indep.peaks <- unique(c(indep.peaks$peak1, indep.peaks$peak2))
      
      ### loop through the indepdendent peaks and check if they are not in the dependent group
      for (k in 1:length(indep.peaks)) {
        if (!indep.peaks[k] %in% dep.peaks) {
          d2 <- data.frame(Gene=g, Peak.name = indep.peaks[k], Peak.family = indep.peaks[k])
          filtered.res <- rbind(filtered.res, d2)
        }
        
      } ## end of independent peaks
      
    }  ## end of ELSE
    
  }    ## end of all genes
  
  ### filter and write results
  annot.pks.final <- peaks[peaks$Peak.name %in% filtered.res$Peak.name, ]
  genes.final <- merge(genes, filtered.res, sort =F)
  
  rt = list()
  rt$peaks = annot.pks.final
  rt$genes = genes.final
  return(rt)
}
##########################################################################################################

##########################################################################################################
###################### Function to compute the wilcoxon test for individual SV types #####################
wilcox_test_SVtypes <- function(geneExp, DUP.pats, DEL.pats, INS.pats,INV.pats, BND.pats) { 
  
  svtypes.pvals = NULL
  
  #### DUP - comp 1-1 (SV samples versus non-SVs samples)
  if (length(DUP.pats)!=0 & length(nonSV.pats)!=0 & (length(DUP.pats)>=5 | length(nonSV.pats)>=5) ) {
    DUP_vs_nonDUP <- s_test(geneExp[geneExp$dup.status=="DUP", 'gene.exp'], geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    ## compute logFC and mean expression
    DUP_mean.exp <- mean(geneExp[geneExp$dup.status=="DUP", 'gene.exp'])
    mean_exp_DUP_nonSV <- mean(geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    DUP_vs_NoSV_LogFC <- signif(log2((DUP_mean.exp+0.00001)/(mean_exp_DUP_nonSV+0.00001)), digits = 4)
    c11 <- data.frame(DUP_vs_NoSV_LogFC, DUP_vs_NoSV_PValue=DUP_vs_nonDUP$p.value, DUP_mean.exp)
  } else {
    c11 <- data.frame(DUP_vs_NoSV_LogFC=NA, DUP_vs_NoSV_PValue=NA, DUP_mean.exp=NA)
  }
  
  ### DUP - copm 1-2 (using gene neutral TD samples only)
  neut.DUP = geneExp[geneExp$gene.cn.status=="neut" & geneExp$dup.status=="DUP",'gene.exp']
  neut.nonDUP = geneExp[geneExp$gene.cn.status=="neut" & geneExp$sample.status=="non-SVs", 'gene.exp']
  if(length(neut.DUP)!=0 & length(neut.nonDUP)!=0 & (length(neut.DUP)>=5 | length(neut.nonDUP)>=5 )) {
    g.neutDUP <- s_test(neut.DUP, neut.nonDUP)
    ## compute logFC and mean expression
    DUP.gneut_mean.exp <- mean(neut.DUP)
    mean_exp_nonSV_neut <- mean(neut.nonDUP)
    DUP.gneut_vs_NoSV.gneut_LogFC <- signif(log2((DUP.gneut_mean.exp+0.00001)/(mean_exp_nonSV_neut+0.00001)), digits = 4)
    c12 <- data.frame(DUP.gneut_vs_NoSV.gneut_LogFC, DUP.gneut_vs_NoSV.gneut_PValue=g.neutDUP$p.value, DUP.gneut_mean.exp)
  } else {
    c12 <- data.frame(DUP.gneut_vs_NoSV.gneut_LogFC=NA, DUP.gneut_vs_NoSV.gneut_PValue=NA, DUP.gneut_mean.exp=NA)
  }
  
  
  #### DEL - comp 2-1 (SV samples versus non-SVs samples)
  if (length(DEL.pats)!=0 & length(nonSV.pats)!=0 & (length(DEL.pats)>=5 | length(nonSV.pats)>=5) ) {
    DEL_vs_nonDEL <- s_test(geneExp[geneExp$del.status=="DEL", 'gene.exp'], geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    ## compute logFC and mean expression
    DEL_mean.exp <- mean(geneExp[geneExp$del.status=="DEL", 'gene.exp'])
    mean_exp_DEL_nonSV <- mean(geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    DEL_vs_NoSV_LogFC <- signif(log2((DEL_mean.exp+0.00001)/(mean_exp_DEL_nonSV+0.00001)), digits = 4)
    c21 <- data.frame(DEL_vs_NoSV_LogFC, DEL_vs_NoSV_PValue=DEL_vs_nonDEL$p.value, DEL_mean.exp)
  } else {
    c21 <- data.frame(DEL_vs_NoSV_LogFC=NA, DEL_vs_NoSV_PValue=NA, DEL_mean.exp=NA)
  }
  
  ### DEL - copm 2-2 (using gene neutral DEL samples only)
  neut.DEL = geneExp[geneExp$gene.cn.status=="neut" & geneExp$del.status=="DEL",'gene.exp']
  neut.nonDEL = geneExp[geneExp$gene.cn.status=="neut" & geneExp$sample.status=="non-SVs", 'gene.exp']
  if(length(neut.DEL)!=0 & length(neut.nonDEL)!=0 & (length(neut.DEL)>=5 | length(neut.nonDEL)>=5 )) {
    g.neutDEL <- s_test(neut.DEL, neut.nonDEL)
    ## compute logFC and mean expression
    DEL.gneut_mean.exp <- mean(neut.DEL)
    mean_exp_nonSV_neut <- mean(neut.nonDEL)
    DEL.gneut_vs_NoSV.gneut_LogFC <- signif(log2((DEL.gneut_mean.exp+0.00001)/(mean_exp_nonSV_neut+0.00001)), digits = 4)
    c22 <- data.frame(DEL.gneut_vs_NoSV.gneut_LogFC, DEL.gneut_vs_NoSV.gneut_PValue=g.neutDEL$p.value, DEL.gneut_mean.exp)
  } else {
    c22 <- data.frame(DEL.gneut_vs_NoSV.gneut_LogFC=NA, DEL.gneut_vs_NoSV.gneut_PValue=NA, DEL.gneut_mean.exp=NA)
  }
  
  
  #### INS - comp 3-1 (SV samples versus non-SVs samples)
  if (length(INS.pats)!=0 & length(nonSV.pats)!=0 & (length(INS.pats)>=5 | length(nonSV.pats)>=5) ) {
    INS_vs_nonINS <- s_test(geneExp[geneExp$ins.status=="INS", 'gene.exp'], geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    ## compute logFC and mean expression
    INS_mean.exp <- mean(geneExp[geneExp$ins.status=="INS", 'gene.exp'])
    mean_exp_INS_nonSV <- mean(geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    INS_vs_NoSV_LogFC <- signif(log2((INS_mean.exp+0.00001)/(mean_exp_INS_nonSV+0.00001)), digits = 4)
    c31 <- data.frame(INS_vs_NoSV_LogFC, INS_vs_NoSV_PValue=INS_vs_nonINS$p.value, INS_mean.exp)
  } else {
    c31 <- data.frame(INS_vs_NoSV_LogFC=NA, INS_vs_NoSV_PValue=NA, INS_mean.exp=NA)
  }
  
  ### INS - copm 3-2 (using gene neutral INS samples only)
  neut.INS = geneExp[geneExp$gene.cn.status=="neut" & geneExp$ins.status=="INS",'gene.exp']
  neut.nonINS = geneExp[geneExp$gene.cn.status=="neut" & geneExp$sample.status=="non-SVs", 'gene.exp']
  if(length(neut.INS)!=0 & length(neut.nonINS)!=0 & (length(neut.INS)>=5 | length(neut.nonINS)>=5 )) {
    g.neutINS <- s_test(neut.INS, neut.nonINS)
    ## compute logFC and mean expression
    INS.gneut_mean.exp <- mean(neut.INS)
    mean_exp_nonSV_neut <- mean(neut.nonINS)
    INS.gneut_vs_NoSV.gneut_LogFC <- signif(log2((INS.gneut_mean.exp+0.00001)/(mean_exp_nonSV_neut+0.00001)), digits = 4)
    c32 <- data.frame(INS.gneut_vs_NoSV.gneut_LogFC, INS.gneut_vs_NoSV.gneut_PValue=g.neutINS$p.value, INS.gneut_mean.exp)
  } else {
    c32 <- data.frame(INS.gneut_vs_NoSV.gneut_LogFC=NA, INS.gneut_vs_NoSV.gneut_PValue=NA, INS.gneut_mean.exp=NA)
  }
  
  
  #### INV - comp 4-1 (SV samples versus non-SVs samples)
  if (length(INV.pats)!=0 & length(nonSV.pats)!=0 & (length(INV.pats)>=5 | length(nonSV.pats)>=5) ) {
    INV_vs_nonINV <- s_test(geneExp[geneExp$inv.status=="INV", 'gene.exp'], geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    ## compute logFC and mean expression
    INV_mean.exp <- mean(geneExp[geneExp$inv.status=="INV", 'gene.exp'])
    mean_exp_INV_nonSV <- mean(geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    INV_vs_NoSV_LogFC <- signif(log2((INV_mean.exp+0.00001)/(mean_exp_INV_nonSV+0.00001)), digits = 4)
    c41 <- data.frame(INV_vs_NoSV_LogFC, INV_vs_NoSV_PValue=INV_vs_nonINV$p.value, INV_mean.exp)
  } else {
    c41 <- data.frame(INV_vs_NoSV_LogFC=NA, INV_vs_NoSV_PValue=NA, INV_mean.exp=NA)
  }
  
  ### INV - copm 4-2 (using gene neutral INV samples only)
  neut.INV = geneExp[geneExp$gene.cn.status=="neut" & geneExp$inv.status=="INV",'gene.exp']
  neut.nonINV = geneExp[geneExp$gene.cn.status=="neut" & geneExp$sample.status=="non-SVs", 'gene.exp']
  if(length(neut.INV)!=0 & length(neut.nonINV)!=0 & (length(neut.INV)>=5 | length(neut.nonINV)>=5 )) {
    g.neutINV <- s_test(neut.INV, neut.nonINV)
    ## compute logFC and mean expression
    INV.gneut_mean.exp <- mean(neut.INV)
    mean_exp_nonSV_neut <- mean(neut.nonINV)
    INV.gneut_vs_NoSV.gneut_LogFC <- signif(log2((INV.gneut_mean.exp+0.00001)/(mean_exp_nonSV_neut+0.00001)), digits = 4)
    c42 <- data.frame(INV.gneut_vs_NoSV.gneut_LogFC, INV.gneut_vs_NoSV.gneut_PValue=g.neutINV$p.value, INV.gneut_mean.exp)
  } else {
    c42 <- data.frame(INV.gneut_vs_NoSV.gneut_LogFC=NA, INV.gneut_vs_NoSV.gneut_PValue=NA, INV.gneut_mean.exp=NA)
  }
  
  #### BND - comp 5-1 (SV samples versus non-SVs samples)
  if (length(BND.pats)!=0 & length(nonSV.pats)!=0 & (length(BND.pats)>=5 | length(nonSV.pats)>=5) ) {
    BND_vs_nonBND <- s_test(geneExp[geneExp$bnd.status=="BND", 'gene.exp'], geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    ## compute logFC and mean expression
    BND_mean.exp <- mean(geneExp[geneExp$bnd.status=="BND", 'gene.exp'])
    mean_exp_BND_nonSV <- mean(geneExp[geneExp$sample.status=="non-SVs", 'gene.exp'])
    BND_vs_NoSV_LogFC <- signif(log2((BND_mean.exp+0.00001)/(mean_exp_BND_nonSV+0.00001)), digits = 4)
    c51 <- data.frame(BND_vs_NoSV_LogFC, BND_vs_NoSV_PValue=BND_vs_nonBND$p.value, BND_mean.exp)
  } else {
    c51 <- data.frame(BND_vs_NoSV_LogFC=NA, BND_vs_NoSV_PValue=NA, BND_mean.exp=NA)
  }
  
  ### BND - copm 5-2 (using gene neutral BND samples only)
  neut.BND = geneExp[geneExp$gene.cn.status=="neut" & geneExp$bnd.status=="BND",'gene.exp']
  neut.nonBND = geneExp[geneExp$gene.cn.status=="neut" & geneExp$sample.status=="non-SVs", 'gene.exp']
  if(length(neut.BND)!=0 & length(neut.nonBND)!=0 & (length(neut.BND)>=5 | length(neut.nonBND)>=5 )) {
    g.neutBND <- s_test(neut.BND, neut.nonBND)
    ## compute logFC and mean expression
    BND.gneut_mean.exp <- mean(neut.BND)
    mean_exp_nonSV_neut <- mean(neut.nonBND)
    BND.gneut_vs_NoSV.gneut_LogFC <- signif(log2((BND.gneut_mean.exp+0.00001)/(mean_exp_nonSV_neut+0.00001)), digits = 4)
    c52 <- data.frame(BND.gneut_vs_NoSV.gneut_LogFC, BND.gneut_vs_NoSV.gneut_PValue=g.neutBND$p.value, BND.gneut_mean.exp)
  } else {
    c52 <- data.frame(BND.gneut_vs_NoSV.gneut_LogFC=NA, BND.gneut_vs_NoSV.gneut_PValue=NA, BND.gneut_mean.exp=NA)
  }
  
  ## combine all 
  svtypes.test.res <- cbind(c11,c12, c21,c22, c31,c32, c41,c42, c51,c52 )
  return (svtypes.test.res)
}
##########################################################################################################

### extract feature column 
annot <- read.table(paste0(out.dir,"/processed_data/genes.bed"), header =T, sep="\t", check.names=F, comment.char = "$")
colnames(annot) <- c('chr', 'start', 'stop', 'gene', 'score', 'strand') 
#feature.col <- colnames(annot)[4]

### read expression file 
cat('Reading expression data...')
if (file.exists((exp.file))) {
  exp <- read.table(exp.file, header =T, check.names=F, sep="\t", stringsAsFactors = F)
  exp.data.cols <- colnames(exp)[-1]
} else {
   stop (paste0("Expression file \"", exp.file, "\" was not found!.\n"))
}
cat('done.\n')

### get raw peak calls
cat('Reading all breakpoints...')
pks = read.table(paste0(out.dir, '/processed_data/peaks_overlap_bp.tsv'), header=F, stringsAsFactors=F, sep='\t')
colnames(pks) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample', 'sv.type')
# read all break points
bp = read.table(paste0(out.dir,'/processed_data/all_bp.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(bp) = c('chr', 'start', 'stop', 'name', 'score', 'strand')
bp$pos = (bp$start+bp$stop)/2
bp$sample = gsub('/.*$', '', bp$name)
bp$svtype = gsub('^.*/', '', bp$name)
### extract total number of samples 
samples.with.sv <- unique(bp$sample)
cat('done.\n')

#### read copy number data and overlap with peaks and genes 
cat('Reading copy number data and overlap with peaks/genes...')
cn_file <- read.table(cn.file, header =T, check.names=F, sep="\t", stringsAsFactors = F, comment.char="")
if (!"cn.call" %in% colnames(cn_file)) {
  cn_file$cn.call <- 'neut'
  cn_file[cn_file$cn > t.amp, 'cn.call'] <- 'amp'
  cn_file[cn_file$cn < t.del, 'cn.call'] <- 'del'
  write.table(cn_file, file=paste0(out.dir,"/processed_data/cn.file.tsv"), quote=F, row.names=F, sep="\t", col.names=F)
  system(paste0("grep -v chrom ", out.dir,"/processed_data/cn.file.tsv | intersectBed -wo -a ", out.dir,"/peaks/all_peaks.bed -b - | sort -k4,4 -k12,12 | sort -rk13 | groupBy -g 4,12 -c 13 -o max -full > ", out.dir,"/processed_data/peaks_with_cn.bed"))
  system(paste0("grep -v chrom ", out.dir,"/processed_data/cn.file.tsv | intersectBed -wo -a ", out.dir, "/processed_data/genes.bed -b - | sort -k4,4 -k10,10 | sort -rk11 | groupBy -g 4,10 -c 11 -o max -full > ", out.dir,"/processed_data/genes_with_cn.bed"))

} else {
  system(paste0("grep -v chrom ", cn.file,"| intersectBed -wo -a ", out.dir,"/peaks/all_peaks.bed -b - | sort -k4,4 -k12,12 | sort -rk13 | groupBy -g 4,12 -c 13 -o max -full > ", out.dir,"/processed_data/peaks_with_cn.bed"))
  system(paste0("grep -v chrom ", cn.file,"| intersectBed -wo -a ", out.dir, "/processed_data/genes.bed -b - | sort -k4,4 -k10,10 | sort -rk11 | groupBy -g 4,10 -c 11 -o max -full > ", out.dir,"/processed_data/genes_with_cn.bed"))
}

#### read peaks with copy number file
pks.cn <-read.table(paste0(out.dir,'/processed_data/peaks_with_cn.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(pks.cn) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'samples',
                       'cn.chr', 'cn.start', 'cn.stop', 'sample', 'seg.cn', 'cn.call','dist', 'cn.value')

#### read genes with copy number file
genes.cn <- read.table(paste0(out.dir,'/processed_data/genes_with_cn.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(genes.cn) <- c('g.chr', 'g.start', 'g.stop', 'gene', 'score', 'strand', 'cn.chr', 'cn.start', 'cn.stop', 'sample', 'seg.cn', 'cn.call','dist', 'cn.value')
genes.cn = genes.cn[, c('gene', 'sample', 'cn.call')]
genes.cn = genes.cn[!duplicated(genes.cn[, c('gene', 'sample')]),]
cat('done.\n')

### read annotated peaks summary file 
cat('Reading annotated peaks...')
res = read.table(paste0(out.dir, '/processed_data/annotated_peaks_summary.tsv'), header=T, sep='\t', stringsAsFactors=F)
cat('done.\n')

cat('Examining all peaks...')
#Initiate the bar
#pb <- txtProgressBar(min = 0, max = nrow(res), style = 3)
final.res <- NULL
all.genes.res <- NULL
for (i in 1:nrow(res)){
  pk <- res$Peak.name[i]
  pk.locus = paste0(res$Chr[i],":",res$Start[i],"-",res$End[i])
  pk.num.samples <- res$Number.SV.samples[i]
  pk.perc.samples <- res$Percentage.SV.samples[i]
  
  genes.in.peak <- c(unlist(strsplit(res$Overlapped.genes[i],  "\\|")), unlist(strsplit(res$Nearby.genes[i],  "\\|")))
  ### remove genes wuthout expression data 
  genes.in.peak <- genes.in.peak[genes.in.peak %in% exp[,1]]
  if (length(genes.in.peak)==0 ) { next }
  
  ### extract peak data
  pp = pks[pks$p.name==pk & pks$sample %in% samples.with.sv, ]
  #pp.cn = pks.cn[pks.cn$p.name==pk & pks.cn$sample %in% samples.with.sv, ]
  #### include samples with no SVs as neutral samples 
  #pp.cn = rbind(pp.cn, data.frame(p.name=pk, sample=samples.with.no.SVs2, cn.value=0, cn.call="neut"))
  
  #### keep samples with expression data only 
  pp = pp[pp$sample %in% exp.data.cols, ]
  
  ### extract sv/other samples 
  sv.pats <- unique(pp$sample)
  nonSV.pats <- unique(samples.with.sv[!samples.with.sv %in% sv.pats])
  
  ### extract sv types  
  DUP.pats <- unique(pp[pp$sv.type=="DUP", 'sample'])
  DEL.pats <- unique(pp[pp$sv.type=="DEL", 'sample'])
  INS.pats <- unique(pp[pp$sv.type=="INS", 'sample'])    
  INV.pats <- unique(pp[pp$sv.type=="INV", 'sample'])    
  BND.pats <- unique(pp[pp$sv.type=="BND", 'sample']) 
  
  #p.genes.res <- NULL
  for (j in 1:length(genes.in.peak)) {
    g = genes.in.peak[j]
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
    
    #### add sample status column
    g.exp$sample.status <- 'non-SVs'
    g.exp[g.exp$sample %in% sv.pats, 'sample.status'] <- 'SVs'
    
    #### add gene copy number status 
    g.exp$gene.cn.status <- "unknown"
    g.exp[g.exp$sample %in% g.neut.samples, "gene.cn.status"] <- "neut"
    g.exp[g.exp$sample %in% g.amp.samples, "gene.cn.status"] <- "amp"
    g.exp[g.exp$sample %in% g.del.samples, "gene.cn.status"] <- "del"

    #### add SV types status 
    g.exp$dup.status <- 'nonDUP'
    g.exp[g.exp$sample %in% DUP.pats, 'dup.status'] <- 'DUP'
    
    g.exp$del.status <- 'nonDEL'
    g.exp[g.exp$sample %in% DEL.pats, 'del.status'] <- 'DEL'
    
    g.exp$ins.status <- 'nonINS'
    g.exp[g.exp$sample %in% INS.pats, 'ins.status'] <- 'INS'
    
    g.exp$inv.status <- 'nonINV'
    g.exp[g.exp$sample %in% INV.pats, 'inv.status'] <- 'INV'
    
    g.exp$bnd.status <- 'nonBND'
    g.exp[g.exp$sample %in% BND.pats, 'bnd.status'] <- 'BND'
    
    ################################ Perform wilcoxon test ##################################################
    ### comp 1 (SV samples versus non-SVs samples)
    if (length(sv.pats)!=0 & length(nonSV.pats)!=0 & (length(sv.pats)>=5 | length(nonSV.pats)>=5) ) {
    	  SVs_vs_non.SVs <- s_test(g.exp[g.exp$sample.status=="SVs", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
        ## computer logFC and mean expresison
        SV_mean.exp <- mean(g.exp[g.exp$sample.status=="SVs", 'gene.exp'])
        NoSV_mean.exp <- mean(g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
        SV_vs_NoSV_LogFC = signif(log2((SV_mean.exp+0.00001)/(NoSV_mean.exp+0.00001)), digits = 4)
        c1 <- data.frame(SV_vs_NoSV_LogFC, SV_vs_NoSV_PValue=SVs_vs_non.SVs$p.value, SV_mean.exp,NoSV_mean.exp)
    } else {
	      c1 <- data.frame(SV_vs_NoSV_LogFC=NA, SV_vs_NoSV_PValue=NA, SV_mean.exp=NA, NoSV_mean.exp=NA)
    } 

    ### copm 2 (using gene neutral samples only)
    neut.SVs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="SVs",'gene.exp']
    neut.nonSVs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(neut.SVs)!=0 & length(neut.nonSVs)!=0 & (length(neut.SVs)>=5 | length(neut.nonSVs)>=5)) {
    	  g.neut.comp <- s_test(neut.SVs, neut.nonSVs)
    	  ## computer logFC and mean expresison
    	  SV.gneut_mean.exp <- mean(neut.SVs)
    	  NoSV.gneut_mean.exp <- mean(neut.nonSVs)
    	  SV.gneut_vs_NoSV.gneut_LogFC = signif(log2((SV.gneut_mean.exp+0.00001)/(NoSV.gneut_mean.exp+0.00001)), digits = 4)
    	  c2 <- data.frame(SV.gneut_vs_NoSV.gneut_LogFC,SV.gneut_vs_NoSV.gneut_PValue=g.neut.comp$p.value, SV.gneut_mean.exp, NoSV.gneut_mean.exp)
    	 
    } else {
        c2 <-data.frame(SV.gneut_vs_NoSV.gneut_LogFC=NA, SV.gneut_vs_NoSV.gneut_PValue=NA, SV.gneut_mean.exp=NA, NoSV.gneut_mean.exp=NA)
    }
    
    ### copm 3 (using gene amplified samples only)
    amp.SVs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="SVs",'gene.exp']
    amp.nonSVs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(amp.SVs)!=0 & length(amp.nonSVs)!=0 & (length(amp.SVs)>=5 | length(amp.nonSVs)>=5 )) {
    	  g.amp.comp <- s_test(amp.SVs, amp.nonSVs)
    	  ## computer logFC and mean expresison
    	  SV.ggain_mean.exp <- mean(amp.SVs)
    	  NoSV.ggain_mean.exp <- mean(amp.nonSVs)
    	  SV.ggain_vs_NoSV.ggain_LogFC = signif(log2((SV.ggain_mean.exp+0.00001)/(NoSV.ggain_mean.exp+0.00001)), digits = 4)
    	  c3 <- data.frame(SV.ggain_vs_NoSV.ggain_LogFC, SV.ggain_vs_NoSV.ggain_PValue=g.amp.comp$p.value, SV.ggain_mean.exp,NoSV.ggain_mean.exp)
    } else {
        c3 <- data.frame(SV.ggain_vs_NoSV.ggain_LogFC=NA, SV.ggain_vs_NoSV.ggain_PValue=NA, SV.ggain_mean.exp=NA, NoSV.ggain_mean.exp=NA)
    }
    
    ### copm 4 (using gene deleted samples only)
    del.SVs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sample.status=="SVs",'gene.exp']
    del.nonSVs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(del.SVs) !=0 & length(del.nonSVs) !=0 & (length(del.SVs)>=5 | length(del.nonSVs)>=5 )) {
         g.del.comp <- s_test(del.SVs, del.nonSVs)
         ## computer logFC and mean expresison
         SV.gloss_mean.exp <- mean(del.SVs)
         NoSV.gloss_mean.exp <- mean(del.nonSVs)
         SV.gloss_vs_NoSV.gloss_LogFC = signif(log2((SV.gloss_mean.exp+0.00001)/(NoSV.gloss_mean.exp+0.00001)), digits = 4)
         c4 <- data.frame(SV.gloss_vs_NoSV.gloss_LogFC, SV.gloss_vs_NoSV.gloss_PValue=g.del.comp$p.value, SV.gloss_mean.exp, NoSV.gloss_mean.exp)
    } else {
      	 c4 <- data.frame(SV.gloss_vs_NoSV.gloss_LogFC=NA,SV.gloss_vs_NoSV.gloss_PValue=NA, SV.gloss_mean.exp=NA, NoSV.gloss_mean.exp=NA) 
    }
    ###########################################################################################################
    
    ################################# Perform wilcoxon test for individual SV types ###############################
    svtype_stat_res <- suppressWarnings (wilcox_test_SVtypes(g.exp, DUP.pats, DEL.pats, INS.pats,INV.pats, BND.pats))
    ## combine all 
    stat_test_res <- cbind(c1,c2,c3,c4, svtype_stat_res)
    
    ### extract the smallest p-value
    pvalues.cols <- colnames(stat_test_res)[grep("PValue", colnames(stat_test_res))]
    g.pval <- min(stat_test_res[, pvalues.cols], na.rm = T)
    if (is.na(g.pval) ) { g.pval = 1 } 
 
    ### copmute mean expression of sv samples and other samples 
    # SVs_mean_exp <- mean(g.exp[g.exp$sample.status=="SVs", 'gene.exp'])
    # nonSVs_mean_exp <- mean(g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
    # log.fc = signif(log2((SVs_mean_exp+0.00001)/(nonSVs_mean_exp+0.00001)), digits = 4)
    # status <- 'nc'
    # if (log.fc > 0) { status <- 'up'}  
    # if (log.fc < 0) { status <- 'dn'}  
    
    #### combine results
    #pvals.res = as.data.frame(t(pvals))
    #rownames(pvals.res) = NULL
    #cols = as.character(pvals$comp)
    #colnames(pvals.res) = cols 
    #pvals.res = pvals.res[-1,]

    d <- data.frame(Gene=g, Peak.name=pk, Peak.locus = pk.locus, Number.SV.samples = pk.num.samples, 
                    Percentage.SV.samples = pk.perc.samples, Min.pval=g.pval, stat_test_res)
    all.genes.res <- rbind(all.genes.res, d)  ### for statistical information about genes 
    
  }  ### end of genes in the current peak 
  
}   ### end of peaks 

#### compute FDR and write all results
pval.cols <- colnames(all.genes.res)[grepl('PValue', colnames(all.genes.res))]
all.genes.res$Fisher.combined.PValue = apply(all.genes.res[, pval.cols], 1, function(r) fisher.p(r))
all.genes.res$Fisher.FDR = p.adjust(all.genes.res$Fisher.combined.PValue, method='BH')

### write all results 
write.table(all.genes.res, file=paste0(out.dir, '/processed_data/de_results_for_all_genes.tsv'), sep="\t", quote=F, row.names=F)

#### extract significant genes using min p-value threshold
sig.genes = all.genes.res[all.genes.res$Fisher.FDR < 0.05 & all.genes.res$Min.pval <= 0.05, ]
sig.genes <- sig.genes[order(sig.genes$Fisher.FDR), ]
#sig.genes = all.genes.res[all.genes.res$Min.pval < pval, ]
#sig.genes <- sig.genes[order(sig.genes$Min.pval), ]

### merge and write resulls
if (nrow(sig.genes) > 0 ) {
   
   final.res <- aggregate(Gene ~ Peak.name, data=sig.genes, FUN=paste, collapse='|')
   colnames(final.res) <- c('Peak.name', 'Associated.genes')
   final.res <- merge(res, final.res, sort =F)
   final.res <- final.res[order(final.res$Percentage.SV.samples, decreasing = T), ]
   
   #### filter results by selecting the top peaks based on the significance of overlap 
   #pickTopPeaks(final.res, sig.genes, length(samples.with.sv))
   #### filter results by selecting the top peaks based on the significance of overlap
   top.peaks = pickTopPeaks(final.res, sig.genes, length(samples.with.sv)) #length(samples.with.sv)) ## length(samples.with.sv)) == 101
   
   #### write results 
   write.table(top.peaks$peaks, file=paste0(out.dir, '/annotated_peaks_summary.tsv'), sep="\t", quote=F, row.names=F)
   
   cols = colnames(top.peaks$genes)[!colnames(top.peaks$genes) %in% c("Peak.name", "Gene", "Peak.family")]
   top.peaks$genes = top.peaks$genes[,c("Peak.name", "Gene", "Peak.family", cols)]
   write.table(top.peaks$genes, file=paste0(out.dir, '/genes.associated.with.SVs.tsv'), sep="\t", quote=F, row.names=F)

} else {
   cat(paste("No associated genes were detected using p-value cutoff of", pval, "\n"))   
}

cat('done.\n')
cat(length(final.res$Peak.name), 'hotspots (peaks) were identified for this analysis.\n')
cat(length(unique(sig.genes$Gene)), 'genes were detected to be associated with identified hotspots.\n')


############### plot circos plot ############
built.in.genomes = c('hg18', 'hg19', 'hg38', 'mm9', 'mm10', 'dm3', 'dm6', 'rn4', 'rn5','rn6')
if (genome %in% built.in.genomes) {
  cat('Generating circos plot for identified peaks...')
  
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
    file.path= paste0(tool_path, "/annotations/cytoband_info/", genome, ".cytoBand.txt.gz")
    if (file.exists(file.path)) {
      cyto.info <- read.table(file.path, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t", stringsAsFactors = FALSE)  
      colnames(cyto.info) = c('Chromosome','chromStart','chromEnd','Name','Stain')
    } else {
      warning(paste0("Circos plot cannot be generated for identified peaks because cytoband information is not available for genome ", genome,".\n"))
    }
    
  }
 
  peaks.with.genes = NULL
  for (i in 1:nrow(final.res)){
    d = data.frame(Peak.name = final.res$Peak.name[i], p.chr = final.res$Chr[i], p.start = final.res$Start[i], p.end = final.res$End[i],
                   gene = unlist(strsplit(final.res$Associated.genes[i],"\\|")), pct.samples= final.res$Percentage.SV.samples[i])
    peaks.with.genes = rbind(peaks.with.genes, d)
  }
  
  #### merge meta columns with gene results 
  peaks.with.genes2 <- merge(peaks.with.genes, annot, by='gene', sort =F)
  peaks.with.genes2 = peaks.with.genes2[order(peaks.with.genes2$pct.samples, decreasing = T),]
  
  ### extract data 
  mydata = unique(peaks.with.genes2[,c('p.chr','p.start','p.end', 'Peak.name', 'pct.samples')])
    
  RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL, tracks.inside=2, tracks.outside=1)  ## for V2
  
  ##### modify plot parameters 
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$text.size = 0.5
  rcircos.params$track.height <- 0.25
  rcircos.params$track.background <- 'white'
  RCircos.Reset.Plot.Parameters(rcircos.params)
  #RCircos.List.Plot.Parameters()
  
  #### initialize the graphic device 
  pdf(paste0(out.dir,"/Circos_plot_of_all_peaks.pdf"),  height=8, width=8)
  RCircos.Set.Plot.Area()
  title("Distribution of identified peaks across the genome")
  
  ###### plot Chromosome Ideogram 
  RCircos.Chromosome.Ideogram.Plot();
  
  ################### plot the histogram ################################
  ## extract top peaks for each chromosome 
  mydata.dt = as.data.table(mydata)
  top.pks = unique(as.character(mydata.dt[mydata.dt[, .I[pct.samples==max(pct.samples)], by=p.chr]$V1]$Peak.name))
  mydata$PlotColor <- "gray60"
  mydata[mydata$Peak.name %in% top.pks, "PlotColor"] = "blue"
  RCircos.Histogram.Plot(mydata, data.col=5, track.num=1, "in", min.value=0, max.value=ceiling(max(mydata$pct.samples))+2)
  ### plot legend 
  legend (-0.45,-1.4, legend=c('All peaks', 'Top peaks'), fill=c("red","blue"),  horiz = TRUE, bty="n",
          border="white", cex=0.7, x.intersp=0.5)
  
  ###################### plot gene names ################################
  #### draw genes of interest names 
  if (GsOI !=0) {
    genes.of.int = scan(GsOI, 'character') 
    genes.of.int <- unique(annot[annot$gene %in% genes.of.int, c('chr', 'start', 'stop', 'gene')])
    RCircos.Gene.Connector.Plot(genes.of.int,  track.num=2, side="in")
    RCircos.Gene.Name.Plot(genes.of.int, name.col=4,track.num=3, side="in")
  }
  
  ################### add max and min values  ##########################
  mtext(paste0('Min. value = ',round(min(mydata$pct.samples)),"%"), side =3, adj=0.51, padj = 47.2, cex=0.7, font=2)
  mtext(paste0('Max. value = ',round(max(mydata$pct.samples)),"%"), side =3, adj=0.51, padj = 48.5, cex=0.7, font=2)
  
  dev.off()
  
}




