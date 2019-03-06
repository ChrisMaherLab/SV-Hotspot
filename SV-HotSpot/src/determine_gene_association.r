#!/usr/bin/Rscript

# Find regions/peaks whose SVs altered expression of the downstream genes
# Written by Abdallah Eteleeb <eteleeb@gmail.com> 

library(plyr)

args = commandArgs(T)

out.dir = args[1]
exp.file = args[2]
cn.file = args[3]
t.amp = as.numeric(args[4])
t.del = as.numeric(args[5])
pval = as.numeric(args[6])
fdr = as.numeric(args[7])
stat_test = args[8]

### determine t-test used 
if (!stat_test %in% c('wilcox.test', 't.test')) {
  stop(paste0("Statistical test ", stat_test, " is not available, please choose either \'wilcox.test'\ or \'t.test\'"))
}
s_test = eval(parse(text = stat_test))

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
#for (i in 1:10){
  #cat (i,'\n')
  pk <- res$Peak.name[i]
  genes.in.peak <- c(unlist(strsplit(res$Overlapped.genes[i],  "\\|")), unlist(strsplit(res$Nearby.genes[i],  "\\|")))
  ### remove genes wuthout expression data 
  genes.in.peak <- genes.in.peak[genes.in.peak %in% exp[,1]]
  if (length(genes.in.peak)==0 ) { next }
  
  ### extract peak data
  pp = pks[pks$p.name==pk & pks$sample %in% samples.with.sv, ]
  #pp.cn = pks.cn[pks.cn$p.name==pk & pks.cn$sample %in% samples.with.sv, ]
  #### include samples with no SVs as neutral samples 
  #pp.cn = rbind(pp.cn, data.frame(p.name=pk, sample=samples.with.no.SVs2, cn.value=0, cn.call="neut"))
  
  ### extract sv/other samples 
  sv.pats <- unique(pp$sample)
  #sv.pats <- unique(strsplit(res[res$Peak.name==pk, 'sample'], ","))
  nonSV.pats <- unique(samples.with.sv[!samples.with.sv %in% sv.pats])
  #nonSV.pats <- c(unique(samples.with.sv[!samples.with.sv %in% unique(pp$sample)]), samples.with.no.SVs) 
  td.pats <- unique(pp[pp$sv.type=="DUP", 'sample'])
  
  #### keep samples with expression data only 
  sv.pats = sv.pats[sv.pats %in% exp.data.cols]
  nonSV.pats = nonSV.pats[nonSV.pats %in% exp.data.cols]
  td.pats  = td.pats[td.pats %in% exp.data.cols]

  p.genes.res <- NULL
  for (j in 1:length(genes.in.peak)) {
    g = genes.in.peak[j]
    ### extract gene copy number samples  
    g.neut.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="neut",'sample'])
    g.amp.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="amp",'sample'])
    g.del.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="del",'sample'])
 
#   ### extract peak copy number samples
#   pk.neut.samples <- unique(pp.cn[pp.cn$cn.call=="neut",'sample'])  
#   pk.amp.samples <-  unique(pp.cn[pp.cn$cn.call=="amp",'sample'])   
#   pk.del.samples <-  unique(pp.cn[pp.cn$cn.call=="del",'sample'])   

    ### extract expression 
    #g.exp <- as.data.frame(t(exp[exp[,feature.col]==g, shared.samples ]))
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

    #### add peak copy number status 
#     g.exp$pk.cn.status <- "unknown"
#     g.exp[g.exp$sample %in% pk.neut.samples, "pk.cn.status"] <- "neut"
#     g.exp[g.exp$sample %in% pk.amp.samples, "pk.cn.status"] <- "amp"
#     g.exp[g.exp$sample %in% pk.del.samples, "pk.cn.status"] <- "del"

    #### add tandem duplication status 
    g.exp$td.status <- 'non-TD'
    g.exp[g.exp$sample %in% td.pats, 'td.status'] <- 'TD'
         
    ################################ Perform wilcoxon test ##################################################
    pvals = NULL
    ### comp 1 (SV samples versus non-SVs samples)
    if (length(sv.pats)!=0 & length(nonSV.pats)!=0 & (length(sv.pats)>=5 | length(nonSV.pats)>=5) ) {
    	  SVs_vs_non.SVs <- s_test(g.exp[g.exp$sample.status=="SVs", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
        #pvals <- c(pvals, SVs_vs_non.SVs$p.value)
        pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-all', pval=SVs_vs_non.SVs$p.value))
    } else {
	      pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-all', pval=NA))    
    } 

    ### copm 2 (using gene neutral samples only)
    neut.SVs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="SVs",'gene.exp']
    neut.nonSVs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(neut.SVs)!=0 & length(neut.nonSVs)!=0 & (length(neut.SVs)>=5 | length(neut.nonSVs)>=5)) {
    	  g.neut.comp <- s_test(neut.SVs, neut.nonSVs)
        #pvals <- c(pvals, g.neut.comp$p.value)
    	  pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-neut', pval=g.neut.comp$p.value))
    } else {
          pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-neut', pval=NA))
    }
    ### copm 3 (using gene amplified samples only)
    amp.SVs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="SVs",'gene.exp']
    amp.nonSVs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(amp.SVs)!=0 & length(amp.nonSVs)!=0 & (length(amp.SVs)>=5 | length(amp.nonSVs)>=5 )) {
    	  g.amp.comp <- s_test(amp.SVs, amp.nonSVs)
        #pvals <- c(pvals, g.amp.comp$p.value)
    	  pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-amp', pval=g.amp.comp$p.value))
    } else {
          pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-amp', pval=NA))
    }
    ### copm 4 (using gene deleted samples only)
    del.SVs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sample.status=="SVs",'gene.exp']
    del.nonSVs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(del.SVs) !=0 & length(del.nonSVs) !=0 & (length(del.SVs)>=5 | length(del.nonSVs)>=5 )) {
         g.del.comp <- s_test(del.SVs, del.nonSVs)
         #pvals <- c(pvals, g.del.comp$p.value)
         pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-del', pval=g.del.comp$p.value))
    } else {
      	 pvals <- rbind(pvals, data.frame(comp='SVs-vs-nonSVs-del', pval=NA)) 
    }
    ###########################################################################################################
    
    ################################ Perform wilcoxon test for TD samples only ###############################
    ### comp 1 (SV samples versus non-SVs samples)
    if (length(td.pats)!=0 & length(nonSV.pats)!=0 & (length(td.pats)>=5 | length(nonSV.pats)>=5) ) {
        TD_vs_nonTD <- s_test(g.exp[g.exp$td.status=="TD", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
        #pvals <- c(pvals, TD_vs_nonTD$p.value)
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-all', pval=TD_vs_nonTD$p.value))
    } else {
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-all', pval=NA))
    }
    ### copm 2 (using gene neutral TD samples only)
    neut.TDs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$td.status=="TD",'gene.exp']
    neut.nonTDs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(neut.TDs)!=0 & length(neut.nonTDs)!=0 & (length(neut.TDs)>=5 | length(neut.nonTDs)>=5 )) {
        g.neutTD.comp <- s_test(neut.TDs, neut.nonTDs)
        #pvals <- c(pvals, g.neutTD.comp$p.value)
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-neut', pval=g.neutTD.comp$p.value))
    } else {
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-neut', pval=NA))
    }
    ### copm 3 (using gene amplified TD samples)
    amp.TDs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$td.status=="TD",'gene.exp']
    amp.nonTDs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(amp.TDs)!=0 & length(amp.nonTDs)!=0 & (length(amp.TDs)>=5 | length(amp.nonTDs)>=5 )) {
        g.ampTD.comp <- s_test(amp.TDs, amp.nonTDs)
        #pvals <- c(pvals, g.ampTD.comp$p.value)
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-amp', pval=g.ampTD.comp$p.value))
    } else {
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-amp', pval=NA))
    }
    ### copm 4 (using gene deleted TD samples only)
    del.TDs = g.exp[g.exp$gene.cn.status=="del" & g.exp$td.status=="TD",'gene.exp']
    del.nonTDs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sdelle.status=="non-SVs", 'gene.exp']
    if(length(del.TDs)!=0 & length(del.nonTDs)!=0 & (length(del.TDs)>=5 | length(del.nonTDs)>=5 )) {
        g.delTD.comp <- s_test(del.TDs, del.nonTDs)
        #pvals <- c(pvals, g.delTD.comp$p.value)
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-del', pval=g.delTD.comp$p.value))
    } else {
        pvals <- rbind(pvals, data.frame(comp='TDs-vs-nonTDs-del', pval=NA))
    }
    ###########################################################################################################

    ### extract the smallest p-value 
    pvals$pval = signif(pvals$pval, digits=3)
    g.pval <- signif(min(pvals[!is.na(pvals$pval), 'pval']), digits=4)
    if (is.na(g.pval) ) { g.pval = 1 } 
 
    ### copmute mean expression of sv samples and other samples 
    SVs_mean_exp <- mean(g.exp[g.exp$sample.status=="SVs", 'gene.exp'])
    nonSVs_mean_exp <- mean(g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
    log.fc = signif(log2((SVs_mean_exp+0.00001)/(nonSVs_mean_exp+0.00001)), digits = 4)
    status <- 'nc'
    if (log.fc > 0 & g.pval < pval ) { status <- 'up'}  
    if (log.fc < 0 & g.pval < pval ) { status <- 'dn'}  
    
    #### combine results
    pvals.res = as.data.frame(t(pvals))
    rownames(pvals.res) = NULL
    cols = as.character(pvals$comp)
    colnames(pvals.res) = cols 
    pvals.res = pvals.res[-1,]

    d <- data.frame(peak_name=pk, gene=g, logFC=log.fc, min.pval=g.pval, pvals.res, SVs.vs.nonSVs.status=status, SVs_mean_exp, nonSVs_mean_exp)
    p.genes.res <- rbind(p.genes.res, d)      ### for peaks final result summary
    all.genes.res <- rbind(all.genes.res, d)  ### for statistical information about genes 
    
  }  ### end of genes in the current peak 
  
  ### extract genes where SVs altered expression 
  pk.with.assoc.genes <- unique(as.character(p.genes.res[p.genes.res$min.pval < pval, 'gene']))
  if (length(pk.with.assoc.genes) == 0) { next }
  d2 <- data.frame(Peak.name=pk, Associated.genes = paste(pk.with.assoc.genes,collapse = "|"))
  
  #### combine results 
  final.res <- rbind(final.res, d2)
  
  # Update the progress bar
  #setTxtProgressBar(pb, i)

}   ### end of peaks 

### merge and write resulls
res.final <- merge(res, final.res, sort =F)
res.final <- res.final[order(res.final$Percentage.SV.samples, decreasing = T), ]

#### extract significant genes using FDR threshold
fdr.res = p.adjust(all.genes.res$min.pval, method = "BH")
all.genes.res$FDR = fdr.res
sig.genes = all.genes.res[all.genes.res$FDR < fdr, ]
sig.genes <- sig.genes[order(sig.genes$FDR), ]

#### write final results 
write.table(res.final, file=paste0(out.dir, '/annotated_peaks_summary_final.tsv'), sep="\t", quote=F, row.names=F)
pval.cols = colnames(sig.genes)[grepl(".vs.", colnames(sig.genes)) & colnames(sig.genes) !="SVs.vs.nonSVs.status"]
colnames(sig.genes) = c("Peak.name","Gene","LogFC","Min.pval",paste0(pval.cols,".pval"),"SVs.vs.nonSVs.status", "SVs.mean.exp","nonSVs.mean.exp","FDR")
write.table(sig.genes, file=paste0(out.dir, '/genes.associated.with.SVs.tsv'), sep="\t", quote=F, row.names=F)

cat('done.')
