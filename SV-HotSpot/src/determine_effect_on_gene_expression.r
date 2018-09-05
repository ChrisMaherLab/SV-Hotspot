#!/gapp/x64linux/opt/R3.1.2/bin/Rscript

# Find regions/peaks whose SVs altered expression of the downstream genes
# Written by Abdallah Eteleeb

library(plyr)

args = commandArgs(T)

out.dir = args[1]
#annot.file = args[2]
exp.file = args[2]
cn.file = args[3]
t.amp = as.numeric(args[4])
t.del = as.numeric(args[5])
pval = as.numeric(args[6])

#if (!"cn.call" %in% colnames(cn_file)) {
#  cn_file$cn.call <- 'neut'
#  cn_file[cn_file$cn < t.del, 'cn.call'] <- 'del'
#  cn_file[cn_file$cn < t.amp, 'cn.call'] <- 'amp'
#  write.table(cn_file, file=paste0(out.dir,"/temp/cn.file.tsv"), quote=F, row.names=F, sep="\t", col.names=F)

### extract feature column 
annot <- read.table(paste0(out.dir,"/temp/genes.bed"), header =T, sep="\t", check.names=F, comment.char = "$")
colnames(annot) <- c('chr', 'start', 'stop', 'gene', 'score', 'strand') 
feature.col <- colnames(annot)[4]

### read expression file 
if (file.exists(exp.file)) {
  cat('Reading expression data...\n\n')
  exp <- read.table(exp.file, header =T, check.names=F, sep="\t", stringsAsFactors = F)
  data.cols <- colnames(exp)[-1]
  
  ## if multiple rows found for a gene, select the top expressed gene
  if (any(duplicated(exp[,feature.col]))) {
    mean.expr = rowMeans(exp[, data.cols])
    exp = exp[order(mean.expr, decreasing=T),]
    exp = exp[!duplicated(exp[[feature.col]]),]
  }
} else {
  stop ("Expression file wasn't found!.\n")
}

### check if the feature names match 
if (!feature.col==colnames(exp)[1]) {
  stop ('feature name in expression file should match annotation\n')
}

#### read copy number data and overlap with peaks and genes 
cat('Reading copy number data...\n\n')
if (file.exists(cn.file)) {
   cn_file <- read.table(cn.file, header =T, check.names=F, sep="\t", stringsAsFactors = F, comment.char="")
   system(paste0("grep -v chrom ", cn.file,"| intersectBed -wo -a ", out.dir,"/peaks/all_peaks.bed -b - | sort -k4,4 -k10,10 | groupBy -g 4,10 -c 11 -o max > ", out.dir,"/temp/peaks_with_cn.bed"))
   system(paste0("grep -v chrom ", cn.file,"| intersectBed -wo -a ", out.dir, "/temp/genes.bed -b - | sort -k4,4 -k10,10 | groupBy -g 4,10 -c 11 -o max > ", out.dir,"/temp/genes_with_cn.bed"))
   #system(paste0("intersectBed -wo -a ",annot.file," -b ",cn.file," | sort -k4,4 -k10,10 | groupBy -g 4,10 -c 11 -o max > ", out.dir,"/temp/genes_with_cn.bed"))

   #### add cn.call column to peaks with copy number  
   pks.cn <-read.table(paste0(out.dir,'/temp/peaks_with_cn.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
   colnames(pks.cn) <- c('p.name', 'sample','cn.value')
   pks.cn$cn.call <- "neut"
   pks.cn[pks.cn$cn.value > t.amp, 'cn.call'] <- "amp"
   pks.cn[pks.cn$cn.value < t.del, 'cn.call'] <- "del"
   write.table(pks.cn, file=paste0(out.dir,'/temp/peaks_with_cn.bed'), sep="\t", quote=F, row.names=F) #### to be used later for visualization 

   #### add cn.call column to genes with copy number 
   genes.cn <- read.table(paste0(out.dir,'/temp/genes_with_cn.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
   colnames(genes.cn) <- c('gene', 'sample','cn.value')
   genes.cn$cn.call <- "neut"
   genes.cn[genes.cn$cn.value > t.amp, 'cn.call'] <- "amp"
   genes.cn[genes.cn$cn.value < t.del, 'cn.call'] <- "del"
   write.table(genes.cn, file=paste0(out.dir,'/temp/genes_with_cn.bed'), sep="\t", quote=F, row.names=F) #### to be used later for visualization 

} else {
  stop ("Copy number file wasn't found!.\n")
}

### get raw peak calls
pks = read.table(paste0(out.dir, '/temp/peaks_overlap_bp.tsv'), header=F, stringsAsFactors=F, sep='\t')
colnames(pks) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'pct.samples', 'sample', 'sv.type')

# read all break points
bp = read.table(paste0(out.dir,'/temp/all_bp.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(bp) = c('chr', 'start', 'stop', 'name', 'score', 'strand')
bp$pos = (bp$start+bp$stop)/2
bp$sample = gsub('/.*$', '', bp$name)
bp$svtype = gsub('^.*/', '', bp$name)
all.samples <- unique(bp$sample)

### keep shared samples only 
shared.samples <- intersect(all.samples, colnames(exp))

if (length(shared.samples)==0) { 
  stop("Samples names do not match. Make sure samples names in all files are the same.")
}

# read annotated peaks summary file 
res = read.table(paste0(out.dir, '/temp/annotated_peaks_summary.tsv'), header=T, sep='\t', stringsAsFactors=F)

#cat("Determine the effect of SVs on gene expression ...\n")
final.res <- NULL
all.genes.res <- NULL
for (i in 1:nrow(res)){
  #cat ('Processing peak:', i,'\n')
  p <- res$p.name[i]
  genes.in.peak <- c(unlist(strsplit(res$overlap.genes[i],  "\\|")), unlist(strsplit(res$nearby.genes[i],  "\\|")))
  genes.in.peak <- genes.in.peak[genes.in.peak %in% exp[,feature.col]]
  if (length(genes.in.peak)==0 ) { next }
  
  ### extract peak data
  pp = pks[pks$p.name==p & pks$sample %in% shared.samples, ]
  pp.cn = pks.cn[pks.cn$p.name==p & pks.cn$sample %in% shared.samples, ]
  
  ### extract sv/other samples 
  sv.pats <- unique(pp$sample)
  other.pats <- unique(shared.samples[!shared.samples %in% unique(pp$sample)])
  td.pats <- unique(pp[pp$sv.type=="DUP", 'sample'])
  
  p.genes.res <- NULL
  for (j in 1:length(genes.in.peak)) {
    g = genes.in.peak[j]
     
    ### extract gene copy number samples  
    g.neut.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="neut",'sample'])
    g.amp.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="amp",'sample'])
    g.del.samples <- unique(genes.cn[genes.cn$gene==g & genes.cn$cn.call =="del",'sample'])
 
    ### extract peak copy number samples
    pk.neut.samples <- unique(pp.cn[pp.cn$cn.call=="neut",'sample'])  
    pk.amp.samples <-  unique(pp.cn[pp.cn$cn.call=="amp",'sample'])   
    pk.del.samples <-  unique(pp.cn[pp.cn$cn.call=="del",'sample'])   

    ### extract expression 
    g.exp <- as.data.frame(t(exp[exp[,feature.col]==g, shared.samples ]))
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
    g.exp$pk.cn.status <- "unknown"
    g.exp[g.exp$sample %in% pk.neut.samples, "pk.cn.status"] <- "neut"
    g.exp[g.exp$sample %in% pk.amp.samples, "pk.cn.status"] <- "amp"
    g.exp[g.exp$sample %in% pk.del.samples, "pk.cn.status"] <- "del"

    #### add tandem duplication status 
    g.exp$td.status <- 'non-TD'
    g.exp[g.exp$sample %in% td.pats, 'td.status'] <- 'TD'
         
    ################################ Perform wilcoxon test ##################################################
    pvals = NULL
    ### comp 1 (SV samples versus non-SVs samples)
    if (length(sv.pats)!=0 & length(other.pats)!=0 & (length(sv.pats)>=5 | length(other.pats)>=5) ) {
    	SVs_vs_non.SVs <- wilcox.test(g.exp[g.exp$sample.status=="SVs", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
        pvals <- c(pvals, SVs_vs_non.SVs$p.value)
    } 
    ### copm 2 (using gene neutral samples only)
    neut.SVs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="SVs",'gene.exp']
    neut.nonSVs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(neut.SVs)!=0 & length(neut.nonSVs)!=0 & (length(neut.SVs)>=5 | length(neut.nonSVs)>=5 )) {
    	g.neut.comp <- wilcox.test(neut.SVs, neut.nonSVs)
        pvals <- c(pvals, g.neut.comp$p.value)
    } 
    ### copm 3 (using gene amplified samples only)
    amp.SVs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="SVs",'gene.exp']
    amp.nonSVs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(amp.SVs)!=0 & length(amp.nonSVs)!=0 & (length(amp.SVs)>=5 | length(amp.nonSVs)>=5 )) {
    	g.amp.comp <- wilcox.test(amp.SVs, amp.nonSVs)
        pvals <- c(pvals, g.amp.comp$p.value)
    }
    ### copm 4 (using gene deleted samples only)
    del.SVs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sample.status=="SVs",'gene.exp']
    del.nonSVs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(del.SVs) !=0 & length(del.nonSVs) !=0 & (length(del.SVs)>=5 | length(del.nonSVs)>=5 )) {
         g.del.comp <- wilcox.test(del.SVs, del.nonSVs)
         pvals <- c(pvals, g.del.comp$p.value)
    } 
    ###########################################################################################################
    
    ################################ Perform wilcoxon test for TD samples only ###############################
    ### comp 1 (SV samples versus non-SVs samples)
    if (length(td.pats)!=0 & length(other.pats)!=0 & (length(td.pats)>=5 | length(other.pats)>=5) ) {
        TD_vs_nonTD <- wilcox.test(g.exp[g.exp$td.status=="TD", 'gene.exp'], g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
        pvals <- c(pvals, TD_vs_nonTD$p.value)
    }
    ### copm 2 (using gene neutral TD samples only)
    neut.TDs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$td.status=="TD",'gene.exp']
    neut.nonTDs = g.exp[g.exp$gene.cn.status=="neut" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(neut.TDs)!=0 & length(neut.nonTDs)!=0 & (length(neut.TDs)>=5 | length(neut.nonTDs)>=5 )) {
        g.neutTD.comp <- wilcox.test(neut.TDs, neut.nonTDs)
        pvals <- c(pvals, g.neutTD.comp$p.value)
    }
    ### copm 3 (using gene amplified TD samples)
    amp.TDs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$td.status=="TD",'gene.exp']
    amp.nonTDs = g.exp[g.exp$gene.cn.status=="amp" & g.exp$sample.status=="non-SVs", 'gene.exp']
    if(length(amp.TDs)!=0 & length(amp.nonTDs)!=0 & (length(amp.TDs)>=5 | length(amp.nonTDs)>=5 )) {
        g.ampTD.comp <- wilcox.test(amp.TDs, amp.nonTDs)
        pvals <- c(pvals, g.ampTD.comp$p.value)
    }
    ### copm 4 (using gene deleted TD samples only)
    del.TDs = g.exp[g.exp$gene.cn.status=="del" & g.exp$td.status=="TD",'gene.exp']
    del.nonTDs = g.exp[g.exp$gene.cn.status=="del" & g.exp$sdelle.status=="non-SVs", 'gene.exp']
    if(length(del.TDs)!=0 & length(del.nonTDs)!=0 & (length(del.TDs)>=5 | length(del.nonTDs)>=5 )) {
        g.delTD.comp <- wilcox.test(del.TDs, del.nonTDs)
        pvals <- c(pvals, g.delTD.comp$p.value)
    }
    ###########################################################################################################

    ### extract the smallest p-value 
    g.pval <- signif(min(pvals), digits=4)
    if (is.na(g.pval) ) { g.pval = 1 } 
 
    ### copmute mean expression of sv samples and other samples 
    SV_mean_exp <- mean(g.exp[g.exp$sample.status=="SVs", 'gene.exp'])
    nonSV_mean_exp <- mean(g.exp[g.exp$sample.status=="non-SVs", 'gene.exp'])
    log.fc = signif(log2((SV_mean_exp+0.00001)/(nonSV_mean_exp+0.00001)), digits = 4)
    status <- 'nc'
    if (log.fc > 0 & g.pval < pval ) { status <- 'up'}  
    if (log.fc < 0 & g.pval < pval ) { status <- 'dn'}  
    
    #### combine results 
    d <- data.frame(peak_name=p, gene=g, logFC=log.fc, p_value=g.pval, status, SV_mean_exp, nonSV_mean_exp)
    p.genes.res <- rbind(p.genes.res, d)      ### for peaks final result summary
    all.genes.res <- rbind(all.genes.res, d)  ### for statistical information about genes 
    
  }  ### end of genes 
  
  ### extract genes where SVs altered expression 
  pk.genes.effected <- unique(as.character(p.genes.res[p.genes.res$p_value < pval, 'gene']))
  
  if (length(pk.genes.effected) != 0) {
    d2 <- data.frame(p.name=p, genes.effected = paste(pk.genes.effected,collapse = "|"))
  } else {
    d2 <- data.frame(p.name=p, genes.effected=NA)
  }
  
  #### combine results 
  final.res <- rbind(final.res, d2)
  
}   ### end of peaks 


### merge and write resulls
res.final <- merge(res, final.res, sort =F)
write.table(res.final, file=paste0(out.dir, '/annotated_peaks_summary_final.tsv'), sep="\t", quote=F, row.names=F)

#### extract significant genes and write results 
sig.genes <- all.genes.res[all.genes.res$p_value < pval, ]
sig.genes <- sig.genes[order(sig.genes$p_value), ]
write.table(sig.genes, file=paste0(out.dir, '/genes.effected.by.SVs.tsv'), sep="\t", quote=F, row.names=F)
