#!/usr/bin/Rscript

library(reshape2)

# Find regions/peaks whose SVs altered expression of the downstream genes

args = commandArgs(T)

# debug
# run here: /gscmnt/gc5105/research/hdang/tmp
# dockit chrismaherlab/sv-hotspot
#args = unlist(strsplit('/gscmnt/gc5105/research/hdang/SV-HotSpot/code/SV-Hotspot/src/determine_gene_association_v2.r ooo/sv-hotspot-output /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_exp_release28/PRAD-CA.exp.res.tsv /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_cna_release28/PRAD-CA.res.cna.tsv 2.99 1.35 0.05 wilcox.test', ' '))[-1]
#args = unlist(strsplit('aaa/sv-hotspot-output /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_exp_release28/PAEN-AU.exp.res.tsv /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_cna_release28/PAEN-AU.res.cna.tsv 2.99 1.35 0.05 wilcox.test', ' '))

#args = unlist(strsplit('/gscmnt/gc5105/research/hdang/SV-HotSpot/code/SV-Hotspot/src/determine_gene_association_v2.r aaa/sv-hotspot-output /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_exp_release28/PAEN-AU.exp.res.tsv /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_cna_release28/PAEN-AU.res.cna.tsv 2.99 1.35 0.05 wilcox.test', ' '))[-1]

#args = unlist(strsplit('/gscmnt/gc5105/research/hdang/SV-HotSpot/code/SV-Hotspot/src/determine_gene_association_v2.r bbb/sv-hotspot-output /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_exp_release28/PRAD-CA.exp.res.tsv /gscmnt/gc6127/research/dqpham/project/SV_pancan/SV_cna_release28/PRAD-CA.res.cna.tsv 2.99 1.35 0.05 wilcox.test', ' '))[-1]

out.dir = args[1]
exp.file = args[2]
cn.file = args[3]
t.amp = as.numeric(args[4])
t.del = as.numeric(args[5])
pval = as.numeric(args[6])
stat_test = args[7]
pre.fitler = args[8]

PF = FALSE
if (!is.na(pre.fitler)) {
  PF = TRUE
  min.exp = as.numeric(unlist(strsplit(pre.fitler, '/'))[1])
  min.n = as.numeric(unlist(strsplit(pre.fitler, '/'))[2])
}

### determine t-test used 
if (!stat_test %in% c('wilcox.test', 't.test')) {
  stop(paste0("Statistical test ", stat_test, " is not available, please choose either \'wilcox.test'\ or \'t.test\'"))
}
s_test = eval(parse(text = stat_test))

############################ Function to select the top peaks for each gene ##############################
pickTopPeaks <- function (peaks, genes, total_samples) {
  ### extract all genes
  assoc.genes <- as.character(unique(genes$Gene))
  
  ### function to extract the top peak
  selectTopPeak <- function (pks, by='pvalue') {
    
    if (by == 'pvalue'){
      gg = genes[genes$Peak.name %in% pks,]
      gg = gg[order(gg$fisher.fdr),]
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

### read all break points
cat('Reading all breakpoints...')
bp = read.table(paste0(out.dir,'/processed_data/all_bp.bed'), header=F, sep='\t', quote='', stringsAsFactors=F)
colnames(bp) = c('chr', 'start', 'stop', 'name', 'score', 'strand')
bp$pos = (bp$start+bp$stop)/2
bp$sample = gsub('/.*$', '', bp$name)
bp$svtype = gsub('^.*/', '', bp$name)
### extract total number of samples 
samples.with.sv <- unique(bp$sample)
cat('done.\n')

### read peak info
cat('Reading peak information...')
x = read.table(file.path(out.dir, 'processed_data/peaks_with_overlap_nearby_genes.tsv'), header=F, stringsAsFactors=F, sep='\t')
pinf = unique(x[, c(4,7,12)]); colnames(pinf) = c('Peak.name', 'Pct.samples', 'Gene')
p2g = unique(pinf[, c('Peak.name', 'Gene')])

# get peak vs. sample vs. svtype
x = read.table(paste0(out.dir, '/processed_data/peaks_overlap_bp.tsv'), header=F, stringsAsFactors=F, sep='\t')
colnames(x) <- c('p.chr', 'p.start', 'p.stop', 'p.name', 'p.id', 'num.samples', 'pct.samples', 'sample', 'sv.type')
x = x[, c('p.name', 'sv.type', 'sample')]
colnames(x)[1] = 'Peak.name'
x = unique(x)
z = merge(x, p2g)
cat('done.\n')


#stop('DEBUG0')

# need to create a dummy data frame for
# samples w/o any SVs within reported peaks
nosv.samples = setdiff(samples.with.sv, z$sample)
if (length(nosv.samples) > 0){
    dummy = data.frame(Peak.name='dummy', sv.type='del', sample=nosv.samples,
        Gene='dummy', stringsAsFactors=F)
}else{dummy = NULL}
#z = rbind(z, dummy)


### read expression file 
cat('Reading expression data...')
if (file.exists((exp.file))) {
  x = read.table(exp.file, header =T, check.names=F, sep="\t", stringsAsFactors = F)
  #x = x[x$gene %in% z$Gene,]
  x = x[x[,1] %in% z$Gene,]
  x = melt(x, id.var='gene')
  colnames(x) = c('Gene', 'sample', 'gene.exp')
  x$sample = as.character(x$sample)
  e = x
  if (!is.null(dummy)){
      dummy$Gene = e$Gene[1] #take a random gene for dummy df
      z = rbind(z, dummy)
  }
  z = merge(z,x)
  if (nrow(z) == 0){
      stop('ERROR: Fail to merge SVs and expression.
This could be due to different sample naming between SV bedpe and expresion matrix/file\n')
  }
  z$pkgid = paste0(z$Peak.name, '/', z$Gene)
  z$Gene = NULL; z$Peak.name = NULL
  colnames(z)[2] = 'status'
  z$status = tolower(z$status)
} else {
  stop (paste0("Expression file \"", exp.file, "\" was not found!.\n"))
}
cat('done.\n')

### create df for sv samples
sv = z[!grepl('dummy', z$pkgid),]
sv$status = 'sv'
sv = unique(sv)

#stop('DEBUG1')

# identify nonSV samples 
a = dcast(pkgid~sample, data=z, value.var='status', fun.aggregate=length)
a = melt(a, id.var='pkgid')
colnames(a) = c('pkgid', 'sample', 'svcnt')
a$sample = as.character(a$sample)
a = a[a$svcnt == 0,]
a$status = 'nosv'
a$Gene = sub('^.*/(.*)$', '\\1', a$pkgid)
a = merge(a, e)
a = a[, colnames(z)]
z = z[!grepl('dummy', z$pkgid),]
z = rbind(z,sv,a)

#stop()

#### read copy number data and overlap with peaks and genes 
cat('Reading copy number data and overlap with genes...')
cn_file <- read.table(cn.file, header =T, check.names=F, sep="\t", stringsAsFactors = F, comment.char="")
if (!"cn.call" %in% colnames(cn_file)) {
  cn_file$cn.call <- 'neut'
  cn_file[cn_file$cn > t.amp, 'cn.call'] <- 'amp'
  cn_file[cn_file$cn < t.del, 'cn.call'] <- 'del'
  write.table(cn_file, file=paste0(out.dir,"/processed_data/cn.file.tsv"), quote=F, row.names=F, sep="\t", col.names=F)
  #system(paste0("grep -v chrom ", out.dir,"/processed_data/cn.file.tsv | intersectBed -wo -a ", out.dir,"/peaks/all_peaks.bed -b - | sort -k4,4 -k12,12 | sort -rk13 | groupBy -g 4,12 -c 13 -o max -full > ", out.dir,"/processed_data/peaks_with_cn.bed"))
  system(paste0("grep -v chrom ", out.dir,"/processed_data/cn.file.tsv | intersectBed -wo -a ", out.dir, "/processed_data/genes.bed -b - | sort -k4,4 -k10,10 | sort -rk11 | groupBy -g 4,10 -c 11 -o max -full > ", out.dir,"/processed_data/genes_with_cn.bed"))

} else {
  #system(paste0("grep -v chrom ", cn.file,"| intersectBed -wo -a ", out.dir,"/peaks/all_peaks.bed -b - | sort -k4,4 -k12,12 | sort -rk13 | groupBy -g 4,12 -c 13 -o max -full > ", out.dir,"/processed_data/peaks_with_cn.bed"))
  system(paste0("grep -v chrom ", cn.file,"| intersectBed -wo -a ", out.dir, "/processed_data/genes.bed -b - | sort -k4,4 -k10,10 | sort -rk11 | groupBy -g 4,10 -c 11 -o max -full > ", out.dir,"/processed_data/genes_with_cn.bed"))
}

## process cn data
cn = read.table(paste0(out.dir, '/processed_data/genes_with_cn.bed'),
                header=F, stringsAsFactors=F, sep='\t')
colnames(cn) <- c('g.chr', 'g.start', 'g.stop', 'gene', 'score', 'strand', 'cn.chr', 'cn.start', 'cn.stop', 'sample', 'seg.cn', 'cn.call','dist', 'cn.value')
cn = cn[, c('gene', 'sample', 'cn.call')]
cn = cn[!duplicated(cn[, c('gene', 'sample')]),]
cat('done.\n')

cat('Preparing data for statistical test...')

# catch a bug if cn neut is not provided
if (!('neut' %in% z$cn.call)){
}

### identify copy neutral
zn = z
zn$gene = sub('^.*/(.*)$', '\\1', zn$pkgid)
zn = merge(zn, cn, all.x=T)
zn = zn[which(zn$cn.call == 'neut'),]



### prepare data matrix for vectorization of the statistical test
### this is a sparse matrix that has twice the number of total samples
prep <- function(x){
    svtypes = setdiff(unique(x$status), 'nosv')
    e = dcast(x, pkgid + status ~ sample, value.var='gene.exp')
    e1 = e[e$status %in% svtypes,]
    e2 = e[e$status == 'nosv',]
    e = NULL
    for (svt in svtypes){
        e2$status = svt
        e = rbind(e, e2)
    }
    colnames(e)[-c(1:2)] = paste0(colnames(e)[-c(1:2)], '.none')
    colnames(e1)[-c(1:2)] = paste0(colnames(e1)[-c(1:2)], '.sv')
    e = merge(e, e1)
    rownames(e) = paste0(e$pkgid, '/', e$status)
    e = as.matrix(e[, -c(1:2)])
    return(e)
}

#stop('DEBUG')

e = prep(z)
en = prep(zn)
rownames(en) = paste0(rownames(en), '_gneut')
#if (!all(colnames(e) == colnames(en))){stop('ERR: cols mismatch!')}
if (ncol(e) != ncol(en) || any(colnames(e) != colnames(en))){
    all.samples = unique(sub('(.none$|.sv$)', '', c(colnames(e), colnames(en))))
    all.cols = c(paste0(all.samples, '.none'), paste0(all.samples, '.sv'))
    for (col in setdiff(all.cols, colnames(e))){
        e = cbind(e,NA); colnames(e)[ncol(e)] = col
    }
    for (col in setdiff(all.cols, colnames(en))){
        en = cbind(en,NA); colnames(en)[ncol(en)] = col
    }
    e = e[, all.cols]
    en = en[, all.cols]
    #message('WARN: cols mismatch!')
}

e = rbind(e,en)

# sample grouping
lbl = ifelse(grepl('none', colnames(e)), 0, 1)
ctrl = lbl==0
trmt = lbl==1

### pre-filtering
if (PF) {
  #min.exp = 10
  #min.n = 5
  enough = rowSums(!is.na(e[, lbl==0])) >= min.n & rowSums(!is.na(e[, lbl==1])) >= min.n
  mean.nosv.exp = rowMeans(e[,ctrl], na.rm=T)
  mean.sv.exp = rowMeans(e[,trmt], na.rm=T)
  expressed = mean.nosv.exp > min.exp | mean.sv.exp > min.exp
  sel = expressed & enough
  e = e[sel,]
}
cat('done.\n')

#### test for association 
cat('Running statistical test...')
pvals = suppressWarnings( apply(e, 1, function(row) s_test(row[trmt], row[ctrl])$p.value)   )
stats = suppressWarnings( apply(e, 1, function(row) s_test(row[trmt], row[ctrl])$statistic) )
test = data.frame(stat=stats, pval=pvals, stringsAsFactors=F)
test$fdr = p.adjust(pvals, method='fdr')
cat('done.\n')

# calculate logFC
test$mean.grp0.expr = rowMeans(e[, lbl==0], na.rm=T)
test$mean.grp1.expr = rowMeans(e[, lbl==1], na.rm=T)
test$logfc = log2((test$mean.grp1.expr+0.001)/(test$mean.grp0.expr+0.001))

if (any(rownames(test) != rownames(e))){stop('ERROR: something wrong!')}
test = cbind(tag=sub('^(.*)/(.*)/(.*)$', '\\3', rownames(test)),
             Gene=sub('^(.*)/(.*)/(.*)$', '\\2', rownames(test)),
             Peak.name=sub('^(.*)/(.*)/(.*)$', '\\1', rownames(test)),
             test, stringsAsFactors=F)

# Fisher's combination
# Fisher p
fisher.p <- function(pvals, max.p=1){
    pvals = pvals[!is.na(pvals)]
    return(pchisq( -2*sum(log(pvals)), df=length(pvals), lower.tail=FALSE))
}

x = test[, c('Peak.name', 'Gene', 'tag', 'pval')]
x$id = rownames(x)
x = dcast(x, Peak.name + Gene  ~ tag, value.var='pval')
nn = ncol(x) # max = 14
colnames(x)[3:nn] = paste0(colnames(x)[3:nn], '.pval')
x$min.pval = apply(x[, 3:nn], 1, function(r) min(r, na.rm=T))
x$min.fdr = p.adjust(x$min.pval, method='BH')

## combine pvalues using Fisher's method
x$fisher.pval = apply(x[, 3:nn], 1, function(r) fisher.p(r))
x$fisher.fdr = p.adjust(x$fisher.pval, method='BH')
z = x

# aggregate logfc
x = dcast(test, Peak.name + Gene  ~ tag, value.var='logfc')
colnames(x)[3:nn] = paste0(colnames(x)[3:nn], '.logfc')
fc = x

# group mean.exp
x = dcast(test, Peak.name + Gene  ~ tag, value.var='mean.grp1.expr')
colnames(x)[3:nn] = paste0(colnames(x)[3:nn], '.mean_expr')
grp1.mean.expr = x

# base group mean.exp
x = dcast(test,#test[test$tag %in% c('sv', 'sv_gneut'),],
          Peak.name + Gene  ~ tag, value.var='mean.grp0.expr')
x$nosv_gneut.mean_expr = rowMeans(x[, grepl('_gneut', colnames(x))], na.rm=T)
x$nosv.mean_expr = rowMeans(x[, colnames(x) %in%
	c('dup', 'del', 'bnd', 'inv', 'ins', 'sv')], na.rm=T)
x = x[, c('Peak.name', 'Gene', 'nosv.mean_expr', 'nosv_gneut.mean_expr')]
grp0.mean.expr = x

x = merge(grp0.mean.expr, grp1.mean.expr)
x = merge(x, fc)
z = merge(z,x,all=T)
fisher.res = z[order(z$fisher.pval),]

## write result for all peaks-genes (before group peak families)
write.table(fisher.res, file=paste0(out.dir, '/processed_data/fisher-results.tsv'), row.names=F, quote=F, sep='\t')


##############################################################################################
##### group peaks into peak families 
##############################################################################################
res = read.table(paste0(out.dir, '/processed_data/annotated_peaks_summary.tsv'), header=T, sep='\t', stringsAsFactors=F)

##### extract significant genes 
sig.genes = fisher.res[fisher.res$fisher.fdr < pval & fisher.res$min.pval <= pval, ]

if (nrow(sig.genes) == 0){
    cat('No peak found to be associated with gene expression\n')
}else{
    sig.genes <- sig.genes[order(sig.genes$fisher.fdr), ]

    final.res <- aggregate(Gene ~ Peak.name, data=sig.genes, FUN=paste, collapse='|')
    colnames(final.res) <- c('Peak.name', 'Associated.genes')
    final.res <- merge(res, final.res, sort =F)
    final.res <- final.res[order(final.res$Percentage.SV.samples, decreasing = T), ]

    #### filter results by selecting the top peaks based on the significance of overlap
    top.peaks = pickTopPeaks(final.res, sig.genes, length(samples.with.sv)) #length(samples.with.sv)) ## length(samples.with.sv)) == 101

    #### write results 
    write.table(top.peaks$peaks, file=paste0(out.dir, '/annotated_peaks_summary.tsv'), sep="\t", quote=F, row.names=F)

    cols = colnames(top.peaks$genes)[!colnames(top.peaks$genes) %in% c("Peak.name", "Gene", "Peak.family")]
    top.peaks$genes = top.peaks$genes[,c("Peak.name", "Gene", "Peak.family", cols)]
    # rename Peak.family col to more appropriate name (ie. Peaks.in.family)
    zzz = top.peaks$genes; colnames(zzz)[colnames(zzz) == 'Peak.family'] = 'Peaks.in.family'
    write.table(zzz, file=paste0(out.dir, '/genes.associated.with.SVs.tsv'), sep="\t", quote=F, row.names=F)
    #write.table(top.peaks$genes, file=paste0(out.dir, '/genes.associated.with.SVs.tsv'), sep="\t", quote=F, row.names=F)

}
