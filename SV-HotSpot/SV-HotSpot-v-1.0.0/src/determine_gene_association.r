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
#tool.path = args[10]

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
  cat (i,'\n')
  pk <- res$Peak.name[i]
  pk.locus = paste0(res$Chr[i],":",res$Start[i],"-",res$End[i])

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

  #p.genes.res <- NULL
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
    if (log.fc > 0) { status <- 'up'}  
    if (log.fc < 0) { status <- 'dn'}  
    
    #### combine results
    pvals.res = as.data.frame(t(pvals))
    rownames(pvals.res) = NULL
    cols = as.character(pvals$comp)
    colnames(pvals.res) = cols 
    pvals.res = pvals.res[-1,]

    d <- data.frame(Peak.name=pk, Peak.locus = pk.locus, gene=g, logFC=log.fc, min.pval=g.pval, pvals.res, SVs.vs.nonSVs.status=status, SVs_mean_exp, nonSVs_mean_exp)
    all.genes.res <- rbind(all.genes.res, d)  ### for statistical information about genes 
    
  }  ### end of genes in the current peak 
  

}   ### end of peaks 



#### compute FDR and write all results 
#fdr.res = p.adjust(all.genes.res$min.pval, method = "BH")
#all.genes.res$FDR = fdr.res
write.table(all.genes.res, file=paste0(out.dir, '/processed_data/de_results_for_all_genes.tsv'), sep="\t", quote=F, row.names=F)

#### extract significant genes using min p-value threshold
sig.genes = all.genes.res[all.genes.res$min.pval < pval, ]
sig.genes <- sig.genes[order(sig.genes$min.pval), ]

### merge and write resulls
if (nrow(sig.genes) > 0 ) {
   final.res <- aggregate(gene ~ Peak.name, data=sig.genes, FUN=paste, collapse='|')
   colnames(final.res) <- c('Peak.name', 'Associated.genes')
   final.res <- merge(res, final.res, sort =F)
   final.res <- final.res[order(final.res$Percentage.SV.samples, decreasing = T), ]
   write.table(final.res, file=paste0(out.dir, '/annotated_peaks_summary.tsv'), sep="\t", quote=F, row.names=F)

   #### write resutls for genes assoicated with SV peaks 
   pval.cols = colnames(sig.genes)[grepl(".vs.", colnames(sig.genes)) & colnames(sig.genes) !="SVs.vs.nonSVs.status"]
   colnames(sig.genes) = c("Gene", "Peak.name", "Peak.locus", "LogFC","Min.pval",paste0(pval.cols,".pval"),"SVs.vs.nonSVs.status", "SVs.mean.exp","nonSVs.mean.exp")
   write.table(sig.genes, file=paste0(out.dir, '/genes.associated.with.SVs.tsv'), sep="\t", quote=F, row.names=F)
} else {
   cat(paste("No associated genes were detected using p-value cutoff of", pval, "\n"))   
}

cat('done.\n')
cat(length(final.res$Peak.name), 'hotspots (peaks) were identified for this analysis.\n')
cat(length(unique(sig.genes$Gene)), 'genes were detected to be associated with identified hotspots.\n')


##### plot circos plot ####
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
    file.path= paste0("annotations/cytoband_info/", genome, ".cytoBand.txt.gz")
    if (file.exists(file.path)) {
      cyto.info <- read.table(file.path, colClasses = c("character", "numeric", "numeric", "character", "character"), sep = "\t", stringsAsFactors = FALSE)  
      colnames(cyto.info) = c('Chromosome','chromStart','chromEnd','Name','Stain')
    } else {
      stop(paste0("SV-Hotspot cannot generate Circos plot for identified peaks because cytoband information file was found for genome ", genome,"!.\n"))
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
  RCircos.Reset.Plot.Parameters(rcircos.params)
  RCircos.List.Plot.Parameters()
  
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
  mydata$PlotColor <- "red"
  mydata[mydata$Peak.name %in% top.pks, "PlotColor"] = "blue"
  RCircos.Histogram.Plot(mydata, data.col=5, track.num=1, "in", min.value=0, max.value=100, )
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




