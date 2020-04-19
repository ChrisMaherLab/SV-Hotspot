# Filtering SV-HotSpot results

usage = '

Post-filtering of SV-HotSpot results. Output will be written
to SV-HotSpot result directory:
    annotated_peaks_summary.filtered.tsv
    genes.associated.with.SVs.filtered.tsv


Usage:  Rscript filter.r
    <SV-HotSpot_result_dir>
    <max.p-value-to-infer-expression-association>
    <min.logfc>
    <min.group.expr>
    <max.number.of.associated.gene.per.peak, peak w/ more genes associated are ignored>
    <max.peak.length, peak wider than this w/o assoc. w/ known genes are ignored>
    <min.peak.percent.samples, 0-100, peaks with percent of samples < this are ignored>
    <tsv_file_of_cosmic_census_genes>
    <csv_other_known_genes_to_keep>

    eg. Rscript filter.r sv-hotspot-output 0.05 0.22 10 9 500000 15 \
            data/cosmic_census_genes.tsv AR,ERG,PTEN,TP53

'

args = commandArgs(T)

if (length(args) < 8){stop(usage)}

in.res.dir = args[1]
max.pval = as.numeric(args[2])
min.logfc = as.numeric(args[3])
min.grp.mean.exp = as.numeric(args[4])
max.num.associated.genes.per.peak = as.numeric(args[5])
max.peak.length = as.numeric(args[6])
min.peak.pct.samples = as.numeric(args[7])
cosmic.census.file = args[8]
known.genes = unlist(strsplit(args[9], ','))

# read col changes
#cc = read.table('cols.tsv', header=T, stringsAsFactors=F, sep='\t')

# read cosmic census genes
cosmic = read.table(cosmic.census.file, header=T, stringsAsFactors=F, sep='\t')
cosmic = cosmic[, c('Gene.Symbol', 'Hallmark', 'Tier', 'Role.in.Cancer', 'Tumour.Types.Somatic.')]
colnames(cosmic) = c('Gene', 'Hallmark', 'Tier', 'Role.In.Cancer', 'Cancer.Types')
colnames(cosmic)[-1] = paste0('COSMIC.', colnames(cosmic)[-1])

# filtering genes/peaks
gg = known.genes

# read peaks
x = read.table(file.path(in.res.dir, 'annotated_peaks_summary.tsv'),
               header=T, stringsAsFactors=F, sep='\t')
x$peak.len = x$End - x$Start + 1
x$peak.density = x$Percentage.SV.samples/x$peak.len*1000

# read genes associated w/ peaks & merge w/ peak info
z = read.table(file.path(in.res.dir, 'genes.associated.with.SVs.tsv'),
               header=T, stringsAsFactors=F, sep='\t')
z = merge(z,x[, c('Peak.name', 'peak.len', 'peak.density',
                  'Percentage.SV.types', 'Dominant.svtype')])

#for (i in 1:nrow(cc)){
#    colnames(z)[colnames(z) == cc$current[i]] = cc$new[i]
#}

# remove comparisons with gene gain/loss
z = z[, !grepl('ggain|gloss', colnames(z))]


exp.cols = grep('mean.exp', colnames(z), value=T)
p.cols = grep('PValue', colnames(z), value=T)
logfc.cols = grep('LogFC', colnames(z), value=T)
comps = sub('_LogFC', '', logfc.cols)

# identify significant comparisons
cmp.keys = c()
for (cmp in comps){
    tt = unlist(strsplit(cmp, '_vs_'))
    ctrl = tt[2]; svgrp = tt[1]
    p.col = paste0(cmp, '_PValue')
    logfc.col = paste0(cmp, '_LogFC')
    ctrl.exp.col = paste0(ctrl, '_mean.exp')
    svgrp.exp.col = paste0(svgrp, '_mean.exp')
    if (!all(c(p.col, logfc.col, ctrl.exp.col, svgrp.exp.col) %in% colnames(z))){
        stop('ERROR: missing some columns\n')
    }
    stat.col = paste0(cmp, '_Status')
    z[[stat.col]] = NA
    z[[stat.col]][!is.na(z[[p.col]])] = 'unchanged'
    z[[stat.col]][which(z[[p.col]] <= max.pval & z[[logfc.col]] >= min.logfc)] = 'upregulated'
    z[[stat.col]][which(z[[p.col]] <= max.pval & z[[logfc.col]] <= -min.logfc)] = 'downregulated'
    z[[stat.col]][z[[ctrl.exp.col]] < min.grp.mean.exp &
                  z[[svgrp.exp.col]] < min.grp.mean.exp] = 'low_expr'

}

# count number of up/down comparisions
z$Num.Comp.Upregulated = rowSums(z[, grep('_Status', colnames(z))] == 'upregulated', na.rm=T)
z$Num.Comp.Downregulated = rowSums(z[, grep('_Status', colnames(z))] == 'downregulated', na.rm=T)

# determine types of SVs associated with expression
comp.keys = sub('_vs.*$', '', comps)
z$Upregulated.With = apply(z[, paste0(comps, '_Status')], 1,
     function(row) paste(comp.keys[which(row == 'upregulated')], collapse=','))
z$Downregulated.With = apply(z[, paste0(comps, '_Status')], 1,
     function(row) paste(comp.keys[which(row == 'downregulated')], collapse=','))
z$Upregulated.With.Uniq = apply(z[, paste0(comps, '_Status')], 1,
     function(row) paste(unique(sub('.gneut|.gloss|.ggain', '',
        comp.keys[which(row == 'upregulated')])), collapse=','))
z$Downregulated.With.Uniq = apply(z[, paste0(comps, '_Status')], 1,
     function(row) paste(unique(sub('.gneut|.gloss|.ggain', '',
        comp.keys[which(row == 'downregulated')])), collapse=','))


z = z[order(z$Min.pval),]

# eval
# z[z$Gene %in% gg, c('Gene', 'Num.Comp.Upregulated', 'Num.Comp.Downregulated',
#                    'Upregulated.With', 'Downregulated.With')]
# z[z$Num.Comp.Upregulated > 3 | z$Num.Comp.Downregulated > 3,
#  c('Gene', 'Num.Comp.Upregulated', 'Num.Comp.Downregulated',
#    'Upregulated.With', 'Downregulated.With')]


# remove genes with ambiguous directions of DE (>=2 in both up/down)
# with(z, table(Num.Comp.Upregulated, Num.Comp.Downregulated))
ambiguous = z$Num.Comp.Upregulated > 1 & z$Num.Comp.Downregulated > 1
DE = z$Num.Comp.Upregulated != 0 | z$Num.Comp.Downregulated != 0
z = z[DE & !ambiguous,]

# Determine dominant DE direction in the presence of SVs
z$Dominant.DE.Direction.WithSV = 'Upregulated'
z$Dominant.DE.Direction.WithSV[z$Num.Comp.Upregulated <
                         z$Num.Comp.Downregulated] = 'Downregulated'
z$Dominant.DE.Direction.WithSV[z$Num.Comp.Upregulated ==
                         z$Num.Comp.Downregulated] = 'Upregulated/Downregulated'
#z[z$Dominant.DE.Direction.WithSV == 'Upregulated/Downregulated',
#  c('Gene', 'Upregulated.With', 'Downregulated.With', 'peak.len')]

# for genes with 1 up and 1 down directions
# use comparison within group whose gene is copy neutral
ambiguous = z$Num.Comp.Upregulated == z$Num.Comp.Downregulated
z$Dominant.DE.Direction.WithSV[ambiguous & grepl('gneut', z$Upregulated.With)
                     & !grepl('gneut', z$Downregulated.With)] = 'Upregulated'
z$Dominant.DE.Direction.WithSV[ambiguous & !grepl('gneut', z$Upregulated.With)
                     & grepl('gneut', z$Downregulated.With)] = 'Downregulated'

# remove contradicting direction after cleaning up
z = z[z$Dominant.DE.Direction.WithSV %in% c('Upregulated', 'Downregulated'),]

# count num of genes per peak
pkg = aggregate(Gene ~ Peak.name, data=z, FUN=length)
colnames(pkg) = c('Peak.name', 'Num.associated.genes')
pkg = pkg[order(pkg$Num.associated.genes),]
z = merge(z, pkg)

# eval
# z[z$Gene %in% gg, c('Gene', 'Num.associated.genes')]
# z[z$Num.associated.genes > 9, c('Gene', 'Num.associated.genes', 'peak.len')]
# table(z$Num.associated.genes > 9)

# annotate with Cosmic
z = merge(z, cosmic, all.x=T)
for (col in grep('COSMIC', colnames(z))){
    z[[col]][is.na(z[[col]])] = ''
}

#z = z[(z$peak.len < 1000000 & z$Num.associated.genes < 8) | z$Gene %in% gg,]
z = z[(z$peak.len <= max.peak.length | z$Gene %in% c(cosmic$Gene, known.genes)) &
      z$Num.associated.genes <= max.num.associated.genes.per.peak &
      z$Percentage.SV.samples >= min.peak.pct.samples, ]
z = z[order(z$Min.pval),]
z = cbind(Rank = 1:nrow(z), z, stringsAsFactors=F)

# eval
# z[z$Gene %in% gg, c('Rank', 'Gene', 'Num.Comp.Upregulated',
#         'Num.Comp.Downregulated', 'Upregulated.With', 'Downregulated.With')]
# z[z$Gene %in% gg, c('Rank', 'Gene', 'Num.Comp.Upregulated',
#          'Num.Comp.Downregulated', 'peak.len', 'Dominant.DE.Direction.WithSV')]

cat('Num. of peaks:', length(unique(z$Peak.name)), '\n')
cat('Num. of genes:', nrow(z), '\n')

# write gene vs peak list
g = z
write.table(g, file=file.path(in.res.dir, 'genes.associated.with.SVs.filtered.tsv'),
                           row.names=F, sep='\t', quote=F)

# select peaks and write
pk = aggregate(Gene ~ Peak.name, data=g, FUN=paste, collapse=',')
colnames(pk) = c('Peak.name', 'Associated.genes')
x$Associated.genes = NULL
pk = merge(x, pk)
pk = merge(pk, pkg)
pk = pk[, c('Peak.name', 'Chr', 'Start', 'End', 'peak.len',
            'Percentage.SV.samples', 'Percentage.SV.types', 'Dominant.svtype',
            'Num.associated.genes', 'Associated.genes', 'Overlapped.GeneHancers')]
pk = pk[order(pk$Percentage.SV.samples, decreasing=T),]
write.table(pk , file=file.path(in.res.dir, 'annotated_peaks_summary.filtered.tsv'),
            row.names=F, sep='\t', quote=F)

# create summarization figures
library(ggplot2)
library(reshape2)
pc = as.data.frame(table(pk$Chr))
colnames(pc) = c('chr', 'num.peaks')
gc = as.data.frame(table(sub(':.*$', '', g$Peak.locus)))
colnames(gc) = c('chr', 'num.associated.genes')
d = merge(gc, pc)
d$chr = sub('chr', '', d$chr)
if (!('Y' %in% d$chr)){
    d = rbind(d, data.frame(chr='Y', num.associated.genes=0, num.peaks=0))
}
matched_chrs = match(c(1:22, 'X', 'Y'), d$chr)
d = d[matched_chrs[!is.na(matched_chrs)],]

d$chr = factor(d$chr, levels=d$chr)
d = melt(d, id.var='chr')
colnames(d) = c('chr', 'grp', 'count')
d$grp = factor(d$grp, levels=c('num.peaks', 'num.associated.genes'))

p = (ggplot(d, aes(x=chr, y=count, fill=grp, color=grp))
     + geom_bar(stat='identity', position='dodge', width=0.8)
     + annotate('text', x=24, y=0.95*max(d$count), label=paste0('Total peaks: ',
              nrow(pk), '\nTotal genes:', nrow(z)), size=2.5, hjust=1, vjust=1)
     + ylab('Count') + xlab('Chromosome')
     + theme_bw()
     + scale_color_manual(values=c('darkgreen', 'cadetblue'))
     + scale_fill_manual(values=c('darkgreen', 'cadetblue'))
     + theme(legend.position='top', panel.grid=element_blank()
     )
     
)
ggsave(p, file=paste0(in.res.dir, '/peaks-genes-per-chr.pdf'), width=5, height=1.75)


