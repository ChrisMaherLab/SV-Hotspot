#!/usr/bin/perl -w

# SV-HSD: strucutral varaint tool for detecting host spots 
# Created by Abdallah Eteleb <eteleeb@gmail.com> and Ha Dang Ha X. Dang <haxdang@gmail.com> 
#
# Version 1.0.0  6/30/2018
#

use strict;
use warnings; 
use POSIX;
use Getopt::Std; 
use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

my $TOOL_PATH='/gscmnt/gc5111/research/eteleeb/projects/SV-HotSpot';

#define input options with default values 
my $sv_file='0';
my $genome='hg38';
my $sliding_w_size = 100000;
my $sliding_w_step = 1000; 
my $output_dir = getcwd();
my $annot_file='0';
my $peakPick_win=100;
my $peakPick_minsd=5;
my $pct_samples_t=5;
my $expr_file='0'; 
my $cn_file='0';
my $t_amp=2.8;
my $t_del= 0.5;
my $pvalue=0.05;
my $genes_of_int='0'; 
my $region_of_int='0'; 
my $chrom="ALL";
my $sv_type="ALL";
my $merge_dist=1000;
my $temp_keep='0';
my $plot_peaks="F";
my $num_nearby_genes=4;
my $chip_cov='0';
my $chip_cov_lbl="";
my $roi_lbl ="";

GetOptions
(
    'sv=s' =>\$sv_file,
    'g|genome=s' =>\$genome,
    'w|sliding-win-size=i' =>\$sliding_w_size,
    's|sliding-win-step=i' =>\$sliding_w_step,
    'o|output=s' =>\$output_dir,
    'a|annot=s' =>\$annot_file,  
    'W|peakPick-window-size=i' =>\$peakPick_win,
    'P|peakPick-min-sd=i' =>\$peakPick_minsd,
    'T|pct-samples=i' =>\$pct_samples_t,
    'e|expr=s' =>\$expr_file,
    'c|cn=s' =>\$cn_file,
    't-amp=f' =>\$t_amp,
    't-del=f' =>\$t_del,
    'p|pval=f' =>\$pvalue,
    'r|region-of-int=s' =>\$region_of_int,
    'G|genes-of-int=s' =>\$genes_of_int,
    'C|chrom=s' => \$chrom,
    'S|sv-type=s' => \$sv_type,
    'd|group-dist-size=i' => \$merge_dist,
    'k|num-nearby-genes=i' => \$num_nearby_genes,
    'chip-cov=s' =>\$chip_cov,
    'plot-peaks=s' => \$plot_peaks,
    'chip-cov-lbl=s' => \$chip_cov_lbl,
    'roi-lbl=s' => \$roi_lbl,
    'keep-temp=s' => \$temp_keep
);

usage() if (!$sv_file | !$genome | !$expr_file | !$cn_file);

#### define the name of the chromosme lengths file and output directroy 
my $chromsize_file = $TOOL_PATH.'/annotations/'.$genome.'/chromsize.tsv';

### Creating results directories";
$output_dir = $output_dir.'/sv-hotspot-output';
system ("rm -rf $output_dir; mkdir $output_dir");
system ("rm -rf $output_dir/temp; mkdir $output_dir/temp");

#### check if annotation file and genes of interest file were provided otherwise use defaults
if (!$annot_file) {
   $annot_file = $TOOL_PATH.'/annotations/'.$genome.'/genes.bed'
}

### prepare genes of interest file 
if (!$genes_of_int) {
   $genes_of_int = $TOOL_PATH.'/annotations/genes-of-interest.txt';
   system("cat $genes_of_int | select-rows.pl 0 $annot_file 3 > $output_dir/temp/genes-of-interest.bed");
   $genes_of_int = $output_dir."/temp/genes-of-interest.bed"
} else {
   system("cat $genes_of_int | select-rows.pl 0 $annot_file 3 > $output_dir/temp/genes-of-interest.bed");
   $genes_of_int = $output_dir."/temp/genes-of-interest.bed"
}

#### show all inputes 
print "########################################################\n";
print "######             SV-HotSpot v1.0.0              ######\n";
print "########################################################\n";
print " Genome: $genome\n";
print " Sliding window size: $sliding_w_size\n";
print " Sliding window step: $sliding_w_step\n"; 
print " Output directory: $output_dir\n";
print " Annotation file: $annot_file\n";
print " peakPick window size: $peakPick_win\n";
print " peakPick minimum SD: $peakPick_minsd\n";
print " Percentage of samples threshold: $pct_samples_t\n";
print " Expression file: $expr_file\n";  
print " Copy number file: $cn_file\n";
print " Amplification threshold: $t_amp\n";
print " Deletion threshold: $t_del\n";
print " P-value threshold: $pvalue\n";
print " Geens of interest: $genes_of_int\n"; 
print " Region of interest: $region_of_int\n"; 
print " Chromosmes to analyze: $chrom\n";
print " SV types to analyze: $sv_type\n"; 
print " Distance for merging peaks: $merge_dist\n"; 
print " Number of nearby genes: $num_nearby_genes\n"; 
print " chip-Seq file: $chip_cov\n"; 
#print " Plot peaks $plot_peaks="";
print " Keep temporary file: $temp_keep\n"; 
print "########################################################\n\n";

#print "\n-------------------------------------------\n";
#print "Preparing breakpoint data ...\n";
#print "-------------------------------------------\n";
system("cat $sv_file | awk -F\"\\t\" '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$7,\$8,\$9; print \$4,\$5,\$6,\$7,\$8,\$10}' | grep -v \"chrom\" > $output_dir/temp/all_bp.bed");
system ("awk '{if(\$2!=0 || \$3!=0) {print \$0}}' $output_dir/temp/all_bp.bed > $output_dir/temp/all2_bp.bed; mv $output_dir/temp/all2_bp.bed $output_dir/temp/all_bp.bed");
system ("for type in `cut -f4 $output_dir/temp/all_bp.bed | cut -f2 -d'/' | sort | uniq`; do cat $output_dir/temp/all_bp.bed | grep \"\/\$type\" > $output_dir/temp/bp.\$type.bed; done");

print "\n--------------------------------------------------\n";
print "STEP 1: Identifying Peaks (hotspot regions) \n";
print "--------------------------------------------------\n";

print "STEP 1: Segmenting the genome into sliding windows\n";
system ("genome-to-sliding-window.r $chromsize_file $sliding_w_size $sliding_w_step $output_dir");

print "Overlapping breakpoints with sliding windows\n";
system ("intersectBed -wao -a $output_dir/temp/genome.segments.bed -b $output_dir/temp/all_bp.bed > $output_dir/temp/genome.segments.with.bps.bed");

print "\nSummarizing sample counts\n";
system ("summarize-sample-count.r $output_dir $plot_peaks"); 

print "\nCall structural variant peaks (hot spots)\n";
system ("find-peaks.r $sv_type $chrom $peakPick_win $peakPick_minsd $pct_samples_t $output_dir $merge_dist $TOOL_PATH $genes_of_int");

print "\n--------------------------------------------------\n";
print "STEP 2: Annotating Peaks \n";
print "--------------------------------------------------\n";

system("grep -v chrom $annot_file | sort -k1,1 -k2,2n > $output_dir/temp/genes_sorted.bed; mv $output_dir/temp/genes_sorted.bed $output_dir/temp/genes.bed");
$annot_file = $output_dir."/temp/genes.bed";
if ($region_of_int) {
   system("grep -v chrom $region_of_int | sort -k1,1 -k2,2n > $output_dir/temp/reg_of_int.bed")
}
system("annotate_peaks.sh $genome $region_of_int $output_dir $num_nearby_genes $TOOL_PATH"); 

print "\n--------------------------------------------------\n";
print "STEP 3: Determining the effect of SVs on gene expression\n";
print "--------------------------------------------------\n";

if ($expr_file & $cn_file) {
  system ("determine_effect_on_gene_expression.r $output_dir $expr_file $cn_file $t_amp $t_del $pvalue");
} else {
  print "To determine the effect of SVs on gene expression, both expression and copy number data are required. \n";
  exit(0); 
}

print "\n--------------------------------------------------\n";
print "STEP 4: Visualizing hotspot regions \n";
print "--------------------------------------------------\n";

if ($expr_file & $cn_file) {
   # check if the user provided chip-seq data, if yes run the script to prcess it (average over a windwo)
   if ($chip_cov) {
      print "Processing chip-seq data, please wait as this may take several minutes\n";
      system("process_chip_data.r $chip_cov $output_dir");
   } 
   system ("plot_sv_region.r $sv_file $output_dir $expr_file $cn_file $chip_cov $t_amp $t_del $chip_cov_lbl $roi_lbl");
} else {
  print "Both expression and copy number data are required to generate visualization\n";
  exit(0);
}

if ($temp_keep eq "T"){
    #system ("rm -rf $output_dir/temp"); 
}


sub usage
{
   #use Term::ANSIColor;
   print "\n";
   print "USAGE:\n      sv-hsd [OPTIONS] -g/--genome <genomeName> --sv <structuralVariants> -e/--expr <expression> -c/--cn <copynumber>\n";  
   print "\n      NOTE:\n\t(1) Genome name must be in the form of abbreviation (e.g. hg18, hg19, mm10, mm9, rn4, dm3)\n";
   print "\t(2) Structutal variants file should be in bedpe format\n";
   print "\t(3) Both expression and copy number data are required to run this tool\n";
   print "\t(4) Region of interest (e.g. promoters, enhancers, UTRs, etc.) file should be in bed format\n";
   print "\nOPTIONS:\n";  

   #print("\t-g/--genome\t\t\tsliding window size\t<int>\t\t[ length of the sliding window. default: 1kb ]\n");
   print("\t-w/--sliding-win-size\t\tsliding window size\t<int>\t\t[ length of the sliding window. default: 100kb ]\n");
   print("\t-s/--sliding-win-step\t\tsliding window step \t<int>\t\t[ step of the sliding window. default: 1kb ]\n");
   print("\t-a/--annot\t\t\tAnnotation file \t<filename>\t[ an annotation file in BED format ]\n");
   print("\t-W/--peakPick-window-size\tpeak calling window \t<int>\t\t[ length of the peak calling window. default: 100bp ]\n");
   print("\t-P/--peakPick-min-sd\t\tpeak calling min. SD\t<int>\t\t[ peak calling minimum standard deviation. default: 5 ]\n");
   print("\t-T/--pct-samples\t\tpercentage of samples\t<int>\t\t[ percentage of samples threshold used to call peaks. default: 5 ]\n");
   print("\t-o/--output\t\t\toutput directory\t<string>\t[ default: ./ ]\n");
   #print("\t-e/--expr-file\t\t\texpression file\t\t<filename>\t[ If expression file is provided, the effect of expression will be tested ]\n");
   #print("\t-c/--cn-file\t\t\tcopy number file\t<filename>\t[ If copy number file is provided, the effect of copy number will be tested ]\n");
   print("\t-p/--pval\t\t\tpvalue cuttoff\t\t<float>\t\t[ pvalue threshold used for differential expression. default: 0.05 ]\n");
   print("\t-G/--genes-of-int\t\tlist of genes\t\t<filename>\t[ If genes-of-interest file is provided, genes will be plotted with peaks ]\n");
   print("\t-r/--region-of-int\t\tregion of interest\t<filename>\t[ If a region of interest file is provided, peaks will be overlapped with this file ]\n");
   print("\t-C/--chrom\t\t\tchromosome name \t<string>\t[ chromosome name used to detect hotspot regions. default: ALL ]\n");
   print("\t-t/--sv-type\t\t\tstructural variant type\t<string>\t[ structural variant type used to detect breakpoints. default: ALL ]\n");
   print("\t-d/--group-dist-size\t\tdistance size\t\t<int>\t\t[ distance used to merge adjacnet peaks to form group peaks. default: 1000bp ]\n");
   print("\t-k/--num-nearby-genes\t\tNumber nearby genes\t<int>\t\t[ Number of nearby genes to the peak. default: 3 ]\n");
   print("\t--t-amp\t\t\t\tamplification threshold\t<float/int>\t[ threshold for copy number amplifications. default: 0.1 ]\n");
   print("\t--t-del\t\t\t\tdeletion threshold\t<float/int>\t[ threshold for copy number deletions. default: -0.1 ]\n");
   print("\t--chip-cov\t\t\tchip-seq coverage\t<filename>\t[ If ChIP-Seq coverage file is provided, peaks will be overlapped with this file ]\n");
   
   print("\t--chip-cov-lbl\t\t\tchip-seq coverage label\t<string>\t[ the chip-seq coverage label used in the plot (e.g. histone name) ]\n");
   print("\t--roi-lbl\t\t\tregion of int. label\t<string>\t[ the region of interest label used in  the plot (e.g. enhancers) ]\n");

   #print("\t--plot-peaks\t\t\tplot peaks\t\t<string>\t[ if this option is enabled (T), all peaks with genes of interest will be plotted. default: F ]\n");
   print("\t--keep-temp\t\t\tkeep intermediate files\t<string>\t[ if this option is enabled (T), all intermediate temporary files will be kept. default: F ]\n");

   print ("\n");
   
   exit 0;
}


