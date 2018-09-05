#!/usr/bin/perl -w

# Plot peak: this scripr is used to plot a structural varaint peak  
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
my $peak="";i
my $res_dir="";
my $sv_file='0';
my $output_dir = getcwd();
my $annot_file='0';
my $expr_file='0'; 
my $cn_file='0';
my $t_amp=2.8;
my $t_del= 0.5;
my $region_of_int='0'; 
my $chip_cov='0';
my $chip_cov_lbl="";
my $roi_lbl ="";

GetOptions
(
    'p|peak=s' => \$peak,
    'd|res-dir=s' =>\$res_dir,
    'sv=s' =>\$sv_file,
    'o|output=s' =>\$output_dir,
    'a|annot=s' =>\$annot_file,  
    'e|expr=s' =>\$expr_file,
    'c|cn=s' =>\$cn_file,
    't-amp=f' =>\$t_amp,
    't-del=f' =>\$t_del,
    'r|region-of-int=s' =>\$region_of_int,
    'chip-cov=s' =>\$chip_cov,
    'chip-cov-lbl=s' => \$chip_cov_lbl,
    'roi-lbl=s' => \$roi_lbl
);

usage() if (!$sv_file | !$peak | !$expr_file | !$cn_file | !$res_dir);

#### check if annotation file and genes of interest file were provided otherwise use defaults
if (!$annot_file) {
   $annot_file = $TOOL_PATH.'/annotations/'.$genome.'/genes.bed'
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
   system ("plot_peak_region.r $sv_file $output_dir $expr_file $cn_file $chip_cov $t_amp $t_del $chip_cov_lbl $roi_lbl");
} else {
  print "Both expression and copy number data are required to generate visualization\n";
  exit(0);
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


