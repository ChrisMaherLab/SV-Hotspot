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
my $peak="";
my $res_dir="";
my $sv_file='0';
my $output_dir = getcwd();
my $expr_file='0'; 
my $cn_file='0';
my $t_amp=2.8;
my $t_del= 0.5;
my $region_of_int='0'; 
my $chip_cov='0';
my $chip_cov_lbl="";
my $roi_lbl ="";
my $left_ext = 0;
my $rigth_ext = 0;

GetOptions
(
    'p|peak=s' => \$peak,
    'd|res-dir=s' =>\$res_dir,
    'sv=s' =>\$sv_file,
    'o|output=s' =>\$output_dir,
    'e|expr=s' =>\$expr_file,
    'c|cn=s' =>\$cn_file,
    't-amp=f' =>\$t_amp,
    't-del=f' =>\$t_del,
    'r|region-of-int=s' =>\$region_of_int,
    'chip-cov=s' =>\$chip_cov,
    'chip-cov-lbl=s' => \$chip_cov_lbl,
    'roi-lbl=s' => \$roi_lbl,
    'left-ext=i' => \$left_ext,
    'rigth-ext=i' => \$rigth_ext
);

usage() if (!$peak | !$sv_file | !$res_dir | !$expr_file | !$cn_file );

print "\n---------------------------------------------\n";
print "            Plotting Peak: $peak\n";
print "---------------------------------------------\n";

if ($expr_file && $cn_file) {
   # check if the user provided chip-seq data, if yes run the script to prcess it (average over a windwo)
   #if ($chip_cov) {
   #   print "Processing chip-seq data, please wait as this may take several minutes\n";
   #   system("process_chip_data.r $chip_cov $output_dir");
   #} 
   system ("plot_peak_region.r $peak $res_dir $sv_file $output_dir $expr_file $cn_file $chip_cov $t_amp $t_del $chip_cov_lbl $roi_lbl $left_ext $rigth_ext");
} else {
  print "Both expression and copy number data are required to generate visualization\n";
  exit(0);
}


sub usage
{
   #use Term::ANSIColor;
   print "\n";
   print "USAGE:\n      plot-peak [OPTIONS] -p <peakName> --sv <structuralVariants> --res-dir <resultsDirectory> -e/--expr <expression> -c/--cn <copynumber>\n";  
   print "\n      NOTE:\n\t(1) Results directory should be the same as output directory used when you run sb-hotspot tool\n";
   #print "\t(2) Structutal variants file should be in bedpe format\n";
   #print "\t(3) Both expression and copy number data are required to run this tool\n";
   #print "\t(4) Region of interest (e.g. promoters, enhancers, UTRs, etc.) file should be in bed format\n";
   print "\nOPTIONS:\n";  

   print("\t-a/--annot\t\t\tAnnotation file \t<filename>\t[ an annotation file in BED format ]\n");
   print("\t-o/--output\t\t\toutput directory\t<string>\t[ default: ./ ]\n");
   print("\t--t-amp\t\t\t\tamplification threshold\t<float/int>\t[ threshold for copy number amplifications. default: 0.1 ]\n");
   print("\t--t-del\t\t\t\tdeletion threshold\t<float/int>\t[ threshold for copy number deletions. default: -0.1 ]\n");
   print("\t--chip-cov\t\t\tchip-seq coverage\t<filename>\t[ If ChIP-Seq coverage file is provided, peaks will be overlapped with this file ]\n");
   print("\t--chip-cov-lbl\t\t\tchip-seq coverage label\t<string>\t[ the chip-seq coverage label used in the plot (e.g. histone name) ]\n");
   print("\t--roi-lbl\t\t\tregion of int. label\t<string>\t[ the region of interest label used in  the plot (e.g. enhancers) ]\n");

   print ("\n");
   
   exit 0;
}


