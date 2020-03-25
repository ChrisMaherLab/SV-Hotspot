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

#define input options with default values 
my $peak="";
my $res_dir="";
my $sv_file=0;
my $output_dir = '/data';
my $expr_file=0; 
my $cn_file=0;
my $t_amp=2.8;
my $t_del= 0.5;
my $region_of_int=0; 
my $chip_cov=0;
#my $chip_cov_lbl="";
my $chip_cov_lbl='chip-seq\ncov.';
#my $roi_lbl ="";
my $left_ext = 0;
my $right_ext = 0;
#my $use_dom = '0';

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
    #'roi-lbl=s' => \$roi_lbl,
    'left-ext=i' => \$left_ext,
    'right-ext=i' => \$right_ext
    #'use-dom=s' => \$use_dom
);

#usage() if (!$peak | !$sv_file | !$res_dir);

print "\n---------------------------------------------\n";
print "            Plotting Peak: $peak\n";
print "---------------------------------------------\n";

### check if chip-cov data was provided 
#if ($chip_cov) { 
#   my $max = `awk '{print \$3-\$2}' $chip_cov | awk 'BEGIN{a=0}{if (\$1>a) a=\$1 fi} END{print a}'`;
#   if ($max < 10000) { 
#      print "\nWarning:\n  it seems the chip coverage file was not averaged using a window approach suggested in the documentation. Visualizing hotspots with the raw chip covergae data results in a very long running time.". 
#            "\n  It is recommended you average chip coverage data using at least 10K window. We have provided a script \"process_chip_cov.r\" for this process. You may consider using it.". 
#            "\n  For more information, please refer to the documentation page on https://github.com/ChrisMaherLab/SV-Hotspot\n\n";
#      exit(0);  
#   } 
#}

if ($expr_file & $cn_file) { 
   system ("Rscript plot_peak_region.r $peak $res_dir $sv_file $output_dir $expr_file $cn_file $region_of_int $chip_cov $t_amp $t_del $chip_cov_lbl $left_ext $right_ext");
} elsif (!$expr_file) {
   system ("Rscript plot_peak_region_with_no_exp.r $peak $res_dir $sv_file $output_dir $cn_file $region_of_int $chip_cov $t_amp $t_del $chip_cov_lbl $left_ext $right_ext");
}

my $output_folder = $output_dir;
$output_folder =~ s/\/data\///;

if ($output_folder ne '/data') {
   print "----------------------------------------------------------------------------------------------------------\n";
   print " Done.\n"; 
   print " Peak(s) plots can be found at \"$output_folder/peaks-plots\" at the local folder you provided with -v option\n";
   print "-----------------------------------------------------------------------------------------------------------\n";
} else {
   print "------------------------------------------------------------------------------------------------------------\n";
   print " Done.\n"; 
   print " Peak(s) plots can be found at \"peaks-plots\" folder at the local folder you provided with -v option\n";
   print "------------------------------------------------------------------------------------------------------------\n";
}


sub usage
{
   #use Term::ANSIColor;
   print "\n";
   print "USAGE:\n      plot-peak.pl [OPTIONS] -p <peakName1,peakName2,...> --sv <structuralVariants> --res-dir <resultsDirectory> -e/--expr <expression> -c/--cn <copynumber>\n";  
   print "\n      NOTE:\n\t(1) Results directory should be the same as the output directory used with sv-hotspot.pl\n";
   #print "\t(2) Structutal variants file should be in bedpe format\n";
   #print "\t(3) Both expression and copy number data are required to run this tool\n";
   #print "\t(4) Region of interest (e.g. promoters, enhancers, UTRs, etc.) file should be in bed format\n";
   print "\nOPTIONS:\n";  

   print("\t-a/--annot\t\t\tAnnotation file \t<filename>\t[ an annotation file in BED format ]\n");
   print("\t-o/--output\t\t\toutput directory\t<string>\t[ default: ./ ]\n");
   print("\t--t-amp\t\t\t\tamplification threshold\t<float/int>\t[ threshold for copy number amplifications. default: 2.8 ]\n");
   print("\t--t-del\t\t\t\tdeletion threshold\t<float/int>\t[ threshold for copy number deletions. default: 0.5 ]\n");
   print("\t--chip-cov\t\t\tchip-seq coverage\t<filename>\t[ If ChIP-Seq coverage file is provided, peaks will be overlapped with this file ]\n");
   print("\t--chip-cov-lbl\t\t\tchip-seq coverage label\t<string>\t[ the chip-seq coverage label used in the plot (e.g. histone name) ]\n");
   #print("\t--roi-lbl\t\t\tregion of int. label\t<string>\t[ the region of interest label used in  the plot (e.g. enhancers) ]\n");
   print("\t--left-ext\t\t\tsize of left extension\t<int>\t\t[ number of extended bases from the left side of the peak. default: 0 ]\n");
   print("\t--right-ext\t\t\tsize of right extension\t<int>\t\t[ number of extended bases from the right side of the peak. default: 0 ]\n");
   #print("\t--use-dom\t\t\tuse dominant SV type\t<T/F>\t\t[ if this option is enabled (T), only the dominant SV type will be ploteed . default: F ]\n");

   print ("\n");
   
   exit 0;
}


