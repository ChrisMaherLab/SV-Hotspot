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

#define input options wi:qth default values 
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
my $plot_top_peaks=10;
my $num_nearby_genes=1;
my $chip_cov='0';
my $chip_cov_lbl="";
my $roi_lbl ="";
my $left_ext = 0;
my $rigth_ext = 0;

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
    'plot-top-peaks=s' => \$plot_top_peaks,
    'chip-cov-lbl=s' => \$chip_cov_lbl,
    'roi-lbl=s' => \$roi_lbl,
    'left-ext=i' => \$left_ext,
    'rigth-ext=i' => \$rigth_ext
);

usage() if (!$sv_file | !$genome | !$expr_file | !$cn_file);

#### define the name of the chromosme lengths file and output directroy 
my $chromsize_file = $TOOL_PATH.'/annotations/'.$genome.'/chromsize.tsv';

#### check if annotation file was provided otherwise use built-in file
if (!$annot_file) {
     $annot_file = $TOOL_PATH.'/annotations/'.$genome.'/genes.bed'
}

### define varaible for checking chip coverage data 
my $max; 

### Creating output directories";
$output_dir = $output_dir.'/sv-hotspot-output';
system ("rm -rf $output_dir; mkdir $output_dir");
system ("rm -rf $output_dir/processed_data; mkdir $output_dir/processed_data");

#### show all inputes 
print "########################################################\n";
print "######             SV-HotSpot v1.0.0              ######\n";
print "########################################################\n";
print " Structural variants file: $sv_file\n";
print " Genome: $genome\n";
print " Sliding window size : $sliding_w_size\n";
print " Sliding window step : $sliding_w_step\n"; 
print " Output directory    : $output_dir\n";
print " Annotation file     : $annot_file\n";
print " peakPick window size: $peakPick_win\n";
print " peakPick minimum SD : $peakPick_minsd\n";
print " Perc. of samples threshold: $pct_samples_t\n";
print " Expression file: $expr_file\n";  
print " Copy number file: $cn_file\n";
print " Ampl:qification threshold: $t_amp\n";
print " Deletion threshold: $t_del\n";
print " P-value threshold: $pvalue\n";
print " Geens of interest: $genes_of_int\n"; 
print " Region of interest: $region_of_int\n"; 
print " Chromosmes to analyze: $chrom\n";
print " SV types to analyze: $sv_type\n"; 
print " Distance for merging peaks: $merge_dist\n"; 
print " Number of nearby genes: $num_nearby_genes\n"; 
print " chip-Seq file: $chip_cov\n"; 
print " Plot top peaks: $plot_top_peaks\n";
print " Chip coverage label: $chip_cov_lbl\n";
print " Region of interest label: $roi_lbl\n";
print " Region left extension: $left_ext\n";
print " Region right extension: $rigth_ext\n";
print "########################################################\n\n";

############################################ RUN TOOL FUNCTIONS #################################################
verify_input();
prepare_annot();
prepare_SVs();
identify_peaks();
annotate_peaks();
determine_effect();
visualize_res();
#################################################################################################################


#################################################################################################################
########################################### Verify input data ###############################################
#################################################################################################################
sub verify_input
{
   
   print "\n--------------------------------------------------\n";
   print "------- Checking the format of input files -------\n"; 
   print "--------------------------------------------------\n";
   
   ##### check the header of the SV file 
   open my $sv_f, '<', $sv_file;
   my $sv_header = <$sv_f>;
   if ($sv_header !~ /chrom1/ || $sv_header!~ /start1/ || $sv_header!~ /end1/ || $sv_header!~ /chrom2/ || $sv_header!~ /start2/ || 
       $sv_header!~ /end2/ || $sv_header!~ /name/ || $sv_header!~ /score/ || $sv_header!~ /strand1/ ||    $sv_header !~ /strand2/) {
      print "The header of structural variants file has format different from what the tool accepts.\n";
      exit(0);
   }
   close $sv_f;

   ##### check the header of annotation file
   if ($annot_file) { 
	open my $annot, '<', $annot_file;
	my $annot_header = <$annot>;
	if ($annot_header !~ /chrom/ || $annot_header!~ /start/ || $annot_header !~ /end/ || $annot_header !~ /gene/ || $annot_header !~ /score/ || $annot_header !~ /strand/) {
      		print "The header of annotation file has a format different from what the tool accepts.\n";
      		exit(0);
   	}
   close $annot; 
   }
   
   ##### check the header of region of interes file
   if ($region_of_int) { 
	open my $roi, '<', $region_of_int;
	my $roi_header = <$roi>;
	if ($roi_header !~ /chrom/ || $roi_header !~ /start/ || $roi_header !~ /end/ || $roi_header !~ /name/) {
      		print "The header of region of interest file has a format different from what the tool accepts.\n";
      		exit(0);
   	}
        close $roi; 
   }
   
   ##### check the header of chip coverage file
   if ($chip_cov) { 
	open my $chipCov, '<', $chip_cov;
	my $chipCov_header = <$chipCov>;
	if ($chipCov_header !~ /chrom/ || $chipCov_header !~ /start/ || $chipCov_header !~ /end/ || $chipCov_header !~ /cov/) {
      		print "The header of chip coverage file has a format different from what the tool accepts.\n";
      		exit(0);
   	}
   	close $chipCov; 

   	### check if the chip coverage file is an averaged file 
   	$max = `awk '{print \$3-\$2}' $chip_cov | awk 'BEGIN{a=0}{if (\$1>a) a=\$1 fi} END{print a}'`;
        $max = $max + 1;
   	if ($max < 1000) { 
        	print "\nWarning:\n  it seems the chip coverage file was not averaged using a window approach suggested in the documentation. Visualizing hotspots with the raw chip covergae data results in a very long running time.". 
            "\n  It is recommended you average chip coverage data using a window size range form 1-10k. We have provided a script \"process_chip_cov.r\" for this process. You may consider using it.\n\n"; 
   	} 
   }
  

   ##### check the header of copy number file 
   if ($cn_file) { 
	open my $cn, '<', $cn_file;
	my $cn_header = <$cn>;
	if ($cn_header !~ /chrom/ || $cn_header !~ /start/ || $cn_header !~ /end/ || $cn_header !~ /sample/ || $cn_header !~ /cn/) {
      		print "The header of copy number file has a format different from what the tool accepts.\n";
      		exit(0);
   	}
   close $cn; 
   }
   
   #### check if feature name in annotation match the one in the expression 
   if ($annot_file && $expr_file) {
      open( my $annot, '<', $annot_file) or die "Can't open file: $!";
      open( my $exp, '<', $expr_file) or die "Can't open file: $!";
      my $annot_header = <$annot>; 
      chomp($annot_header);
      my @annot_cols = split( "\t", $annot_header ); 
      my $exp_header = <$exp>;
      my @feature = split( "\t", $exp_header );  
      my $Found = "0"; 
      foreach (@annot_cols) {
         if ($_ eq $feature[0]) { $Found = "1" } 
      }
      if (!$Found) {
       	 print "Feature name in expression file doesn't match feature name in the annotation file!.\n";
      	 exit(0);
      }
   close $annot; close $exp;
   }
 
   #### check if the expression file has no duplicated rows 
   if ($expr_file) {
      open( my $exp, '<', $expr_file) or die "Can't open file: $!";
      my @genes; 
      my %count;
      while (my $line = <$exp>) {
        next if $. == 1; 
        my @row = split( "\t", $line );
        push (@genes, $row[0]);
      }
      close $exp;
      for(@genes) {
         $count{$_}++;
         if ($count{$_} > 1) { 
            print "Expression file has duplicated rows!.\n";
            exit(0); 
         }
      }
   }

   print "\n--------------------------------------------------\n";
   print "\tALL INPUT DATA LOOKS GOOD.\n";
   print "--------------------------------------------------\n";
}


#################################################################################################################
############################################ Preapre Annotation #################################################
#################################################################################################################
sub prepare_annot 
{
   ### sort gene annotation file 
   system("grep -v chrom $annot_file | sort -k1,1 -k2,2n > $output_dir/processed_data/genes_sorted.bed; mv $output_dir/processed_data/genes_sorted.bed $output_dir/processed_data/genes.bed");
   $annot_file = $output_dir."/processed_data/genes.bed";
 
   ### prepare genes of interest file 
   if ($genes_of_int) {
   #   $genes_of_int = $TOOL_PATH.'/annotations/genes-of-interest.txt';
   #   system("cat $genes_of_int | select-rows.pl 0 $annot_file 3 > $output_dir/processed_data/genes-of-interest.bed");
   #   $genes_of_int = $output_dir."/processed_data/genes-of-interest.bed"
   #} else {
      system("cat $genes_of_int | select-rows.pl 0 $annot_file 3 > $output_dir/processed_data/genes-of-interest.bed");
      $genes_of_int = $output_dir."/processed_data/genes-of-interest.bed"
   }

}

#################################################################################################################
############################################# Preapre SVs data ##################################################
#################################################################################################################
sub prepare_SVs
{
   system("cat $sv_file | awk -F\"\\t\" '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$7,\$8,\$9; print \$4,\$5,\$6,\$7,\$8,\$10}' | grep -v \"chrom\" > $output_dir/processed_data/all_bp.bed");
   system ("awk '{if(\$2!=0 || \$3!=0) {print \$0}}' $output_dir/processed_data/all_bp.bed > $output_dir/processed_data/all2_bp.bed; mv $output_dir/processed_data/all2_bp.bed $output_dir/processed_data/all_bp.bed");
   #system ("for type in `cut -f4 $output_dir/processed_data/all_bp.bed | cut -f2 -d'/' | sort | uniq`; do cat $output_dir/processed_data/all_bp.bed | grep \"\/\$type\" > $output_dir/processed_data/bp.\$type.bed; done");

 ### generate bedpe file for DEL and DUP events
   system ("generate-bedpe-for-DEL-DUP.r $sv_file $output_dir");
   system ("cat $output_dir/processed_data/all_bp.bed $output_dir/processed_data/del_dup_sv.bed > $output_dir/processed_data/all_bp_tmp.bed; mv $output_dir/processed_data/all_bp_tmp.bed $output_dir/processed_data/all_bp.bed");

}


#################################################################################################################
############################################# Identify Peaks ####################################################
#################################################################################################################
sub identify_peaks
{
   print "\n--------------------------------------------------\n";
   print "STEP 1: Identifying Peaks (hotspot regions) \n";
   print "--------------------------------------------------\n";
   print "\nSegmenting the genome into sliding windows\n";
   system ("genome-to-sliding-window.r $chromsize_file $sliding_w_size $sliding_w_step $output_dir");

   print "Overlapping breakpoints with sliding windows\n";
   system ("intersectBed -wao -a $output_dir/processed_data/genome.segments.bed -b $output_dir/processed_data/all_bp.bed > $output_dir/processed_data/genome.segments.with.bps.bed");
   
   ### split the result by chromosme 
   print "\nSplitting the file by chromosome\n";
   system ("rm -rf $output_dir/processed_data/segments_with_bps_per_chr; mkdir $output_dir/processed_data/segments_with_bps_per_chr");
   system("awk '{print \$0 >> \"$output_dir/processed_data/segments_with_bps_per_chr/\"\$1\".bed\"}' $output_dir/processed_data/genome.segments.with.bps.bed");

   ### extract chromosme file names and summarize counts 
   system ("rm -rf $output_dir/processed_data/counts; mkdir $output_dir/processed_data/counts");

   opendir chrsDir, "$output_dir/processed_data/segments_with_bps_per_chr" || die "$!"; 
   my @chr_files = grep {/^chr*.*?\.bed?/} readdir chrsDir; 
   close chrsDir;
   @chr_files = sort @chr_files;
   foreach (@chr_files) { 
   	my $chr = $_;
	system ("summarize-sample-count.r $chr $output_dir"); 
   }
   
   ### combine all counts for all chromosomes
   system("combine-counts-files.r $output_dir");

   ### remove folders 
   system ("rm -rf $output_dir/processed_data/counts");
   system ("rm -rf $output_dir/processed_data/segments_with_bps_per_chr");
   
   print "\nCall structural variant peaks (hotspots)\n";
   system ("find-peaks.r $sv_type $chrom $peakPick_win $peakPick_minsd $pct_samples_t $output_dir $merge_dist $TOOL_PATH $genes_of_int");
}


#################################################################################################################
############################################# Annotate Peaks ####################################################
#################################################################################################################
sub annotate_peaks
{
   print "\n--------------------------------------------------\n";
   print "STEP 2: Annotating Peaks \n";
   print "--------------------------------------------------\n";
  
   system("annotate_peaks.sh $genome $region_of_int $output_dir $num_nearby_genes $TOOL_PATH"); 
}


#################################################################################################################
############################################# Determine effect ##################################################
#################################################################################################################
sub determine_effect
{
   print "\n--------------------------------------------------\n";
   print "STEP 3: Determining the effect of SVs on gene expression\n";
   print "--------------------------------------------------\n";

   if ($expr_file & $cn_file) {
     system ("determine_gene_association.r $output_dir $expr_file $cn_file $t_amp $t_del $pvalue");
   } else {
     print "To determine the effect of SVs on gene expression, both expression and copy number data are required. \n";
   exit(0); 
   }

}


#################################################################################################################
############################################# Visualize results #################################################
#################################################################################################################
sub visualize_res
{
   print "\n--------------------------------------------------\n";
   print "STEP 4: Visualizing hotspot regions \n";
   print "--------------------------------------------------\n";

   ### check if chip-cov data was provided 
   if ($chip_cov) { 
     if ($max < 1000) { 
      print "\nWarning:\n  it seems the chip coverage file was not averaged using a window approach suggested in the documentation. Visualizing hotspots with the raw chip covergae data results in a very long running time.". 
            "\n  It is recommended you average chip coverage data using at least 10K window. We have provided a script \"process_chip_cov.r\" for this process. You may consider using it.". 
            "\n  For more information, please refer to the documentation page on https://github.com/ChrisMaherLab/SV-Hotspot\n\n";
      exit(0);  
     } 
   }

   if ($expr_file && $cn_file) {
      ### extract top peaks 
      my $pks = `cut -f1 $output_dir/annotated_peaks_summary_final.tsv | grep -v "p.name" | head -$plot_top_peaks |  paste -sd, -`;
      chomp($pks);
      #system ("plot_peaks.r $sv_file $output_dir $expr_file $cn_file $chip_cov $t_amp $t_del $chip_cov_lbl $roi_lbl $plot_top_peaks $left_ext $rigth_ext");
      system ("plot_peak_region.r $pks $output_dir $sv_file $output_dir $expr_file $cn_file $region_of_int $chip_cov $t_amp $t_del $chip_cov_lbl $roi_lbl $left_ext $rigth_ext");
   } else {
      print "Both expression and copy number data are required to generate visualization\n";
      exit(0);
   }

}


#################################################################################################################
################################################# TOOL USAGE ####################################################
#################################################################################################################
sub usage
{
   #use Term::ANSIColor;
   print "\n";
   print "USAGE:\n      sv-hotspot.pl [OPTIONS] -g/--genome <genomeName> --sv <structuralVariants> -e/--expr <expression> -c/--cn <copynumber>\n";  
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

   print("\t--plot-top-peaks\t\tplot top # peaks\t<int>\t\t[ Plot the top # number of peaks. default: 10 ]\n");
   print("\t--left-ext\t\t\tsize of left extension\t<int>\t\t[ number of extended bases from the left side of the peak. default: 0 ]\n");
   print("\t--right-ext\t\t\tsize of right extension\t<int>\t\t[ number of extended bases from the right side of the peak. default: 0 ]\n");

   #print("\t--keep-processed_data\t\t\tkeep intermediate files\t<string>\t[ if this option is enabled (T), all intermediate processed_dataorary files will be kept. default: F ]\n");

   print ("\n");
   
   exit 0;
}


