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

####################################################################################
## You need to change this path to the location where you installed the tool  
my $TOOL_PATH='/gscmnt/gc5111/research/eteleeb/projects/SV-HotSpot';
####################################################################################


#define input options wi:qth default values 
my $sv_file=0;
my $genome='hg38';
my $sliding_w_size = 100000;
my $sliding_w_step = 1000; 
my $output_dir = getcwd();
my $annot_file=0;
my $peakPick_win=100;
my $peakPick_minsd=5;
my $pct_samples_t=10;
my $expr_file=0; 
my $cn_file=0;
my $t_amp=2.8;
my $t_del= 0.5;
my $pvalue=0.05;
my $fdr=0.05;
my $genes_of_int=0; 
my $region_of_int=0; 
my $chrom="ALL";
my $sv_type="ALL";
my $merge_dist=10000;
my $plot_top_peaks=10;
my $num_nearby_genes=1;
my $t_stat="wilcox.test"; 
my $chip_cov=0;
my $chip_cov_lbl='chip-seq\ncov.';
my $roi_lbl = '0';
my $left_ext = 0;
my $right_ext = 0;

GetOptions
(
    'sv=s' =>\$sv_file,
    'g|genome=s' =>\$genome,
    'w|sliding-win-size=i' =>\$sliding_w_size,
    's|sliding-win-step=i' =>\$sliding_w_step,
    'o|output=s' =>\$output_dir,
    'a|annot=s' =>\$annot_file,  
    'W|peakPick-window-size=i' =>\$peakPick_win,
    'm|peakPick-min-sd=i' =>\$peakPick_minsd,
    't|pct-samples=i' =>\$pct_samples_t,
    'e|expr=s' =>\$expr_file,
    'c|cn=s' =>\$cn_file,
    't-amp=f' =>\$t_amp,
    't-del=f' =>\$t_del,
    'p|pval=f' =>\$pvalue,
    'q|FDR=f' =>\$fdr,
    'r|region-of-int=s' =>\$region_of_int,
    'G|genes-of-int=s' =>\$genes_of_int,
    'C|chrom=s' => \$chrom,
    'S|sv-type=s' => \$sv_type,
    'd|merge-dist-size=i' => \$merge_dist,
    'k|num-nearby-genes=i' => \$num_nearby_genes,
    'stat-test=s' =>\ $t_stat,
    'chip-cov=s' =>\$chip_cov,
    'plot-top-peaks=i' => \$plot_top_peaks,
    'chip-cov-lbl=s' => \$chip_cov_lbl,
    'roi-lbl=s' => \$roi_lbl,
    'left-ext=i' => \$left_ext,
    'right-ext=i' => \$right_ext
);

usage() if (!$sv_file);

my $max; ### for checking chip-seq coverage data

#################################################################################################################
######################### check if the genome name is in the list of bulit-in genomes ###########################
#################################################################################################################
my $chromsize_file; 
my @genomes = ('hg38', 'hg19', 'hg18', 'mm10', 'mm9', 'rn5','rn6','dm6','dm3');
my $Found=0; 
foreach my $g (@genomes) {
   if ($genome eq $g) {
        $Found = 1;
        last;
   }
}

if (!$Found) { 
    print "\nError:\nThe genome name \"$genome\" is not in the list of built-in genomes. If your genome size is available, please prepare a file with chromosome names and sizes and place it in the annotations folder.\n".
          "For more information on how to extract this file, please refer to the documentation page on https://github.com/ChrisMaherLab/SV-Hotspot\n\n";
    exit(0); 
} else {
  $chromsize_file = $TOOL_PATH.'/annotations/'.$genome.'/chromsize.tsv';
}

#################################################################################################################
##################### check if annotation file was provided otherwise use built-in file #########################
#################################################################################################################
if (!$annot_file) {
     $annot_file = $TOOL_PATH.'/annotations/'.$genome.'/genes.bed'
}

#################################################################################################################
######################################## Creat output directories ###############################################
#################################################################################################################
$output_dir = $output_dir.'/sv-hotspot-output';
system ("rm -rf $output_dir; mkdir $output_dir");
system ("rm -rf $output_dir/processed_data; mkdir $output_dir/processed_data");

#################################################################################################################
########################################## Show all inputs ######################################################
#################################################################################################################
print "########################################################\n";
print "######             SV-HotSpot v1.0.0              ######\n";
print "########################################################\n";
print " SVs file              : $sv_file\n";
print " Genome name           : $genome\n";
print " Sliding window size   : $sliding_w_size\n";
print " Sliding window step   : $sliding_w_step\n"; 
print " Output directory      : $output_dir\n";
print " Annotation file       : $annot_file\n";
print " peakPick window size  : $peakPick_win\n";
print " peakPick minimum SD   : $peakPick_minsd\n";
print " % of samples cutoff   : $pct_samples_t\n";
print " Expression file       : $expr_file\n";  
print " Copy number file      : $cn_file\n";
print " Amplification cutoff  : $t_amp\n";
print " Deletion cutoff       : $t_del\n";
print " FDR threshold         : $fdr\n";
print " Geens of interest     : $genes_of_int\n"; 
print " Region of interest    : $region_of_int\n"; 
print " Chromosomes to analyze: $chrom\n";
print " SV types to analyze   : $sv_type\n"; 
print " Peak merge distance   : $merge_dist\n"; 
print " Number of nearby genes: $num_nearby_genes\n"; 
print " ChIP-Seq cov. file    : $chip_cov\n"; 
print " Statistical test      : $t_stat\n";
print " Plot top peaks        : $plot_top_peaks\n";
print " ChIP-seq coverage     : $chip_cov_lbl\n";
#print " Region of interest    : $roi_lbl\n";
print " Left extension size   : $left_ext\n";
print " Right extension size  : $right_ext\n";
print "########################################################\n\n";

#################################################################################################################
####################################### RUN TOOL FUNCTIONS ######################################################
#################################################################################################################
### set start tim 
my $start = localtime();

verify_input();
prepare_annot();
prepare_SVs();
identify_peaks();
annotate_peaks();
determine_association();
visualize_res();

### set end time 
my $end = localtime();
print "\n";
print "----------------------------------------------------------------------------------\n";
print " SV-Hotspot pipeline finished successfully. Results can be found at $output_dir\n";
print " Started at $start\n";
print " Results reported on $end\n";
print "----------------------------------------------------------------------------------\n";

#################################################################################################################
########################################### Verify input data ###################################################
#################################################################################################################
sub verify_input
{
   
   print "\n--------------------------------------------------\n";
   print "Input Verification Step\n"; 
   print "--------------------------------------------------\n";
   
   ##### check the header of the SV file 
   print " Checking structural varaints file format...";
   open my $sv_f, '<', $sv_file;
   my $sv_header = <$sv_f>;
   if ($sv_header !~ /\bchrom1\b/ || $sv_header!~ /\bstart1\b/ || $sv_header!~ /\bend1\b/ || $sv_header!~ /\bchrom2\b/ || $sv_header!~ /\bstart2\b/ || 
       $sv_header!~ /\bend2\b/ || $sv_header !~ /\bname\b/ || $sv_header!~ /\bscore\b/ || $sv_header!~ /\bstrand1\b/ ||    $sv_header !~ /\bstrand2\b/) {
      print "The header of structural variants file has format different from what the tool accepts.\n";
      exit(0);
   }
   close $sv_f;
   print "PASS.\n";

   ##### check the header of annotation file
   if ($annot_file) { 
        print " Checking annotation file format...";
	open my $annot, '<', $annot_file;
	my $annot_header = <$annot>;
	if ($annot_header !~ /\bchrom\b/ || $annot_header!~ /\bstart\b/ || $annot_header !~ /\bend\b/ || $annot_header !~ /\bgene\b/ || $annot_header !~ /\bscore\b/ || $annot_header !~ /\bstrand\b/) {
      		print "The header of annotation file has a format different from what the tool accepts.\n";
      		exit(0);
   	}
   close $annot; 
   }
   print "PASS.\n";

   ##### check the header of region of interes file
   if ($region_of_int) { 
   print " Checking region of interest file(s) format...";
      my @roi_files = split(",", $region_of_int);
      foreach (@roi_files) {
          open my $roi, '<', $_;
	  my $roi_header = <$roi>;
	  if ($roi_header !~ /\bchrom\b/ || $roi_header !~ /\bstart\b/ || $roi_header !~ /\bend\b/ || $roi_header !~ /\bname\b/) {
      		print "The header of region of interest file has a format different from what the tool accepts.\n";
      		exit(0);
   	  }
          close $roi; 
      }
	
   }
   print "PASS.\n";

   ##### check the header of chip coverage file
   if ($chip_cov) { 
        print " Checking ChIP-Seq coverage file format...";
	open my $chipCov, '<', $chip_cov;
	my $chipCov_header = <$chipCov>;
	if ($chipCov_header !~ /\bchrom\b/ || $chipCov_header !~ /\bstart\b/ || $chipCov_header !~ /\bend\b/ || $chipCov_header !~ /\bcov\b/) {
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
   print "PASS.\n";

   ##### check the header of copy number file 
   if ($cn_file) { 
   print " Checking copy number file format...";
	open my $cn, '<', $cn_file;
	my $cn_header = <$cn>;
	if ($cn_header !~ /\bchrom\b/ || $cn_header !~ /\bstart\b/ || $cn_header !~ /\bend\b/ || $cn_header !~ /\bsample\b/ || $cn_header !~ /\bcn\b/) {
      		print "The header of copy number file has a format different from what the tool accepts.\n";
      		exit(0);
   	}
   close $cn; 
   }
   print "PASS.\n";
   
   #### check if feature name in annotation match the one in the expression 
   if ($annot_file && $expr_file) {
      print " Checking if feature name in the exprsssion matches the one in the annotation file...";
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
   print "PASS.\n";

   #### check if the expression file has no duplicated rows 
   if ($expr_file) {
      print " Checking if the exprsssion file has no duplicated rows...";
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
   print "PASS.\n";

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
   
   ### remove del_dup_sv.bed 
   unlink("$output_dir/processed_data/del_dup_sv.bed");
    
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
   system ("genome-to-sliding-window.r $chromsize_file $sliding_w_size $sliding_w_step $output_dir $chrom");

   print "Overlapping breakpoints with sliding windows\n";
   system ("intersectBed -wao -a $output_dir/processed_data/genome.segments.bed -b $output_dir/processed_data/all_bp.bed > $output_dir/processed_data/genome.segments.with.bps.bed");
   
   ### split the result by chromosme 
   print "\nSplitting the overlapped file by chromosome\n";
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
   system ("find-peaks.r $sv_type $chrom $peakPick_win $peakPick_minsd $pct_samples_t $output_dir $merge_dist $TOOL_PATH $genes_of_int $chromsize_file");

   ### remove intermediate files 
   unlink("$output_dir/processed_data/genome.segments.bed"); 
   unlink("$output_dir/processed_data/genome.segments.with.bps.bed");
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
   print ("Summarizing annotated peaks ...\n");
   system ("summarize_peaks_results.r $output_dir $region_of_int $roi_lbl"); 
}


#################################################################################################################
######################################### Determine association #################################################
#################################################################################################################
sub determine_association
{
   print "\n------------------------------------------------------------------------------------\n";
   print "STEP 3: Determining the association between SVs and the expressssion of neearby genes\n";
   print "-------------------------------------------------------------------------------------\n";

   if ($expr_file & $cn_file) {
     system ("determine_gene_association.r $output_dir $expr_file $cn_file $t_amp $t_del $pvalue $fdr $t_stat");
   } else {
     print "To determine the association between SVs and gene expression, both expression and copy number data are required. \n";
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
      my $pks = `cut -f1 $output_dir/annotated_peaks_summary_final.tsv | grep -v "Peak.name" | head -$plot_top_peaks |  paste -sd, -`;
      chomp($pks);
      #system ("plot_peaks.r $sv_file $output_dir $expr_file $cn_file $chip_cov $t_amp $t_del $chip_cov_lbl $roi_lbl $plot_top_peaks $left_ext $rigth_ext");

      system ("plot_peak_region.r $pks $output_dir $sv_file $output_dir $expr_file $cn_file $region_of_int $chip_cov $t_amp $t_del $chip_cov_lbl $left_ext $right_ext ");
   } else {
      print "Expression and copy number data are required to visualize hotspot reagions.\n";
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
   print "USAGE:\n      sv-hotspot.pl [OPTIONS] -g/--genome <genomeName> --sv <structuralVariants>\n";  
   print "\n      NOTE:\n\t(1) Genome name should be one of the UCSC genome releases (https://genome.ucsc.edu/FAQ/FAQreleases.html#release1). 
            - Built-in Genomes: hg18, hg19, hg38, mm9, mm10, dm3, dm6, rn5, rn6.
            - Please refer to the documentation in case your genome is not listed above. \n";
   print "\t(2) Structutal variants file should be in \"BEDPE\" format.\n";
   print "\t(3) Gene expression data and copy number segments are required to visualize hotspot regions.\n";
   print "\t(4) Region of interest file(s) (e.g. promoters, enhancers, chip-seq, etc.) should be in \"BED\" format\n";
   print "\nOPTIONS:\n";  

   print("\t-w/--sliding-win-size\t\tsliding window size\t<int>\t\t[ sliding window size. default: 100kb ]\n");
   print("\t-s/--sliding-win-step\t\tsliding window step \t<int>\t\t[ sliding window step. default: 1kb ]\n");
   print("\t-a/--annot\t\t\tannotation file \t<filename>\t[ an annotation file in \"BED\" format ]\n");
   print("\t-e/--exp\t\t\texpression file \t<filename>\t[ expression file in a \"matrix\" format ]\n");
   print("\t-c/--cn\t\t\t\tCopy number file \t<filename>\t[ copy number segments in BED format ]\n");
   print("\t-W/--peakPick-window-size\tpeak calling window \t<int>\t\t[ peakPick window size. default: 100bp ]\n");
   print("\t-m/--peakPick-min-sd\t\tpeak calling min. SD\t<int>\t\t[ peakPick standard deviation cutoff. default: 5 ]\n");
   print("\t-t/--pct-samples\t\tpercentage of samples\t<int>\t\t[ percentage of samples cutoff to call peaks. default: 10 ]\n");
   print("\t-o/--output\t\t\toutput directory\t<string>\t[ results output directory. default: ./ ]\n");
   print("\t-p/--pval\t\t\tpvalue cuttoff\t\t<float>\t\t[ pvalue threshold used for significance. default: 0.05 ]\n");
   print("\t-q/--FDR\t\t\tqvalue cuttoff\t\t<float>\t\t[ qvalue (minimum FDR) cutoff to call significant genes. default: 0.05 ]\n");
   print("\t-G/--genes-of-int\t\tlist of genes\t\t<filename>\t[ list of genes of interest to be used for visualization ]\n");
   print("\t-r/--region-of-int\t\tregion(s) of interest\t<filename>\t[ region of interest file(s) in \"BED\" foramt separated by comma ]\n");
   print("\t-C/--chrom\t\t\tchromosome name \t<string>\t[ chromosome name used to detect hotspots. default: ALL ]\n");
   print("\t-t/--sv-type\t\t\tstructural variant type\t<string>\t[ SV type used to detect hotspots. default: ALL ]\n");
   print("\t-d/--merge-dist-size\t\tdistance size\t\t<int>\t\t[ distance cutoff used to merge adjacnet peaks. default: 10kb ]\n");
   print("\t-k/--num-nearby-genes\t\tNumber nearby genes\t<int>\t\t[ number of up/downstream genes to the peak. default: 1 ]\n");
   print("\t--t-amp\t\t\t\tamplification threshold\t<float/int>\t[ amplification cutoff. default: 2.8 ]\n");
   print("\t--t-del\t\t\t\tdeletion threshold\t<float/int>\t[ deletion cutoff. default: 0.5 ]\n");
   print("\t--stat-test\t\t\tstatistical test\t<string>\t[ statistical test used for comparison (wilcox.test or t.test). default: wilcox.test ]\n");
   print("\t--chip-cov\t\t\tchip-seq coverage\t<filename>\t[ chip-seq coverage file in \"BED\" foramt ]\n");
   print("\t--chip-cov-lbl\t\t\tchip-seq coverage label\t<string>\t[ label used for chip-seq coverage ]\n");
   #print("\t--roi-lbl\t\t\tregion of int. label(s)\t<string>\t[ labels used for region(s) of interest separated by comma  ]\n");
   print("\t--plot-top-peaks\t\tplot top # peaks\t<int>\t\t[ number of top peaks to plot. default: top 10 ]\n");
   print("\t--left-ext\t\t\tsize of left extension\t<int>\t\t[ size of the left extension of the peak. default: 0bp ]\n");
   print("\t--right-ext\t\t\tsize of right extension\t<int>\t\t[ size of the right extension of the peak. default: 0bp ]\n");

   #print("\t--keep-processed_data\t\t\tkeep intermediate files\t<string>\t[ if this option is enabled (T), all intermediate processed_dataorary files will be kept. default: F ]\n");

   print ("\n");
   
   exit 0;
}


