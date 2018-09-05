#!/bin/bash

## Annotate detected peaks 
# Created by: Abdallah Eteleeb <eteleeb@gmial.com>

args=("$G")

genome=$1
#annot=$2
region_of_interest=$2
output_dir=$3
num_nearby_genes=$4
sv_path=$5
#r_path=$6

################# 6. Annotate grouped peaks by overlapping with genes flanked up/downstream by 2kb
#combine all peak groups
cat $output_dir/peaks/*.peak.group.bed | awk '{if($0 !~"pct.samples" && $1 !="") { print $0}}' > $output_dir/peaks/all_peaks.bed    

################# add peak.name column 
awk -v OFS="\t" '{print $0, "p"substr($1,4)"."$4 }' $output_dir/peaks/all_peaks.bed | awk -v OFS="\t" '{print $1,$2,$3,$6,$4,$5}' > $output_dir/peaks/all_peaks2.bed; mv $output_dir/peaks/all_peaks2.bed $output_dir/peaks/all_peaks.bed  

#################################################################################
################# overlap with annotation and extract nearest genes #############
#################################################################################
intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b $output_dir/temp/genes.bed > $output_dir/temp/genes_overlap_with_peaks.tsv   
awk -v OFS="\t" '{print $0, "overlap"}' $output_dir/temp/genes_overlap_with_peaks.tsv > $output_dir/temp/genes_overlap_with_peaks2.tsv; mv $output_dir/temp/genes_overlap_with_peaks2.tsv $output_dir/temp/genes_overlap_with_peaks.tsv

closestBed -D a -io -k $num_nearby_genes -a $output_dir/peaks/all_peaks.bed -b $output_dir/temp/genes.bed > $output_dir/temp/peaks_with_nearest_genes.tsv
awk -v OFS="\t" '{print $0, "nearby"}' $output_dir/temp/peaks_with_nearest_genes.tsv > $output_dir/temp/peaks_with_nearest_genes2.tsv;mv $output_dir/temp/peaks_with_nearest_genes2.tsv $output_dir/temp/peaks_with_nearest_genes.tsv

###### cobmine all peaks in one file 
cat $output_dir/temp/genes_overlap_with_peaks.tsv $output_dir/temp/peaks_with_nearest_genes.tsv > $output_dir/peaks_with_overlap_nearby_genes.tsv
awk -F"\t" '{if($7 ~ "^chr") {print $0}}' $output_dir/peaks_with_overlap_nearby_genes.tsv > $output_dir/peaks_with_overlap_nearby_genes2.tsv
mv $output_dir/peaks_with_overlap_nearby_genes2.tsv $output_dir/peaks_with_overlap_nearby_genes.tsv

#################################################################################
############# overlap peaks with region if interest if provided #################
#################################################################################
if [ "$region_of_interest" != "0" ]; then

  intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b $output_dir/temp/reg_of_int.bed > $output_dir/temp/peaks_overlap_with_region_of_interest.tsv  
  awk -v OFS="\t" '{print $0, "overlap"}' $output_dir/temp/peaks_overlap_with_region_of_interest.tsv > $output_dir/temp/peaks_overlap_with_region_of_interest2.tsv; mv $output_dir/temp/peaks_overlap_with_region_of_interest2.tsv $output_dir/temp/peaks_overlap_with_region_of_interest.tsv 

  closestBed -D a -io -k $num_nearby_genes -a $output_dir/peaks/all_peaks.bed -b $output_dir/temp/reg_of_int.bed > $output_dir/temp/peaks_nearby_region_of_interest.tsv
  awk -v OFS="\t" '{print $0, "nearby"}' $output_dir/temp/peaks_nearby_region_of_interest.tsv > $output_dir/temp/peaks_nearby_region_of_interest2.tsv; mv $output_dir/temp/peaks_nearby_region_of_interest2.tsv $output_dir/temp/peaks_nearby_region_of_interest.tsv 

 ###### cobmine all peaks in one file 
 cat $output_dir/temp/peaks_overlap_with_region_of_interest.tsv $output_dir/temp/peaks_nearby_region_of_interest.tsv > $output_dir/peaks_overlap_nearby_region_of_interest.tsv
 awk -F"\t" '{if($7 ~ "^chr") {print $0}}' $output_dir/peaks_overlap_nearby_region_of_interest.tsv > $output_dir/peaks_overlap_nearby_region_of_interest2.tsv
 mv $output_dir/peaks_overlap_nearby_region_of_interest2.tsv $output_dir/peaks_overlap_nearby_region_of_interest.tsv

fi

##### summarize results for peaks (compute the percentage of sv types falling in each peak )
### overlap preakpoints with all peaks 
intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b $output_dir/temp/all_bp.bed | cut -f1-6,10 | sed 's/\//\t/' > $output_dir/temp/peaks_overlap_bp.tsv

### summarize results 
summarize_peaks_results.r $output_dir $region_of_interest

### end of scripts 










