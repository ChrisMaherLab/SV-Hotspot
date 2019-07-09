#!/bin/bash

## Annotate detected peaks 
# Created by: Abdallah Eteleeb <eteleeb@wustl.com> and Ha X. Dang <hdang@wustl.edu>

args=("$G")

genome=$1
#annot=$2
region_of_interest=$2
output_dir=$3
num_nearby_genes=$4

################# Annotate grouped peaks by overlapping with genes flanked up/downstream by 2kb
#combine all peak groups
cat $output_dir/peaks/*.peak.group.bed | awk '{if($0 !~"pct.samples" && $1 !="") { print $0}}' > $output_dir/peaks/all_peaks.bed    

################# add peak.name column 
awk -v OFS="\t" '{print $0, "p"substr($1,4)"."$4 }' $output_dir/peaks/all_peaks.bed | awk -v OFS="\t" '{print $1,$2,$3,$8,$4,$5,$6,$7}' > $output_dir/peaks/all_peaks2.bed; mv $output_dir/peaks/all_peaks2.bed $output_dir/peaks/all_peaks.bed  

#################################################################################
################# overlap with annotation and extract nearest genes #############
#################################################################################
### get overlap genes
intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b $output_dir/processed_data/genes.bed > $output_dir/processed_data/genes_overlap_with_peaks.tsv   
awk -v OFS="\t" '{print $0, "overlap"}' $output_dir/processed_data/genes_overlap_with_peaks.tsv > $output_dir/processed_data/genes_overlap_with_peaks2.tsv; mv $output_dir/processed_data/genes_overlap_with_peaks2.tsv $output_dir/processed_data/genes_overlap_with_peaks.tsv

### get up/downstream genes
# produce two peaks with different strand from a non-stranded peak
cat $output_dir/peaks/all_peaks.bed | awk -v OFS="\t" '{print $1,$2,$3,$4,".","+",$5,$6,$7,$8;print $1,$2,$3,$4,".","-",$5,$6,$7,$8}' > $output_dir/processed_data/all_peaks_with_strands.bed

#closestBed -D a -io -k $num_nearby_genes -a $output_dir/peaks/all_peaks.bed -b $output_dir/processed_data/genes.bed > $output_dir/processed_data/peaks_with_nearest_genes.tsv
closestBed -D a -io -iu -s -k $num_nearby_genes -a $output_dir/processed_data/all_peaks_with_strands.bed -b $output_dir/processed_data/genes.bed | cut -f1-4,7- > $output_dir/processed_data/peaks_with_nearest_genes.tsv
awk -v OFS="\t" '{print $0, "nearby"}' $output_dir/processed_data/peaks_with_nearest_genes.tsv > $output_dir/processed_data/peaks_with_nearest_genes2.tsv;mv $output_dir/processed_data/peaks_with_nearest_genes2.tsv $output_dir/processed_data/peaks_with_nearest_genes.tsv

###### cobmine all peaks in one file 
cat $output_dir/processed_data/genes_overlap_with_peaks.tsv $output_dir/processed_data/peaks_with_nearest_genes.tsv > $output_dir/processed_data/peaks_with_overlap_nearby_genes.tsv
awk -F"\t" '{if($9 ~ "^chr") {print $0}}' $output_dir/processed_data/peaks_with_overlap_nearby_genes.tsv > $output_dir/processed_data/peaks_with_overlap_nearby_genes2.tsv
mv $output_dir/processed_data/peaks_with_overlap_nearby_genes2.tsv $output_dir/processed_data/peaks_with_overlap_nearby_genes.tsv

rm $output_dir/processed_data/genes_overlap_with_peaks.tsv 
rm $output_dir/processed_data/peaks_with_nearest_genes.tsv
rm $output_dir/processed_data/all_peaks_with_strands.bed

#################################################################################
############# overlap peaks with region if interest if provided #################
#################################################################################

### extract the region of interest files 
IFS=',' read -r -a roi <<< "$region_of_interest"

if [ "${#roi[@]}" != "0" ]; then
   rm -rf $output_dir/processed_data/ROIs; mkdir $output_dir/processed_data/ROIs
   for item in "${roi[@]}"
   do
     roi_name=$(basename ${item%.*})
     grep -v "chrom" $item | intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b - > $output_dir/processed_data/ROIs/peaks_overlap_with_$roi_name.tsv
   done 

  #grep -v "chrom" $item | intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b - > $output_dir/processed_data/ROIs/peaks_overlap_with_region_of_interest.tsv

  ### combine regions of interest files 
  combine_roi_files.r $output_dir



#  intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b $output_dir/processed_data/reg_of_int.bed > $output_dir/processed_data/peaks_overlap_with_region_of_interest.tsv  
#  awk -v OFS="\t" '{print $0, "overlap"}' $output_dir/processed_data/peaks_overlap_with_region_of_interest.tsv > $output_dir/processed_data/peaks_overlap_with_region_of_interest2.tsv
#  mv $output_dir/processed_data/peaks_overlap_with_region_of_interest2.tsv $output_dir/processed_data/peaks_overlap_with_region_of_interest.tsv

#  closestBed -D a -io -k $num_nearby_genes -a $output_dir/peaks/all_peaks.bed -b $output_dir/processed_data/reg_of_int.bed > $output_dir/processed_data/peaks_nearby_region_of_interest.tsv
#  awk -v OFS="\t" '{print $0, "nearby"}' $output_dir/processed_data/peaks_nearby_region_of_interest.tsv > $output_dir/processed_data/peaks_nearby_region_of_interest2.tsv; mv $output_dir/processed_data/peaks_nearby_region_of_interest2.tsv $output_dir/processed_data/peaks_nearby_region_of_interest.tsv 

 ###### cobmine all peaks in one file 
# cat $output_dir/processed_data/peaks_overlap_with_region_of_interest.tsv $output_dir/processed_data/peaks_nearby_region_of_interest.tsv > $output_dir/peaks_with_overlap_nearby_region_of_interest.tsv
# awk -F"\t" '{if($9 ~ "^chr") {print $0}}' $output_dir/peaks_with_overlap_nearby_region_of_interest.tsv > $output_dir/peaks_with_overlap_nearby_region_of_interest2.tsv
# mv $output_dir/peaks_with_overlap_nearby_region_of_interest2.tsv $output_dir/peaks_with_overlap_nearby_region_of_interest.tsv

 #rm $output_dir/processed_data/peaks_overlap_with_region_of_interest.tsv 
 #rm $output_dir/processed_data/peaks_nearby_region_of_interest.tsv

fi

##### summarize results for peaks (compute the percentage of sv types falling in each peak )
### overlap preakpoints with all peaks 
intersectBed -wo -a $output_dir/peaks/all_peaks.bed -b $output_dir/processed_data/all_bp.bed | cut -f1-7,12 | sed 's/\//\t/' > $output_dir/processed_data/peaks_overlap_bp.tsv

### summarize results 
#summarize_peaks_results.r $output_dir $region_of_interest $roi_lbl

### end of scripts 










