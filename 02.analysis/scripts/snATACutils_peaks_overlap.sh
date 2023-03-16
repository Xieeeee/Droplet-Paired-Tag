#!/bin/bash
### merge peaks using Yang Li's method
### https://github.com/yal054/snATACutils/tree/master/03.peak_calling

current=/home/y2xie/scratch/28.FC_16k_Droplet_PT_221210
peak_dir=${current}/09.macs2
script_dir=${current}/scripts

### set names
REP1_PEAK=${peak_dir}/DNA_cluster/corrected/FC_H3K27ac_valid_YX794_peaks.narrowPeak
REP2_PEAK=${peak_dir}/DNA_cluster/corrected/FC_H3K27ac_valid_YX835_peaks.narrowPeak
REP3_PEAK=${peak_dir}/DNA_cluster/corrected/FC_H3K27ac_valid_YX837_peaks.narrowPeak
REP1_SUMMIT=${peak_dir}/DNA_cluster/corrected/FC_H3K27ac_valid_YX794_summits.bed
REP2_SUMMIT=${peak_dir}/DNA_cluster/corrected/FC_H3K27ac_valid_YX835_summits.bed
REP3_SUMMIT=${peak_dir}/DNA_cluster/corrected/FC_H3K27ac_valid_YX837_summits.bed
POOLED_SUMMIT=${peak_dir}/DNA_cluster/corrected/rep_merged/FC_H3K27ac_replicates_merged_summits_total.bed
naiveSummitList=${peak_dir}/DNA_cluster/corrected/rep_merged/FC_H3K27ac_replicates_merged_summits.bed
REP1_REP2_PEAK=${peak_dir}/DNA_cluster/corrected/rep_merged/YX794_YX835_valid_peaks.narrowPeak
REP1_REP3_PEAK=${peak_dir}/DNA_cluster/corrected/rep_merged/YX794_YX837_valid_peaks.narrowPeak
REP2_REP3_PEAK=${peak_dir}/DNA_cluster/corrected/rep_merged/YX835_YX837_valid_peaks.narrowPeak

mkdir ${peak_dir}/DNA_cluster/corrected/rep_merged
cat ${REP1_SUMMIT} ${REP2_SUMMIT} ${REP3_SUMMIT} > ${POOLED_SUMMIT}

intersectBed -wo -a ${REP1_PEAK} -b ${REP2_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.3) || ($21/s2 >= 0.3)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq > ${REP1_REP2_PEAK}
intersectBed -wo -a ${REP1_PEAK} -b ${REP3_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.3) || ($21/s2 >= 0.3)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq > ${REP1_REP3_PEAK}
intersectBed -wo -a ${REP2_PEAK} -b ${REP3_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.3) || ($21/s2 >= 0.3)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq > ${REP2_REP3_PEAK}

join -1 1 -2 4 <(cat ${REP1_REP2_PEAK} ${REP1_REP3_PEAK} ${REP2_REP3_PEAK} | cut -f 4 - | sort) <(sort -k4,4 ${POOLED_SUMMIT}) -t$'\t' | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1,$5}' | sort -k1,1 -k2,2n > ${naiveSummitList} 

echo -e "FC_H3K27ac_2rep\t${naiveSummitList}" > ${peak_dir}/DNA_cluster/corrected/rep_merged/FC_H3K27ac_replicates_merged_summits.list.txt
/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/Rscript ${script_dir}/iterative_overlap_peak_merging.R --size 500 -i ${peak_dir}/DNA_cluster/corrected/rep_merged/FC_H3K27ac_replicates_merged_summits.list.txt -d ${peak_dir}/DNA_cluster/corrected/rep_merged/ -o FC_H3K27ac_clusters_replicates_merged
awk 'NR>1{print $1"\t"$2"\t"$3"\t"$7"\t"$6"\t"$11}' ${peak_dir}/DNA_cluster/corrected/rep_merged/FC_H3K27ac_clusters_replicates_merged.filteredNfixed.union.peakSet > ${peak_dir}/DNA_cluster/corrected/rep_merged/FC_H3K27ac_clusters_replicates_merged.filteredNfixed.union.bed
