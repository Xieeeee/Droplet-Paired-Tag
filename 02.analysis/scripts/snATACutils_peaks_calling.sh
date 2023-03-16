#!/bin/bash
#PBS -q hotel
#PBS -N snATACutils_macs2_${s}
#PBS -l nodes=1:ppn=2
#PBS -l walltime=8:00:00

### call peaks using Yang Li's method
### https://github.com/yal054/snATACutils/tree/master/03.peak_calling

#####
current=/home/y2xie/scratch/28.FC_16k_Droplet_PT_221210/
macs2_dir=${current}/09.macs2/DNA_cluster
bam_dir=${current}/03.mapping/10X/DNA_cluster
mm10_chromsize=/projects/ps-renlab/y2xie/projects/genome_ref/mm10.main.chrom.sizes
hg38_chromsize=/projects/ps-renlab/y2xie/projects/genome_ref/hg38.main.chrom.sizes
#####

# s=`basename $bam .bam`
genome=mm10
if [ $genome == "mm10" ]; then chromsize=${mm10_chromsize}; macs2_genome=mm; fi
if [ $genome == "hg38" ]; then chromsize=${hg38_chromsize}; macs2_genome=hs; fi

mkdir -p ${bam_dir}/bedpe ${macs2_dir}/tagAlign ${macs2_dir}/tnf5 ${macs2_dir}/corrected ${macs2_dir}/spm

### reads pairs to tagAlign
samtools sort -n ${bam_dir}/${s}.bam | bedtools bamtobed -bedpe -i stdin | gzip -c > ${bam_dir}/bedpe/${s}.bedpe.gz

### split into reads
zcat ${bam_dir}/bedpe/${s}.bedpe.gz | awk '{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' - | gzip -c > ${macs2_dir}/tagAlign/${s}.tagAlign.gz

### strand shift
zcat ${macs2_dir}/tagAlign/${s}.tagAlign.gz | awk 'BEGIN{OFS=FS="\t"}{if($1=="."){next}if($6=="+"){$2=$2+4}else if($6=="-"){$3=$3-5} print $0}' - | gzip -nc > ${macs2_dir}/tnf5/${s}.tnf5_tag.gz

### peak calling
(macs2 callpeak -t ${macs2_dir}/tnf5/${s}.tnf5_tag.gz -f BED -q 0.01 -g ${macs2_genome} -n ${s} --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR --call-summits --outdir ${macs2_dir}/corrected) 2>&1> ${macs2_dir}/corrected/${s}.macs2.log

# grep 'chr' ${macs2_dir}/corrected/${s}_treat_pileup.bdg | egrep -v "chrM|chrUn|random" | sort -k1,1 -k2,2n > ${macs2_dir}/corrected/${s}_treat_pileup.srt.bdg
# bedGraphToBigWig ${macs2_dir}/corrected/${s}_treat_pileup.srt.bdg ${chromsize} ${macs2_dir}/corrected/${s}_treat_pileup.srt.bw
# rm ${macs2_dir}/corrected/${s}_treat_pileup.srt.bdg ${macs2_dir}/corrected/${s}_treat_pileup.bdg

### SPM corrected (narrowPeak only)
/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/Rscript /projects/ps-renlab/y2xie/scripts/scifi/snATACutils_spm_corrected.R ${macs2_dir}/corrected/${s}_peaks.narrowPeak ${macs2_dir}/spm/${s}_peaks.spm.narrowPeak
