#!/bin/bash
#PBS -q home
#PBS -N cellranger_RNA
#PBS -l nodes=1:ppn=4
#PBS -l walltime=72:00:00

### Author: Yang Xie (y2xie@health.ucsd.edu)

###############################
bcl2="221212_VH00454_115_AAAV7J7HV"
current="/home/y2xie/scratch/28.FC_16k_Droplet_PT_221210"
bcl2_dir="${current}/${bcl2}/"
fastq_dir="${current}/RNA_reseq_221219/"
trim_dir="${current}/02.trimmed/"
map_dir="${current}/03.mapping/"
mtx_dir="${current}/04.matrices/"
script_dir="${current}/scripts/"
cellranger_mm10="/projects/ps-renlab/y2xie/projects/genome_ref/mm10"
cellranger_hg38="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
################################

# cd ${current} 
# cellranger-arc mkfastq --run ${bcl2_dir} --csv ${script_dir}/SampleSheet.rna.csv --output-dir ${fastq_dir}

cd ${fastq_dir}

while read s genome c
do
	if [[ "${genome}" == "mm10" ]]
	then 
		cellranger_ref=${cellranger_mm10}
	elif [[ "${genome}" == "hg38" ]]
	then
		cellranger_ref=${cellranger_hg38}
	fi			
	# for f in ${fastq_dir}/${bcl2##*_}/*
	# do 
	# 	fname=`basename $f`
	# 	s=${fname%%_*}
	cellranger count --id=${s} --project=${bcl2##*_} --transcriptome=${cellranger_ref} --fastqs=${fastq_dir}/ --sample=${s} --include-introns --chemistry=ARC-v1
done <${script_dir}/Multiome_RNA_ref2.txt
