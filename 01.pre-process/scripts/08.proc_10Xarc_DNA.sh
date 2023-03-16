#!/bin/bash
#PBS -q home
#PBS -N cellranger_arc_DNA
#PBS -l nodes=1:ppn=8
#PBS -l walltime=72:00:00

### Author: Yang Xie (y2xie@health.ucsd.edu)

###############################
bcl2="221208_VH00454_113_AAAVYGGHV"
current="/home/y2xie/scratch/28.FC_16k_Droplet_PT_221210"
bcl2_dir="${current}/${bcl2}/"
fastq_dir="${current}/DNA_reseq_221221/"
map_dir="${current}/03.mapping/10X/"
bw_dir="${current}/06.bw/"
macs2_dir="${current}/09.macs2/"
script_dir="${current}/scripts/"
mm10_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"
hg38_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
################################

### demultiplex
# cd ${current}
# cellranger-arc mkfastq --run ${bcl2_dir} --csv ${script_dir}/SampleSheet.atac.csv --output-dir ${fastq_dir}

### run atac processing
cd ${fastq_dir}/
# for f in ${fastq_dir}/${bcl2##*_}/YX*
while read s genome type c #name, genome, peak type, note
do
	if [[ "${genome}" == "mm10" ]]
	then 
		cellranger_ref=${mm10_ref}
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/mm10-blacklist.v2.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vM25.annotation.gtf
		macs2_genome="mm"
	elif [[ "${genome}" == "hg38" ]]
	then
		cellranger_ref=${hg38_ref}
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/hg38-blacklist.v2.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vH35.annotation.gtf
		macs2_genome="hs"
	fi	

	### cellranger	
	# cellranger-atac count --id=${s} --project=${bcl2##*_} --reference=${cellranger_ref} --fastqs=${fastq_dir}/ --sample=${s} --chemistry=ARC-v1

	### rmdup
	python ${script_dir}/scifi.CB_to_BB.py --in ${current}/${s}/outs/possorted_bam.bam
	java -Xmx8G -XX:ParallelGCThreads=16 -jar /projects/ps-renlab/y2xie/packages/picard.jar MarkDuplicates INPUT=${current}/${s}/outs/possorted_bam.bam.BB.bam TMP_DIR=${map_dir} METRICS_FILE=${map_dir}/${s}_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE OUTPUT=${map_dir}/${s}_rmdup.bam BARCODE_TAG=BB REMOVE_DUPLICATES=TRUE

	### bamCoverage
	samtools index ${map_dir}/${s}_rmdup.bam
	bamCoverage -b ${map_dir}/${s}_rmdup.bam -o ${bw_dir}/${s}_rmdup.bw -p max --normalizeUsing RPKM -bl ${bl}
	computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R ${ref_peak} -S ${bw_dir}/${s}_rmdup.bw --skipZeros -o ${bw_dir}/${s}_rmdup.mtx.gz -p max
	plotProfile -m ${bw_dir}/${s}_rmdup.mtx.gz --refPointLabel "TSS" -o ${bw_dir}/${s}_rmdup_profile.pdf --outFileNameData ${bw_dir}/${s}_rmdup_profile.txt
	plotHeatmap -m ${bw_dir}/${s}_rmdup.mtx.gz --refPointLabel "TSS" -o ${bw_dir}/${s}_rmdup_heatmap.pdf --outFileNameMatrix {bw_dir}/${s}_rmdup_heatmap.mtx.gz

	### peak calling
	if [[ "${type}" == "narrow" ]]
	then
		(macs2 callpeak -t ${map_dir}/${s}_rmdup.bam -g ${macs2_genome} -n ${s} -f BAMPE --outdir ${macs2_dir} -q 0.05 --bdg) 2>&1> ${macs2_dir}/log/${s}_rmdup.macs2.log
	elif [[ "${type}" == "broad" ]]
	then
		### estimated extend size
		size=$(/projects/ps-renlab/y2xie/anaconda3/bin/python /projects/ps-renlab/y2xie/scripts/random/getSize.py --bam ${map_dir}/${s}_rmdup.bam)
		size=${size##* } ### median size
		(macs2 callpeak -t ${map_dir}/${s}_rmdup.bam -g ${macs2_genome} -n ${s} -f BAMPE --outdir ${macs2_dir} -q 0.05 --nomodel --extsize ${size} --nolambda --broad-cutoff 0.1 --broad --bdg) 2>&1> ${macs2_dir}/log/${s}_rmdup.macs2.log
	fi

	### FRiP
	zcat ${current}/${s}/outs/fragments.tsv.gz | awk '{if(substr($1,1,1)=="#")next; print}' - | sort -k1,1 -k2,2n - > ${macs2_dir}/${s}_clean_fragments.tsv
	bgzip ${macs2_dir}/${s}_clean_fragments.tsv
	tabix -p bed ${macs2_dir}/${s}_clean_fragments.tsv.gz
	zcat ${macs2_dir}/${s}_clean_fragments.tsv.gz | bedtools intersect -wa -u -a stdin -b ${macs2_dir}/${s}_peaks.${type}Peak > ${macs2_dir}/${s}_FRiP_fragments.tsv
	perl /projects/ps-renlab/y2xie/scripts/scifi/count_10X_fragment.pl ${macs2_dir}/${s}_clean_fragments.tsv.gz
	perl /projects/ps-renlab/y2xie/scripts/scifi/count_10X_fragment.pl ${macs2_dir}/${s}_FRiP_fragments.tsv
	
	rm ${macs2_dir}/${s}_peaks.${type}Peak.tmp ${macs2_dir}/${s}_clean_fragments.tsv.gz
done <${script_dir}/Multiome_DNA_ref2.txt


### FRiP calculation
