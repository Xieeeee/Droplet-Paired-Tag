#!/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/Rscript

### usage: Rscript 10X_countDup.R [narrowpeak] [outfile] 

args <- commandArgs(trailingOnly = TRUE)
peaks <- args[1] #### 10 columns narrowpeak from macs2, 5th column is log10 q-value used to calculated "score per million"  
prefix <- args[2] ### output name?

tmp <- read.table(peaks, header = F, sep = "\t")
mlogq <- tmp[,5]
tmp$spm <-  10^6 * mlogq / sum(mlogq) ### spm
write.table(tmp, file = prefix, row.names = F, col.names = F, sep = "\t", quote = F) ### export narrowpeak. column 11 is spm score!