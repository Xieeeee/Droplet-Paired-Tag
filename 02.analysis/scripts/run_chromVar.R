#!/usr/bin/env Rscript
### Author: Yang Xie (y2xie@health.ucsd.edu)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--input", required=TRUE, help="cell-by-peak seurat object")
# parser$add_argument("-g", "--genome", default = "mm10", help="genome?")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
args <- parser$parse_args()
# args <- commandArgs(trailingOnly = TRUE)

inF = args$input
output = args$output
# genoF = args$genome

####################################
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(JASPAR2020))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressPackageStartupMessages(library(patchwork))
set.seed(921)

obj <- readRDS(inF)
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
obj <- RunChromVAR(obj, genome = BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(obj) <- 'chromvar'
Idents(obj) <- "rna_Anno0.8"
valid <- rownames(obj@meta.data[obj$DNA_keep_cells == "True" & obj$RNA_keep_cells != "False", ]) 
obj_valid <- subset(obj, cells = valid)
objMotif <- FindAllMarkers(obj_valid, only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff")

objMotif <- objMotif %>% dplyr::filter(p_val_adj < 0.05)
write.table(objMotif, file = paste0(output, "_enriched_motif.xls"), 
    row.names = F, col.names = T, sep = "\t", quote = F)
saveRDS(obj, inF)

### if error threw out during chromVAR needs to run from scratch
# library(chromVAR)
# library(motifmatchr)
# library(SummarizedExperiment)
# library(BiocParallel)
# register(SerialParam())

# obj <- readRDS(inF)
# mtx <- obj[["DNA"]]@counts[-which(rowSums(obj[["DNA"]]@counts) == 0),]
# peak <- rownames(mtx) %>% stringr::str_split_fixed(pattern = "-", n = 3) %>% as.data.frame() %>% setNames(c("chr", "start", "end"))
# write.table(peak, file = "FC_H3K27me3_signac_1kbin.bed", row.names = F, col.names = F, sep = "\t", quote = F)
# peaks <- getPeaks("FC_H3K27me3_signac_1kbin.bed", sort_peaks = TRUE)
# fragment_counts <- SummarizedExperiment(assays = list(mtx), rowRanges = peaks)
# fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Mmusculus.UCSC.mm10)
# ### or will mess up the getBackGround
# row.data <- data.frame(rowData(fragment_counts))
# row.data[is.na(row.data)] <- 0
# rowData(fragment_counts) <- row.data
# motif_ix <- matchMotifs(pfm, fragment_counts, genome=BSgenome.Mmusculus.UCSC.mm10)
# assayNames(fragment_counts) <- "counts"
# dev <- computeDeviations(object = fragment_counts, annotations = motif_ix) ### !!!

# saveRDS(dev, paste0(output, "_chromvar.rds"))

