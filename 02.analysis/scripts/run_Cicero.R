#!/usr/bin/env Rscript
### Author: Yang Xie (y2xie@health.ucsd.edu)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--input", required=TRUE, help="cell-by-peak matrix rds")
parser$add_argument("-c", "--cell", default = "None", help="PF cells for filtering")
parser$add_argument("-m", "--meta", required=TRUE, help="metadata for input matrix. Should contain umap_1 and umap_2.")
parser$add_argument("-g", "--genome", default="/projects/ps-renlab/y2xie/projects/genome_ref/mm10.main.chrom.sizes", help="chromInfo for specified genome")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
args <- parser$parse_args()
# args <- commandArgs(trailingOnly = TRUE)

inF = args$input
pfF = args$cell
meta = args$meta
output = args$output
genoF = args$genome

####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(cicero))

dna.meta <- read.table(meta, header = T, row.names = 1, sep = "\t")
indata <- readRDS(inF)
if(pfF != "None"){
  PF_cells <- read.table(pfF, header = F, sep = "\t")
  ### filtering cells!
  indata <- indata[, intersect(colnames(indata), PF_cells$V1)]
  dna.meta <- dna.meta[intersect(rownames(dna.meta), PF_cells$V1), ]
  ### doubel secure
  if(ncol(indata) != nrow(dna.meta)){
    real_valid <- intersect(rownames(dna.meta), colnames(indata))
    indata <- indata[, real_valid]
    dna.meta <- dna.meta[real_valid, ]
  }
}

# indata@x[indata@x > 0] <- 1
cellinfo <- colnames(indata) %>% as.data.frame() %>% setNames("V1")
rownames(cellinfo) <- cellinfo$V1

peakinfo <- rownames(indata) %>% as.data.frame() %>% setNames("cood")
rownames(peakinfo) <- peakinfo$cood
peakinfo[, c("chr", "bp1", "bp2")] <- stringr::str_split_fixed(peakinfo$cood, pattern = "[:-]", n = 3)
peakinfo <- peakinfo %>% select(c("chr", "bp1", "bp2"))

### to cds
input_cds <- suppressWarnings(new_cell_data_set(indata, cell_metadata = cellinfo, gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

### consider providing umap
if(("umap_1" %in% colnames(dna.meta))&("umap_2" %in% colnames(dna.meta))){
  print("extract umap coordinates from metadata...")
  umap_coords <- dna.meta[,c("umap_1", "umap_2")]
  }else{
  print("run umap on cds...")
  set.seed(921)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
  umap_coords <- reducedDims(input_cds)$UMAP
} 
umap_coords <- as.matrix(umap_coords)
### cicero
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
saveRDS(cicero_cds, paste0(output, "_cicero_cds.rds"))

genome <- read.csv(genoF, sep = "\t", header = F)
conns <- run_cicero(cicero_cds, genome)#, window = 1e+06, distance_constraint = 5e+05) ### default window 500k, distance 250k
write.table(conns, file = paste0(output, "_cicero_conns_unfilt.xls"), row.names = F, col.names = T)

### filter coaccess
CCAN_assigns <- generate_ccans(conns)
write.table(CCAN_assigns, file = paste0(output, "_cicero_conns_CCAN.xls"), row.names = F, col.names = T)


