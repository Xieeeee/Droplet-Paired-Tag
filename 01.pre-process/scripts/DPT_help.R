
# ## 2022-12-19
# ## Yang Xie (y2xie@ucsd.edu)
# ## Functions for Droplet Paired-Tag analysis in R

ImportArcFRiP <- function(raw_count, frip_count){
	### read in fragment counts to calculate frip ###
    tmp1 <- read.table(paste0(raw_count), header = F, row.names = 1)
    tmp2 <- read.table(paste0(frip_count), header = F, row.names = 1)
    frip <- merge(x = tmp1, y = tmp2, by = 0) %>% setNames(c("bc", "counts_total", "counts_in_peaks"))
    frip$frip <- frip$counts_in_peaks / frip$counts_total
    return(frip)
}

PlotArcFRiP <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){
	### Plot calculated frip ###
    valid <- frip[frip$frip > ycut_low & frip$frip < ycut_high & frip$counts_total > xcut_low, ]
    label <- data.frame(anno = paste0("Total reads cutoff: ", xcut_low, "\n",
                                      "FRiP cutoff: ", ycut_low, "-", ycut_high, "\n", 
                                      "PF cells: ", nrow(valid), "\n")) 
                                      # "Median total reads: ", as.integer(median(valid$counts_total))))
    t1 <- frip %>%
    ggplot(aes(x = counts_total, y = 100*frip)) + 
    geom_point(size = 1, alpha = 0.5, color = "grey") + 
    geom_vline(xintercept = xcut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_high, color = colfunc2(1), linetype="dashed") +
    geom_point(data = valid, aes(x = counts_total, y = 100*frip),
               size = 1, alpha = 0.5, color = colfunc2(3)[[2]]) + 
    geom_text(data = label, aes(x = max(frip$counts_total), y = max(100*as.numeric(frip$frip)), 
                                hjust = 1, vjust = 1, label = anno), size = 2.5) +  
    theme_classic() + xlab("log10(Fragments)") + ylab("% FRiP") + 
    scale_x_log10()
    
    write.table(valid, file = paste0(prefix, "_PF_cells.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    ggsave(t1, filename = paste0(prefix, "_FRiP_macs2.pdf"), height = 6, width = 6, dpi = 300)
    return(valid)
}

PlotArcFragment <- function(frag_path, prefix){
	### Plot CUT&TAG library fragment length pattern ###
    suppressMessages(library(data.table))
    frag <- fread(frag_path)
    frag$len <- frag$V3 - frag$V2
    frags <- frag[sample(nrow(frag), 10000, replace = F), ]
    frags <- frags[, c("len", "V5")]
    frags_pv <- frags[rep(1:nrow(frags), frags[["V5"]]), ]
    t6 <- frags_pv %>% dplyr::filter(len < 1000) %>%
    ggplot(aes(x = len)) + 
    geom_histogram(color = "white", bins = 100) +
    theme_classic() + 
    scale_x_continuous(breaks = seq(0,1000,100)) + 
    xlab("Fragment length")
    ggsave(t6, filename = paste0(prefix, "_fragments.pdf"), height = 6, width = 6, dpi = 300)
    return(t6)
}

PairArc <- function(dmat, rmat, names = "atac"){
    ### in house processing of 10X Multiome. Merge both modalities for plotting ###
    translate <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    # rownames(translate) <- paste0(translate$rna, "-1")
    if (names == "atac"){ ### cells name in dna
        idx <- match(colnames(rmat), translate$rna)
        colnames(rmat) <- translate[idx, names]
        colnames(dmat) <- colnames(dmat)
    }else if (names == "rna"){ ### cells name in rna
        idx <- match(colnames(dmat), translate$atac)
        colnames(dmat) <- translate[idx, names]
        colnames(rmat) <- colnames(rmat)
    }
    dna <- colSums(dmat) %>% as.data.frame() %>% setNames("dna") %>% tibble::rownames_to_column("bc")
    rna <- colSums(rmat) %>% as.data.frame() %>% setNames("rna") %>% tibble::rownames_to_column("bc")
    merge1 <- merge(x = dna, y = rna, by = "bc")
    return(merge1)
}

PlotArcPair <- function(pair, dcutoff = 100, rcutoff = 100, prefix, names = "atac"){
    ### Plot reads count for each pair of barcode. Save plot and segmented cells with name prefix ###
    if(!identical(colnames(pair), c("bc", "dna", "rna"))){
        stop('colnames of input dataframe should be c("bc", "dna", "rna")')
    }
    translate <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    
    pair_pf <- pair[pair$dna > dcutoff & pair$rna > rcutoff, ]
    label <- data.frame(anno = paste0("DNA cutoff: ", min(pair_pf$dna), "\n", "DNA PF cells: ", nrow(pair[pair$dna > dcutoff, ]), "\n", 
                                      "RNA cutoff: ", min(pair_pf$rna), "\n",  "RNA PF cells: ", nrow(pair[pair$rna > rcutoff, ]), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))
    t1 <- pair %>% dplyr::filter(dna > 5 & rna > 5) %>%
    ggplot(aes(x = dna, y = rna)) + 
    geom_point(size = 0.5, color = "grey") + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_point(data = pair_pf, aes(x = dna, y = rna), size = 0.5, color = colfunc2(3)[[2]]) + 
    theme_classic() + 
    geom_text(data = label, aes(x = max(pair$dna), y = min(pair$rna), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("RNA reads") + xlab("DNA reads") + 
    scale_x_log10() + 
    scale_y_log10()
    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    idx <- match(pair_pf$bc, translate[, names])
    valid_bc <- translate[idx, c("atac", "rna")]
    valid_bc <- merge(x = valid_bc, y = pair_pf, by.x = names, by.y = "bc") %>% setNames(c("dna_bc", "rna_bc", "dna_counts", "rna_counts"))
    write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}


OP2 <- function (x) {
    ### dgT to dgC ###
    nms <- colnames(x)
    uniquenms <- unique(nms)
    sparseMatrix(i = x@i + 1, j = match(nms, uniquenms)[x@j + 
        1], x = x@x, dimnames = list(rownames(x), uniquenms), 
        repr = "C")
}

ArcXPM <- function (obj_mtx, meta, method = "CPM", group.by, gname = "gene_id", gene_length = "/projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed") {
    ### calculate cpm/rpkm by cluster (aggreagated cells for calculation) ###
    ### for histone, nCount_histone is used. for RNA, nCount_RNA is used. ###
    colnames(obj_mtx) <- meta[colnames(obj_mtx), group.by]
    obj_mtx <- as(obj_mtx, "dgTMatrix")
    obj_mtx_collapse <- OP2(obj_mtx)
    spars <- length(obj_mtx_collapse@x)/obj_mtx_collapse@Dim[1]/obj_mtx_collapse@Dim[2]
    cat(paste0("sparsity: ", spars, "\n"))
    if (spars > 0.2) {
        obj_mtx_collapse <- as(obj_mtx_collapse, "matrix")
        cat("coarse dgTMatrix into Matrix.\n")
    }
    if (method == "CPM") {
        readSums <- aggregate(meta$nCount_histone, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^6/readSums)
    }
    else if (method == "RPKM") {
        ### mm10: /projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed
        ### hg38: 
        ### gname needs to be either ensembl or gene_id
        length <- read.table(file = gene_length, header = F)
        length <- length %>% setNames(c("chr", "start", "end", "strand", "ensembl", "gene_id"))
        length$length <- length$end - length$start
        len_mtx <- as.data.frame(length[match(rownames(obj_mtx_collapse), length[,gname]), ])
        readSums <- aggregate(meta$nCount_RNA, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^9/readSums)/len_mtx$length
    }
    return(obj_collapse_XPM)
}

### calculate big cor:
### propagate::bigcor
bigcor <- function(x, y, chunks = 2000, method = "pearson"){
    cor <- list()
    nchunks <- ceiling(seq_along(1:ncol(x))/chunks)
    range_chunks <- split(seq_len(ncol(x)), nchunks)
    for(chunk in names(range_chunks)){
        x_qry <- x[,range_chunks[[chunk]]]
        cor[[chunk]] <- cor(x_qry, y, method = method)
    }
    cor <- do.call(rbind, cor)
    return(cor)
}

### from Yang Li
### calculate overlap score after integration
cal_ovlpScore <- function(t1, t2){ 
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- dplyr::full_join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}

