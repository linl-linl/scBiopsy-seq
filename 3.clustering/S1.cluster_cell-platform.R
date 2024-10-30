library(scater)
library(SingleCellExperiment)
library(SC3)
library(tidyr)
library(ggplot2)

data <- read.table("DGE_count_clean.txt", header=T, row.names=1)
data2 <- read.table("anno.txt", header=T)
rownames(data2) <- data2$Sample
data3 <- data[which(rowSums(data)>0), rownames(data2)]

sc3_plot_expression2 <- function(object, k, show_pdata) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset(object)
    if (!is.null(metadata(object)$sc3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- data.frame(colData(object)[show_pdata])
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            # make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    if(nrow(dataset) > 100) {
        heatmap <- do.call( pheatmap::pheatmap, c(list(dataset, cluster_cols=hc, kmeans_k=100, cutree_cols=k, show_rownames=FALSE, show_colnames=TRUE, angle_col=c("45"), fontsize_col=4), list(annotation_col=ann)[add_ann_col]) )
    } else {
        heatmap <- do.call( pheatmap::pheatmap, c(list(dataset, cluster_cols=hc, cutree_cols=k, show_rownames=FALSE, show_colnames=TRUE, angle_col=c("45"), fontsize_col=4), list(annotation_col=ann)[add_ann_col]) )
    }
    return(heatmap)
}

sce <- SingleCellExperiment(assays=list(counts=as.matrix(data3),logcounts=log2(as.matrix(data3)+1)),colData=data2)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- sc3_estimate_k(sce)
sce@metadata$sc3$k_estimation

sce <- sc3(sce, ks=2:5, biology=TRUE, gene_filter=T, n_cores=10, 
           d_region_min=0.06, d_region_max=0.25, 
           pct_dropout_min=20, pct_dropout_max=80)
pdf('heatmap.pdf',width=5,height=4)
sc3_plot_expression2(sce, k=4, show_pdata=c("Cell", "Platform", "Time"))
dev.off()

sce2 <- sce[sce@rowRanges@elementMetadata@listData[["sc3_gene_filter"]]=='TRUE', ]
sce2 <- runPCA(sce2)
plotPCA(sce2, colour_by = "Cell")
sce2 <- runTSNE(sce2)
pdf('tSNE.pdf',width=3.8,height=3)
plotTSNE(sce2, colour_by = "Cell")
dev.off()
sce2 <- runUMAP(sce2)
pdf('UMAP.pdf',width=3.8,height=3)
plotUMAP(sce2, colour_by = "Cell")
dev.off()

tt=sce2@int_colData$reducedDims$TSNE
rownames(tt)=colnames(data3)
write.table(tt, 'tSNE.txt', sep='\t')
uu=sce2@int_colData$reducedDims$UMAP
rownames(uu)=colnames(data3)
write.table(uu, 'UMAP.txt', sep='\t')
