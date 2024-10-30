library(scater)
library(SingleCellExperiment)
library(SC3)
library(tidyr)
library(ggplot2)

data <- read.table("DGE_count_clean.txt", header=T, row.names=1)
anno <- read.table("anno.txt", header=T)
anno$TIME[anno$TIME=='t0']='0h'
anno$TIME[anno$TIME=='t1']='1h'
anno$TIME[anno$TIME=='t4']='4h'
anno$TIME[anno$TIME=='t6']='6h'
data <- data[,as.character(anno$Sample)]
data2 <- data[which(rowSums(data)>0), ]

e=data2
e.tmp=e
e.tmp[e.tmp>0]=1
data2=e[which(rowSums(e.tmp)==dim(e.tmp)[2]),]

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


sce <- SingleCellExperiment(assays=list(counts=as.matrix(data2),logcounts=log2(as.matrix(data2)+1)),colData=anno)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- sc3_estimate_k(sce)
sce@metadata$sc3$k_estimation
sce <- sc3(sce, ks=4:5, biology=TRUE, gene_filter=T, n_cores=10, 
           d_region_min=0.06, d_region_max=0.25, 
           pct_dropout_min=20, pct_dropout_max=80)
pdf('heatmap.pdf',width=8,height=5)
sc3_plot_expression2(sce, k=4, show_pdata=c("TIME"))
dev.off()

sce2 <- sce[sce@rowRanges@elementMetadata@listData[["sc3_gene_filter"]]=='TRUE', ]
sce2 <- runPCA(sce2)
plotPCA(sce2, colour_by = "TIME")
sce2 <- runTSNE(sce2)
pdf('tSNE.pdf',width=3.5,height=3)
plotTSNE(sce2, colour_by = "TIME")
dev.off()
sce2 <- runUMAP(sce2)
pdf('UMAP.pdf',width=3.5,height=3)
plotUMAP(sce2, colour_by = "TIME")
dev.off()

tt=sce2@int_colData$reducedDims$TSNE
rownames(tt)=colnames(data3)
write.table(tt, 'tSNE.txt', sep='\t')
uu=sce2@int_colData$reducedDims$UMAP
rownames(uu)=colnames(data3)
write.table(uu, 'UMAP.txt', sep='\t')

