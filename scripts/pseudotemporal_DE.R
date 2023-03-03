###############################################
# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(tradeSeq)

###############################################
# params
input_adata = 'data/PMCA_Jak2_Ery_raw_count.h5ad'
outdir = "data/"
outname_prefix = "PMCA_Jak2"

###############################################
# Load data
sce <- readH5AD(input_adata)

###############################################
# tradeSeq args
cw <- rep(1,ncol(sce)) # cell weights
pseudotime <- colData(sce)$dpt_pseudotime # pseudotime
batch <- as.factor(colData(sce)$dataset)

###############################################
# Fit GAM using tradeSeq
set.seed(5)
if (length(unique(batch))==1){
  sceGAM <- fitGAM(counts = as.matrix(assays(sce)$X),
                   conditions = factor(colData(sce)$Condition),
                   pseudotime=pseudotime,
                   cellWeights=cw,
                   nknots=6,
                   verbose=T)
} else {
  # create a model matrix
  U <- model.matrix(~batch)
  
  sceGAM <- fitGAM(counts = as.matrix(assays(sce)$X),
                   conditions = factor(colData(sce)$Condition),
                   U = U,
                   pseudotime=pseudotime,
                   cellWeights=cw,
                   nknots=6,
                   verbose=T)
}

saveRDS(sceGAM,file = paste0(outdir,outname_prefix,"_sceGAM.rds"))

###############################################
# Differential expression between conditions
condRes <- conditionTest(sceGAM, l2fc = 0.5)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")

condRes = condRes[order(condRes$waldStat,decreasing=T),]
condRes = condRes[order(condRes$padj,decreasing=F),]

write.table(condRes,
            paste0(outdir,outname_prefix,"_pseudotemporal_DE.txt"),
            sep='\t',quote=F,row.names=T,col.names=NA)

###############################################
# plot smoothers for top 5 genes
gene.ls = rownames(condRes)[1:5]

for (gene_to_plot in gene.ls) {
  pdf(paste0(outdir,outname_prefix,"_",gene_to_plot,".pdf"))
  print(plotSmoothers(sceGAM, as.matrix(assays(sce)$X), gene_to_plot) + ggplot2::ggtitle(gene_to_plot))
  dev.off()
}

cat("Analysis completed\n")
