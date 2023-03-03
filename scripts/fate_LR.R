###############################################
# Load packages
options(warn=-1)
library(Seurat)
library(SeuratDisk)

###############################################
# params
input_probs = 'data/Fate_prob_for_LR.csv'
input_metadata = 'data/metadata_for_LR.csv'
outdir = "data/"
test_condition = "Jak2_Homo"
ctrl_condition = "Jak2_WT"
outname_prefix = "PMCA_Jak2"

###############################################
# Load data
probs = read.csv(input_probs,row.names=1)
df = as.data.frame(t(probs))
cell_meta = read.csv(input_metadata,row.names=1)

###############################################
# Create object
obj = CreateSeuratObject(counts = df, meta.data = cell_meta)

Idents(obj) = 'celltype'
for (ct in setdiff(unique(obj$celltype),"Unassigned")){
  obj.tmp = subset(obj, idents = ct)
  
  test.median = apply(GetAssayData(obj.tmp[,obj.tmp$Condition==test_condition]),1,median)
  ctrl.median = apply(GetAssayData(obj.tmp[,obj.tmp$Condition==ctrl_condition]),1,median)
  
  median.diff = test.median - ctrl.median
  median.diff = as.data.frame(median.diff)
  
  Idents(obj.tmp)='Condition'
  
  # LR test
  if (length(unique(obj$dataset))==1){
    markers = FindMarkers(obj.tmp, ident.1=test_condition, ident.2=ctrl_condition,
                          test.use = "LR",
                          min.cells.group = 1,
                          min.cells.feature = 1,
                          min.pct = 0,
                          logfc.threshold = 0)
  } else {
    markers = FindMarkers(obj.tmp, ident.1=test_condition, ident.2=ctrl_condition,
                          test.use = "LR",
                          min.cells.group = 1,
                          min.cells.feature = 1,
                          min.pct = 0,
                          logfc.threshold = 0,
                          latent.vars = "dataset")
  }
  
  markers$median.diff = median.diff[rownames(markers),]
  markers$avg_log2FC = NULL
  markers$pct.1 = NULL
  markers$pct.2 = NULL
  
  write.table(markers,paste0(outdir,outname_prefix,"_",ct,"_fate_",test_condition,"_vs_",ctrl_condition,".csv"),sep=",",quote=F,row.names=T,col.names=NA)
  
}

cat("Analysis completed\n")
