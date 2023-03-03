###############################################
# Load packages
options(warn=-1)
library(Seurat)
library(SeuratDisk)

###############################################
# params
input_flux = 'data/Flux.csv'
input_metadata = 'data/metadata_for_LR.csv'
outdir = "data/"
test_condition = "Jak2_Homo"
ctrl_condition = "Jak2_WT"
outname_prefix = "PMCA_Jak2"

###############################################
# Load data
flux = read.csv(input_flux,row.names=1)
df = as.data.frame(t(flux))
width = apply(df,1,function(x){max(x)-min(x)})
df.filt = df[which(width>=1e-4),]
cell_meta = read.csv(input_metadata,row.names=1)

###############################################
# Create object
obj = CreateSeuratObject(counts = df.filt, meta.data = cell_meta)

Idents(obj) = 'celltype'
for (ct in setdiff(unique(obj$celltype),"Unassigned")){
  obj.tmp = subset(obj, idents = ct)
  
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
  
  markers$avg_log2FC = NULL
  markers$pct.1 = NULL
  markers$pct.2 = NULL
  
  write.table(markers,paste0(outdir,outname_prefix,"_",ct,"_metab_",test_condition,"_vs_",ctrl_condition,".csv"),sep=",",quote=F,row.names=T,col.names=NA)
  
}

cat("Analysis completed\n")
