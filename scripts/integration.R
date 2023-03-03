###########################
# Load packages
library(Seurat)
library(SeuratDisk)
library(tidyverse)

###########################
# params
vargenes_file = "scripts/integration_genes.txt"
ref_file = "data/Dahlin_atlas.h5ad"
query_file = "data/PMCA_Jak2_SLX11516_SIGAC2.h5ad"
outdir = "data/"
outname_prefix = "PMCA_Jak2_SLX11516_SIGAC2"

###########################
# Use atlas HVGs for integration
vargenes = read.delim(vargenes_file,header = F)
vargenes = as.character(vargenes$V1)

###########################
# Convert h5ad to h5seurat
ref_file_h5seurat = str_replace(ref_file,".h5ad",".h5seurat")
query_file_h5seurat = str_replace(query_file,".h5ad",".h5seurat")

if (!file.exists(ref_file_h5seurat)) {
  Convert(ref_file, dest = "h5seurat")
}

if (!file.exists(query_file_h5seurat)) {
  Convert(query_file, dest = "h5seurat")
}

###########################
# Load h5seurat files to integrate
ref.obj = LoadH5Seurat(ref_file_h5seurat)
sample.obj = LoadH5Seurat(query_file_h5seurat)

###########################
# pre-calc CCA
ref.obj.scale <- ScaleData(object = ref.obj, features = vargenes, do.scale = T, verbose = FALSE)
sample.obj.scale <- ScaleData(object = sample.obj, features = vargenes, do.scale = T, verbose = FALSE)

combined.ob <- RunCCA(
  object1 = ref.obj.scale,
  object2 = sample.obj.scale,
  features = vargenes,
  num.cc = 30,
  renormalize = FALSE,
  rescale = FALSE,
  verbose = T)

###########################
# store CCA results in the original obj
ref.obj[["cca"]] <- CreateDimReducObject(embeddings = Embeddings(combined.ob[["cca"]])[colnames(ref.obj),],
                                         key = "CC_", assay = DefaultAssay(ref.obj))

sample.obj[["cca"]] <- CreateDimReducObject(embeddings = Embeddings(combined.ob[["cca"]])[colnames(sample.obj),],
                                            key = "CC_", assay = DefaultAssay(sample.obj))

###########################
# Split datasets
ref.list = list(REF=ref.obj)
sample.list = list(QUERY=sample.obj)

###########################
# Make lists
comb.list = c(ref.list,sample.list)
reference = which(names(comb.list) == "REF")

###########################
# Check lists
comb.list

###########################
# Find anchors
anchors = FindIntegrationAnchors(object.list = comb.list, dims = 1:30,
                                 reference = reference, anchor.features = vargenes,
                                 k.anchor = 5, k.filter = 200, k.score = 30)

###########################
# Integrate
combined = IntegrateData(anchorset = anchors,dims = 1:30,
                         weight.reduction = lapply(comb.list,function(x){return(x[['cca']])}),
                         features.to.integrate = row.names(ref.obj))

###########################
# Save integration results
integrated.sample = subset(x = combined, subset = library != "REF")
SaveH5Seurat(integrated.sample, filename = paste0(outdir, outname_prefix, "_Integrated.h5seurat"))
Convert(paste0(outdir, outname_prefix, "_Integrated.h5seurat"), dest = "h5ad")

###########################
# Proceed to label transfer steps
# Find transfer anchors
transfer.anchors = FindTransferAnchors(reference = ref.obj, query = sample.obj,
                                       reduction = "cca",features = vargenes,
                                       dims = 1:30, k.anchor = 5, k.filter = 200, k.score = 30)

###########################
# Transfer labels
predictions <- TransferData(anchorset = transfer.anchors, weight.reduction = "cca",
                            refdata = ref.obj$celltype, dims = 1:30)

###########################
# Save results
write.table(predictions, paste0(outdir, outname_prefix, "_label_transfer.csv"),
            sep=",",quote=F,row.names = T,col.names = NA)

###########################
# Deleting intermediate files
file.remove(ref_file_h5seurat)
file.remove(query_file_h5seurat)
file.remove(paste0(outdir, outname_prefix, "_Integrated.h5seurat"))

cat("Analysis completed\n")
