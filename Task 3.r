#########################################
### Single-cell analysis using Seurat ###
#########################################
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
# Load the single-cell dataset
csf.data = read.table(
  'counts_TPM_ALL.csv',
  header = T,
  row.names = 1,
  sep = '\t'
)
csf.data = log(csf.data + 1)
celltypes = unlist(lapply(colnames(csf.data), ExtractField, 1))
table(celltypes)
# Create Seurat object. Keep all genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
csf <-
  CreateSeuratObject(
    raw.data = csf.data,
    min.cells = 3,
    min.genes = 200,
    project = 'CSF_2018'
  )
# AddMetaData adds columns to object@meta.data:
sc_cell_info <-
  read.table('sc_cell_info.txt', header = T, row.names = 1)
csf <-
  AddMetaData(
    object = csf,
    metadata = sc_cell_info,
    col.name = c('Twin', 'Case', 'Sample',
                 'index.sort', 'Clones')
  )
# QC and selecting cells for further analysis. We calculate the percentage of mitochondrial genes here
# and store it in percent.mito using AddMetaData. We use object@raw.data since this represents nontransformed and non -
#  log - normalized counts. The % of counts mapping to MT - genes is a common
# scRNAseq QC metric.
Sys.setlocale('LC_ALL', 'C')
mito.genes <-
  grep(pattern = "^MT-",
       x = rownames(x = csf@data),
       value = TRUE)
percent.mito <-
  Matrix::colSums(csf@raw.data[mito.genes,]) / Matrix::colSums(csf@raw.data)
csf <-
  AddMetaData(object = csf,
              metadata = percent.mito,
              col.name = 'percent.mito')
VlnPlot(
  object = csf,
  features.plot = c('nGene', 'percent.mito'),
  nCol = 2,
  x.lab.rot = TRUE
)
csf <-
  FilterCells(
    object = csf,
    subset.names = c('nGene', 'percent.mito'),
    low.thresholds = c(200,-Inf),
    high.thresholds = c(6000, 0.05)
  )
# Normalizing the data
## Further normalization is performed since the dataset used to create the Seurat object is a merge of
# normalized datasets into one file. That means that, in the final merged file, the data from different
# dataset are not normalized to the same sequencing depth. Therefore, a further normalization is
# required.
csf <-
  NormalizeData(
    object = csf,
    normalization.method = 'LogNormalize',
    scale.factor = 10000
  )
# Detection of variable genes across the single cells
csf <-
  FindVariableGenes(
    object = csf,
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    x.low.cutoff = 0.0125,
    x.high.cutoff = 3,
    y.cutoff = 0.5
  )
length(x = csf@var.genes)
# Scaling the data and removing unwanted sources of variation
csf <- ScaleData(object = csf, vars.to.regress = 'percent.mito')
# Perform linear dimensional reduction
csf <-
  RunPCA(
    object = csf,
    pc.genes = csf@var.genes,
    do.print = TRUE,
    pcs.print = 1:5,
    genes.print = 5
  )
VizPCA(object = csf, pcs.use = 1:2)
PCAPlot(object = csf,
        dim.1 = 1,
        dim.2 = 2)
# ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their
# correlation with the calculated components
csf <- ProjectPCA(object = csf, do.print = FALSE)
# Heatmaps based on the PCA
PCHeatmap(
  object = csf,
  pc.use = 1,
  cells.use = 500,
  do.balanced = TRUE,
  label.columns = FALSE
)
PCHeatmap(
  object = csf,
  pc.use = 1:20,
  cells.use = 500,
  do.balanced = TRUE,
  label.columns = FALSE,
  use.full = FALSE
)
PrintPCA(
  object = csf,
  pcs.print = 1:20,
  genes.print = 5,
  use.full = FALSE
)
# Determine statistically significant principal components
csf <- JackStraw(object = csf, num.replicate = 100)
JackStrawPlot(object = csf, PCs = 1:20)
PCElbowPlot(object = csf)
# Cluster the cells
csf <-
  FindClusters(
    object = csf,
    reduction.type = 'pca',
    dims.use = 1:10,
    resolution = 0.6,
    print.output = 0,
    save.SNN = TRUE
  )
PrintFindClustersParams(object = csf)
# Run Non-linear dimensional reduction (tSNE)
csf <-
  RunTSNE(
    object = csf,
    dims.use = 1:10,
    do.fast = TRUE,
    check_duplicates = FALSE
  )
TSNEPlot(object = csf, do.label = T)
# QC of each cluster
VlnPlot(
  object = csf,
  features.plot = c('nGene', 'percent.mito'),
  nCol = 2,
  x.lab.rot = TRUE
)
## Cells in Cluster #3 that were characterized by low number of detected genes and low mitochondrial
# gene transcripts (indicating low - quality cells) were removed
new.subset <- SubsetData(object = csf, ident.remove = "3")
csf_v1 <- csf
save(csf_v1, file = "./csf_v1_2018.Rda")
csf <- new.subset
csf <-
  FilterCells(
    object = csf,
    subset.names = c('nGene', 'percent.mito'),
    low.thresholds = c(200,-Inf),
    high.thresholds = c(6000, 0.025)
  )
## Additional round of clustering after SubsetData
# re-running FindVariableGenes() and ScaleData()
csf <- FindVariableGenes(
  object = csf,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.0125,
  x.high.cutoff = 3,
  y.cutoff = 0.5
)
length(x = csf@var.genes)
csf <- ScaleData(object = csf, vars.to.regress = 'percent.mito')
csf <-
  RunPCA(
    object = csf,
    pc.genes = csf@var.genes,
    do.print = TRUE,
    pcs.print = 1:5,
    genes.print = 5
  )
VizPCA(object = csf, pcs.use = 1:2)
PCAPlot(object = csf,
        dim.1 = 1,
        dim.2 = 2)
csf <- ProjectPCA(object = csf, do.print = FALSE)
PCHeatmap(
  object = csf,
  pc.use = 1,
  cells.use = 500,
  do.balanced = TRUE,
  label.columns = FALSE
)
PCHeatmap(
  object = csf,
  pc.use = 1:20,
  cells.use = 500,
  do.balanced = TRUE,
  label.columns = FALSE,
  use.full = FALSE
)
PrintPCA(
  object = csf,
  pcs.print = 1:20,
  genes.print = 5,
  use.full = FALSE
)
csf <- JackStraw(object = csf, num.replicate = 100)
JackStrawPlot(object = csf, PCs = 1:20)
PCElbowPlot(object = csf)
csf <-
  FindClusters(
    object = csf,
    reduction.type = 'pca',
    dims.use = 1:11,
    resolution = 0.6,
    print.output = 0,
    save.SNN = TRUE,
    force.recalc = TRUE
  )
PrintFindClustersParams(object = csf)
csf <-
  RunTSNE(
    object = csf,
    dims.use = 1:11,
    do.fast = TRUE,
    check_duplicates = FALSE
  )
TSNEPlot(object = csf, do.label = T)
# QC of each cluster
VlnPlot(
  object = csf,
  features.plot = c('nGene', 'percent.mito'),
  nCol = 2,
  x.lab.rot = TRUE
)
## Cells in Cluster #6 that were characterized primarily by mitochondrial and ribosomal gene transcripts
# (indicating low - quality cells) were removed
new.subset <- SubsetData(object = csf, ident.remove = "6")
csf_v2 <- csf
csf_v2
save(csf_v2, file = "./csf_v2_2018.Rda")
csf <- new.subset
## Additional round of clustering after SubsetData
# re-running FindVariableGenes() and ScaleData()
csf <- FindVariableGenes(
  object = csf,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.0125,
  x.high.cutoff = 3,
  y.cutoff = 0.5
)
length(x = csf@var.genes)
csf <- ScaleData(object = csf, vars.to.regress = 'percent.mito')
csf <-
  RunPCA(
    object = csf,
    pc.genes = csf@var.genes,
    do.print = TRUE,
    pcs.print = 1:5,
    genes.print = 5
  )
VizPCA(object = csf, pcs.use = 1:2)
PCAPlot(object = csf,
        dim.1 = 1,
        dim.2 = 2)
csf <- ProjectPCA(object = csf, do.print = FALSE)
PCHeatmap(
  object = csf,
  pc.use = 1,
  cells.use = 500,
  do.balanced = TRUE,
  label.columns = FALSE
)
PCHeatmap(
  object = csf,
  pc.use = 1:20,
  cells.use = 500,
  do.balanced = TRUE,
  label.columns = FALSE,
  use.full = FALSE
)
PrintPCA(
  object = csf,
  pcs.print = 1:20,
  genes.print = 5,
  use.full = FALSE
)
csf <- JackStraw(object = csf, num.replicate = 100)
JackStrawPlot(object = csf, PCs = 1:20)
PCElbowPlot(object = csf)
csf <-
  FindClusters(
    object = csf,
    reduction.type = 'pca',
    dims.use = 1:11,
    resolution = 0.6,
    print.output = 0,
    save.SNN = TRUE,
    force.recalc = TRUE
  )
PrintFindClustersParams(object = csf)
csf <-
  RunTSNE(
    object = csf,
    dims.use = 1:11,
    do.fast = TRUE,
    check_duplicates = FALSE
  )
TSNEPlot(object = csf, do.label = T)
csf.markers <-
  FindAllMarkers(
    object = csf,
    only.pos = TRUE,
    min.pct = 0.25,
    thresh.use = 0.25
  )
csf.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(csf.markers %>% group_by(cluster) %>% top_n(10, avg_logFC),
            'csf.markers.tsv',
            sep = '\t')
top10 <- csf.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(
  object = csf,
  genes.use = top10$gene,
  slim.col.label = TRUE,
  remove.key = TRUE
)
## Heatmap shows that cells in clusters #0 and #3 share same top markers, and therefore were clustered
# together. And the same holds for clusters #1 and #2. Tiny cluster between clusters #5 and #9 which
# contains platelet - like cells coming from the blood samples, was removed for the figure in the publication.
TSNEPlot(
  object = csf,
  do.return = TRUE,
  group.by = 'index.sort',
  do.label = TRUE
)
TSNEPlot(
  object = csf,
  do.return = TRUE,
  group.by = 'Clones',
  do.label = TRUE
)
VlnPlot(
  object = csf,
  features.plot = c('nGene', 'percent.mito'),
  nCol = 2,
  x.lab.rot = TRUE
)
write.table(csf@meta.data, 'meta.data.tsv', sep = '\t')
save(csf, file = "./csf_2018.Rda")
###################################
sessionInfo()
# R version 3.4.4 (2018-03-15)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.1 LTS
# Matrix products: default
# BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# locale:
# [1] LC_CTYPE=C LC_NUMERIC=C
# [3] LC_TIME=C LC_COLLATE=C
# [5] LC_MONETARY=C LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=de_DE.UTF-8 LC_NAME=C
# [9] LC_ADDRESS=C LC_TELEPHONE=C
# [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C
# attached base packages:
# [1] stats graphics grDevices utils datasets methods # base
# other attached packages:
# [1] bindrcpp_0.2.2 dplyr_0.7.8 Seurat_2.3.4 Matrix_1.2-15 # cowplot_0.9.3
# [6] ggplot2_3.1.0 usethis_1.4.0 devtools_2.0.1