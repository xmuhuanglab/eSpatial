
# eSpatial Tutorial for human melanoma dataset


## Library packages and source functions

```{r}
# omics
library(Signac)
library(Seurat)
library(SeuratDisk)
library(SummarizedExperiment)
library(GenomeInfoDb)

# plots
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)

# other
library(dplyr)
library(purrr)
library(Matrix)
library(data.table)
library(future)
library(pbapply)

# eSpatial functions
source('./R/functions.R')
```

## Create Seurat Object
The data downloaded from the Broad Institute Single Cell Portal under the following accession numbers SCP2176. The object created using the Seurat standard pipeline including atac assays, rna assays, cell coordinate information, cell type annotation

```{r}
cre.mat <- readRDS("./Data/cre.mat.rds")
exp.mat <- readRDS("./Data/exp.mat.rds")
meta.data <- readRDS("./Data/meta.data.rds")
cells.coordinates <- readRDS("./Data/cells.coordinates.rds")
chrom.assays = CreateChromatinAssay(counts = cre.mat,
                              assays = "peaks")
obj = CreateSeuratObject(chrom.assays,
                         assay = "peaks",
                         meta.data =meta.data)
obj[["RNA"]] = CreateAssayObject(counts = exp.mat)
obj@reductions[["cod"]] = cells.coordinates
obj
```
## Calculate  specificty scores
For each gene/cis-regulatory element in a specific spatial domain, eSpatial calculated its fold change of the average log-normalized counts in a specific spatial domain versus the rest spatial domains.

```{r}
clusters <- as.character(unique(obj$finalType))
Idents(obj) <- obj$finalType
exp.avg <- SpecificityScore(object = obj,
                            assays = "RNA",
                            clusters = clusters)

peaks.avg <- SpecificityScore(object = obj,
                              assays = "peaks",
                              clusters = clusters)
head(exp.avg)
head(peaks.avg)                            
```

## Define enhancer clusters of genes
eSpatial defines enhancer clusters of genes based on the correlation between enhancer chromatin accessibility and gene expression across spots/spatial neighbor networks.

THIS PROCESS NEED AROUND 1 HOUR!!!

```{r}
exp.mat = obj@assays$RNA@data
cre.mat = obj@assays$peaks@data
GPTab = GPCor(cre.mat = cre.mat,
              exp.mat = exp.mat,
              genome = "hg38")
```

```{r}
GPTabFilt <- subset(GPTab, !is.na(estimate) & estimate > 0 & class == "corr")
# Remove multi-mapping peaks (force 1-1 mapping)
cat("Keeping max correlation for multi-mapping peaks ..\n")
GPTabFilt <- GPTabFilt %>% group_by(Peak) %>% filter(estimate==max(estimate))

```
##  Depict spatial patterns of enhancers for genes expressed in multiple spatial domains
eSpatial first identified gene expression modules using k-means clustering. This clustering grouped genes showing similar expression patterns across spatial domains.

```{r}
exp.mat <- as.matrix(exp.avg)
genes <- intersect(rownames(exp.mat), GPTabFilt$Gene)
exp.mat <- exp.mat[genes,]
clusters = c("tumour_1","tumour_2" ,
             "Tcell_1","Tcell_2")
gene.module = geneModules(exp.mat = exp.mat,
                         k = 6,
                         clusters = clusters)
```
Visualize the gene modules as heatmap

```{r}
matrix <- exp.mat[gene.module$genes, clusters]
plotHeatmap(matrix = matrix,
            row_split = gene.module$kmeans)
```

For each gene module, eSpatial ordered their enhancers based on the specificity scores and visualize.

```{r}
cre.mat = as.matrix(peaks.avg)
colMeans(cre.mat)
cre.mat = scale(cre.mat,center = F)
colMeans(cre.mat)
enh.order = enhancerOrder(exp.mat = exp.mat,
                          cre.mat = cre.mat,
                          genemodule = gene.module,
                          gppairs = GPTabFilt,
                          clusters = clusters,
                          cutoff = 0.02)
matrix <- cre.mat[enh.order$peak, clusters]
plotHeatmap(matrix = matrix,
            row_split = enh.order$genemodule,
            col = c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))            
```

For each gene module, eSpatial defined the spatial patterns of those enhancers regulating the genes within that module.
```{r}
enh.pattern = enhancerPattern(exp.mat = exp.mat,
                              cre.mat = cre.mat,
                              genemodule = gene.module,
                              gppairs = GPTabFilt,
                              clusters = clusters,
                              cutoff.cluster = 0.02,
                              cutoff.bin = 0)


```
##  Decode the spatially divergent combinations of enhancers
eSpatial quantified the diversity of spatial patterns of enhancers within a specific enhancer cluster

THIS PROCESS NEED AROUND 15 MINUTES!
```{r}
enh.combination = enhancerCombination(enh.pattern = enh.pattern,
                                      gppairs = GPTabFilt)
head(enh.combination)
```

