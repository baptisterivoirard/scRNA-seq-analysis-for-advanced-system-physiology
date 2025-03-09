rm(list = ls())


# fix the random seed for reproducibility of the results

the.seed <-42
set.seed(the.seed)

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)


#########load the data (in the sparse format) ##########
getwd()

# pré exercice
subject1pre <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject1pre')
subject2pre <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject2pre')
subject3pre <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject3pre')
subject3pre2 <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject3pre2')

# post exercice
subject1post <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject1post')
subject2post <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject2post')
subject3post <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject3post')
subject3post2 <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Free project/Subject3post2')

#### Création objet seurat pré exercice ####
subject1pre <- CreateSeuratObject(counts = subject1pre, min.cells = 5, min.features = 200,)
subject2pre <- CreateSeuratObject(counts = subject2pre, min.cells = 5, min.features = 200,)
subject3pre <- CreateSeuratObject(counts = subject3pre, min.cells = 5, min.features = 200,)
subject3pre2 <- CreateSeuratObject(counts = subject3pre2, min.cells = 5, min.features = 200,)

# Barcodes names 
head(colnames(subject1pre))
# number of barcodes 
ncol(subject1pre)
# features names (genes) 
head(rownames(subject1pre))
# numebr of features (genes) 
nrow(subject1pre)#
# Affichage de l'objet seurat
subject3pre2






###### Anchoring et integration ###

# Ajouter les métadonnées pour identifier les sujets et conditions
subject1pre$subject <- "Subject1"
subject1pre$condition <- "Pre"
subject2pre$subject <- "Subject2"
subject2pre$condition <- "Pre"
subject3pre$subject <- "Subject3"
subject3pre$condition <- "Pre"
subject3pre2$subject <- "Subject3"
subject3pre2$condition <- "Pre"

# Sélectionner les caractéristiques communes pour l'intégration
features <- SelectIntegrationFeatures(object.list = list(subject1pre, subject2pre, subject3pre), nfeatures = 2000)

subject1pre <- NormalizeData(subject1pre)
subject2pre <- NormalizeData(subject2pre)
subject3pre <- NormalizeData(subject3pre)
subject3pre2 <- NormalizeData(subject3pre2)

# Normalisation et réduction en dimensions pour chaque objet
subject1pre <- ScaleData(subject1pre, features = features)
subject1pre <- RunPCA(subject1pre, features = features)

subject2pre <- ScaleData(subject2pre, features = features)
subject2pre <- RunPCA(subject2pre, features = features)

subject3pre <- ScaleData(subject3pre, features = features)
subject3pre <- RunPCA(subject3pre, features = features)

subject3pre2 <- ScaleData(subject3pre2, features = features)
subject3pre2 <- RunPCA(subject3pre2, features = features)


# Trouver les ancres entre les objets
anchors <- FindIntegrationAnchors(object.list = list(subject1pre, subject2pre, subject3pre), anchor.features = features)
# Intégrer les données
integrated_pre <- IntegrateData(anchorset = anchors)

head(colnames(integrated_pre))
head(rownames(integrated_pre))
integrated_pre


## Descriptive analysis 


# Meta.data 

integrated_pre@meta.data$nCount_RNA
integrated_pre$nCount_RNA
integrated_pre$nFeature_RNA

summary(integrated_pre$nCount_RNA)
summary(integrated_pre$nFeature_RNA)


### Basic préprocessing #####


VlnPlot(integrated_pre, features = c("nFeature_RNA"),log = TRUE,alpha = 0.05)+ geom_hline(yintercept = c(270,1000))


VlnPlot(integrated_pre, features = c("nCount_RNA"),log= TRUE ,alpha=0.05)+ geom_hline(yintercept = c(400,1500))

# DEFINE THE CUTOFF VALUES FOR EACH VARIABLE AND VISUALIZE THE RESULTS.
minGene    = 270
minUMI     = 400
maxGene    = 1000
maxUMI     = 1500

# we create a new seurat object containing the filtred cells
pre_filtrd <- subset(integrated_pre, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & nCount_RNA > minUMI & nCount_RNA < maxUMI )

VlnPlot(pre_filtrd, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,log=TRUE, alpha=0.1)


# Vérifie les assays disponibles dans l'objet
Assays(pre_filtrd)

# Vérifie l'assay actif
DefaultAssay(pre_filtrd)
DefaultAssay(pre_filtrd) <- "RNA"
mito_patern <- "^MT-"
pre_filtrd$percent_mt <- PercentageFeatureSet(pre_filtrd, pattern = mito_patern)
summary(pre_filtrd$percent_mt)


qplot(pre_filtrd$percent_mt, geom = "histogram", bins = 100, main = "percent_mt", xlab = "percent_mt") 
maxpct_mt  = 15
# we create a new seurat object containing the filtred cells
pre_filtrd <- subset(pre_filtrd, subset =  percent_mt < maxpct_mt)

rm(anchors,subject3pre2,subject1pre,subject2pre,subject3pre)



DefaultAssay(pre_filtrd) <- "integrated"
simplified <- CreateSeuratObject(
  counts = GetAssayData(pre_filtrd, layer = "data", assay = DefaultAssay(pre_filtrd)),
  meta.data = pre_filtrd@meta.data
)
simplified
Layers(simplified)  # Devrait être vide ou "counts"
sceobj <- as.SingleCellExperiment(simplified)




counts_data <- GetAssayData(pre_filtrd, layer  = "data", assay = "integrated")

# Métadonnées des cellules
cell_metadata <- pre_filtrd@meta.data

# Métadonnées des gènes
feature_metadata <- data.frame(row.names = rownames(pre_filtrd))

# Créer l'objet SingleCellExperiment
library(SingleCellExperiment)

sceobj <- SingleCellExperiment(
  assays = list(logcounts = counts_data),  # Utilise les données intégrées
  colData = cell_metadata,
  rowData = feature_metadata
)

# Prédiction
doublet_res <- scDblFinder(sceobj)


# Sauvegarde des résultats dans sobj
pre_filtrd$doublets.class <- doublet_res$scDblFinder.class

#  % of doublets 
table(pre_filtrd$doublets.class)

pre_filtrd <- subset(pre_filtrd, subset = doublets.class == "singlet")
Layers(pre_filtrd)
Assays(pre_filtrd)

sceobj <- as.SingleCellExperiment(pre_filtrd)

# Prédiction
doublet_res <- scDblFinder(sceobj)


# Sauvegarde des résultats dans sobj
sobj_filtrd$doublets.class <- doublet_res$scDblFinder.class

#  % of doublets 
table(sobj_filtrd$doublets.class)

sobj_filtrd <- subset(sobj_filtrd, subset = doublets.class == "singlet")