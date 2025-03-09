rm(list = ls())


# fix the random seed for reproducibility of the results

the.seed <-42
set.seed(the.seed)

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)


#load the data (in the sparse format)
getwd()
scmat <- Read10X('C:/Users/rivoi/Documents/M2MAD/AdSP/Guidedproject/rep2')

sample.name <- "rep2"

# sparse matrix
class(scmat)
# extract some data
scmat[1:5,1:4]
print(dim(scmat))

## Création objet seurat 
sobj <- CreateSeuratObject(counts = scmat, 
                           min.cells = 5, 
                           min.features = 200,
                           project = sample.name)
# Classe de l'objet seurat
class(sobj)

# Barcodes names 
head(colnames(sobj))
# number of barcodes 
ncol(sobj)
# features names (genes) 
head(rownames(sobj))
# numebr of features (genes) 
nrow(sobj)#
# Affichage de l'objet seurat
sobj


######### CONVERT GENE ID TO GENE SYMBOL with biomaRt###############

library(biomaRt)
# Connect to Ensembl biomart for Drosophila melanogaster gene annotation
mart <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

# Retrieve mapping from Ensembl gene IDs to gene symbols
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = rownames(sobj),
                      mart = mart)

# Replace gene IDs with gene symbols in Seurat object
# Create a lookup table to map gene IDs to gene symbols
gene_lookup <- setNames(gene_mapping$external_gene_name, gene_mapping$ensembl_gene_id)

# Update Seurat object's gene names
rownames(sobj) <- ifelse(rownames(sobj) %in% names(gene_lookup),
                         gene_lookup[rownames(sobj)],
                         rownames(sobj))

# View updated Seurat object with gene symbols
sobj
head(rownames(sobj))

rm(scmat)


## Descriptive analysis 
# number total of UMI dans la matrice : 

sum(sobj@assays$RNA@layers$counts)

# Meta.data 

sobj@meta.data$nCount_RNA
sobj$nCount_RNA
sobj$nFeature_RNA

summary(sobj$nCount_RNA)
summary(sobj$nFeature_RNA)


### Basic préprocessing #####

## Low quality cells : 

summary(sobj$nFeature_RNA)
hist( sobj$nFeature_RNA,
      breaks = 200,
      xlab = "n genes",
      main = "Number of genes detected / cell"
)
#violin plot + horizontal line to help set threshold
VlnPlot(sobj, features = c("nFeature_RNA"),alpha = 0.05)+ geom_hline(yintercept = c(350,1500))
VlnPlot(sobj, features = c("nFeature_RNA"),log = TRUE,alpha = 0.05)+ geom_hline(yintercept = c(1500,10000))


summary(sobj$nCount_RNA)
#violin plot + horizontal line to help set threshold
VlnPlot(sobj, features = c("nCount_RNA"),alpha=0.1)+ geom_hline(yintercept = c(1000,20000))
# in a log scale
VlnPlot(sobj, features = c("nCount_RNA"),log= TRUE ,alpha=0.05)+ geom_hline(yintercept = c(10000,150000))





## Courbe count / gène 
FeatureScatter(sobj, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") +
  geom_hline(yintercept = 50000, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 5000, linetype = "dashed", color = "blue")


# DEFINE THE CUTOFF VALUES FOR EACH VARIABLE AND VISUALIZE THE RESULTS.
minGene    = 1500
minUMI     = 10000
maxGene    = 15000
maxUMI     = 150000

# we create a new seurat object containing the filtred cells
sobj_filtrd <- subset(sobj, subset = 
                        nFeature_RNA > minGene & nFeature_RNA < maxGene &    # Filtrage par nFeature_RNA
                        nCount_RNA > minUMI & nCount_RNA < maxUMI )


VlnPlot(sobj_filtrd, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,log=TRUE, alpha=0.1)


###############LIST MITOCHONDRIAL GENES ####################
# Query mitochondrial genes in biomart
mito_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "chromosome_name",
  values = "mitochondrion_genome",
  mart = mart
)
# from this list can you guess the pattern of mitochondrial genes in drosophila ? 
head(mito_genes)
## Je dirait : mt: donc 
mito_patern <- "^mt:"


# Mito
sobj_filtrd$percent_mt <- PercentageFeatureSet(sobj, pattern = mito_patern)
summary(sobj_filtrd$percent_mt)


# Violinplot
VlnPlot(sobj_filtrd, features = c("percent_mt"))
# Histogrammes
qplot(sobj_filtrd$percent_mt, geom = "histogram", bins = 100, main = "percent_mt", xlab = "percent_mt") 

maxpct_mt  = 10
# we create a new seurat object containing the filtred cells
sobj_filtrd <- subset(sobj_filtrd, subset =  percent_mt < maxpct_mt)
VlnPlot(sobj_filtrd, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
qplot(sobj_filtrd$percent_mt, geom = "histogram", bins = 100, main = "percent_mt", xlab = "percent_mt") 


rm(sobj)

######### Doublet elimination ########

# convert seurat object to a SingleCellExperiment object
sceobj <- as.SingleCellExperiment(sobj_filtrd)

# Prédiction
doublet_res <- scDblFinder(sceobj)


# Sauvegarde des résultats dans sobj
sobj_filtrd$doublets.class <- doublet_res$scDblFinder.class

#  % of doublets 
table(sobj_filtrd$doublets.class)

sobj_filtrd <- subset(sobj_filtrd, subset = doublets.class == "singlet")
VlnPlot(sobj_filtrd, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

###### Normalisation ##### 


## Méthode SCTransform 
sobj_filtrd <- SCTransform(sobj_filtrd, vars.to.regress = c("nCount_RNA", "percent_mt"))

### PCA #####
DefaultAssay(sobj_filtrd) <- "SCT"
VariableFeatures(sobj_filtrd)



# Vérifie le nombre de gènes et de cellules dans l'objet
dim(sobj_filtrd)

sobj_filtrd <- RunPCA(sobj_filtrd, features = VariableFeatures(sobj_filtrd))



  # Examine and visualize PCA results :
print(sobj_filtrd[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sobj_filtrd, dims = 1:3, reduction = "pca")
# Plot expression of gene driving principal component 
DimHeatmap(sobj_filtrd, dims = 1:1)

DimPlot(sobj_filtrd, reduction = "pca")
ElbowPlot(sobj_filtrd)


#### GO enrichement des 3 première PC####

library(clusterProfiler)
library(org.Dm.eg.db)


top_genes_PC1 <- TopFeatures(sobj_filtrd[["pca"]], dim = 1, nfeatures = 50)
top_genes_PC2 <- TopFeatures(sobj_filtrd[["pca"]], dim = 2, nfeatures = 50)
top_genes_PC3 <- TopFeatures(sobj_filtrd[["pca"]], dim = 3, nfeatures = 50)
# Analyse GO
go_enrichmentPC1 <- enrichGO(gene          = top_genes_PC1,
                          OrgDb         = org.Dm.eg.db,
                          keyType       = "SYMBOL",  # Type d'identifiants (SYMBOL, ENSEMBL, etc.)
                          ont           = "BP",     # BP = Biological Process
                          pAdjustMethod = "BH",     # Correction pour tests multiples
                          pvalueCutoff  = 0.05,
                          readable      = TRUE)

go_enrichmentPC2 <- enrichGO(gene          = top_genes_PC2,
                          OrgDb         = org.Dm.eg.db,
                          keyType       = "SYMBOL",  # Type d'identifiants (SYMBOL, ENSEMBL, etc.)
                          ont           = "BP",     # BP = Biological Process
                          pAdjustMethod = "BH",     # Correction pour tests multiples
                          pvalueCutoff  = 0.05,
                          readable      = TRUE)


go_enrichmentPC3 <- enrichGO(gene          = top_genes_PC3,
                          OrgDb         = org.Dm.eg.db,
                          keyType       = "SYMBOL",  # Type d'identifiants (SYMBOL, ENSEMBL, etc.)
                          ont           = "BP",     # BP = Biological Process
                          pAdjustMethod = "BH",     # Correction pour tests multiples
                          pvalueCutoff  = 0.05,
                          readable      = TRUE)

head(go_enrichmentPC1)
barplot(go_enrichmentPC1, showCategory = 10)
dotplot(go_enrichment, showCategory = 10)

head(go_enrichmentPC2)
barplot(go_enrichmentPC2, showCategory = 10)

head(go_enrichmentPC3)
barplot(go_enrichmentPC3, showCategory = 10)

### Clustering #####


# number of PC kepts for the analysis
nPC = 15
# s-nn graph
sobj_filtrd <- FindNeighbors(sobj_filtrd, dims = 1:nPC) 
# make the clusters (here we try 3 resolutions)
sobj_filtrd <- FindClusters(sobj_filtrd, resolution = c(0.2,0.5,1)) 
# the clusters are stored in  `metadata`:
head(sobj_filtrd@meta.data)


sobj_filtrd <- RunUMAP(sobj_filtrd, dims = 1:nPC)
DimPlot(sobj_filtrd, reduction = "umap", group.by = "SCT_snn_res.0.2")  # Pour résolution 0.2
DimPlot(sobj_filtrd, reduction = "umap", group.by = "SCT_snn_res.0.5")  # Pour résolution 0.5
DimPlot(sobj_filtrd, reduction = "umap", group.by = "RNA_snn_res.1")    # Pour résolution 1


VlnPlot(sobj_filtrd, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), group.by = "seurat_clusters")

Idents(sobj_filtrd)

Idents(sobj_filtrd) <- "SCT_snn_res.0.2"

###############LIST MARKER GENES FOR CLUSTER IDENTIFICATION###
#Here is a list provided by biologists.
# You should check first that the genes are expressed in the dataset

# genes majoritarily expressed in muscles
list_muscle  = c('blow','Grip','Him','htl','Pdp1', 'sls')

# genes majoritarily expressed in  muscle precursors
list_pre_muscle  = c('sns','Strn-Mlck','TpnC73F','sis') 

# genes majoritarily expressed in tendons
list_tendon     = c('Alk','Hand','kon','vkg')

# genes majoritarily expressed in SOP (sop = précurseurs des organes sensoriels externes)
list_sop        = c( 'ase','cpo', 'sca','sens')

# genes majoritarily expressed in hemocytes 
list_hemo       = c('eater','et','NimC4', 'srp')

# genes majoritarily expressed in epihelium
list_epi        = c('arm','baz','grh','hth','shg','zip','eyg')


# Vérifier si les gènes sont dans le jeu de données
all_genes <- rownames(sobj_filtrd)
list_muscle[!list_muscle %in% all_genes]  # Affiche les gènes absents
list_pre_muscle[!list_pre_muscle %in% all_genes]  
list_hemo[!list_hemo %in% all_genes] 
list_sop[!list_sop %in% all_genes]  
list_tendon[!list_tendon %in% all_genes]  
list_epi[!list_epi %in% all_genes]  




VlnPlot(sobj_filtrd, features = list_muscle,pt.size = 0.1)
FeaturePlot(sobj_filtrd, features = list_muscle)
for (gene in list_muscle) {
  print(FeaturePlot(sobj_filtrd, features = gene) + ggtitle(paste("UMAP for", gene)))
}

VlnPlot(sobj_filtrd, features = list_pre_muscle,pt.size = 0.1)
FeaturePlot(sobj_filtrd, features = list_pre_muscle)
for (gene in list_pre_muscle) {
  print(FeaturePlot(sobj_filtrd, features = gene) + ggtitle(paste("UMAP for", gene)))
}

VlnPlot(sobj_filtrd, features = list_hemo,pt.size = 0.1)
FeaturePlot(sobj_filtrd, features = list_pre_muscle)
for (gene in list_pre_muscle) {
  print(FeaturePlot(sobj_filtrd, features = gene) + ggtitle(paste("UMAP for", gene)))
}

VlnPlot(sobj_filtrd, features = list_sop,pt.size = 0.5)
FeaturePlot(sobj_filtrd, features = list_pre_muscle)
for (gene in list_pre_muscle) {
  print(FeaturePlot(sobj_filtrd, features = gene) + ggtitle(paste("UMAP for", gene)))
}

VlnPlot(sobj_filtrd, features = list_tendon,pt.size = 0.1)
FeaturePlot(sobj_filtrd, features = list_pre_muscle)
for (gene in list_pre_muscle) {
  print(FeaturePlot(sobj_filtrd, features = gene) + ggtitle(paste("UMAP for", gene)))
}

VlnPlot(sobj_filtrd, features = list_epi,pt.size = 0.5)
FeaturePlot(sobj_filtrd, features = list_pre_muscle)
for (gene in list_pre_muscle) {
  print(FeaturePlot(sobj_filtrd, features = gene) + ggtitle(paste("UMAP for", gene)))
}



### Add module score

library(RColorBrewer)

marker_genes <- list(
  muscle  = c('blow','Grip','Him','htl','Pdp1', 'sls'),
  pre_muscle  = c('sns','Strn-Mlck','TpnC73F','sis') ,
  tendon     = c('Alk','Hand','kon','vkg'),
  sop        = c('ase','cpo','sca','sens'),
  hemo       = c('eater','et','NimC4', 'srp'),
  epi        = c('arm','baz','grh','hth','shg','zip','eyg'))


#Missing : "sis" "eater" "zip"


indices_to_remove <- c()
j <- 1

for (i in seq(1, length(marker_genes))) {
  marker_genes[[i]] <- marker_genes[[i]][which(marker_genes[[i]] %in% rownames(sobj_filtrd))]
  if(length(marker_genes[[i]]) > 1) {
    sobj_filtrd <- AddModuleScore(sobj_filtrd,
                                  features = marker_genes[i],
                                  pool = NULL,
                                  nbin = 5,
                                  seed = 1,
                                  ctrl = length(marker_genes[i]),
                                  k = FALSE,
                                  name = names(marker_genes[i]))
    col_name <- names(marker_genes)[[i]]
    col_name_val <- which(colnames(sobj_filtrd[[]]) == paste0(col_name, 1))
    colnames(sobj_filtrd@meta.data)[col_name_val] <- col_name
  } else {
    indices_to_remove[j] <- i
    j <- j + 1
  }
}
marker_genes[indices_to_remove] <- NULL
print(marker_genes)

#Scores for UMAP
for (f in names(marker_genes)){
  print(FeaturePlot(sobj_filtrd,
                    features = f,
                    dims = c(1, 2),
                    cols = c("grey90", brewer.pal(9,"YlGnBu")),
                    pt.size = 0.2,
                    ncol = 1)+NoLegend())}

for (f in names(marker_genes)){
  print(VlnPlot(sobj_filtrd, features = f))

table(Idents(sobj_filtrd))

sobj_filtrd <- RenameIdents(
  object = sobj_filtrd,
  '0' = 'Epi',
  '1' = 'Epi',
  '2' = 'Epi',
  '3' = 'Muscle',
  '4' = 'Hemo',
  '5'='SOP',
  '6'='Pre_Muscle',
  '7'='Tendon',
  '8'='Epi'
)

DimPlot(sobj_filtrd, reduction = "umap", label = TRUE) + NoLegend()


cluster_markers <- FindMarkers(
  sobj_filtrd,
  ident.1 = "0",            # Identifier les marqueurs pour le cluster 0
  min.pct = 0.25,
  logfc.threshold = 0.25
)
markers <- FindAllMarkers(
  sobj_filtrd,
  only.pos = TRUE,       # Ne conserver que les gènes exprimés positivement
  min.pct = 0.25,        # Gènes exprimés dans au moins 25% des cellules d'un cluster
  logfc.threshold = 0.25 # Seuil pour le Log Fold Change
)

# Afficher les 10 premiers marqueurs identifiés
head(markers, 10)


top_genes_all <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3)  # Top 3 gènes avec le plus grand log2FC
print (n=28,top_genes_all)
write.csv(top_genes_all, "Top3_FindAllMarkers.csv", row.names = FALSE)


all_clusters <- levels(Idents(sobj_filtrd))  # Liste de tous les clusters
markers_list <- lapply(all_clusters, function(cluster) {
  FindMarkers(sobj_filtrd, ident.1 = cluster, min.pct = 0.25, logfc.threshold = 0.25)
})
names(markers_list) <- all_clusters

write.csv(markers, "FindAllMarkers_results.csv", row.names = FALSE)


all_markers_combined <- do.call(rbind, lapply(names(markers_list), function(cluster) {
  markers <- markers_list[[cluster]]
  markers$cluster <- cluster  # Ajouter une colonne avec le nom du cluster
  markers
}))
write.csv(all_markers_combined, "FindMarkers_combined_results.csv", row.names = FALSE)

top_genes_list <- lapply(names(markers_list), function(cluster) {
  markers <- markers_list[[cluster]]
  top_genes <- markers %>%
    arrange(desc(avg_log2FC)) %>%
    head(3)  # Prendre les 3 premiers gènes avec le plus grand log2FC
  top_genes$cluster <- cluster
  top_genes
})
top_genes_combined <- do.call(rbind, top_genes_list)
print(( top_genes_combined))
write.csv(top_genes_combined, "Top3_FindMarkers.csv", row.names = FALSE)


  
}
