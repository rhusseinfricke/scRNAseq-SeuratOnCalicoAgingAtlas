# March 2024

# Roux, A. E., Yuan, H., Podshivalova, K., Hendrickson, D., Kerr, R., Kenyon, C., & Kelley, D. (2023).
#Individual cell types in C. elegans age differently and activate distinct cell-protective responses. Cell reports, 42(8), 112902.
#https://doi.org/10.1016/j.celrep.2023.112902

#single-cell RNA-seq from published data

# Reproducing Calico data set using Seurat package

# single-cell RNA-seq analysis - Quality Control (QC)

#Exploring Calico Data Set

#In this study the authors characterise a substanitally large amount of C. elegans cell types isolated from aging worm samples.
#The data in this study identifies CELL-TYPE-SPECIFIC changes in gene expression in aging worms. 
#The significance of this work is that the changes in gene expression could mean that the aging cells are activating mechanisms that
#increase cellular resilience through proteostasis, repair and other mechanisms.
#While the discovery of aging-related genes in C. elegans has facilitated discovery of conserved aging cell regulators, genome-wide analysis
#has been limited to whole animal studies. As a consequence, there is a dearth of longevity data at the cellular level and in turn in the
#various tissue of multicellular organisms.

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(SeuratData)
library(SparseArray)
library(sparseMatrixStats)
library(SparseM)

#Adjusted names of the features, barcodes and matrix files in the folder of each technical replicate so that they are named
#features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz respectively so that they can be read with the Read10X command (do not give the files a prefix to identify the sample, only
#name folders with the sample name, otherwise error message = file missing)

#setwd("~/Desktop/Calico/Data"): set this as your working directory

#day_1_T1_counts <- Read10X(data.dir = "day_1_T1")
#day_1_T2_counts <- Read10X(data.dir = "day_1_T2")
#day_3_T1_counts <- Read10X(data.dir = "day_3_T1")
#day_3_T2_counts <- Read10X(data.dir = "day_3_T2")
#day_5_T1_counts <- Read10X(data.dir = "day_5_T1")
#day_5_T2_counts <- Read10X(data.dir = "day_5_T2")
#day_8_T1_counts <- Read10X(data.dir = "day_8_T1")
#day_8_T2_counts <- Read10X(data.dir = "day_8_T2")
#day_11_T1_counts <- Read10X(data.dir = "day_11_T1")
#day_11_T2_counts <- Read10X(data.dir = "day_11_T2")
#day_15_counts <- Read10X(data.dir = "day_15")

#Create Seurat object for all the technical replicate. Seurat automatically creates meta data for each cell when Read10X is used
#The meta data is stored in meta.data
#Read10X part of loop below does not work

for (file in c("day_1_T1", "day_1_T2", "day_3_T1", "day_3_T2", "day_5_T1", "day_5_T2", "day_8_T1", "day_8_T2", "day_11_T1", "day_11_T2", "day_15")){
  seurat_data_ALLT <- Read10X(data.dir = paste0("~/Desktop/Calico/Data/", file))
  seurat_obj_ALLT <- CreateSeuratObject(counts = seurat_data_ALLT, min.features = 100, project = file)
  assign(file, seurat_obj_ALLT)
}

#Explore the meta data created by CreateSeuratObject
head(seurat_obj_ALLT@meta.data)
tail(seurat_obj_ALLT@meta.data)

#orig. ident. column contains the SAMPLE identity if known (this defaults to "SeuratProject")
#nCount_RNA is the number of UMIs per cell
#nFeature_RNA is the number of genes per cell

#Merge the different layers of the Seurat object into a single Seurat object so that it's easier
#to compare data quality for all the samples and run QC for all the samples at the same time.

merged_seuratObjALLT <- merge(x = day_15, y = c(day_1_T1, day_1_T2, day_3_T1, day_3_T2, day_5_T1, day_5_T2, day_8_T1, day_8_T2, day_11_T1, day_11_T2), add.cell.id = c("D15", "D1T1", "D1T2", "D3T1", "D3T2", "D5T1", "D5T2", "D8T1", "D8T2", "D11T1", "D11T2"))


#Since the same cell can stem from a different sample, the same CELL ID can be used so, a SAMPLE-SPECIFIC PREFIX is added to 
#each cell id. Prefixes in the row names can be seen in merged Seurat metadata

head(merged_seuratObjALLT@meta.data)
tail(merged_seuratObjALLT@meta.data)

#A look at the Seurat object metadata (merged_seuratObjALLT) shows columns labeled:
#orig.ident: this is the sample name (was argument added to 'project =' when data was loaded to create Seurat object)
#nCount_RNA: number of Unique Molecular Identifiers per cell (UMI per cell)
#nFeature_RNA: number of genes expressed per cell

#View(merged_seuratObjALLT)

#NOVELTY SCORE
#This is number of genes expressed per UMI and is calculated by taking log10 of nCount_RNA and dividing it by 
#log10 of nFeature_RNA

merged_seuratObjALLT$log10GenePerUMI <- log10(merged_seuratObjALLT$nFeature_RNA)/log10(merged_seuratObjALLT$nCount_RNA)


#MITOCHONDRIAL RATIO
#Use PercentageFeatureSet() function which computes the percentage/proportion of transcripts that map to mitochondrial genes. 
#Function take sum count of features (genes) that are mitochondrial (identified by MTCE- in the case of C. elegans) and divides by sum count of all features (genes)
#then multiplies by 100. Function uses 'pattern=' argument in code to identify mitochondrial gene identifiers with the MTCE- nomenclature pattern.
#Since we are searching C. elegans mitochondrial genes we must use the pattern 'MTCE-' instead of 'MT-" used for human mitochondrial gene identifiers.

merged_seuratObjALLT$MitoRatio <- PercentageFeatureSet(object = merged_seuratObjALLT, pattern = "^MTCE.")

#Will use ratio instead of percentage in analysis so divide mitochondrail output in merged_seuratObjALLT$MitoRatio by 100

merged_seuratObjALLT$MitoRatio <- merged_seuratObjALLT@meta.data$MitoRatio/100

#Preparation of quality metrics for data assessment is now complete. 
#Additional metadata information can be added to metadata slots within Seurat object using $ operator (like we did above).
#To avoid impacting metadata in Seurat object, this metadata will be extracted and kept in a separate data frame (by
#assigning it to variable) and then we will continue to add to the metadata slots.
#The following extracts the metadata slot from the Seurat object.

metadataCalico <- merged_seuratObjALLT@meta.data

#The metadata dataframe rownames are cell identifiers. Will add a cell identifier column to the metadata dataframe by duplicating the rownames and making it a metadata column called cells.

metadataCalico$cells <- rownames(metadataCalico)

#Using the cell ID prefixes showing the time points the samples were collected, an additional column called 'sample'
#is added to the metadata with the values in that column indicating the time points.

metadataCalico$sample <- NA
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D15_"))] <- "D15"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D11T2_"))] <- "D11T2"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D11T1_"))] <- "D11T1"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D8T2_"))] <- "D8T2"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D8T1_"))] <- "D8T1"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D5T2_"))] <- "D5T2"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D5T1_"))] <- "D5T1"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D3T2_"))] <- "D3T2"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D3T1_"))] <- "D3T1"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D1T2_"))] <- "D1T2"
metadataCalico$sample[which(str_detect(metadataCalico$cells, "^D1T1_"))] <- "D1T1"

#Rename metadatq columns that were automatically created by Seurat Read10X(), namely: orig.ident (sample identity)
#, nCount_RNA (number of UMI per cell) and nFeature_RNA (number of genes per cell).

metadataCalico <- metadataCalico %>% dplyr::rename(SeqFolder = orig.ident, nUMI = nCount_RNA, nGenes = nFeature_RNA)

#Add metadata data frame (metadataCalico) back to Seurat object

merged_seuratObjALLT@meta.data <- metadataCalico

#Create .RData object. The object can be loaded at any time.

save(merged_seuratObjALLT, file = "~/Desktop/Calico/Data/mergedFiltered_seuratObjALLT.RData")

#Sometimes there can be a higher number of cellular barcodes than cells (numbers in data can be higher than what was loaded)
#Ex.: a hydrogel cn have more than one barcode in it.
#Visual showing number of cell counts per sample

metadataCalico %>% ggplot(aes(x=sample, fill=sample)) + 
                            geom_bar() +
                            theme_classic() +
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                            theme(plot.title = element_text(hjust=0.5, face="bold")) +
                            ggtitle("Number of cell counts per time point")

#Visualize the number of UMI counts (number of transcripts) per cell

metadataCalico %>% 
  ggplot(aes(color=sample, x=SeqFolder, y=nUMI, fill= sample)) + 
  geom_boxplot(alpha = 0.2) + 
  scale_y_log10() +
  theme_classic() +
  xlab("day and technical replicate") +
  ylab("log10(UMI count)") +
  theme(axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 14)) +
  geom_hline(yintercept = 500)

#Visualize the gene distribution per cell (expected to be lower than UMI counts per cell) but similar
#expecation.

metadataCalico %>% 
  ggplot(aes(color=sample, x=SeqFolder, y=nGenes, fill= sample)) + 
  geom_violin(alpha = 0.2) + 
  theme_classic() +
  scale_y_log10() +
  xlab("day and technical replicate") +
  ylab("Number of Genes per Cell") +
  theme(axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 14)) +
  geom_hline(yintercept = 300)


#Visualize the novelty score to determine the complexity of the data. The novelty score gives the 
#number of genes per nUMI. If we have a large UMI number and a small number of genes then there might be
#a contamination or something that led to capturing a small number of genes then sequencing transcripts
#repeatedly from that small number of genes.

metadataCalico %>% ggplot(aes(x=log10GenePerUMI, color=sample, fill=sample))+geom_density(alpha=0.2)+theme_classic()+geom_vline(xintercept = 0.8)

#Since density plot makes it hard to discern data points, used boxplot

metadataCalico %>% ggplot(aes(x=SeqFolder, y=log10GenePerUMI, color=sample, fill=sample))+
  geom_boxplot(alpha=0.2)+
  theme_classic()+geom_hline(yintercept = 0.8)+
  theme(axis.text.x = element_text(angle=90, size=14), axis.text.y = element_text(size = 14))

#Visualize the distribution of mitochondrial gene expression per cell to get an idea of the number of dead or dying cells.
#Poor quality samples have a mitochondrial counts ratio that is higher than 0.2.

metadataCalico %>% ggplot(aes(color=sample, x=MitoRatio, fill=sample))+
  geom_density(alpha=0.2)+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

#To better visualize the distribution, violin plot

metadataCalico %>% ggplot(aes(color=sample, x=SeqFolder, y=MitoRatio, fill=sample))+
  geom_violin(alpha=0.2)+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 14)) +
  xlab("day and technical replicate") +
  ylab("MitoRatio") +
  geom_hline(yintercept = 0.2)
  
#Visualize the number of genes vs. the number of UMIs to understand the quality per cell
#Is there a high presence of cells with low numbers of genes and UMIs?
#Visualize correlation between the number of genes and number of UMIs
#Tried adjusting the xintercept and yintercept values to 100

metadataCalico %>% ggplot(aes(x=nUMI, y=nGenes, color=MitoRatio))+
  geom_point()+
  scale_color_gradient(low="gray90", high="black")+
  stat_smooth(method=lm)+
  scale_x_log10()+
  scale_y_log10()+
  theme_classic()+
  geom_vline(xintercept = 500)+
  geom_hline(yintercept = 250)+
  xlab("Number of UMI") +
  ylab("Number of Genes") +
  facet_wrap(~SeqFolder)

#Cell-level filtering
#Filter out low quality cells using thresholds selected by user (the thresholds vary with experiment)

filtered_seuratObjALLT <- subset(x=merged_seuratObjALLT, subset=(nUMI>=100)&
  (nGenes >=100)&
  (log10GenePerUMI>0.8)&
  (MitoRatio<0.2))

#Gene-level filtering:
#Use GetAssayData() to extract counts slot from cell-level filtered object

countCalico <- GetAssayData(object = filtered_seuratObjALLT, slot = "counts") 


#received following error:
#The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
#ℹ Please use the `layer` argument instead.
#This warning is displayed once every 8 hours.
#Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 

#GetAssayData has been updated for Seurat v.5. Use LayerData().

countCalico <- LayerData(object = filtered_seuratObjALLT, layer = "counts")

#Matrix with logical output that specifies whether or not there is a result of greater than 0 for each gene
nonzero <- countCalico > 0

#Genes expressed in a small number of cells are filtered out because they will skew down the average for all cells.
#By keeping genes that are expressed in only 10 or more cells, genes that have a zero count expression in all cells will be removed

keep_genesCalico <- Matrix::rowSums(nonzero) >= 10

#Genes expressed in 10 or more cells are kept

filtered_seuratObjALLT_Counts <- countCalico[keep_genesCalico, ]

#Use the above filtered counts to create a new Seurat object by reassigning the counts objects to the filtered Seurat object.

filteredObjALLT <- CreateSeuratObject(counts = filtered_seuratObjALLT_Counts, meta.data = filtered_seuratObjALLT@meta.data)

#Create .RData object for filtered cell object. This can be loaded at any time and this is what is used for clustering and for identifying markers.












#OR

#Read each technical replicate file and create Seurat object for each technical replicate

for (file in c("day_1_T1", "day_1_T2")){
  seurat_dataD1 <- Read10X(data.dir = paste0("Data/", file))
  seurat_objD1 <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)}

for (file in c("day_3_T1", "day_3_T2")){
  seurat_dataD3 <- Read10X(data.dir = paste0("Data/", file))
  seurat_objD3 <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)}
  
for (file in c("day_5_T1", "day_5_T2")){
  seurat_dataD5 <- Read10X(data.dir = paste0("Data/", file))
  seurat_objD5 <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)}
    
for (file in c("day_8_T1", "day_8_T2")){
    seurat_dataD8 <- Read10X(data.dir = paste0("Data/", file))
    seurat_objD8 <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
    assign(file, seurat_obj)}

for (file in c("day_11_T1", "day_11_T2")){
  seurat_dataD11 <- Read10X(data.dir = paste0("Data/", file))
  seurat_objD11 <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)}

for (file in c("day_15")){
  seurat_dataD15 <- Read10X(data.dir = paste0("Data/", file))
  seurat_objD15 <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)}

#Convert technical replicates Day1 to dataframe and look for duplicated rows (gene names)
#"In many cases gene names are repeated because there isn't consensus over which coding sequence represents the common name"

seurDataDF <- as.data.frame(seuratDataT1)
DPLseuratDataT1 <- which(duplicated(seurDataDF))

#Test: name the first
rownames(seuratDataT1) <- d1_T1T2_gene_ids




#Concatenate Read10X matrix output to combine technical replicates. Technical replicates have the same row names but different 
#column names and numbers so cbind() command was used

c10mtx_day_1_counts <- cbind(day_1_T1_counts,day_1_T2_counts)
c10mtx_day_3_counts <- cbind(day_3_T1_counts,day_3_T2_counts)
c10mtx_day_5_counts <- cbind(day_5_T1_counts,day_5_T2_counts)
c10mtx_day_8_counts <- cbind(day_8_T1_counts,day_8_T2_counts)
c10mtx_day_11_counts <- cbind(day_11_T1_counts,day_11_T2_counts)
c10mtx_day_15_counts <- day_15_counts #give day 15 matrix same c10mtx prefix in variable for consistency even though no technical replicates in day 15

#Convert combined technical replicate matrices into a dataframe

#DFc10mtx_day_1_counts <- as.data.frame(c10mtx_day_1_counts)

#Remove duplicated rows from combined matrices of DAY 1 technical replicates

#ROWREMc10mtx_day_1_countsDF <- !duplicated(DFc10mtx_day_1_counts)

#Remove duplicated columns from combined matrices of DAY 1 technical replicates

#DPLc10mtx_day_1_DF <- duplicated(t(DFc10mtx_day_1_counts)) #transpose matrix columns to rows and remove duplicate rows

#COLREMc10mtx_day_1_DF <- DFc10mtx_day_1_counts[, !DPLc10mtx_day_1_counts] #Not duplicated values are placed in the columns of the new matrix






#If Read10X command does not work, read in each file type with the read_tsv command:

# Read in day 1, day 3, day 5, day 8 and day 11 technical replicate 1 (T1) sample matrix.mtx files
d1_T1_mtx_counts <- readMM("~/Desktop/Calico/Data/day_1_T1/GSM6338211_matrix_TC2_d1_1.mtx.gz")

d3_T1_mtx_counts <- readMM("~/Desktop/Calico/Data/day_3_T1/GSM6338212_matrix_TC2_d3_1.mtx.gz")

d5_T1_mtx_counts <- readMM("~/Desktop/Calico/Data/day_5_T1/GSM6338213_matrix_TC2_d5_1.mtx.gz")

d8_T1_mtx_counts <- readMM("~/Desktop/Calico/Data/day_8_T1/GSM6338214_matrix_TC2_d8_1.mtx.gz")

d11_T1_mtx_counts <- readMM("~/Desktop/Calico/Data/day_11_T1/GSM6338215_matrix_TC2_d11_1.mtx.gz")


# Read in day 1, day 3, day 5, day 8 and day 11 technical replicate 2 (T2) sample matrix.mtx files
d1_T2_mtx_counts <- readMM("~/Desktop/Calico/Data/day_1_T2/GSM6338216_matrix_TC2_d1_2.mtx.gz")

d3_T2_mtx_counts <- readMM("~/Desktop/Calico/Data/day_3_T2/GSM6338217_matrix_TC2_d3_2.mtx.gz")

d5_T2_mtx_counts <- readMM("~/Desktop/Calico/Data/day_5_T2/GSM6338218_matrix_TC2_d5_2.mtx.gz")

d8_T2_mtx_counts <- readMM("~/Desktop/Calico/Data/day_8_T2/GSM6338219_matrix_TC2_d8_2.mtx.gz")

d11_T2_mtx_counts <- readMM("~/Desktop/Calico/Data/day_11_T2/GSM6338220_matrix_TC2_d11_2.mtx.gz")

# Read in day 15 sample matrix.mtx
d15_mtx_counts <- readMM("~/Desktop/Calico/Data/day_15/GSM6338221_matrix_TC2_d15_1.mtx.gz")


# Read in day 1, day 3, day 5, day 8 and day 11 technical replicate 1 (T1) sample genes.tsv files
d1_T1_genes <- read_tsv("~/Desktop/Calico/Data/day_1_T1/GSM6338211_features_TC2_d1_1.tsv.gz", col_names = FALSE)
d1_T1_gene_ids <- d1_T1_genes$X1

d3_T1_genes <- read_tsv("~/Desktop/Calico/Data/day_3_T1/GSM6338212_features_TC2_d3_1.tsv.gz", col_names = FALSE)
d3_T1_gene_ids <- d3_T1_genes$X1

d5_T1_genes <- read_tsv("~/Desktop/Calico/Data/day_5_T1/GSM6338213_features_TC2_d5_1.tsv.gz", col_names = FALSE)
d5_T1_gene_ids <- d5_T1_genes$X1

d8_T1_genes <- read_tsv("~/Desktop/Calico/Data/day_8_T1/GSM6338214_features_TC2_d8_1.tsv.gz", col_names = FALSE)
d8_T1_gene_ids <- d8_T1_genes$X1

d11_T1_genes <- read_tsv("~/Desktop/Calico/Data/day_11_T1/GSM6338215_features_TC2_d11_1.tsv.gz", col_names = FALSE)
d11_T1_gene_ids <- d11_T1_genes$X1

# Read in day 1, day 3, day 5, day 8 and day 11 technical replicate 2 (T2) sample genes.tsv files
d1_T2_genes <- read_tsv("/Users/reemhussein-fricke/Desktop/Calico/Data/day_1_T2/GSM6338216_features_TC2_d1_2.tsv.gz", col_names = FALSE)
d1_T2_gene_ids <- d1_T2_genes$X1

d3_T2_genes <- read_tsv("~/Desktop/Calico/Data/day_3_T2/GSM6338217_features_TC2_d3_2.tsv.gz", col_names = FALSE)
d3_T2_gene_ids <- d3_T2_genes$X1

d5_T2_genes <- read_tsv("~/Desktop/Calico/Data/day_5_T2/GSM6338218_features_TC2_d5_2.tsv.gz", col_names = FALSE)
d5_T2_gene_ids <- d5_T2_genes$X1

d8_T2_genes <- read_tsv("~/Desktop/Calico/Data/day_8_T2/GSM6338219_features_TC2_d8_2.tsv.gz", col_names = FALSE)
d8_T2_gene_ids <- d8_T2_genes$X1

d11_T2_genes <- read_tsv("~/Desktop/Calico/Data/day_11_T2/GSM6338220_features_TC2_d11_2.tsv.gz", col_names = FALSE)
d11_T2_gene_ids <- d11_T2_genes$X1

# Read in day 15 sample genes.tsv files
d15_genes <- read_tsv("~/Desktop/Calico/Data/day_15/GSM6338221_features_TC2_d15_1.tsv.gz", col_names = FALSE)
d15_gene_ids <- d15_genes$X1

# Read in day 1, day 3, day 5, day 8 and day 11 technical replicate 1 (T1) sample barcodes.tsv cell id files
d1_T1_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_1_T1/GSM6338211_barcodes_TC2_d1_1.tsv.gz", col_names = FALSE)$X1

d3_T1_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_3_T1/GSM6338212_barcodes_TC2_d3_1.tsv.gz", col_names = FALSE)$X1

d5_T1_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_5_T1/GSM6338213_barcodes_TC2_d5_1.tsv.gz", col_names = FALSE)$X1

d8_T1_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_8_T1/GSM6338214_barcodes_TC2_d8_1.tsv.gz", col_names = FALSE)$X1

d11_T1_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_11_T1/GSM6338215_barcodes_TC2_d11_1.tsv.gz", col_names = FALSE)$X1

# Read in day 1, day 3, day 5, day 8 and day 11 technical replicate 2 (T2) sample barcodes.tsv cell id files
d1_T2_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_1_T2/GSM6338216_barcodes_TC2_d1_2.tsv.gz", col_names = FALSE)$X1

d3_T2_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_3_T2/GSM6338217_barcodes_TC2_d3_2.tsv.gz", col_names = FALSE)$X1

d5_T2_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_5_T2/GSM6338218_barcodes_TC2_d5_2.tsv.gz", col_names = FALSE)$X1

d8_T2_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_8_T2/GSM6338219_barcodes_TC2_d8_2.tsv.gz", col_names = FALSE)$X1

d11_T2_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_11_T2/GSM6338220_barcodes_TC2_d11_2.tsv.gz", col_names = FALSE)$X1

# Read in day 15 sample barcodes.tsv cell id files
d15_cell_ids <- read_tsv("~/Desktop/Calico/Data/day_15/GSM6338221_barcodes_TC2_d15_1.tsv.gz", col_names = FALSE)$X1

# Create a sparse matrix for more efficient computation (*this step is only necessary if Read10X command from Seurat package was not used).

d1_T1_mtx_counts <- as(d1_T1_mtx_counts,"sparseMatrix")
d3_T1_mtx_counts <- as(d3_T1_mtx_counts,"sparseMatrix")
d5_T1_mtx_counts <- as(d5_T1_mtx_counts,"sparseMatrix")
d8_T1_mtx_counts <- as(d8_T1_mtx_counts,"sparseMatrix")
d11_T1_mtx_counts <- as(d11_T1_mtx_counts,"sparseMatrix")
d1_T2_mtx_counts <- as(d1_T2_mtx_counts,"sparseMatrix")
d3_T2_mtx_counts <- as(d3_T2_mtx_counts,"sparseMatrix")
d5_T2_mtx_counts <- as(d5_T2_mtx_counts,"sparseMatrix")
d8_T2_mtx_counts <- as(d8_T2_mtx_counts,"sparseMatrix")
d11_T2_mtx_counts <- as(d11_T2_mtx_counts,"sparseMatrix")
d15_mtx_counts <- as(d15_mtx_counts,"sparseMatrix")

# Add row names to the technical replicate 1 (T1) counts matrix to be the gene ids 
#(*this step is only necessary if Read10X command from Seurat package was not used).
rownames(d1_T1_mtx_counts) <- d1_T1_gene_ids

rownames(d3_T1_mtx_counts) <- d3_T1_gene_ids

rownames(d5_T1_mtx_counts) <- d5_T1_gene_ids

rownames(d8_T1_mtx_counts) <- d8_T1_gene_ids

rownames(d11_T1_mtx_counts) <- d11_T1_gene_ids


# Add row names to the technical replicate 2 (T2) counts matrix to be the gene ids
#(*this step is only necessary if Read10X command from Seurat package was not used).
rownames(d1_T2_mtx_counts) <- d1_T2_gene_ids

rownames(d3_T2_mtx_counts) <- d3_T2_gene_ids

rownames(d5_T2_mtx_counts) <- d5_T2_gene_ids

rownames(d8_T2_mtx_counts) <- d8_T2_gene_ids

rownames(d11_T2_mtx_counts) <- d11_T2_gene_ids

# Add row names to day 15 sample counts matrix to be the gene ids
#(*this step is only necessary if Read10X command from Seurat package was not used).
rownames(d15_mtx_counts) <- d15_gene_ids

# Add the column names to the technical replicate 1 (T1) count matrix to be the cell ids
#(*this step is only necessary if Read10X command from Seurat package was not used).

colnames(d1_T1_mtx_counts) <- d1_T1_cell_ids

colnames(d3_T1_mtx_counts) <- d3_T1_cell_ids

colnames(d5_T1_mtx_counts) <- d5_T1_cell_ids

colnames(d8_T1_mtx_counts) <- d8_T1_cell_ids

colnames(d11_T1_mtx_counts) <- d11_T1_cell_ids


# Add the column names to the technical replicate 2 (T2) count matrix to be the cell ids
#(*this step is only necessary if Read10X command from Seurat package was not used).

colnames(d1_T2_mtx_counts) <- d1_T2_cell_ids

colnames(d3_T2_mtx_counts) <- d3_T2_cell_ids

colnames(d5_T2_mtx_counts) <- d5_T2_cell_ids

colnames(d8_T2_mtx_counts) <- d8_T2_cell_ids

colnames(d11_T2_mtx_counts) <- d11_T2_cell_ids

# Add the column names to day 15 sample count matrix to be the cell ids
#(*this step is only necessary if Read10X command from Seurat package was not used).

colnames(d15_mtx_counts) <- d15_cell_ids

#Check number of rows and number of columns of each technical replicate count matrix
#(*this step is only necessary if Read10X command from Seurat package was not used).

dim(d1_T1_mtx_counts)
dim(d1_T2_mtx_counts)
dim(d3_T1_mtx_counts)
dim(d3_T2_mtx_counts)
dim(d5_T1_mtx_counts)
dim(d5_T2_mtx_counts)
dim(d8_T1_mtx_counts)
dim(d8_T2_mtx_counts)
dim(d11_T1_mtx_counts)
dim(d11_T2_mtx_counts)
dim(d15_mtx_counts)

#Check that the row names of the technical replicates from each time point match
#(*this step is only necessary if Read10X command from Seurat package was not used).

d1_rowNameCheck <- rownames(d1_T1_mtx_counts)==rownames(d1_T2_mtx_counts)
d3_rowNameCheck <- rownames(d3_T1_mtx_counts)==rownames(d3_T2_mtx_counts)
d5_rowNameCheck <- rownames(d5_T1_mtx_counts)==rownames(d5_T2_mtx_counts)
d8_rowNameCheck <- rownames(d8_T1_mtx_counts)==rownames(d8_T2_mtx_counts)
d11_rowNameCheck <- rownames(d11_T1_mtx_counts)==rownames(d11_T2_mtx_counts)
d15_rowNameCheck <- rownames(d15_mtx_counts)==rownames(d15_mtx_counts)


#Check that the column names of the technical replicates from each time point are unique
#(*this step is only necessary if Read10X command from Seurat package was not used).

d1_colNameCheck <- colnames(d1_T1_mtx_counts)==colnames(d1_T2_mtx_counts)
d3_colNameCheck <- colnames(d3_T1_mtx_counts)==colnames(d3_T2_mtx_counts)
d5_colNameCheck <- colnames(d5_T1_mtx_counts)==colnames(d5_T2_mtx_counts)
d8_colNameCheck <- colnames(d8_T1_mtx_counts)==colnames(d8_T2_mtx_counts)
d11_colNameCheck <- colnames(d11_T1_mtx_counts)==colnames(d11_T2_mtx_counts)
d15_colNameCheck <- colnames(d15_mtx_counts)==colnames(d15_mtx_counts) #day 15 has no technical replicates so should be TRUE

#Concatenate the two technical replicate matrices for each time point to generate one count matrix
#(*this step is only necessary if Read10X command from Seurat package was not used).

cmtx_d1_mtx_counts <- cbind(d1_T1_mtx_counts,d1_T2_mtx_counts)
cmtx_d3_mtx_counts <- cbind(d3_T1_mtx_counts,d3_T2_mtx_counts)
cmtx_d5_mtx_counts <- cbind(d5_T1_mtx_counts,d5_T2_mtx_counts)
cmtx_d8_mtx_counts <- cbind(d8_T1_mtx_counts,d8_T2_mtx_counts)
cmtx_d11_mtx_counts <- cbind(d11_T1_mtx_counts,d11_T2_mtx_counts)
cmtx_d15_mtx_counts <- d15_mtx_counts #give day 15 matrix same cmtx prefix in variable even though no technical replicates in day 15

#Check dimensions of concatenated matrices
#(*this step is only necessary if Read10X command from Seurat package was not used).

dim(cmtx_d1_mtx_counts)
dim(cmtx_d3_mtx_counts)
dim(cmtx_d5_mtx_counts)
dim(cmtx_d8_mtx_counts)
dim(cmtx_d11_mtx_counts)
dim(cmtx_d15_mtx_counts)

#Convert concatenated matrices into data frame and view data frame of concatenated matrices
#(*this step is only necessary if Read10X command from Seurat package was not used).

DF_cmtx_d1_mtx_counts <- as.data.frame(cmtx_d1_mtx_counts)
DF_cmtx_d3_mtx_counts <- as.data.frame(cmtx_d3_mtx_counts)
DF_cmtx_d5_mtx_counts <- as.data.frame(cmtx_d5_mtx_counts)
DF_cmtx_d8_mtx_counts <- as.data.frame(cmtx_d8_mtx_counts)
DF_cmtx_d11_mtx_counts <- as.data.frame(cmtx_d11_mtx_counts)
DF_cmtx_d15_mtx_counts <- as.data.frame(cmtx_d15_mtx_counts)






