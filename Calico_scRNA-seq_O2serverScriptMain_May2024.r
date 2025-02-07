'''
##create a folder for package installations
mkdir -p ~/R/4.4.0/library
##modify the environment to redirect the installations to the above folder
export R_LIB_USER="~/R/4.4.0/library"
##check the contents of the environment variable R_LIB_USER
echo $R_LIB_USER
'''

#!/usr/bin/env Rscript

#Work started on RStudio (lines 10 - 319) will be used to begin R Script for O2 server. 

# March 2024

# Roux, A. E., Yuan, H., Podshivalova, K., Hendrickson, D., Kerr, R., Kenyon, C., & Kelley, D. (2023).
#Individual cell types in C. elegans age differently and activate distinct cell-protective responses. Cell reports, 42(8), 112902.
#https://doi.org/10.1016/j.celrep.2023.112902

#single-cell RNA-seq from published data

# Reproducing Calico data set

# single-cell RNA-seq analysis - Quality Control (QC)

#Install required packages

#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("SingleCellExperiment")

#install.packages('Seurat'). Note following Seurat package installation: The downloaded source packages are in
 #   ‘/tmp/RtmpHRk1gP/downloaded_packages’

#install.packages("tidyverse")

#install.packages("Matrix")

#install.packages("scales")

#install.packages("cowplot")

#install.packages("RCurl")

#BiocManager::install("AnnotationHub")

#BiocManager::install("ensembldb")

#install.packages("ggplot2_3.5.1") Note: Received this message 'A version of this package for your version of R might be available elsewhere'

#install.packages("devtools_2.4.5") Note: Received this message 'A version of this package for your version of R might be available elsewhere'

#BiocManager::install("SparseArray")

#BiocManager::install("sparseMatrixStats")

#install.packages("SparseM_1.81") Note: Received this message 'A version of this package for your version of R might be available elsewhere'

#install.packages("devtools") Note: This package was installed

#devtools::install_github('satijalab/seurat-data')


# Load libraries
library(SingleCellExperiment)
library(devtools)
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

#Metadata:
#The somatic cells were processed into libraries after flow-sorting. The somatic cells were processed into libraries using 10X Genomics Chromium Single Cell 3' V3 Chemistry
#Libraries were processed on Illumina cBot and HiSeq4000 sequencer
#This study seeks to profile C. elegans cell types during aging
#A time series scRNA-seq is performed and six time points are collected (day 1, day 3, day 5, day 8, day 11 and day 15) over the
#animals adult lifespan
#To enrich somatic cells, germ cell number was reduced by using temperature-sensitive gon-2(q388ts) mutants and by sorting somatic and germ cells based on DNA ploidy in silico
#Ambient RNA corrected for via CellBender and manually
#Final data set contained 47,423 cells and 20,305 genes(2,175 UMI vs. 156 UMI per cell, 644 genes vs. 52 genes per cell which is higher reads per cell compared to previous C. elegans adult scRNA-seq data [see Preston et al.])
#Note: Last time point (day 15) recovered less cells and less UMI per cell
#Raw scRNA-seq data is available in NCBI GEO under ID GSE208154
#Since these are somatic C. elegans cells, the cell types expected are:
#Nerve cells
#Hypodermal cells
#Muscle cells
#Intestinal cells
#glial cells
#epithelial cells
#seam cells
#Since the gon-2 mutation is not completely penetrant, some germline cells are expected:
#vulva

#Unzip .mtx.gz files from data set using a for loop:
#for i in *.gz; do gzip -d $i; done

#Adjusted names of the features, barcodes and matrix files in the folder of each technical replicate so that they are named
#features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz respectively so that they can be read with the Read10X command (do not give the files a prefix to identify the sample, only
#name folders with the sample name, otherwise error message = file missing)

day_1_T1_counts <- Read10X(data.dir = "day_1_T1")
day_1_T2_counts <- Read10X(data.dir = "day_1_T2")
day_3_T1_counts <- Read10X(data.dir = "day_3_T1")
day_3_T2_counts <- Read10X(data.dir = "day_3_T2")
day_5_T1_counts <- Read10X(data.dir = "day_5_T1")
day_5_T2_counts <- Read10X(data.dir = "day_5_T2")
day_8_T1_counts <- Read10X(data.dir = "day_8_T1")
day_8_T2_counts <- Read10X(data.dir = "day_8_T2")
day_11_T1_counts <- Read10X(data.dir = "day_11_T1")
day_11_T2_counts <- Read10X(data.dir = "day_11_T2")
day_15_counts <- Read10X(data.dir = "day_15")

#Concatenate Read10X matrix output to combine technical replicates. Technical replicates have the same row names but different 
#column names and numbers so cbind() command was used

c10mtx_day_1_counts <- cbind(day_1_T1_counts,day_1_T2_counts)
c10mtx_day_3_counts <- cbind(day_3_T1_counts,day_3_T2_counts)
c10mtx_day_5_counts <- cbind(day_5_T1_counts,day_5_T2_counts)
c10mtx_day_8_counts <- cbind(day_8_T1_counts,day_8_T2_counts)
c10mtx_day_11_counts <- cbind(day_11_T1_counts,day_11_T2_counts)
c10mtx_day_15_counts <- day_15_counts #give day 15 matrix same c10mtx prefix in variable for consistency even though no technical replicates in day 15

# Turn count matrix into a Seurat object (output is a Seurat object)

combined_mtx_day_1 <- CreateSeuratObject(count = c10mtx_day_1_counts, min.features = 5)
combined_mtx_day_3 <- CreateSeuratObject(count = c10mtx_day_1_counts, min.features = 5)
combined_mtx_day_5 <- CreateSeuratObject(count = c10mtx_day_1_counts, min.features = 5)
combined_mtx_day_8 <- CreateSeuratObject(count = c10mtx_day_1_counts, min.features = 5)
combined_mtx_day_11 <- CreateSeuratObject(count = c10mtx_day_1_counts, min.features = 5)
combined_mtx_day_15 <- CreateSeuratObject(count = c10mtx_day_1_counts, min.features = 5)


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






