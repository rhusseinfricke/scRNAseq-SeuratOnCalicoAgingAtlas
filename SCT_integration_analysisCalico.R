#Calico single-cell RNA-seq: NORMALIZATION

#Exploring "uninteresting" sources of variation.

#Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

#SIMPLE NORMALIZATION

#This simple normalization explores sources of variation in the data 
#and allows us to check if data correction is needed
#Use filtered Seurat object that was updated with filtered cell data in the QC analysis.

simpleNorm_seuratCalico <- NormalizeData(filtered_seuratObjALLT)

#The most common source of biological ("uninteresting") variation 
#that requires data correction is from the effect of the cell cycle.
#Use function CellCycleScoring() to assign each cell a score based on its 
#expression of G2/S and M markers


#C. elegans cell cycle markers file was generated using human cell cycle marker orthologs and
#contributed to hbc/tinyatlas.

#Download C. elegans cell cycle markers file from hbc/tinyatlas 
#Import cell cycle marker data set from downloaded .csv file

Caenorhabditis_elegans <- read_csv("cell_cycle/Caenorhabditis_elegans.csv")

#Save cell cycle data set as .rda file to data folder

save(Caenorhabditis_elegans, file = "~/Desktop/Calico/Data/Caenorhabditis_elegans.rda")

#Load cell cycle markers

load("~/Desktop/Calico/Data/Caenorhabditis_elegans.rda")

#Use CellCycleScoring() function to assign scores to cells based on the G2/M and S phase markers
#they express

CalicoSeurat_phase <- CellCycleScoring(simpleNorm_seuratCalico, g2m.features = )





