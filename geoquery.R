library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)

data <- getGEO(GEO = "GSE58831", GSEMatrix = T, getGPL = T, AnnotGPL = T)
assay <- data$GSE58831_series_matrix.txt.gz@assayData$exprs
pheno <- data$GSE58831_series_matrix.txt.gz@phenoData@data
feature <- data$GSE58831_series_matrix.txt.gz@featureData@data

#for exchanging prob ids into gene symbols of feature file of GEO
#removing repetitive genes in data
exdat <- assay[!duplicated(feature$`Gene symbol`[match(rownames(assay), rownames(feature))]),]
#exchanging the unique rownames into gene symbols
rownames(exdat) <- feature$`Gene symbol`[!duplicated(feature$`Gene symbol`)]






