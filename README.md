# GEO microarray_analysis

# **GEOquery package**

![](https://img.shields.io/badge/version-2.27.0-green?style=plastic)
![](https://img.shields.io/badge/source-GEO-%23007fff?style=plastic)
![](https://img.shields.io/badge/platform-all-red?style=plastic)
![](https://img.shields.io/badge/data%20availability-NCBI-aqua?style=plastic)

The **NCBI** **Gene Expression Omnibus (GEO)** serves as a public repository for a wide range of high-throughput experimental data. These data include *single* and *dual* channel microarray-based experiments measuring mRNA, genomic DNA, and protein abundance. 

[](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html)

## Getting Started using GEOquery 

The function `getGEO` interprets its input and determines how to get the data from GEO.

#### Installation

Based on *bioconductor* and `BiocManager`

```
{r}
BiocManager::install("GEOquery")
```

```{r}
library(GEOquery)
```

##### GEO accession series

A Series record defines a set of related Samples considered to be part of a group, how the Samples are related, and if and how they are ordered. Each Series record is assigned a unique and stable **GEO accession number (GSE...)**. The matrix file is accessible and downloaded by identifying *GSE number* of a particular datasets. 

* Here I choose a GSE for an example

+ A network connection is required to access the GSE dataset

```{r}
data <- getGEO(GEO = "GSE30812", GSEMatrix = T, getGPL = T, AnnotGPL = T)
```

* *GSEMatrix* argument is defined as providing a matrix includes information of data

* *getGPL* argument is defined as presenting features of samples and clinical information

* *AnnotGPL* argument is defined as providing annotation of probes to know probe IDs and genes


