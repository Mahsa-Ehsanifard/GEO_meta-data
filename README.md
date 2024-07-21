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

```{r}
BiocManager::install("GEOquery")
```

```{r}
library(GEOquery)
```

##### GEO accession series

A Series record defines a set of related Samples considered to be part of a group, how the Samples are related, and if and how they are ordered. Each Series record is assigned a unique and stable **GEO accession number (GSE...)**. The matrix file is accessible and downloaded by identifying *GSE number* of a particular datasets. 

* Here I choose a GSE related to MDS data for an example

+ A network connection is required to access the GSE dataset

```{r}
data <- getGEO(GEO = "GSE58831", GSEMatrix = T, getGPL = T, AnnotGPL = T)
```

* *GSEMatrix* argument is defined as providing a matrix includes information of data

* *getGPL* argument is defined as presenting features of samples and clinical information

* *AnnotGPL* argument is defined as providing annotation of probes to know probe IDs and genes


Now, we can obtain to expression matrix as an *assaydata* inside the *data* 

```{r}
assayData <- data$GSE30812_series_matrix.txt.gz@assayData$exprs
```

```{r}
head (assayData)
```

```
           GSM1420393 GSM1420394 GSM1420395 GSM1420396 GSM1420397 GSM1420398 GSM1420399
1007_s_at   4.440616   4.308647   4.308647   4.307359   4.307359   4.300252   4.308647
1053_at     8.905845   9.073703   9.600805   6.967785   8.379297   8.216474   7.740933
117_at      2.739938   3.449308   2.734385   2.831919   2.983352   3.287753   2.664737
121_at      4.218233   4.218233   4.218233   4.218233   4.218233   4.149675   4.218233
1255_g_at   2.254394   2.254394   2.254394   2.254394   2.254394   2.254394   2.254394
           GSM1420400 GSM1420401 GSM1420402 GSM1420403 GSM1420404 GSM1420405 GSM1420406
1007_s_at   4.335793   4.939411   6.172363   4.308647   4.151949   4.308647   4.230704
1053_at     8.392995   7.511040   7.715622   8.112302   8.577034   8.298663   7.856461
117_at      2.958662   3.348214  11.884425   2.734385   7.764825   2.665857   5.207729
121_at      4.508521   4.218233   4.218233   4.218233   4.218233   4.218233   4.316007
1255_g_at   2.254394   2.254394   2.254394   2.254394   2.254394   2.254394   2.254394
```

* GSM columns are defined as sample IDs

The *phenotypes* or *clinical* information can be released from *phenoData* in *data*.

This command provides all information about each sample, and the types of samples are shown.

```{r}
pheno <- data$GSE30812_series_matrix.txt.gz@phenoData@data
```

There is another column called *featureData* in *data* which is defined as the gene annotation or probe IDs

```{r}
names(feature)
```

```
"ID"                    "Gene title"            "Gene symbol"          
 [4] "Gene ID"               "UniGene title"         "UniGene symbol"       
 [7] "UniGene ID"            "Nucleotide Title"      "GI"                   
[10] "GenBank Accession"     "Platform_CLONEID"      "Platform_ORF"         
[13] "Platform_SPOTID"       "Chromosome location"   "Chromosome annotation"
[16] "GO:Function"           "GO:Process"            "GO:Component"         
[19] "GO:Function ID"        "GO:Process ID"         "GO:Component ID" 
```

## DEG analysis

**Differential Expressed Genes** are employed for detecting hub genes based on significant differential expression between two groups.

* Using `edgeR` and `limma` package, we can analyze differential expressions of *microarray* datasets with **TMM (trimmed mean of M-values)** method of normalization for microarray. Besides, *Adj.P.Value* is added identifying the significance of differential expression, and **logFC (log Fold Change)** will be added to describe the counts of differential exprssion by *positive or negative* values. 

Here, I normalize and filter the raw data first using limma and edgeR package based on value distribution by TMM normalization method.

```{r}
dge <- DGEList(EXdata)
keep <- filterByExpr(dge, design = design)
filt <- dge[keep,,keep.lib.sizes=F]
norm <- calcNormFactors.DGEList(filt, method = "TMM")
```

* `voom` plot is provided for checking normalization and distribution of values using `limma` package

```{r}
v <- voom(norm,design = design, plot = T)
```

* this voom plot shows good distribution of normalized expression values.

* Now, we can analyze DEG result

* *Ex* expression matrix is released from `voom` function in *logarithm* scale normalized values

```
                     GSM1420410        GSM1420411       GSM1420412       GSM1420413       GSM1420414 
 ENSG00000000003 3.49401269701033 5.28960394224448 7.31855304943793 6.79353365184211 4.34146623558799
 ENSG00000000005 2.48196276733284  2.4564806549044 2.47516837620647 2.47052677624909 2.48200884881274
 ENSG00000000419 12.0297749782973 11.6021575137356 11.8016795438964 11.5994084078562 11.9931222883087
```

```{r}
fit <- lmFit(Ex, design = design)
contrast <- makeContrasts(non_V600E-v600e, levels = design)
fit2 <- contrasts.fit(fit = fit, contrasts = contrast)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, number = Inf)
```

```
                  logFC  AveExpr         t      P.Value    adj.P.Val        B 
ENSG00000114302  2.059983 9.016187  9.425437 5.230565e-18 7.353333e-14 30.26257    
ENSG00000135709 -1.866375 2.860313 -9.361123 8.118053e-18 7.353333e-14 29.83572   
ENSG00000106733 -1.558516 2.832475 -9.152146 3.356666e-17 1.843978e-13 28.45755      
ENSG00000157954 -2.057048 2.989810 -9.123564 4.071489e-17 1.843978e-13 28.27013      
ENSG00000164944  1.905614 7.507759  9.077754 5.544899e-17 2.009028e-13 27.97029      
ENSG00000187098 -1.353065 2.716279 -9.004784 9.056078e-17 2.533801e-13 27.49412
```














