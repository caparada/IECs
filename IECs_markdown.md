Intracranial Epidermoid Cysts
================
Carolina Parada
2025-01-09

# About this project

Intracranial Epidermoid Cysts (IECs) are ultra-rare tumors. They usually
have a benign prognosis, however, some IECs can adhere to critical brain
structures resulting in impairment and morbidity.

The primary treatment is surgery. Cyst adherence often complicates
complete removal, resulting in notably high progression rates after
subtotal resection.

Due to the rarity of IECs and small cohort of patients, there has been
limited research focused on understanding the molecular mechanisms of
the disease to advance therapeutic options.

As a result, there are currently no effective drug therapies available
for treating patients with IECs.

Targeted therapy is an effective form of treatment for tumors. It
focuses on mutations that turn healthy cells into tumor cells. The
mutation profile of IECs has not been investigated.

Here we applied Whole Exome Sequencing (WES) to determine the somatic
landscape of six IECs resected during surgery between 1995 and 2021 at
the University of Washington hospitals (Seattle, WA, USA).

The data analysis focused on extracting insights into tumor biology and
identifying mutations that could potentially serve as targets for
drug-therapies.

``` r
# Clean up global environment
rm(list = ls())

# Set up your working dir
setwd("F:\\Epidermoids\\epidermoids\\Git\\")
```

# Load required packages

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(maftools)
library(clusterProfiler)
```

    ## 

    ## clusterProfiler v4.12.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(org.Hs.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
    ##     tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     rename

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     slice

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

``` r
library(ggplot2)
library(ReactomePA)
```

    ## ReactomePA v1.48.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use ReactomePA in published research, please cite:
    ## Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479

``` r
library(STRINGdb)
library(sigminer)
```

    ## Registered S3 method overwritten by 'sigminer':
    ##   method      from
    ##   print.bytes Rcpp

    ## sigminer version 2.3.1
    ## - Star me at https://github.com/ShixiangWang/sigminer
    ## - Run hello() to see usage and citation.

    ## 
    ## Attaching package: 'sigminer'

    ## The following object is masked from 'package:maftools':
    ## 
    ##     MAF

# Load data set

Data: WES of six IECs

This data is publicly available on
<https://www.mdpi.com/2072-6694/16/20/3487>

File name: Table S3_EpidermoidVariantsMultimerge.csv

``` r
data <- read.csv("Review_SupplTableS3_EpidermoidVariantsMultimerge.csv", header = T)
# The data contains 1221 obs (variants) and 61 variables (columns)
```

# Overview of the somatic mutational signature of IECs

Using the package maftools for Exploratory Data Analysis.

``` r
# Inspect variant classification 
unique(data$Variant_Classification)
```

    ## [1] "Missense_Mutation"      "Frame_Shift_Del"        "In_Frame_Del"          
    ## [4] "Nonsense_Mutation"      "In_Frame_Ins"           "Nonstop_Mutation"      
    ## [7] "Frame_Shift_Ins"        "Translation_Start_Site"

``` r
# Converting the CSV file into a MAF object
maf <- read.maf(data,vc_nonSyn = data$Variant_Classification) 
```

    ## -Validating
    ## --Removed 12 duplicated variants
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.060s elapsed (0.060s cpu)

``` r
# Maftools only includes non-synonymous variants for analysis since these
# variants are more likely to have a direct effect on gene function. 
# Here I wanted to include all variant types identified.
# (e.g: Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, etc..)
# Then, the argument vc_nonSyn = data$Variant_Classification was set to include
# all mutation types in addition to non-synonymous mutations.
```

# Summarizing Mutational Data

``` r
# Getting sample summary
SampleSummary <- getSampleSummary(maf)
```

``` r
# Getting gene summary
GeneSummary <- getGeneSummary(maf)
```

``` r
# Getting clinical summary
getClinicalData(maf) # it shows sample names
```

    ##    Tumor_Sample_Barcode
    ##                  <char>
    ## 1:                 1777
    ## 2:                 2939
    ## 3:                 4151
    ## 4:                 4316
    ## 5:                  810
    ## 6:                 H766

``` r
# Getting MAF columns
getFields(maf) # it shows all cols in maf
```

    ##  [1] "Tumor_Sample_Barcode"        "Chromosome"                 
    ##  [3] "Start_Position"              "End_Position"               
    ##  [5] "Reference_Allele"            "Tumor_Seq_Allele2"          
    ##  [7] "Hugo_Symbol"                 "Variant_Classification"     
    ##  [9] "tx"                          "exon"                       
    ## [11] "txChange"                    "aaChange"                   
    ## [13] "Variant_Type"                "Func.refGene"               
    ## [15] "Gene.refGene"                "GeneDetail.refGene"         
    ## [17] "ExonicFunc.refGene"          "AAChange.refGene"           
    ## [19] "cytoBand"                    "SIFT_score"                 
    ## [21] "SIFT_pred"                   "Polyphen2_HDIV_score"       
    ## [23] "Polyphen2_HDIV_pred"         "Polyphen2_HVAR_score"       
    ## [25] "Polyphen2_HVAR_pred"         "LRT_score"                  
    ## [27] "LRT_pred"                    "MutationTaster_score"       
    ## [29] "MutationTaster_pred"         "MutationAssessor_score"     
    ## [31] "MutationAssessor_pred"       "FATHMM_score"               
    ## [33] "FATHMM_pred"                 "PROVEAN_score"              
    ## [35] "PROVEAN_pred"                "VEST3_score"                
    ## [37] "CADD_raw"                    "CADD_phred"                 
    ## [39] "DANN_score"                  "fathmm.MKL_coding_score"    
    ## [41] "fathmm.MKL_coding_pred"      "MetaSVM_score"              
    ## [43] "MetaSVM_pred"                "MetaLR_score"               
    ## [45] "MetaLR_pred"                 "integrated_fitCons_score"   
    ## [47] "integrated_confidence_value" "GERP.._RS"                  
    ## [49] "phyloP7way_vertebrate"       "phyloP20way_mammalian"      
    ## [51] "phastCons7way_vertebrate"    "phastCons20way_mammalian"   
    ## [53] "SiPhy_29way_logOdds"         "cosmic70"                   
    ## [55] "CLINSIG"                     "CLNDBN"                     
    ## [57] "CLNACC"                      "CLNDSDB"                    
    ## [59] "CLNDSDBID"                   "Otherinfo12"                
    ## [61] "Otherinfo13"

``` r
# Writing maf summaries to an output file with basename "IECs" in your 
# working directory
write.mafSummary(maf = maf,basename = "IECs")
```

# Plotting sumaries

The bar plots summarize the variant classification, variant type, SNV
class in the dataset, and the number of variants identified in each
sample.

The median number of variants per sample is displayed.

Most variants identified were missense mutations (1065/1221, 87.22%) and
SNPs (1092/1221, 89.43%).

Predominant SNV alterations included C \> T (n = 408) and T \> C (n =
231).

``` r
plotmafSummary(maf = maf, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE, 
               titvRaw = FALSE,
               fs = 1, # base size for text
               textSize = 0.8, # font size if showBarcodes is TRUE.
               showBarcodes = F, # include sample names in the top bar plot.
               top = 10) # include n top genes. 
```

![](IECs_markdown_files/figure-gfm/plot%20summaries-1.png)<!-- -->

# Pairwise Comparison: identification of mutually exclusive/co-occurring genes

The triangular matrix displays relationships among the top mutated genes
(p \< 0.1).

Green indicates a tendency toward co-occurrence, whereas beige indicates
a tendency toward exclusiveness.

``` r
#exclusive/co-occurrence event analysis on top 50 mutated genes. 
somaticInteractions(maf = maf,
                    pvalue = c(0.05, 0.1), 
                    countsFontSize = 0.8,
                    sigSymbolsFontSize = 0.9,
                    sigSymbolsSize = 1,
                    fontSize = 0.5,
                    showCounts = F,
                    showSum = F,
                    countStats = 'sig',
                    returnAll = T,
                    top=50)
```

![](IECs_markdown_files/figure-gfm/identify%20somatic%20interactions-1.png)<!-- -->

    ##          gene1  gene2     pValue oddsRatio    00    01    11    10      pAdj
    ##         <char> <char>      <num>     <num> <int> <int> <int> <int>     <num>
    ##    1: SLC38A10   GGT2 0.06666667         0     0     4     0     2 0.9393451
    ##    2:   KCNJ12  GSTT4 0.06666667       Inf     2     0     4     0 0.9393451
    ##    3:  KIR2DL1  GSTT4 0.06666667       Inf     2     0     4     0 0.9393451
    ##    4:  KIR2DL3  GSTT4 0.06666667       Inf     2     0     4     0 0.9393451
    ##    5:    OR9G1  GSTT4 0.06666667       Inf     2     0     4     0 0.9393451
    ##   ---                                                                       
    ## 1221:   NBPF26 ZNF717 1.00000000         1     2     2     1     1 1.0000000
    ## 1222:   PABPC1 ZNF717 1.00000000         1     2     2     1     1 1.0000000
    ## 1223: SLC38A10 ZNF717 1.00000000         1     2     2     1     1 1.0000000
    ## 1224:     SVIL ZNF717 1.00000000         1     2     2     1     1 1.0000000
    ## 1225:  TBC1D3D ZNF717 1.00000000         1     2     2     1     1 1.0000000
    ##                    Event             pair event_ratio
    ##                   <char>           <char>      <char>
    ##    1: Mutually_Exclusive   GGT2, SLC38A10         0/6
    ##    2:       Co_Occurence    GSTT4, KCNJ12         4/0
    ##    3:       Co_Occurence   GSTT4, KIR2DL1         4/0
    ##    4:       Co_Occurence   GSTT4, KIR2DL3         4/0
    ##    5:       Co_Occurence     GSTT4, OR9G1         4/0
    ##   ---                                                
    ## 1221: Mutually_Exclusive   NBPF26, ZNF717         1/3
    ## 1222: Mutually_Exclusive   PABPC1, ZNF717         1/3
    ## 1223: Mutually_Exclusive SLC38A10, ZNF717         1/3
    ## 1224: Mutually_Exclusive     SVIL, ZNF717         1/3
    ## 1225: Mutually_Exclusive  TBC1D3D, ZNF717         1/3

# Identifying hypermutated genomic regions

The rainfall plot of log10(inter-event distance) and chromosome number
shows regions where potential changes in inter-event distances are
located, predominantly in sample 2939. These regions are indicated by a
black arrow.

Each point of the rainfall plot is color-coded according to the SNV
class. Changes in inter-event distances (arrow) are located on
chromosomes 7, 13, 17.

``` r
hyperReg = rainfallPlot(maf = maf,
                        detectChangePoints = TRUE,
                        fontSize = 0.8,
                        pointSize = 1)
```

    ## Processing 2939..

    ## Kataegis detected at:

    ##    Chromosome Start_Position End_Position nMuts Avg_intermutation_dist  Size
    ##         <num>          <num>        <num> <int>                  <num> <num>
    ## 1:          7      142772207    142773976     6              353.80000  1769
    ## 2:         13       25096713     25097231     7               86.33333   518
    ## 3:         17       21415555     21416556     8              143.00000  1001
    ##    Tumor_Sample_Barcode   C>A   C>T   T>A   T>C   T>G
    ##                  <fctr> <int> <int> <int> <int> <int>
    ## 1:                 2939     1     3     1     1    NA
    ## 2:                 2939     3     2    NA     1     1
    ## 3:                 2939     1     7    NA    NA    NA

![](IECs_markdown_files/figure-gfm/hypermutated%20regions-1.png)<!-- -->

# Prediction of Disease-Associated Driver Genes.

Prediction was performed by using the algorithm Oncodrive.

The bubble plot of -log10(fdr) and number of variants within clusters
shows predicted driver candidates (marked in red).

The number of closely spaced mutational clusters is highlighted within
brackets.

``` r
maf.sig = oncodrive(maf = maf,
                    AACol = 'aaChange',
                    minMut = 5,
                    pvalMethod = 'zscore')
```

    ## Warning in oncodrive(maf = maf, AACol = "aaChange", minMut = 5, pvalMethod =
    ## "zscore"): Oncodrive has been superseeded by OncodriveCLUSTL. See
    ## http://bg.upf.edu/group/projects/oncodrive-clust.php

    ## No syn mutations found! Skipping background estimation. Using predefined values. (Mean = 0.279; SD = 0.13)

    ## Estimating cluster scores from non-syn variants..

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Comapring with background model and estimating p-values..

    ## Done !

``` r
plotOncodrive(res = maf.sig,
              fdrCutOff = 0.05,
              useFraction = F,
              labelSize = 0.7,
              bubbleSize = 1)
```

![](IECs_markdown_files/figure-gfm/oncodrive-1.png)<!-- -->

# Gene Ontology (GO) over-representation analysis

For GO classification I used the cluterProfiler package, which
implements enrichGO() for GO over-representation.

Over-representation analysis is a statistical method that determines
whether genes from pre-defined sets are over-represented in a subset of
the data.

In this case, the subset is the set of altered genes identified in IECs.

GO was performed with the inclusion of Molecular Function, Biological
Process, and Molecular Component databases.

The bar plots show the top enriched terms.

Bars correspond to terms with significant adjusted p-values \< 0.05.

## Creating input file

``` r
# Creating a gene list and converting gene symbol to EMSEMBL and ENTREZID
gene.df <- bitr(maf@data$Hugo_Symbol,
                fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## Warning in bitr(maf@data$Hugo_Symbol, fromType = "SYMBOL", toType =
    ## c("ENSEMBL", : 1.79% of input gene IDs are fail to map...

## GO Cellular Component (CC)

The GO CC describes the location of a gene product within a cell.

``` r
CC <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'ENTREZID',
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable = T) # print gene names

# The output in a large enrichResult
head(CC, 5)
```

    ##                    ID                                    Description GeneRatio
    ## GO:0042611 GO:0042611                            MHC protein complex     9/672
    ## GO:0098553 GO:0098553 lumenal side of endoplasmic reticulum membrane     9/672
    ## GO:0098576 GO:0098576                       lumenal side of membrane     9/672
    ## GO:0062023 GO:0062023       collagen-containing extracellular matrix    33/672
    ## GO:0042613 GO:0042613                   MHC class II protein complex     6/672
    ##              BgRatio       pvalue     p.adjust       qvalue
    ## GO:0042611  25/19894 6.827230e-08 3.393133e-05 3.205205e-05
    ## GO:0098553  29/19894 2.964621e-07 7.367082e-05 6.959057e-05
    ## GO:0098576  39/19894 4.635574e-06 7.679601e-04 7.254266e-04
    ## GO:0062023 428/19894 1.009420e-05 1.254204e-03 1.184740e-03
    ## GO:0042613  17/19894 1.307226e-05 1.299383e-03 1.227417e-03
    ##                                                                                                                                                                                               geneID
    ## GO:0042611                                                                                                                   HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5
    ## GO:0098553                                                                                                                   HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5
    ## GO:0098576                                                                                                                   HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5
    ## GO:0062023 ACAN/ADAMTSL4/AEBP1/COL11A2/COL4A5/COL5A3/COL6A2/ECM2/EGFL7/FBN2/FGFR2/FLG/FN1/FRAS1/FREM1/HRNR/KRT1/LAMA1/LAMA2/LAMA3/LAMA5/MMP9/P3H2/PLG/PRG4/PRSS1/SSC5D/THBS3/TNC/TPSAB1/VCAN/VTN/ZP3
    ## GO:0042613                                                                                                                                     HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5
    ##            Count
    ## GO:0042611     9
    ## GO:0098553     9
    ## GO:0098576     9
    ## GO:0062023    33
    ## GO:0042613     6

# GO Molecular Function (MF)

GO MF is the activity that a gene product performs at the molecular
level.

``` r
MF <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'ENTREZID',
                    ont = 'MF',
                    readable = T)

# The output in a large enrichResult
head(MF, 5)
```

    ##                    ID                                 Description GeneRatio
    ## GO:0005201 GO:0005201 extracellular matrix structural constituent    20/631
    ## GO:0030246 GO:0030246                        carbohydrate binding    27/631
    ## GO:0032396 GO:0032396    inhibitory MHC class I receptor activity     5/631
    ## GO:0042605 GO:0042605                     peptide antigen binding     9/631
    ## GO:0023023 GO:0023023                 MHC protein complex binding     8/631
    ##              BgRatio       pvalue     p.adjust       qvalue
    ## GO:0005201 166/18522 9.914606e-07 0.0004205789 0.0003921012
    ## GO:0030246 279/18522 1.158620e-06 0.0004205789 0.0003921012
    ## GO:0032396  10/18522 9.873902e-06 0.0023894843 0.0022276909
    ## GO:0042605  47/18522 2.496896e-05 0.0040564753 0.0037818089
    ## GO:0023023  37/18522 2.793716e-05 0.0040564753 0.0037818089
    ##                                                                                                                                                                            geneID
    ## GO:0005201                                                         ACAN/AEBP1/COL11A2/COL4A5/COL5A3/COL6A2/DSPP/FBN2/FBN3/FN1/LAMA1/LAMA2/LAMA3/LAMA5/PRG4/THBS3/TNC/VCAN/VTN/ZP3
    ## GO:0030246 ACAN/ALPK1/CD22/CD33/CEMIP2/CLEC18A/CLEC18B/CLEC18C/CLEC2D/CRYBG2/DBH/FREM1/GAL3ST3/GALNT15/HLA-DRB1/IGF2R/KRT1/MRC1/NOMO1/P3H2/PKD1L2/PRG4/SIGLEC11/SORD/VCAN/VTN/ZP3
    ## GO:0032396                                                                                                                                    KIR3DL1/LILRA2/LILRA4/LILRB1/LILRB2
    ## GO:0042605                                                                                                HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5
    ## GO:0023023                                                                                                    HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5/LILRB1/LILRB2
    ##            Count
    ## GO:0005201    20
    ## GO:0030246    27
    ## GO:0032396     5
    ## GO:0042605     9
    ## GO:0023023     8

# GO Biological Process (BP)

GO BPs are the larger processes or ‘biological programs’ accomplished by
the concerted action of multiple molecular activities.

``` r
BP <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'ENTREZID',
                    ont = 'BP',
                    readable = T)

head(BP, 5)
```

    ##                    ID
    ## GO:0002396 GO:0002396
    ## GO:0002501 GO:0002501
    ## GO:0002399 GO:0002399
    ## GO:0002503 GO:0002503
    ##                                                           Description GeneRatio
    ## GO:0002396                               MHC protein complex assembly     7/634
    ## GO:0002501          peptide antigen assembly with MHC protein complex     7/634
    ## GO:0002399                      MHC class II protein complex assembly     6/634
    ## GO:0002503 peptide antigen assembly with MHC class II protein complex     6/634
    ##             BgRatio       pvalue    p.adjust      qvalue
    ## GO:0002396 21/18888 3.581849e-06 0.008041251 0.007765072
    ## GO:0002501 21/18888 3.581849e-06 0.008041251 0.007765072
    ## GO:0002399 16/18888 8.388920e-06 0.009416563 0.009093148
    ## GO:0002503 16/18888 8.388920e-06 0.009416563 0.009093148
    ##                                                                 geneID Count
    ## GO:0002396 HLA-A/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5     7
    ## GO:0002501 HLA-A/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5     7
    ## GO:0002399       HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5     6
    ## GO:0002503       HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5     6

# Saving all files as CSV files

``` r
# List of dataframes
df_list <- list(CC@result, MF@result, BP@result) 

# Vector with corresponding filenames for each data frame
df_names <- c("CC.csv", "MF.csv", "BP.csv")

# Write each data frame to a separate CSV file
for (i in 1:length(df_list)) {
  write.csv(df_list[[i]], file = df_names[i], row.names = FALSE)
}
```

# Visualizing enrichment results

``` r
## Creating and saving enrichment plots

# List of enrichResults
enrich_list <- list(CC, MF, BP) 

# List of plot filenames
plot_names <- c("CC_barplot.tiff", "MF_barplot.tiff", "BP_barplot.tiff")

# Function to create and save bar plots for each enrichResult
create_enrichment_plots <- function(enrich_list, filename) {
  # Generate the bar plot using clusterProfiler's built-in barplot function
  p <- barplot(enrich_list, showCategory = 5)  # Display top 5 categories
    # Save the plot as a PNG file
  ggsave(filename,
         plot = p,
         width = 8,
         height = 6,
         dpi = 300)
}
  
# Generate and save the plots in batch
plots <- mapply(create_enrichment_plots,
                enrich_list,
                plot_names)

# Outputs will be saved in your working directory
```

## Printing enrichment plots in terminal

``` r
# Function to create and return the bar plots for each enrichResult
print_enrichment_plots <- function(enrich_list) {
  # Generate the bar plot using clusterProfiler's barplot function
  p <- barplot(enrich_list, showCategory = 10)  # Show top 10 categories
  
  # Return the plot object
  return(p)
}

# Generate and store the plots in a list
plots <- lapply(enrich_list, print_enrichment_plots)

# Display each plot in terminal
for (p in plots) {
  print(p)
}
```

![](IECs_markdown_files/figure-gfm/batch%20print_enrichment_plots-1.png)<!-- -->![](IECs_markdown_files/figure-gfm/batch%20print_enrichment_plots-2.png)<!-- -->![](IECs_markdown_files/figure-gfm/batch%20print_enrichment_plots-3.png)<!-- -->

# Pathway Enrichment Analysis (PEA)

Method to identify biological pathways that are more enriched in a gene
list than would be expected by chance.

Helpful to identify potential disease-associated mechanisms.

I explored the KEGG, WikiPathways and Reactome databases.

# KEGG over-representation analysis

``` r
kegg <- enrichKEGG(gene = gene.df$ENTREZID)
```

    ## Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...

    ## Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...

``` r
kegg2 <- setReadable(kegg, 
                     OrgDb = org.Hs.eg.db,
                     keyType='ENTREZID')
head(kegg2, 5)
```

    ##                                      category
    ## hsa04612                   Organismal Systems
    ## hsa05332                       Human Diseases
    ## hsa04514 Environmental Information Processing
    ## hsa04512 Environmental Information Processing
    ## hsa05165                       Human Diseases
    ##                                  subcategory       ID
    ## hsa04612                       Immune system hsa04612
    ## hsa05332                      Immune disease hsa05332
    ## hsa04514 Signaling molecules and interaction hsa04514
    ## hsa04512 Signaling molecules and interaction hsa04512
    ## hsa05165           Infectious disease: viral hsa05165
    ##                                  Description GeneRatio  BgRatio       pvalue
    ## hsa04612 Antigen processing and presentation    19/337  81/8538 2.094300e-10
    ## hsa05332           Graft-versus-host disease    13/337  45/8538 1.043454e-08
    ## hsa04514             Cell adhesion molecules    23/337 158/8538 5.060270e-08
    ## hsa04512            ECM-receptor interaction    16/337  89/8538 3.112514e-07
    ## hsa05165      Human papillomavirus infection    32/337 333/8538 2.504124e-06
    ##              p.adjust       qvalue
    ## hsa04612 6.282900e-08 5.334953e-08
    ## hsa05332 1.565181e-06 1.329031e-06
    ## hsa04514 5.060270e-06 4.296790e-06
    ## hsa04512 2.334386e-05 1.982180e-05
    ## hsa05165 1.387121e-04 1.177836e-04
    ##                                                                                                                                                                                               geneID
    ## hsa04612                                          CD8B/HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5/HSPA6/KIR2DL1/KIR2DL3/KIR2DL4/KIR2DS1/KIR2DS2/KIR2DS4/KIR3DL1/KIR3DL2
    ## hsa05332                                                                                     HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5/KIR2DL1/KIR2DL3/KIR3DL1/KIR3DL2
    ## hsa04514                              CD22/CD8B/CDH15/CDH3/CLDN22/CLDN24/HLA-A/HLA-B/HLA-C/HLA-DPB1/HLA-DQA1/HLA-DQA2/HLA-DQB1/HLA-DRB1/HLA-DRB5/ICAM3/ITGAL/ITGB7/LRRC4B/MADCAM1/PTPRD/SLITRK5/VCAN
    ## hsa04512                                                                                                   COL4A5/COL6A2/DSPP/FN1/FRAS1/FREM1/ITGA7/ITGB7/LAMA1/LAMA2/LAMA3/LAMA5/RELN/THBS3/TNC/VTN
    ## hsa05165 COL4A5/COL6A2/CREB3L1/EGF/EIF4EBP1/FN1/GNAS/HLA-A/HLA-B/HLA-C/IFNA17/IKBKE/ITGA7/ITGB7/JAG1/LAMA1/LAMA2/LAMA3/LAMA5/MAGI1/NOTCH1/NOTCH2/NOTCH3/PATJ/RBL2/RELN/TBP/TCF7L2/THBS3/TNC/TYK2/VTN
    ##          Count
    ## hsa04612    19
    ## hsa05332    13
    ## hsa04514    23
    ## hsa04512    16
    ## hsa05165    32

# Visualizing most significant KEGG enriched terms

``` r
KEGG_dotplot <- dotplot(kegg2, showCategory = 5)
KEGG_dotplot
```

![](IECs_markdown_files/figure-gfm/dotplot%20KEGG-1.png)<!-- -->

# WikiPathways over-representation analysis

``` r
wiki <- enrichWP(gene.df$ENTREZID, organism = "Homo sapiens")
wiki2 <- setReadable(wiki,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID")
head(wiki2, 5)
```

    ##            ID                                               Description
    ## WP3932 WP3932                    Focal adhesion PI3K Akt mTOR signaling
    ## WP4172 WP4172                                        PI3K Akt signaling
    ## WP306   WP306                                            Focal adhesion
    ## WP244   WP244                                  Alpha 6 beta 4 signaling
    ## WP4239 WP4239 Epithelial to mesenchymal transition in colorectal cancer
    ##        GeneRatio  BgRatio       pvalue   p.adjust      qvalue
    ## WP3932    27/318 303/8776 1.336581e-05 0.00721754 0.006781392
    ## WP4172    28/318 341/8776 4.143476e-05 0.01118739 0.010511345
    ## WP306     19/318 199/8776 1.053842e-04 0.01896915 0.017822866
    ## WP244      7/318  33/8776 1.454472e-04 0.01963538 0.018448832
    ## WP4239    16/318 162/8776 2.502206e-04 0.02702382 0.025390806
    ##                                                                                                                                                                  geneID
    ## WP3932      COL11A2/COL5A3/COL6A2/CREB3L1/EGF/EIF4EBP1/FGF14/FGFR2/FN1/HIF1A/IL4R/IL7R/IRS1/ITGA7/ITGAL/ITGB7/KDR/LAMA1/LAMA2/LAMA3/LAMA5/NOS1/RELN/RPTOR/THBS3/TNC/VTN
    ## WP4172 COL4A5/COL6A2/CREB3L1/EGF/EIF4EBP1/FGF14/FGFR2/FLT3/FN1/IFNA17/IL4R/IL7R/IRS1/ITGA7/ITGB7/KDR/LAMA1/LAMA2/LAMA3/LAMA5/MCL1/PIK3AP1/RBL2/RELN/RPTOR/THBS3/TNC/VTN
    ## WP306                                                            COL5A3/COL6A2/DOCK1/EGF/ERBB2/FN1/ITGA7/ITGB7/KDR/LAMA1/LAMA2/LAMA3/LAMA5/PAK3/RELN/SHC1/THBS3/TNC/VTN
    ## WP244                                                                                                                        EIF4EBP1/IRS1/LAMA1/LAMA2/LAMA3/LAMA5/SHC1
    ## WP4239                                                                    CLDN22/CLDN24/COL4A5/FN1/HIF1A/JAG1/JAG2/LRP5/MAP2K3/MMP9/NOTCH1/NOTCH2/NOTCH3/PDCD6/SHC1/VTN
    ##        Count
    ## WP3932    27
    ## WP4172    28
    ## WP306     19
    ## WP244      7
    ## WP4239    16

# Visualizing most significant WikiPathways enriched terms

``` r
wiki_dotplot <- dotplot(wiki2, showCategory = 5)
wiki_dotplot
```

![](IECs_markdown_files/figure-gfm/dotplot%20wiki-1.png)<!-- -->

# Reactome over-representation analysis

``` r
reactome <- enrichPathway(gene=gene.df$ENTREZID,
                          pvalueCutoff = 0.05,
                          readable=TRUE)
head(reactome, 5)
```

    ##                          ID
    ## R-HSA-198933   R-HSA-198933
    ## R-HSA-1474244 R-HSA-1474244
    ## R-HSA-3000178 R-HSA-3000178
    ## R-HSA-1474228 R-HSA-1474228
    ## R-HSA-211935   R-HSA-211935
    ##                                                                            Description
    ## R-HSA-198933  Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell
    ## R-HSA-1474244                                        Extracellular matrix organization
    ## R-HSA-3000178                                                        ECM proteoglycans
    ## R-HSA-1474228                                  Degradation of the extracellular matrix
    ## R-HSA-211935                                                               Fatty acids
    ##               GeneRatio   BgRatio       pvalue     p.adjust       qvalue
    ## R-HSA-198933     25/409 133/11091 1.320532e-11 1.391841e-08 1.365013e-08
    ## R-HSA-1474244    34/409 300/11091 4.877419e-09 2.570400e-06 2.520856e-06
    ## R-HSA-3000178    15/409  76/11091 8.894725e-08 3.125013e-05 3.064779e-05
    ## R-HSA-1474228    20/409 140/11091 1.912708e-07 5.039984e-05 4.942839e-05
    ## R-HSA-211935      6/409  15/11091 9.142635e-06 1.927267e-03 1.890119e-03
    ##                                                                                                                                                                                                          geneID
    ## R-HSA-198933                        CD22/CD33/CD8B/CLEC2D/HLA-A/HLA-B/HLA-C/ICAM3/ITGAL/ITGB7/KIR2DL1/KIR2DL3/KIR2DL4/KIR2DS1/KIR2DS2/KIR3DL1/KIR3DL2/LILRA2/LILRA4/LILRB1/LILRB2/LILRB5/MADCAM1/PILRB/SIGLEC11
    ## R-HSA-1474244 ACAN/CAPN1/CAPN8/CAPN9/COL11A2/COL4A5/COL5A3/COL6A2/DSPP/FBN2/FBN3/FN1/ICAM3/ITGA7/ITGAL/ITGB7/KDR/LAMA1/LAMA2/LAMA3/LAMA5/MADCAM1/MMP9/MUSK/P3H2/PLG/PRSS1/PRSS2/SPOCK3/TLL1/TNC/TPSAB1/VCAN/VTN
    ## R-HSA-3000178                                                                                                                ACAN/COL4A5/COL5A3/COL6A2/DSPP/FN1/ITGA7/LAMA1/LAMA2/LAMA3/LAMA5/MUSK/TNC/VCAN/VTN
    ## R-HSA-1474228                                                                             ACAN/CAPN1/CAPN8/CAPN9/COL11A2/COL4A5/COL5A3/COL6A2/FBN2/FBN3/FN1/LAMA3/LAMA5/MMP9/PLG/PRSS1/PRSS2/SPOCK3/TLL1/TPSAB1
    ## R-HSA-211935                                                                                                                                                         CYP2A7/CYP2D6/CYP2F1/CYP4A11/CYP4B1/CYP4F2
    ##               Count
    ## R-HSA-198933     25
    ## R-HSA-1474244    34
    ## R-HSA-3000178    15
    ## R-HSA-1474228    20
    ## R-HSA-211935      6

## Visualizing most significant Reactome enriched terms

``` r
reactome_dotplot <- dotplot(reactome, showCategory = 5)
reactome_dotplot
```

![](IECs_markdown_files/figure-gfm/dotplot%20reactome-1.png)<!-- -->

## Saving all PEA results as CSV

``` r
# List of PEA dfs
PEA_list <- list(kegg2@result, wiki2@result, reactome@result) 

# Vector with corresponding filenames for each data frame
PEA_filenames <- c("kegg.csv", "wikipathway.csv", "reactome.csv")

# Write each data frame to a separate CSV file
for (i in 1:length(PEA_list)) {
  write.csv(PEA_list[[i]], file = PEA_filenames[i], row.names = FALSE)
}

# Outputs will be saved in your working directory
```

## Saving all dotplots of PEA enrichment analysis

``` r
# List of PEA dotplots
PEA_dotplot <- list(KEGG_dotplot, wiki_dotplot, reactome_dotplot)

# Vector with corresponding plot names
dotplot_names <- c("KEGG_dotplot.tiff",
                   "wikipathway_dotplot.tiff",
                   "reactome_dotplot.tiff")

# Save each plot separately
for (i in seq_along(PEA_dotplot)) {
  filename <- paste0(dotplot_names, i, ".tiff") # Create filename for each plot
  ggsave(filename[[i]],
         plot = PEA_dotplot[[i]],
         width = 8,
         height = 6,
         dpi = 300)  # Save plot to file
}
```

# Gene-Concept Network Analysis

The dotplots only displayed most significant selected enriched terms.

Now, I wanted to visualize which genes are involved in these significant
terms. In order to consider the potentially biological complexities in
which a gene may belong to multiple annotation categories and provide
information of numeric changes if available, we applied the cnetplot()
function to extract the complex association.

The cnetplot() function of clusterProfiler depicts the linkages of genes
and biological concepts (e.g. GO terms or KEGG pathways) as a network.

Here, I applied the cnetplot() to visualize color-coded pathway
categories and the genes associated with each significant term.

``` r
# KEGG cnetplot
kegg_cnetplot <- cnetplot(kegg2,
                       foldChange = NULL,
                       circular = T,
                       layout = "kk",
                       colorEdge = T,
                       node_label = "all",
                       cex_category = 0.5,
                       cex_gene = 0.7,
                       cex_label_gene = 0.9)
```

    ## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(foldChange = your_value)' instead of 'foldChange'.
    ##  The foldChange parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(edge = your_value)' instead of 'colorEdge'.
    ##  The colorEdge parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(category_node = your_value)' instead of 'cex_category'.
    ##  The cex_category parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(gene_node = your_value)' instead of 'cex_gene'.
    ##  The cex_gene parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(gene_label = your_value)' instead of 'cex_label_gene'.
    ##  The cex_label_gene parameter will be removed in the next version.

``` r
kegg_cnetplot
```

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](IECs_markdown_files/figure-gfm/KEGG%20gene%20concept%20network-1.png)<!-- -->

``` r
# KEGG cnetplot
wiki_cnetplot <- cnetplot(wiki2,
                       foldChange = NULL,
                       circular = T,
                       layout = "kk",
                       colorEdge = T,
                       node_label = "all",
                       cex_category = 0.5,
                       cex_gene = 0.7,
                       cex_label_gene = 0.9)
```

    ## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(foldChange = your_value)' instead of 'foldChange'.
    ##  The foldChange parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(edge = your_value)' instead of 'colorEdge'.
    ##  The colorEdge parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(category_node = your_value)' instead of 'cex_category'.
    ##  The cex_category parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(gene_node = your_value)' instead of 'cex_gene'.
    ##  The cex_gene parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(gene_label = your_value)' instead of 'cex_label_gene'.
    ##  The cex_label_gene parameter will be removed in the next version.

``` r
wiki_cnetplot
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](IECs_markdown_files/figure-gfm/wikipathway%20gene%20concept%20network-1.png)<!-- -->

``` r
# KEGG cnetplot
reactome_cnetplot <- cnetplot(reactome,
                       foldChange = NULL,
                       circular = T,
                       layout = "kk",
                       colorEdge = T,
                       node_label = "all",
                       cex_category = 0.5,
                       cex_gene = 0.7,
                       cex_label_gene = 0.9)
```

    ## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(foldChange = your_value)' instead of 'foldChange'.
    ##  The foldChange parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(edge = your_value)' instead of 'colorEdge'.
    ##  The colorEdge parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(category_node = your_value)' instead of 'cex_category'.
    ##  The cex_category parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(gene_node = your_value)' instead of 'cex_gene'.
    ##  The cex_gene parameter will be removed in the next version.

    ## Warning in cnetplot.enrichResult(x, ...): Use 'cex.params = list(gene_label = your_value)' instead of 'cex_label_gene'.
    ##  The cex_label_gene parameter will be removed in the next version.

``` r
reactome_cnetplot
```

    ## Warning: ggrepel: 20 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](IECs_markdown_files/figure-gfm/reactome%20gene%20concept%20network-1.png)<!-- -->

# Saving all cnetplots

``` r
# List of PEA dotplots
cnetplots <- list(kegg_cnetplot, wiki_cnetplot, reactome_cnetplot)

# Vector with corresponding plot names
cnetplot_names <- c("KEGG_cnetplot",
                    "wikipathway_cnetplot",
                    "reactome_cnetplot")

# Saving each plot separately
for (i in seq_along(cnetplots)) {
  filename <- paste0(cnetplot_names, ".tiff") # Create a filename for each plot
  ggsave(filename[[i]], 
         plot = cnetplots[[i]],
         width = 10, 
         height = 10, 
         dpi = 300)  # Save plot to file
}
```

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 4 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
# Outputs saved in working directory
```

GO and PEA analysis demonstrated that many altered genes in IECs are
closely associated with immunity and Extracellular Matrix (ECM),
suggesting the involvement of the immune system and microenvironment as
a mechanism of cyst formation.

To better visualize alterations in the immune repertoire, I assessed
genetic changes in immune-related genes using an oncostrip plot.

Remarkably, somatic mutations in these genes were present in 100% of the
cohort.

The most frequently altered immune-related genes in IECs were Killer
cell immunoglobulin-like receptors 2DL1 (KIR2DL1) and 2DL3 (KIR2DL3),
which harbored multiple missense and multi-hit mutations, affecting
approximately 70% of IECs.

Additionally, I observed DNA variations in other Killer cell
immunoglobulin-like receptors (KIR2DS1/2/4, KIR3DL1/2, KIR2DL4), HLA
genes (HLA-A/B/C, HLA-DRB1/5, HLA-DQB1, HLA-DPB1, HLA-DQA1/2), and other
immune-associated genes, albeit at lower frequencies within the study
cohort.

These alterations in the immune repertoire associated with IECs strongly
suggest mechanisms of tumor immune evasion.

# Oncostrip plot of over-represented gene sets

# Immune-associated genes: highlighted by the KEGG database

``` r
# Gene list: immune genes
immune_genes <- c('CD8B', 'CRLF1','FLT3','IL4R','IL6ST','IL7R', 'HLA-A',
                  'HLA-B', 'HLA-C', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2',
                  'HLA-DQB1', 'HLA-DRB1','HLA-DRB5', 'HLA-DQB1',
                  'HLA-DRB1', 'HLA-DRB5', 'HSPA6', 'KIR2DL1', 'KIR2DL3',
                  'KIR2DL4', 'KIR2DS1', 'KIR2DS2', 'KIR2DS4', 'KIR3DL1',
                  'KIR3DL2', 'LILRA2','LILRA4','LILRB1','LILRB2', 'IFNA17',
                  'IL4R','IL7R', 'JAG1','JAG2','NOTCH1','NOTCH2','NOTCH3',
                  'TYK2', 'CD22','CD33','CD5','CD8B', 'MADCAM1', 'ITGB7',
                  'SHC1','ULBP2', 'ITGAL')


oncostrip(maf = maf, 
          showTumorSampleBarcodes = T,
          SampleNamefontSize =1.0,
          fontSize = 0.45,
          genes = immune_genes)
```

![](IECs_markdown_files/figure-gfm/Immune-associated%20genes-1.png)<!-- -->
\# Focal Adhesion + PI3K-AKT-mTOR: highlighted by The WikiPathway
database.

Next, I assessed the mutation burden in genes associated with
PI3K-Akt-mTOR signaling and ECM-mediated cellular interactions as
enriched by WikiPathways and Reactome.

Somatic mutations on these genes also affected 100% of the cohort.

``` r
# Gene list: focal adhesion + PI3K-AKT-mTOR
PI3K <-c("EGF", "FLT3", "IL4R", "LAMA1", "RELN", "COL4A5", "COL6A2", 
"CREB3L1", "EIF4EBP1", "FGF14", "FGFR2", "FN1", "IFNA17", "IL7R",
"IRS1", "ITGA7", "ITGB7", "KDR", "LAMA2", "LAMA3", "LAMA5",
"MCL1", "PIK3AP1", "RBL2", "RPTOR", "THBS3", "TNC", "VTN",
"COL11A2", "COL5A3", "HIF1A", "ITGAL", "NOS1", "DOCK1", "ERBB2",
"PAK3", "SHC1")

oncostrip(maf = maf, 
          showTumorSampleBarcodes = T,
          SampleNamefontSize =1.0,
          fontSize = 0.45,
          genes = PI3K)
```

![](IECs_markdown_files/figure-gfm/focal%20adhesion%20+%20PI3K-AKT-mTOR-1.png)<!-- -->

## Extracellular Matrix (ECM): highlighted by the Reactome database.

``` r
# Gene list
ECM <- c('ACAN', 'CAPN1', 'CAPN8', 'CAPN9', 'COL11A2', 'COL4A5', 'COL5A3',
         'COL6A2', 'DSPP', 'FN1', 'FBN2', 'FBN3','FLG', 'FN1', 'ICAM3',
         'ITGA7', 'ITGAL', 'ITGB7','MADCAM1', 'MMP9', 'LAMA1', 'LAMA2',
         'LAMA3', 'LAMA5', 'PRG4', 'P3H2', 'PLG', 'PRSS1', 'SSC5D',
         'PRSS2', 'RELN', 'TPSAB1', 'THBS3', 'TNC', 'TLL1', 'VCAN', 'VTN',
         'MMP9', 'KDR', 'SPOCK3', 'MUSK')

oncostrip(maf = maf, 
          showTumorSampleBarcodes = T,
          SampleNamefontSize =1.0,
          fontSize = 0.5,
          genes = ECM)
```

![](IECs_markdown_files/figure-gfm/ecm-1.png)<!-- -->

# Identifying the most altered genes in the cohort

``` r
# Counting the occurrences of mutated genes
gene_counts <- data %>%
  group_by(Hugo_Symbol) %>%
  mutate(gene_counts = n_distinct(Tumor_Sample_Barcode)) %>%
  ungroup()

head(gene_counts)
```

    ## # A tibble: 6 × 62
    ##   Tumor_Sample_Barcode Chromosome Start_Position End_Position Reference_Allele
    ##   <chr>                <chr>               <int>        <int> <chr>           
    ## 1 2939                 chr17            81122300     81122300 C               
    ## 2 1777                 chr17            69149050     69149051 GA              
    ## 3 4316                 chr15            88856755     88856755 C               
    ## 4 4316                 chr7            150220267    150220267 T               
    ## 5 2939                 chr14            70457980     70457980 A               
    ## 6 810                  chr15           100341846    100341848 CAG             
    ## # ℹ 57 more variables: Tumor_Seq_Allele2 <chr>, Hugo_Symbol <chr>,
    ## #   Variant_Classification <chr>, tx <chr>, exon <chr>, txChange <chr>,
    ## #   aaChange <chr>, Variant_Type <chr>, Func.refGene <chr>, Gene.refGene <chr>,
    ## #   GeneDetail.refGene <chr>, ExonicFunc.refGene <chr>, AAChange.refGene <chr>,
    ## #   cytoBand <chr>, SIFT_score <chr>, SIFT_pred <chr>,
    ## #   Polyphen2_HDIV_score <chr>, Polyphen2_HDIV_pred <chr>,
    ## #   Polyphen2_HVAR_score <chr>, Polyphen2_HVAR_pred <chr>, LRT_score <chr>, …

``` r
# Filtering genes that are mutated in >= 50% of the samples
filtered_top_genes <- gene_counts %>%
  filter(gene_counts > 2)

# Save df with top altered genes
write.csv(filtered_top_genes, "top_genes_50.csv", row.names = F)

# Optional: visualize it
vc_nonSyn <- unique(filtered_top_genes$Variant_Classification)
oncostrip(read.maf(filtered_top_genes, vc_nonSyn = vc_nonSyn))
```

    ## -Validating
    ## --Removed 2 duplicated variants
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.030s elapsed (0.020s cpu)

![](IECs_markdown_files/figure-gfm/top%20altered%20genes-1.png)<!-- -->

# Identifying the most recurrent mutations in the cohort

``` r
# Counting the occurrences of mutations 
aa_change_counts <- data %>%
  group_by(AAChange.refGene) %>%
  mutate(count = n_distinct(Tumor_Sample_Barcode)) %>%
  ungroup()

head(aa_change_counts)
```

    ## # A tibble: 6 × 62
    ##   Tumor_Sample_Barcode Chromosome Start_Position End_Position Reference_Allele
    ##   <chr>                <chr>               <int>        <int> <chr>           
    ## 1 2939                 chr17            81122300     81122300 C               
    ## 2 1777                 chr17            69149050     69149051 GA              
    ## 3 4316                 chr15            88856755     88856755 C               
    ## 4 4316                 chr7            150220267    150220267 T               
    ## 5 2939                 chr14            70457980     70457980 A               
    ## 6 810                  chr15           100341846    100341848 CAG             
    ## # ℹ 57 more variables: Tumor_Seq_Allele2 <chr>, Hugo_Symbol <chr>,
    ## #   Variant_Classification <chr>, tx <chr>, exon <chr>, txChange <chr>,
    ## #   aaChange <chr>, Variant_Type <chr>, Func.refGene <chr>, Gene.refGene <chr>,
    ## #   GeneDetail.refGene <chr>, ExonicFunc.refGene <chr>, AAChange.refGene <chr>,
    ## #   cytoBand <chr>, SIFT_score <chr>, SIFT_pred <chr>,
    ## #   Polyphen2_HDIV_score <chr>, Polyphen2_HDIV_pred <chr>,
    ## #   Polyphen2_HVAR_score <chr>, Polyphen2_HVAR_pred <chr>, LRT_score <chr>, …

``` r
# Filtering mutations that occur in >= 50% of samples
filtered_aa_changes <- aa_change_counts %>%
  filter(count > 2) 

# Save df with top recurrrent mutations
write.csv(filtered_aa_changes, "top_mutations_50.csv", row.names = F)

# Optional: Visualizing it
oncostrip(read.maf(filtered_aa_changes))
```

    ## -Validating
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.030s elapsed (0.030s cpu)

![](IECs_markdown_files/figure-gfm/filtering%20top%20recurrent%20mutations-1.png)<!-- -->

``` r
# Filtering mutations that occur in 2 or more samples AND
# are predicted to be deleterious according to the CADD score
recurrent_deleterious <- aa_change_counts %>%
  filter(count > 1 & CADD_phred >= 10)

# Visualizing it in a oncostrip plot
recurrent_deleterious %>% 
  read.maf() %>% # converting to MAF
  oncostrip() # plotting oncostrip
```

    ## -Validating
    ## --Removed 2 duplicated variants
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.030s elapsed (0.040s cpu)

![](IECs_markdown_files/figure-gfm/recurrent%20and%20deleterious%20mutations-1.png)<!-- -->

# Protein-Protein Interaction (PPI)

Protein-protein interaction between altered genes was assessed by using
the STRING database.

The functional protein association network analysis revealed three main
interaction clusters associated with the MHC protein complex, and
antigen processing and presentation.

``` r
# Loading String Database
string_db <- STRINGdb$new(version = "12.0",
                          species = 9606,
                          score_threshold = 200)

# Converting filtered_top_genes into a data frame
STRING_df <- data.frame(filtered_top_genes)

# Mapping gene names to STRING database identifiers
STRING_df_mapped <- string_db$map(STRING_df,
                                  "Hugo_Symbol",
                                  removeUnmappedRows = TRUE )
```

    ## Warning:  we couldn't map to STRING 9% of your identifiers

``` r
# Visualization of string interaction network
hits <- STRING_df_mapped$STRING_id[1:2000]
string_db$plot_network(hits)
```

![](IECs_markdown_files/figure-gfm/visualization%20of%20string%20interaction%20network-1.png)<!-- -->

# Lolliplots

I also utilized lollipop plots to visualize recurrent and predicted
deleterious variants on protein domains.

``` r
# The input here is a list of frequently altered genes that have been associated
# to cancer and display CADD deleterious scores

# Convert recurrent_deleterious into MAF
recurrent_deleterious_maf <- read.maf(recurrent_deleterious)
```

    ## -Validating
    ## --Removed 2 duplicated variants
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.010s elapsed (0.030s cpu)

``` r
# Genes of interest
# My list of genes of interest includes:
# tumor suppressor or oncogenes that are frequently mutated with CADD 
# deleterious variants, and previously 
genes_of_interest <- c("USP8", "HLA-DRB1", "PDPR", "UBXN11", "NOTCH2", "PABPC3")

# Plotting lolliplots iterating over the list of genes of interest
for (gene in genes_of_interest) {
  plot <- lollipopPlot(maf = recurrent_deleterious_maf,
               gene = gene,
               AACol = 'aaChange',
               labelPos = "all",
               labPosAngle = 90,
               labPosSize = 1.0,
               collapsePosLabel = T,
               domainLabelSize = 0.6)
  print(plot)
}
```

    ## 3 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##      HGNC    refseq.ID   protein.ID aa.length
    ##    <char>       <char>       <char>     <num>
    ## 1:   USP8 NM_001128610 NP_001122082      1118
    ## 2:   USP8 NM_001128611 NP_001122083      1118
    ## 3:   USP8    NM_005154    NP_005145      1118

    ## Using longer transcript NM_001128610 for now.

    ## NULL

    ## 2 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##        HGNC    refseq.ID   protein.ID aa.length
    ##      <char>       <char>       <char>     <num>
    ## 1: HLA-DRB1 NM_001243965 NP_001230894       266
    ## 2: HLA-DRB1    NM_002124    NP_002115       266

    ## Using longer transcript NM_001243965 for now.

![](IECs_markdown_files/figure-gfm/lolliplots-1.png)<!-- -->![](IECs_markdown_files/figure-gfm/lolliplots-2.png)<!-- -->

    ## NULL

    ## NULL

    ## 3 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##      HGNC    refseq.ID   protein.ID aa.length
    ##    <char>       <char>       <char>     <num>
    ## 1: UBXN11 NM_001077262 NP_001070730       400
    ## 2: UBXN11    NM_145345    NP_663320       487
    ## 3: UBXN11    NM_183008    NP_892120       520

    ## Using longer transcript NM_183008 for now.

![](IECs_markdown_files/figure-gfm/lolliplots-3.png)<!-- -->

    ## NULL

    ## 2 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##      HGNC    refseq.ID   protein.ID aa.length
    ##    <char>       <char>       <char>     <num>
    ## 1: NOTCH2 NM_001200001 NP_001186930      1235
    ## 2: NOTCH2    NM_024408    NP_077719      2471

    ## Using longer transcript NM_024408 for now.

![](IECs_markdown_files/figure-gfm/lolliplots-4.png)<!-- -->![](IECs_markdown_files/figure-gfm/lolliplots-5.png)<!-- -->

    ## NULL

![](IECs_markdown_files/figure-gfm/lolliplots-6.png)<!-- -->

    ## NULL

# Analysis of somatic copy-number alterations

Copy number inference was performed by using CNVKit (v0.9.9).

The merged .cns output files is publicly available on
<https://www.mdpi.com/2072-6694/16/20/3487>

File name: Review_SupplTableS10_CNVkit.csv

# CNVkit Analysis

``` r
# Loading the merged CNVkit .cns files.
CNV_data <- read.csv("Review_SupplTableS10_CNVkit.csv", header = T)
```

# Identifying COSMIC signatures from records of variant calling data.

The package signminer was used to identify COSMIC somatic signatures and
also generate multiple visualizations of copy-number data.

The COSMIC signatures include three type of signatures: SBS, DBS and ID
(INDEL).

The signature identification procedure has been divided into 4 steps:

Step 1: Read mutation data. Step 2: Tally components: for SBS, it means
classifying SBS records into 96 components (the most common case) and
generate sample matrix. Step 3: Extract signatures: estimate signature
number and identify signatures. Step 4: Visualize signatures

# Step 1: Setting up sex and preparing input

## Setting up sex

``` r
## For cohort contains both males and females,
## set a data.frame with two columns: female and male
sex_df = data.frame(sample = c("H766", "4151",
                               "1777", "2939",
                               "810", "4316"), 
                               sex = c("female", "female",
                                       "male", "female",
                                       "female", "female"))

# The optionsigminer.sex is used to control the processing of sex. 
# If you don’t care the sex chromosomes (X and Y), you can ignore this
# setting after removing the X/Y segments, 
# otherwise the summary in the result "cn" and tally process (which will be
# carried out next) may be biased.
options(sigminer.sex = sex_df)
```

## Preparing input (loading the CNV data)

The input is a data frame with the result of the CNVkit (merged .cns
files) and requires the following information:

- Segment chromosome.
- Segment start.
- Segment end.
- Absolute copy number value for this segment: must be integer.
- Sample ID.

The read_copynumber() function loads the input file.

``` r
# Load copy number data containing cn integers
cn <- read_copynumber(CNV_data,
                      seg_cols = c("chromosome", "start", "end", "p_ttest"),
                      genome_build = "hg38",
                      complement = FALSE,
                      verbose = TRUE)
```

    ## ℹ [2025-01-09 12:02:48.625254]: Started.

    ## ℹ [2025-01-09 12:02:48.78979]: Genome build  : hg38.

    ## ℹ [2025-01-09 12:02:48.806908]: Genome measure: called.

    ## ✔ [2025-01-09 12:02:48.826624]: Chromosome size database for build obtained.

    ## ℹ [2025-01-09 12:02:48.841727]: Reading input.

    ## ✔ [2025-01-09 12:02:48.855879]: A data frame as input detected.

    ## ✔ [2025-01-09 12:02:48.870661]: Column names checked.

    ## ✔ [2025-01-09 12:02:48.885814]: Column order set.

    ## ✔ [2025-01-09 12:02:48.914064]: Chromosomes unified.

    ## ✔ [2025-01-09 12:02:48.940259]: Data imported.

    ## ℹ [2025-01-09 12:02:48.955411]: Segments info:

    ## ℹ [2025-01-09 12:02:48.970274]:     Keep - 3391

    ## ℹ [2025-01-09 12:02:48.984643]:     Drop - 0

    ## ✔ [2025-01-09 12:02:49.000425]: Segments sorted.

    ## ℹ [2025-01-09 12:02:49.015151]: Joining adjacent segments with same copy number value. Be patient...

    ## ✔ [2025-01-09 12:02:49.149321]: 584 segments left after joining.

    ## ✔ [2025-01-09 12:02:49.164421]: Segmental table cleaned.

    ## ℹ [2025-01-09 12:02:49.178168]: Annotating.

    ## ✔ [2025-01-09 12:02:49.203612]: Annotation done.

    ## ℹ [2025-01-09 12:02:49.218324]: Summarizing per sample.

    ## ✔ [2025-01-09 12:02:49.24404]: Summarized.

    ## ℹ [2025-01-09 12:02:49.259028]: Generating CopyNumber object.

    ## ✔ [2025-01-09 12:02:49.273833]: Generated.

    ## ℹ [2025-01-09 12:02:49.2877]: Validating object.

    ## ✔ [2025-01-09 12:02:49.302192]: Done.

    ## ℹ [2025-01-09 12:02:49.316231]: 0.691 secs elapsed.

# Step 2: Tally components

Tally components: for SBS, it means classifying SBS records into 96
components (the most common case) and generating a sample matrix.

``` r
# Option sigminer.copynumber.max is used to control the processing of max copy
# number values. 
options(sigminer.copynumber.max = 20)
## Even you set max_copynumber = 20 in read_copynumber(),
## the segmental copy number may be greater than 20
## because for male samples, the X/Y segmental copy number
## values will be doubled in tally process.
## This setting will make copy number values of all segments
## not greater than 20.

# Using the method designed by Wang, Shixiang et al.
cn_tally <- sig_tally(cn, method = "W")
```

    ## ℹ [2025-01-09 12:02:49.345326]: Started.

    ## ℹ [2025-01-09 12:02:49.361913]: Step: getting copy number features.

    ## 
    ## Attaching package: 'purrr'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     reduce

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     simplify

    ## ℹ [2025-01-09 12:02:49.449316]: Getting breakpoint count per 10 Mb...

    ## ℹ [2025-01-09 12:02:49.565083]: Getting breakpoint count per chromosome arm...

    ## ℹ [2025-01-09 12:02:49.642606]: Getting copy number...

    ## ℹ [2025-01-09 12:02:49.666424]: Getting change-point copy number change...

    ## ℹ [2025-01-09 12:02:49.741945]: Getting length of chains of oscillating copy number...

    ## ℹ [2025-01-09 12:02:49.817288]: Getting (log10 based) segment size...

    ## ℹ [2025-01-09 12:02:49.835718]: Getting the minimal number of chromosome with 50% CNV...

    ## ℹ [2025-01-09 12:02:49.887806]: Getting burden of chromosome...

    ## ✔ [2025-01-09 12:02:49.927652]: Gotten.

    ## ℹ [2025-01-09 12:02:49.943023]: Step: generating copy number components.

    ## ✔ [2025-01-09 12:02:49.956859]: `feature_setting` checked.

    ## ℹ [2025-01-09 12:02:49.97257]: Step: counting components.

    ## ✔ [2025-01-09 12:02:50.203751]: Counted.

    ## ℹ [2025-01-09 12:02:50.218921]: Step: generating components by sample matrix.

    ## ✔ [2025-01-09 12:02:50.233448]: Matrix generated.

    ## ℹ [2025-01-09 12:02:50.248072]: 0.903 secs elapsed.

``` r
# Of note, the sigminer.copynumber.max option only has effect on sig_tally()
# with method “W,” 
# the sigminer.sex option has effects on read_copynumber() and sig_tally() 
# with method “W.”
```

# Step 3: Extract signatures

Once I generated the matrix, I stimated signature number and
identified  
signatures.

``` r
# Extract signatures
cn_tally$nmf_matrix[1:5, 1:5]
```

    ##      BP10MB[0] BP10MB[1] BP10MB[2] BP10MB[3] BP10MB[4]
    ## 1777       302        13         6         1         0
    ## 2939       290        17        11         3         0
    ## 4151       307         5         6         4         0
    ## 4316       272        19        27         1         3
    ## 810        292        15        12         1         2

``` r
# library(NMF)
sig_w <- sig_extract(cn_tally$nmf_matrix, n_sig = 2)
```

    ## NMF algorithm: 'brunet'

    ## Multiple runs: 10

    ## Mode: sequential [foreach:doParallelSNOW]

    ## Runs: |                                                        Runs: |                                                  |   0%Runs: |                                                        Runs: |==================================================| 100%
    ## System time:
    ##    user  system elapsed 
    ##    1.23    0.06    2.59

``` r
# This step returns a list containing information about copy number features, 
# components and matrix, NMF, etc.
```

# Step 4: Visualizing signatures

``` r
# Copy Number Signature Profile
show_sig_profile(sig_w,
                 mode = "copynumber",
                 normalize = "feature",
                 method = "W",
                 style = "cosmic",
                 font_scale = 0.8)
```

![](IECs_markdown_files/figure-gfm/visualizing%20signatues-1.png)<!-- -->

# Additional plots

Histogram of copy-number distribution by segment length of segment and
percentage

``` r
# Plotting a histogram of copy-number distribution by segment length of copy 
# number alterations and percentage (as fraction of all copy number alterations). 
show_cn_distribution(cn, mode = "ld")
```

    ## Warning: The dot-dot notation (`..density..`) was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `after_stat(density)` instead.
    ## ℹ The deprecated feature was likely used in the sigminer package.
    ##   Please report the issue at <https://github.com/ShixiangWang/sigminer/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](IECs_markdown_files/figure-gfm/visualize%20percentage%20of%20copy%20number%20distribution%20per%20length-1.png)<!-- -->

``` r
# Bar plot of copy number distribution per chromosome
show_cn_distribution(cn, mode = 'cd')
```

![](IECs_markdown_files/figure-gfm/visualize%20copy%20number%20distribution%20per%20chromosome-1.png)<!-- -->

``` r
# Circos plot of the copy number variation frequency profile in each sample. 
layout(matrix(1:6, 2, 3))
show_cn_circos(cn, samples = 6)
```

![](IECs_markdown_files/figure-gfm/circus%20plot-1.png)<!-- -->

# GISTIC Analysis

To identify somatic copy number regions that are significantly amplified
or deleted across all six samples, the log2 ratios and segments
generated with CNVkit were used as inputs for GISTIC.

I summarized output files generated by GISTIC using the maftools
package.

For this analysis, you need four GISTIC outputs:

all_lesions.conf_XX.txt,

amp_genes.conf_XX.txt,

del_genes.conf_XX.txt and

scores.gistic,

where XX is the confidence level. See GISTIC documentation for details.

The cohort analysis of the distributions of amplification and deletion
events revealed seven cytobands significantly altered in the cohort,
among them, amplifications of 7q35, 4p16.1, 17q12 (AP_5 and AP_6),
8p23.1, and 14q11.2, and the deletion of 2p11.1.

``` r
# folder with gistic outputs
dir <- "R:\\epidermoid_project\\CNVkit\\CNVfiles\\Gistic\\594894\\594894\\"

# Reading GISTIC directory
IEC.gistic <- readGistic(gisticDir = dir, isTCGA = TRUE)
```

    ## -Processing Gistic files..
    ## --Processing amp_genes.conf_90.txt
    ## --Processing del_genes.conf_90.txt
    ## --Processing scores.gistic
    ## --Summarizing by samples

``` r
IEC.gistic
```

    ## An object of class  GISTIC 
    ##           ID summary
    ##       <char>   <num>
    ## 1:   Samples       6
    ## 2:    nGenes      75
    ## 3: cytoBands       7
    ## 4:       Amp     213
    ## 5:       Del      12
    ## 6:     total     225

For copy-number exploratory data analysis, I called the methods
available to access slots of GISTIC object: getSampleSummary,
getGeneSummary and getCytoBandSummary.

``` r
# Sample summary
samples_sum <- getSampleSummary(IEC.gistic)
```

``` r
# Gene summary
getGeneSummary(IEC.gistic)
```

    ##      Hugo_Symbol   Amp   Del total AlteredSamples
    ##           <char> <int> <int> <num>          <int>
    ##  1:     ACTR3BP2     0     4     4              4
    ##  2:    ARHGEF34P     4     0     4              4
    ##  3:     ARHGEF35     4     0     4              4
    ##  4:      ARHGEF5     4     0     4              4
    ##  5:       CTAGE4     4     0     4              4
    ##  6:       CTAGE8     4     0     4              4
    ##  7:        GGT8P     0     4     4              4
    ##  8: LOC101060389     4     0     4              2
    ##  9: LOC101928605     4     0     4              4
    ## 10:    LOC440434     4     0     4              2
    ## 11:    LOC654342     0     4     4              4
    ## 12:        OR2A1     4     0     4              4
    ## 13:    OR2A1-AS1     4     0     4              4
    ## 14:      OR2A20P     4     0     4              4
    ## 15:        OR2A7     4     0     4              4
    ## 16:       OR2A9P     4     0     4              4
    ## 17:      TBC1D3G     4     0     4              2
    ## 18:      TBC1D3H     4     0     4              2
    ## 19:      TBC1D3K     4     0     4              2
    ## 20:      TBC1D3L     4     0     4              2
    ## 21:     DEFB103A     3     0     3              3
    ## 22:     DEFB104A     3     0     3              3
    ## 23:     DEFB105A     3     0     3              3
    ## 24:     DEFB106A     3     0     3              3
    ## 25:     DEFB107A     3     0     3              3
    ## 26:      DEFB131     3     0     3              3
    ## 27:       DEFB4B     3     0     3              3
    ## 28:       FAM66B     3     0     3              3
    ## 29:    FAM90A10P     3     0     3              3
    ## 30:     FAM90A7P     3     0     3              3
    ## 31:      PRR23D1     3     0     3              3
    ## 32:      PRR23D2     3     0     3              3
    ## 33:      SPAG11A     3     0     3              3
    ## 34:      SPAG11B     3     0     3              3
    ## 35:      USP17L1     3     0     3              3
    ## 36:     USP17L10     3     0     3              3
    ## 37:     USP17L11     3     0     3              3
    ## 38:     USP17L12     3     0     3              3
    ## 39:     USP17L13     3     0     3              3
    ## 40:     USP17L15     3     0     3              3
    ## 41:     USP17L17     3     0     3              3
    ## 42:     USP17L18     3     0     3              3
    ## 43:     USP17L19     3     0     3              3
    ## 44:     USP17L20     3     0     3              3
    ## 45:     USP17L24     3     0     3              3
    ## 46:     USP17L25     3     0     3              3
    ## 47:     USP17L26     3     0     3              3
    ## 48:     USP17L27     3     0     3              3
    ## 49:     USP17L28     3     0     3              3
    ## 50:     USP17L29     3     0     3              3
    ## 51:     USP17L30     3     0     3              3
    ## 52:      USP17L4     3     0     3              3
    ## 53:     USP17L6P     3     0     3              3
    ## 54:     USP17L9P     3     0     3              3
    ## 55:      ZNF705G     3     0     3              3
    ## 56:      BMS1P17     2     0     2              2
    ## 57:      BMS1P22     2     0     2              2
    ## 58:       CCL3L1     2     0     2              2
    ## 59:       CCL3L3     2     0     2              2
    ## 60:       CCL4L1     2     0     2              2
    ## 61:       CCL4L2     2     0     2              2
    ## 62:      DUXAP10     2     0     2              2
    ## 63:    LINC01296     2     0     2              2
    ## 64: LOC100508046     2     0     2              2
    ## 65: LOC101929572     2     0     2              2
    ## 66:       OR11H2     2     0     2              2
    ## 67:        POTEG     2     0     2              2
    ## 68:    POTEH-AS1     2     0     2              2
    ## 69:        POTEM     2     0     2              2
    ## 70:       TBC1D3     2     0     2              2
    ## 71:      TBC1D3B     2     0     2              2
    ## 72:      TBC1D3C     2     0     2              2
    ## 73:      TBC1D3E     2     0     2              2
    ## 74:      TBC1D3F     2     0     2              2
    ## 75:      TBC1D3I     2     0     2              2
    ##      Hugo_Symbol   Amp   Del total AlteredSamples

``` r
# Cytoband summary
getCytobandSummary(IEC.gistic)
```

    ##     Unique_Name nGenes nSamples Variant_Classification Cytoband
    ##          <char>  <int>    <int>                 <char>   <char>
    ## 1:    AP_2:7q35     11        4                    Amp     7q35
    ## 2:  AP_1:4p16.1     19        3                    Amp   4p16.1
    ## 3:   AP_5:17q12     11        2                    Amp    17q12
    ## 4:   AP_6:17q12     11        2                    Amp    17q12
    ## 5:  AP_3:8p23.1     16        3                    Amp   8p23.1
    ## 6: AP_4:14q11.2     10        2                    Amp  14q11.2
    ## 7:  DP_2:2p11.1      3        4                    Del   2p11.1
    ##            Wide_Peak_Limits    qvalues
    ##                      <char>      <num>
    ## 1: chr7:144184283-144388755 1.1384e-06
    ## 2:     chr4:9210658-9450513 1.3311e-06
    ## 3:  chr17:36155548-36397351 2.6841e-04
    ## 4:  chr17:37915697-38260294 1.2234e-03
    ## 5:     chr8:7322050-7858796 1.2838e-01
    ## 6:  chr14:18970677-19723907 1.2838e-01
    ## 7:   chr2:89637308-94581097 1.5796e-01

``` r
# Saving summarized results 
write.GisticSummary(IEC.gistic, basename = "GISTIC_IECs")
```

Summaries are publicly available on
<https://www.mdpi.com/2072-6694/16/20/3487>

File name: Review_SupplTableS11_GisticSummaries.xlsx

``` r
# Gistic chromPlot of top 5 alterations with the lowest q-values.
gisticChromPlot(gistic = IEC.gistic,
                ref.build = "hg38",
                color = c("red", "blue"),
                cytobandTxtSize = 0.8,
                markBands = "all")
```

![](IECs_markdown_files/figure-gfm/GISTIC%20chromPlot-1.png)<!-- -->

``` r
# bubble plot with the number of samples affected by frequently occurring 
# copy-number events.
gisticBubblePlot(gistic = IEC.gistic, color = "red", txtSize = 2)
```

![](IECs_markdown_files/figure-gfm/GISTIC%20bubble%20plot-1.png)<!-- -->

``` r
# Setting up the color labels
CNVcolors=c("red", "blue")
names(CNVcolors)=c("Amp","Del")

# GISTIC oncoplot showing the samples affected by frequently occurring 
# copy-number events. 
gisticOncoPlot(gistic = IEC.gistic,
               SampleNamefontSize = 0.6,
               fontSize = 0.8,
               showTumorSampleBarcodes = T,
               legendFontSize = 1.2,
               colors = CNVcolors)
```

![](IECs_markdown_files/figure-gfm/GISTIC%20oncoplot-1.png)<!-- -->

The cohort analysis of the distributions of amplification and deletion
events revealed seven cytobands significantly altered in the cohort,
among them, amplifications of 7q35, 4p16.1, 17q12 (AP_5 and AP_6),
8p23.1, and 14q11.2, and the deletion of 2p11.1. Of these, the events
with the lowest False Discovery Rate qvalues calculated for aberrant
regions included the amplification of the cytobands 7q35 (AP_5 and
AP_6), 4p16.1, and 17q12, affecting 67% (n = 4/6), 50% (n = 3/6), and
33% (n = 2/6) of the cohort, respectively.

Interestingly, the pattern of alterations in deubiquitinases and
immune-related genes is also present at copy number levels (Output from
GISTIC). Recurrent amplifications were predominantly noted in
DUB/ubiquitin-specific protease 17 (USP17) family, TBC1 domain protein
family, and beta-defensins genes (DEFB4B, DEFB103A, DEFB104A, DEFB105A,
DEFB106A, DEFB107A, and DEFB131) within the cohort, which overexpression
plays a major role in the immune response to cancer.

# Druggable categories

To generate hypotheses about how the potential driver candidates of IECs
might be targeted therapeutically or prioritized for drug development, I
matched the list of the most altered tumor suppressors and oncogenes
(gene driver candidates) against the Drug-Gene Interaction database
(DGIdb), a compendium of drug-gene interactions and potentially
druggable genes, to prioritize drug-gene interactions.

``` r
dgi = drugInteractions(maf = read.maf(filtered_top_genes), # input is a MAF file
                       fontSize = 1.2,
                       plotType = 'bar',
                       top = 100)
```

    ## -Validating
    ## --Removed 2 duplicated variants
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.030s elapsed (0.030s cpu)

![](IECs_markdown_files/figure-gfm/druggable%20categories-1.png)<!-- -->

Above plot shows potential druggable gene categories along with top
genes involved in them. The categorical analysis revealed the PTEN
family which includes the TPTE gene.

Variants in TPTE (p.L332P, p.R57Q, p.R84del, p.R84W) were observed in
70% of the cohort. Given that the pathway enrichment analysis indicated
the significant enrichment by PI3K-Akt-mTOR and ECM pathways, current
PI3K, AKT, mTORC1/mTORC2, and PDK1 inhibitors may be interesting drug
candidates for therapeutic intervention of IECs.

I also extracted information on drug-gene interactions, as showed below:

``` r
# Obtaining claimed drugs for given genes 
dgi.genes <- drugInteractions(genes = filtered_top_genes$Hugo_Symbol,
                              drugs = T)
```

    ## Number of claimed drugs for given genes:
    ##        Gene     N
    ##      <char> <int>
    ## 1:      SMO    44
    ## 2:    HLA-B     9
    ## 3:   NOTCH2     6
    ## 4:   KCNJ12     3
    ## 5:  KIR2DL3     3
    ## 6: HLA-DRB1     2
    ## 7:  KIR2DL1     2
    ## 8: HLA-DRB5     1
    ## 9:    PRSS2     1

``` r
# Matrix with known/reported drugs to interact with the given genes
df.dgi.genes <- dgi.genes[,.(Gene,
                             interaction_types,
                             drug_name,
                             drug_claim_name,
                             drug_chembl_id)]
head(df.dgi.genes)
```

    ##      Gene interaction_types           drug_name drug_claim_name drug_chembl_id
    ##    <char>            <char>              <char>          <char>         <char>
    ## 1:  PRSS2                           CHEMBL27885     DERMOLASTIN    CHEMBL27885
    ## 2:    SMO        antagonist           SARIDEGIB            8198   CHEMBL538867
    ## 3:  HLA-B                              ABACAVIR        Abacavir     CHEMBL1380
    ## 4:    SMO         inhibitor SONIDEGIB PHOSPHATE   CHEMBL3137317  CHEMBL3137317
    ## 5:    SMO         inhibitor          VISMODEGIB    CHEMBL473417   CHEMBL473417
    ## 6:    SMO         inhibitor          VISMODEGIB      VISMODEGIB   CHEMBL473417

This drug-gene interaction data frame indicated inhibitors and antibody
candidates for potential targeted therapies. Among these, the NOTCH2
inhibitor Nirogacestat (PF-03084014) is an FDA-approved chemotherapeutic
agent marketed as Ogsiveo® (SpringWorks Therapeutics). It is indicated
for adult patients with progressing desmoid tumors requiring systemic
treatment.

Additionally, NOTCH2 inhibitors RO4929097 and MK0752 are currently in
trial for pancreatic cancer. NOTCH2 antibodies OMP-59R5 (Tarextumab),
and REGN-421 are under investigation for solid tumors.

# Oncogenic mechanisms driving intracranial epidermoid cysts

The previous pathway enrichment analysis identified Focal Adhesion (FA),
PI3K-Akt-mTOR, immune response cascades, and Extracellular Matrix (ECM)
remodeling as the most significantly altered cascades.

I grouped the altered gene sets that are representative of these
pathways, and after combining them, I explored the KEGG database to
uncover potential relationships among them.

The results highlight multiple mechanisms and pathway cross talk,
providing insight into the mechanisms underlying the oncogenesis of
IECs.

``` r
# Representative gene sets of the most affected pathways:
FA_PI3K <- c("COL11A2", "COL5A3","COL6A2", "CREB3L1", "EGF",
             "EIF4EBP1" ,"FGF14", "FGFR2", "FN1", "HIF1A",
             "IL4R", "IL7R", "IRS1", "ITGA7", "ITGAL",
             "ITGB7", "KDR", "LAMA1", "LAMA2", "LAMA3",
             "LAMA5", "NOS1", "RELN", "RPTOR", "THBS3",
             "TNC", "VTN")

PI3k <- c("COL4A5", "COL6A2", "CREB3L1", "EGF", "EIF4EBP1",
          "FGF14", "FGFR2", "FLT3", "FN1", "IFNA17",
          "IL4R", "IL7R", "IRS1", "ITGA7", "ITGB7",
          "KDR", "LAMA1", "LAMA2", "LAMA3", "LAMA5",
          "MCL1", "PIK3AP1", "RBL2", "RELN", "RPTOR", 
          "THBS3", "TNC", "VTN")

FA <- c("COL5A3", "COL6A2", "DOCK1", "EGF", "ERBB2",
        "FN1", "ITGA7", "ITGB7", "KDR", "LAMA1",
        "LAMA2", "LAMA3", "LAMA5", "PAK3", "RELN",
        "SHC1", "THBS3", "TNC", "VTN")
          
immune <- c("KIR2DL1", "KIR2DL3", "HLA-B", "HLA-DRB1", "HLA-DRB5",
            "KIR2DS1", "KIR2DS2", "KIR3DL1", "LILRB1", "NOTCH2",
            "FLT3", "IL4R", "HLA-A", "HLA-DQB1", "KIR2DL4",
            "LILRB2","MADCAM1", "CD8B", "CRLF1", "IL6ST",
            "IL7R", "HLA-C", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2",
            "HSPA6", "KIR2DS4", "KIR3DL2", "LILRA2", "LILRA4",
            "IFNA17", "JAG1", "JAG2", "NOTCH1", "NOTCH3",
            "TYK2", "CD22", "CD33", "CD5", "ITGB7",
            "SHC1", "ULBP2", "ITGAL")

ECM <- c("ACAN", "CAPN1", "CAPN8", "CAPN9", "COL11A2",
         "COL4A5", "COL5A3", "COL6A2", "DSPP", "FBN2",
         "FBN3", "FN1", "ICAM3", "ITGA7", "ITGAL",
         "ITGB7", "KDR", "LAMA1", "LAMA2", "LAMA3", 
         "LAMA5", "MADCAM1", "MMP9", "MUSK", "P3H2",
         "PLG", "PRSS1", "PRSS2", "SPOCK3", "TLL1",
         "TNC", "TPSAB1", "VCAN", "VTN")

# Combining all gene sets and adding ENTREZIDs
hh = bitr(c(immune, FA, PI3K, FA_PI3K, ECM),
          fromType="SYMBOL",
          toType="ENTREZID",
          OrgDb="org.Hs.eg.db")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
# KEGG pathway enrichment of the gene list
hh2 <- enrichKEGG(gene = hh$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

# Plotting KEGG pathway
browseKEGG(hh2, 'hsa05165')
```

The present work has been published on Cancer (Basels) MDPI - Oct/2024.

The publication is available on
<https://www.mdpi.com/2072-6694/16/20/3487>

*AUTHOR: Carolina Parada*

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] doParallel_1.0.17      iterators_1.0.14       foreach_1.5.2         
    ##  [4] NMF_0.27               cluster_2.1.6          rngtools_1.5.2        
    ##  [7] registry_0.5-1         purrr_1.0.2            sigminer_2.3.1        
    ## [10] STRINGdb_2.16.4        ReactomePA_1.48.0      ggplot2_3.5.1         
    ## [13] org.Hs.eg.db_3.19.1    AnnotationDbi_1.66.0   IRanges_2.38.1        
    ## [16] S4Vectors_0.42.1       Biobase_2.64.0         BiocGenerics_0.50.0   
    ## [19] clusterProfiler_4.12.0 maftools_2.20.0        dplyr_1.1.4           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      shape_1.4.6.1           rstudioapi_0.16.0      
    ##   [4] jsonlite_1.8.8          magrittr_2.0.3          farver_2.1.2           
    ##   [7] rmarkdown_2.27          GlobalOptions_0.1.2     ragg_1.3.2             
    ##  [10] fs_1.6.4                zlibbioc_1.50.0         vctrs_0.6.5            
    ##  [13] memoise_2.0.1           ggtree_3.12.0           htmltools_0.5.8.1      
    ##  [16] curl_5.2.1              plotrix_3.8-4           gridGraphics_0.5-1     
    ##  [19] parallelly_1.38.0       KernSmooth_2.23-24      gsubfn_0.7             
    ##  [22] plyr_1.8.9              cachem_1.1.0            igraph_2.0.3           
    ##  [25] lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.7-0           
    ##  [28] R6_2.5.1                fastmap_1.2.0           gson_0.1.0             
    ##  [31] future_1.34.0           GenomeInfoDbData_1.2.12 digest_0.6.36          
    ##  [34] aplot_0.2.3             enrichplot_1.24.0       colorspace_2.1-0       
    ##  [37] furrr_0.3.1             patchwork_1.2.0         chron_2.3-61           
    ##  [40] textshaping_0.4.0       RSQLite_2.3.7           labeling_0.4.3         
    ##  [43] fansi_1.0.6             httr_1.4.7              polyclip_1.10-6        
    ##  [46] compiler_4.4.1          bit64_4.0.5             withr_3.0.0            
    ##  [49] graphite_1.50.0         BiocParallel_1.38.0     viridis_0.6.5          
    ##  [52] DBI_1.2.3               highr_0.11              R.utils_2.12.3         
    ##  [55] ggforce_0.4.2           gplots_3.2.0            MASS_7.3-61            
    ##  [58] rappdirs_0.3.3          HDO.db_0.99.1           gtools_3.9.5           
    ##  [61] DNAcopy_1.78.0          caTools_1.18.3          tools_4.4.1            
    ##  [64] ape_5.8                 scatterpie_0.2.3        R.oo_1.26.0            
    ##  [67] glue_1.7.0              nlme_3.1-165            GOSemSim_2.30.0        
    ##  [70] grid_4.4.1              shadowtext_0.1.4        gridBase_0.4-7         
    ##  [73] reshape2_1.4.4          fgsea_1.30.0            generics_0.1.3         
    ##  [76] gtable_0.3.5            R.methodsS3_1.8.2       tidyr_1.3.1            
    ##  [79] data.table_1.15.4       tidygraph_1.3.1         utf8_1.2.4             
    ##  [82] XVector_0.44.0          ggrepel_0.9.5           pillar_1.9.0           
    ##  [85] stringr_1.5.1           yulab.utils_0.1.5       circlize_0.4.16        
    ##  [88] splines_4.4.1           tweenr_2.0.3            treeio_1.28.0          
    ##  [91] lattice_0.22-6          survival_3.7-0          bit_4.0.5              
    ##  [94] tidyselect_1.2.1        GO.db_3.19.1            Biostrings_2.72.1      
    ##  [97] knitr_1.48              reactome.db_1.88.0      gridExtra_2.3          
    ## [100] sqldf_0.4-11            xfun_0.45               graphlayouts_1.1.1     
    ## [103] proto_1.0.0             stringi_1.8.4           UCSC.utils_1.0.0       
    ## [106] lazyeval_0.2.2          ggfun_0.1.5             yaml_2.3.9             
    ## [109] evaluate_0.24.0         codetools_0.2-20        ggraph_2.2.1           
    ## [112] tibble_3.2.1            qvalue_2.36.0           hash_2.2.6.3           
    ## [115] graph_1.82.0            ggplotify_0.1.2         cli_3.6.3              
    ## [118] systemfonts_1.1.0       munsell_0.5.1           Rcpp_1.0.12            
    ## [121] GenomeInfoDb_1.40.1     globals_0.16.3          png_0.1-8              
    ## [124] blob_1.2.4              DOSE_3.30.1             bitops_1.0-7           
    ## [127] listenv_0.9.1           viridisLite_0.4.2       tidytree_0.4.6         
    ## [130] scales_1.3.0            crayon_1.5.3            rlang_1.1.4            
    ## [133] cowplot_1.1.3           fastmatch_1.1-4         KEGGREST_1.44.1
