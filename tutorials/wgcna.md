wgcna
================
Jennifer Chang
11/24/2020

<!--
---
title: "WGCNA Gene Correlation Network Analysis"
layout: single
author: Jennifer Chang
author1: Siva Chudalayandi
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---
-->

**Last Update**: 14 Dec 2020 <br/> **R Markdown**:
[WGCNA.Rmd](https://bioinformaticsworkbook.org/tutorials/WGCNA.Rmd)

# Network analysis with WGCNA

There are many gene correlation network builders but we shall provide an
example of the WGCNA R Package.

The **WGCNA R package** builds “weighted gene correlation networks for
analysis” from expression data. It was originally published in 2008 and
cited as the following:

-   Langfelder, P. and Horvath, S., 2008. [WGCNA: an R package for
    weighted correlation network
    analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559).
    BMC bioinformatics, 9(1), p.559.
-   Zhang, B. and Horvath, S., 2005. [A general framework for weighted
    gene co-expression network
    analysis](https://pubmed.ncbi.nlm.nih.gov/16646834/). Statistical
    applications in genetics and molecular biology, 4(1).

**More information**

-   [Recent PubMed
    Papers](https://pubmed.ncbi.nlm.nih.gov/?term=wgcna&sort=date)
-   [Original WGCNA tutorials - Last updated
    2016](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)
-   [Video: ISCB Workshop 2016 - Co-expression network analysis using
    RNA-Seq data (Keith Hughitt)](https://youtu.be/OdqDE5EJSlA)

## Installing WGCNA

We will assume you have a working R environment. If not, please see the
following tutorial:

-   [Seting up an R and RStudio
    Environment](https://bioinformaticsworkbook.org/dataWrangling/R/r-setup.html)

Since WGCNA is an R package, we will need to start an R environment and
install from R’s package manager, CRAN.

``` r
install.packages("WGCNA")   # WGCNA is available on CRAN
library(WGCNA)
```

## Overview

The **WGCNA** pipeline is expecting an input matrix of RNA Sequence
counts. Usually we need to rotate (transpose) the input data so `rows` =
`treatments` and `columns` = `gene probes`.

![WGCNA Overview](Assets/wgcna_overview.png)

The output of **WGCNA** is a list of clustered genes, and weighted gene
correlation network files.

# Example Dataset

We shall start with an example dataset about ER (Endocrine Reticulum)
cell death response. \[Add description of data and maybe link to paper
here\]

-   [All\_Counts\_ER.txt](data/All_Counts_ER.txt)

``` bash
wget https://bioinformaticsworkbook.org/tutorials/data/All_Counts_ER.txt
```

## Load R Libraries

This analysis requires the following R libraries. You might need to
install the library if it’s not already on your system.

We’re going to conform to the [tidyverse]() ecosystem. For a discussion
on its benefits see [“Welcome to the Tidyverse” (Wickham et al,
2019)](https://tidyverse.tidyverse.org/articles/paper.html). This allows
us to organize the pipeline in the following framework ([Fig. from “R
for Data Science” (Wickham and Grolemund,
2017)](https://r4ds.had.co.nz/)):

<img src="https://tidyverse.tidyverse.org/articles/data-science.png" width="400" />

``` r
# Uncomment and modify the following to install any missing packages
# install.packages(c("tidyverse", "magrittr", "WGCNA))
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)        
```

## Tidy the Dataset, and using exploratory graphics

Load and look at the data

``` r
# ==== Load and clean data
data <- readr::read_delim("data/All_Counts_ER.txt", delim="\t")
#> 
#> ── Column specification ────────────────────────────────────────────────────────
#> cols(
#>   .default = col_double(),
#>   Geneid = col_character()
#> )
#> ℹ Use `spec()` for the full column specifications.

data[1:5,1:10]        # Look at first 5 rows and 10 columns
#> # A tibble: 5 x 10
#>   Geneid `0I_S3_L006` `0I_S3_L007` `0II_S11_L006` `0II_S11_L007` `0III_S19_L006`
#>   <chr>         <dbl>        <dbl>          <dbl>          <dbl>           <dbl>
#> 1 Zm000…          103          135             98            103             187
#> 2 Zm000…         1280         1370            876            970            1602
#> 3 Zm000…            0            2              0              0               0
#> 4 Zm000…            0            6              3              2              18
#> 5 Zm000…            0            0              0              0               0
#> # … with 4 more variables: `0III_S19_L007` <dbl>, `3I_S4_L006` <dbl>,
#> #   `3I_S4_L007` <dbl>, `3II_S12_L006` <dbl>

# str(data)           # str = structure of data, useful for debugging data type mismatch errors

names(data)           # Look at the column names
#>  [1] "Geneid"           "0I_S3_L006"       "0I_S3_L007"       "0II_S11_L006"    
#>  [5] "0II_S11_L007"     "0III_S19_L006"    "0III_S19_L007"    "3I_S4_L006"      
#>  [9] "3I_S4_L007"       "3II_S12_L006"     "3II_S12_L007"     "3III_S20_L006"   
#> [13] "3III_S20_L007"    "6I_S5_L006"       "6I_S5_L007"       "6II_S13_L006"    
#> [17] "6II_S13_L007"     "6III_S21_L006"    "6III_S21_L007"    "12I_S6_L006"     
#> [21] "12I_S6_L007"      "12II_S14_L006"    "12II_S14_L007"    "12III_S22_L006"  
#> [25] "12III_S22_L007"   "24I_S7_L006"      "24I_S7_L007"      "24II_S15_L006"   
#> [29] "24II_S15_L007"    "24III_S23_L006"   "24III_S23_L007"   "36I_S8_L006"     
#> [33] "36I_S8_L007"      "36II_S16_L006"    "36II_S16_L007"    "36III_S24_L006"  
#> [37] "36III_S24_L007"   "48I_S9_L006"      "48I_S9_L007"      "48II_S17_L006"   
#> [41] "48II_S17_L007"    "48III_S25_L006"   "48III_S25_L007"   "48mkI_S27_L007"  
#> [45] "48Imk_S10_L006"   "48mkII_S28_L007"  "48IImk_S18_L006"  "48mkIII_S29_L007"
#> [49] "48IIImk_S26_L006"
```

If you are in RStudio, you can also click on the `data` object in the
Environment tab to see an Excel-like view of the data.

![RStudio View](Assets/RStudio_data.png)

From looking at the data, we can come to the following insights:

-   We see that `rows` = `gene probes` which probably means `columns` =
    `treatments` which is the opposite of what’s needed in WGCNA (`rows`
    = `treatment`, `columns` = `gene probes`). This dataset will need to
    be rotated (transposed) before sending to WGCNA.
-   This is also wide data, we will convert this to tidy data before
    visualization. For [Hadley Wickham](http://hadley.nz/)’s tutorial on
    the what and why to convert to tidy data, see
    <https://r4ds.had.co.nz/tidy-data.html>.
-   The column names are prefixed with the hour of treatment
    (e.g. `0I_S3_L006` is 0 hours, `3I_S4_L006` is 3 hours,
    `48mkII_S28_L007` is 48 hours mock.)

The following R commands clean and tidy the dataset for exploratory
graphics.

``` r
col_sel = names(data)[-1]                                  # Get all but first column name
mdata <- tidyr::pivot_longer( data,                        # Convert to Tidy Data
                              col = all_of(col_sel)
                             ) %>%  
  mutate(
    group = gsub("I.*", "", name) %>% gsub("_.*", "", .),  # Add a "Group" column for the hour
  )

mdata$group[grepl("mk", mdata$name)] = "48mk"              # Deal with columns where it's 48mk or 48IIImk,

# This sets the order of the hours in the plot... otherwise 48mk will be between "3" and "6".
mdata$group = factor(mdata$group,
                     levels = c("0", "3", "6", "12", "24","36","48","48mk"))
```

Think through what kinds of plots may tell you something about the
dataset. This example plots the data to identify any outliers.

``` r
# ==== Plot groups (Sample Groups vs RNA Seq Counts) to identify outliers
p <- mdata %>%
    ggplot(., aes(x=name, y=value)) +     # x = treatment, y = RNA Seq count
    geom_violin() +                       # violin plot, show distribution
    geom_point(alpha=0.2) +               # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=90)                          # Rotate treatment text
    ) +
    facet_grid(cols = vars(group), drop=TRUE, scales="free_x")      # Facet by hour

p + labs(x="Treatment Groups", y = "RNA Seq Counts")
```

![](Assets/wgcna_group_by_hour-1.png)<!-- -->

From here, we can see there’s something strange in some of the hour 24
samples. One has very high RNASeq values `24II_S15_L006` with maybe a
wide range, while another has very low range of RNASeq values
`24_S15_L007`. We should follow up with the wet lab folks on an
explanation of those samples, but for now, we’ll remove the 24 hour
group and maybe the 48 hour group.

``` r
keep_cols = names(data) %>% grep("24", .,  invert = T, value = T) %>% grep("48I+_", ., invert=T, value=T)
cdata = data %>% select(all_of(keep_cols))

temp <- cdata[rowSums(cdata[,-1]) > 0.1, ]      # Remove genes with all 0 values
row_var <- apply(temp[,-1], 1, var)             # Remove genes with variance below 10
cdata <- temp[row_var > 10, ]
#cdata[1:5, 1:10]
```

You can look at the `cdata` object (click on item in `environment` or
use `names(cdata)`) to convince yourself that the “24 hour” group is
gone. The original dataset had 46,430 genes (too many to explore),
subsetting by variance and other strange artifacts reduced it down to
28,128 genes. Let’s continue and determine the correlation networks for
these 28,128 genes.

## WGCNA

Now let’s transpose the data and prepare the dataset for WGCNA.

``` r
input_mat = t(as.matrix(cdata[,-1]))
colnames(input_mat) = cdata$Geneid

input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns
#>               Zm00001d027230 Zm00001d027231 Zm00001d027233 Zm00001d027236
#> 0I_S3_L006               103           1280              0              6
#> 0I_S3_L007               135           1370              6             25
#> 0II_S11_L006              98            876              3              3
#> 0II_S11_L007             103            970              2             28
#> 0III_S19_L006            187           1602             18              6
#>               Zm00001d027239 Zm00001d027242 Zm00001d027248 Zm00001d027250
#> 0I_S3_L006               562            432              8              4
#> 0I_S3_L007               411            548             27             18
#> 0II_S11_L006             434            337             16              4
#> 0II_S11_L007             351            404              3              9
#> 0III_S19_L006           1999           1770             98              2
#>               Zm00001d027254 Zm00001d027256
#> 0I_S3_L006                 6             19
#> 0I_S3_L007                 5             14
#> 0II_S11_L006               2             14
#> 0II_S11_L007               4              0
#> 0III_S19_L006             16             16
```

We can see now that the `rows` = `treatments` and `columns` =
`gene probes`. We’re ready to start WGCNA. A correlation network will be
a complete network (all genes are connected to all other genes). Ergo we
will need to pick a threshhold value (if correlation is below threshold,
remove the edge). We assume the true biological network follows a
scale-free structure (see papers by [Albert
Barabasi](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model)).

To do that, WGCNA will try a range of soft thresholds and create a
diagnostic plot. This step will take several minutes so feel free to run
and get coffee.

``` r
#library(WGCNA)
allowWGCNAThreads()          # allow multi-threading (optional)
#> Allowing multi-threading with up to 4 threads.

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(input_mat,             # <= Input data
                        #blockSize = 30,
                        powerVector = powers,
                        verbose = 5
                        )
#> pickSoftThreshold: will use block size 1590.
#>  pickSoftThreshold: calculating connectivity for given powers...
#>    ..working on genes 1 through 1590 of 28128
#>    ..working on genes 1591 through 3180 of 28128
#>    ..working on genes 3181 through 4770 of 28128
#>    ..working on genes 4771 through 6360 of 28128
#>    ..working on genes 6361 through 7950 of 28128
#>    ..working on genes 7951 through 9540 of 28128
#>    ..working on genes 9541 through 11130 of 28128
#>    ..working on genes 11131 through 12720 of 28128
#>    ..working on genes 12721 through 14310 of 28128
#>    ..working on genes 14311 through 15900 of 28128
#>    ..working on genes 15901 through 17490 of 28128
#>    ..working on genes 17491 through 19080 of 28128
#>    ..working on genes 19081 through 20670 of 28128
#>    ..working on genes 20671 through 22260 of 28128
#>    ..working on genes 22261 through 23850 of 28128
#>    ..working on genes 23851 through 25440 of 28128
#>    ..working on genes 25441 through 27030 of 28128
#>    ..working on genes 27031 through 28128 of 28128
#>    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
#> 1      1   0.9700  1.8900         0.9690   15900     17400  20900
#> 2      2   0.8390  0.6230         0.8640   10700     11900  16700
#> 3      3   0.1920  0.1620         0.0549    7760      8550  14000
#> 4      4   0.0547 -0.0926        -0.2100    5910      6340  12000
#> 5      5   0.2610 -0.2700         0.0604    4650      4800  10400
#> 6      6   0.3790 -0.4050         0.2430    3750      3700   9210
#> 7      7   0.4470 -0.5090         0.3530    3090      2890   8210
#> 8      8   0.4950 -0.5910         0.4420    2580      2280   7370
#> 9      9   0.5270 -0.6590         0.5020    2180      1820   6660
#> 10    10   0.5430 -0.7290         0.5430    1860      1470   6060
#> 11    12   0.5750 -0.8320         0.6200    1390       979   5090
#> 12    14   0.5820 -0.9360         0.6600    1070       667   4340
#> 13    16   0.5970 -1.0100         0.7040     843       466   3740
#> 14    18   0.6160 -1.0600         0.7440     674       330   3250
#> 15    20   0.6360 -1.1200         0.7770     546       239   2840

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, col = "red")
```

![](Assets/wgcna_soft_threshold-1.png)<!-- -->

Pick a soft threshold power near the curve of the plot, so here we could
pick 7 or 8. We’ll pick 7 but feel free to experiment with other powers
to see how it affects your results.

``` r
netwk <- blockwiseModules(input_mat,                # <= input here
                          power = 7,                # <= power here
                          minModuleSize = 30,
                          reassignThreshold = 0, 
                          mergeCutHeight = 0.25, 
                          numericLabels = T,
                          pamRespectsDendro = F, 
                          saveTOMs = T, 
                          saveTOMFileBase = "ER",
                          verbose = 3, 
                          maxBlockSize = 40000, 
                          deepSplit = 2, 
                          detectCutHeight = 0.5,
                          networkType = "signed")
```
