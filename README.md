Omnipath to hipathia
================
Martin Garrido Rodriguez-Cordoba

## Abstract

The main goal of this parser is to transform an
[Omnipath](http://omnipathdb.org/) formatted set of interactions into an
object usable by the mechanistic modelling tool
[hipathia](http://hipathia.babelomics.org/). Hipathia uses a signal
propagation algorithm to estimate the activity of receptor-to-effector
signaling circuits. Its current version uses a pathway-centric approach
fed by a mixture of context-specific and generic signaling pathways as
defined by the [KEGG
classification](https://www.genome.jp/kegg/pathway.html). On the other
hand, Omnipath is a meta-resource that contains biological information
from different databases in a ready-to-use tabular format. It comprises
several types of relationships between biological entities as
protein-protein interactions or TF-target relationships, as well as gene
and protein functional annotations.

## Case of use

Load required packages

``` r
library(OmnipathR)
library(hipathia)
library(igraph)
library(dplyr)
library(Rgraphviz)
library(org.Hs.eg.db)
source("R/omnipath-to-hipathia.R")
```

Import omnipath interactions and subset to directed interactions which
are not “auto” interactions.

``` r
# import all interactions in omnipath and intercell information
interactions <- OmnipathR::import_Omnipath_Interactions()
```

    ## Downloaded 43898 interactions

    ## removed 0 interactions during database filtering.

``` r
# filter only to directed interactions (consensus) and remove self interactions
intInteractions <- subset(interactions, 
                          (consensus_stimulation == 1 & consensus_inhibition == 0) | 
                            (consensus_stimulation == 0 & consensus_inhibition == 1)) %>%
  subset(source_genesymbol != target_genesymbol)
# print head of interesting interactions
head(intInteractions) %>% OmnipathR::print_interactions()
```

    ##            source interaction         target nsources nrefs
    ## 5   CAV1 (Q03135)  ==( + )==> TRPC1 (P48995)        2     8
    ## 2  CALM1 (P0DP23)  ==( - )==> TRPC1 (P48995)        1     3
    ## 3  CALM3 (P0DP25)  ==( - )==> TRPC1 (P48995)        1     3
    ## 4  CALM2 (P0DP24)  ==( - )==> TRPC1 (P48995)        1     3
    ## 10 FKBP4 (Q02790)  ==( - )==> TRPC1 (P48995)        1     3
    ## 7   DRD2 (P14416)  ==( + )==> TRPC1 (P48995)        1     1

To test the parser, we will create a network using functional
annotations from KEGG. We will select the interactions between genes
annotated within the “Retrograde endocannabinoid signaling” pathway.

``` r
# import kegg annotations
keggAnnotations <- OmnipathR::import_Omnipath_annotations(filter_databases = "KEGG")
```

    ## Downloaded 16735 annotations.

``` r
# subset to pathway of interest
intGenes <- subset(keggAnnotations, value == "Retrograde endocannabinoid signaling") %>% 
  pull(genesymbol)
selectedInteractions <- subset(intInteractions, source_genesymbol %in% intGenes & target_genesymbol %in% intGenes)
# print head of interesting interactions
head(selectedInteractions) %>% OmnipathR::print_interactions()
```

    ##                source interaction           target nsources nrefs
    ## 10894   GNAQ (P50148)  ==( + )==>   PLCB1 (Q9NQ66)        5     6
    ## 18234 PRKACA (P17612)  ==( - )==>   KCNJ3 (P48549)        5     3
    ## 10865   GRM1 (Q13255)  ==( + )==>    GNAQ (P50148)        3     3
    ## 11505 PRKACA (P17612)  ==( - )==>   PLCB1 (Q9NQ66)        3     2
    ## 11647 PRKACA (P17612)  ==( - )==> CACNA1A (O00555)        3     1
    ## 43619   GNAQ (P50148)  ==( + )==>  MAPK14 (Q16539)        3     0

Use the **omnipathToHipathia** function to create the metaginfo object.

``` r
mgi <- omnipathToHipathia(selectedInteractions)
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Parsing a graph with 7 nodes and 6 edges.

    ## Writing on /tmp/RtmpmpbiZw

    ## Starting hipathia mgi creation from sif...

    ## Loading graphs...

    ## Creating MGI...

    ## Created MGI with 1 pathway(s)

    ## Done!

Plot the resulting graph object, which has been divided into circuits
usable by the hipathia function.

``` r
g <- mgi$pathigraphs$hsa00$graph
plot(g, layout = graphVizLayout(g), vertex.shape = "none")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Note

This example is carried out with an extremely simple network. **The
creation of the MGI object takes a really long time on complex graphs**,
specially when there are lots of edges between nodes which are not
receptors nor effectors.

## Session info

``` r
sessionInfo()
```

    ## R version 4.0.0 (2020-04-24)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] org.Hs.eg.db_3.11.0         AnnotationDbi_1.50.0       
    ##  [3] Rgraphviz_2.32.0            graph_1.66.0               
    ##  [5] dplyr_0.8.5                 hipathia_2.4.0             
    ##  [7] MultiAssayExperiment_1.14.0 SummarizedExperiment_1.18.1
    ##  [9] DelayedArray_0.14.0         matrixStats_0.56.0         
    ## [11] Biobase_2.48.0              GenomicRanges_1.40.0       
    ## [13] GenomeInfoDb_1.24.0         IRanges_2.22.1             
    ## [15] S4Vectors_0.26.0            AnnotationHub_2.20.0       
    ## [17] BiocFileCache_1.12.0        dbplyr_1.4.3               
    ## [19] BiocGenerics_0.34.0         OmnipathR_1.2.0            
    ## [21] igraph_1.2.5               
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6                  lattice_0.20-41              
    ##  [3] tidyr_1.0.2                   assertthat_0.2.1             
    ##  [5] digest_0.6.25                 mime_0.9                     
    ##  [7] R6_2.4.1                      RSQLite_2.2.0                
    ##  [9] evaluate_0.14                 httr_1.4.1                   
    ## [11] pillar_1.4.4                  zlibbioc_1.34.0              
    ## [13] rlang_0.4.6                   curl_4.3                     
    ## [15] blob_1.2.1                    Matrix_1.2-18                
    ## [17] preprocessCore_1.50.0         rmarkdown_2.1                
    ## [19] servr_0.16                    stringr_1.4.0                
    ## [21] RCurl_1.98-1.2                bit_1.1-15.2                 
    ## [23] shiny_1.4.0.2                 compiler_4.0.0               
    ## [25] httpuv_1.5.2                  xfun_0.13                    
    ## [27] pkgconfig_2.0.3               htmltools_0.4.0              
    ## [29] tidyselect_1.0.0              tibble_3.0.1                 
    ## [31] GenomeInfoDbData_1.2.3        interactiveDisplayBase_1.26.0
    ## [33] crayon_1.3.4                  later_1.0.0                  
    ## [35] bitops_1.0-6                  rappdirs_0.3.1               
    ## [37] jsonlite_1.6.1                xtable_1.8-4                 
    ## [39] lifecycle_0.2.0               DBI_1.1.0                    
    ## [41] magrittr_1.5                  stringi_1.4.6                
    ## [43] XVector_0.28.0                promises_1.1.0               
    ## [45] ellipsis_0.3.0                vctrs_0.2.4                  
    ## [47] tools_4.0.0                   bit64_0.9-7                  
    ## [49] glue_1.4.0                    purrr_0.3.4                  
    ## [51] BiocVersion_3.11.1            fastmap_1.0.1                
    ## [53] yaml_2.2.1                    BiocManager_1.30.10          
    ## [55] memoise_1.1.0                 knitr_1.28
