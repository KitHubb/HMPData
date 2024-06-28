This package contains data analyzed using the DADA2 pipeline, derived from HMP data originally generated through pyrosequencing and OTU-based methods.

## Install
```
library(devtools)
devtools::install_github("KitHubb/HMPData")
```


## Dataset
`HMP13`
source: https://qiita.ucsd.edu/study/description/1927
samples: 3,530

`HMP35`
source: https://qiita.ucsd.edu/study/description/1928
samples: 6,346


## Usage 
#### HMP13
```
data('HMP13')
HMP13
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6237 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 6237 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6237 tips and 6204 internal nodes ]
```


#### HMP35(not yes)


## Preprocessing
#### 1) Demultiplexing
```

```
#### 2) Analysis in QIIME2 Env

```


```
#### 3) QIIME2 to Phyloseq object

#### 4) BLAST sequences

