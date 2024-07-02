This package contains data analyzed using the DADA2 pipeline, derived from HMP data originally generated through pyrosequencing and OTU-based methods.

  
## Install
```
library(devtools)
devtools::install_github("KitHubb/HMPData")
```
  
  
## Input data
1. HMP V1V3
- Source: https://qiita.ucsd.edu/study/description/1927
- Samples: 3,530

2. HMP V3V5
- Source: https://qiita.ucsd.edu/study/description/1928
- Samples: 6,346


## Dataset
1. `V13p5`
- Using dada2-pyro plugin
- Truncate the read length to 500 and remove 20 bp from the forward
2. `V13p4`
- Using dada2-pyro plugin
- Truncate the read length to 450 and remove 20 bp from the forward
3. `V13s5`
- Using dada2-single plugin
- Truncate the read length to 500 and remove 20 bp from the forward
4. `V13s4`
- Using dada2-single plugin
- Truncate the read length to 450 and remove 20 bp from the forward

  
## Usage
```
library(phyloseq)
library(HMPData)

data(‘V13p5’) # dada2-pyro, trunc-length 500
data(‘V13p4’) # dada2-pyro, trunc-length 450
data(‘V13s5’) # dada2-single, trunc-length 500
data(‘V13s4’) # dada2-single, trunc-length 450
```

```
V13p5
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 52283 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 52283 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 52283 tips and 51732 internal nodes ]

V13p4
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 26035 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 26035 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 26035 tips and 25774 internal nodes ]

V13s5
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 26754 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 26754 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 26754 tips and 26612 internal nodes ]

V13s4
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 53762 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 53762 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 53762 tips and 53073 internal nodes ]


```

#### HMP35(not yes)


## Statistics
```
summary(sample_sums(V13p5))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   101.0   916.5  1385.1  1945.8 29071.0

summary(sample_sums(V13p4))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   945.5  3578.0  3779.7  5330.8 70684.0 

summary(sample_sums(V13s5))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   100.2   905.0  1362.7  1921.5 28891.0 

summary(sample_sums(V13s4))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     933    3549    3742    5285   70158 
```

## Mathods
- Tools
  - QIIME2(v2024.02)
  - DADA2(v1.32.0)
  - phyloseq(v1.48.0)
  - SILVA Database(v138)

## Preprocessing
#### 1) Prepare dataset
(1) rearrange mapping file  
- In the V1V3 dataset, there are 10 multiplexed groups in 3,530 samples
- mapping file save in `/QIIME_preprocessing/Mapping_files/`
  
(2) change `.fna` to `.fasta`
```
for file in *.fna; do # using root
    cp -- "$file" "${file%.fna}.fasta"
done 
```
  
(3) Merge `.fna` + `.qual` to `.fastq`
```
# instapp biopython
python -m pip install --upgrade pip
pip install biopython
pip install biopython –-upgrade

```
  
reference: https://gist.github.com/necrolyte2/b45a82fb4ecb0ffd70ab#file-fastaqual_too_fastq-py-L1
```
#!/usr/bin/env python

import sys

from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file')
parser.add_argument('qual', help='Qual file')
args = parser.parse_args()

records = PairedFastaQualIterator(
    open(args.fasta),
    open(args.qual)
)
for rec in records:
    sys.stdout.write(rec.format('fastq'))
```
  
```
chmod +x fastaqual_too_fastq.py

for file in $(ls | sed -E 's/\.[^/.]+$//' | sort | uniq); do  ../fastaqual_too_fastq.py \
./${file}.fna ./${file}.qual >  ../FASTQ/${file}.fastq ; done

```
  
check the number of samples
```
ls -al | grep ^- | wc -l
```
  
  
#### 2) Demultiplexing
Make bash script for demultiplexing using QIIME2
- script file save in `/QIIME_preprocessing/script/`
- reference: https://forum.qiime2.org/t/analyzing-454-data-in-qiime-2/25055
```
#!/bin/bash

# Check if at least one argument is passed
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <list_of_SRR_identifiers>"
  exit 1
fi

# List of SRR identifiers from command-line arguments
list="$@"

# Loop through each identifier
for multi in $list; do
  # Create directory for each identifier
  mkdir -p ${multi}/

  # Copy fasta and qual files to the new directory
  cp ./Rawdata/fastaqual/${multi}.fasta ${multi}/
  cp ./Rawdata/fastaqual/${multi}.qual ${multi}/

  # Change to the new directory
  cd ${multi}

  # Rename fasta and qual files
  mv ${multi}.fasta reads.fasta
  mv ${multi}.qual reads.qual

  # Import qiime artifact
  qiime tools import \
    --type MultiplexedSingleEndBarcodeInSequence \
    --input-format MultiplexedFastaQualDirFmt \
    --input-path ./ \
    --output-path ${multi}_seqs.qza

  # Demultiplex the reads
  qiime cutadapt demux-single \
    --i-seqs ${multi}_seqs.qza \
    --m-barcodes-file ../mapping_files/HMPV13_qiime2_mapping_FINAL.txt \
    --m-barcodes-column barcode \
    --o-per-sample-sequences ${multi}_demultiplexed-seqs.qza \
    --o-untrimmed-sequences ${multi}_untrimmed.qza

  # Visualization
  qiime demux summarize \
    --i-data ${multi}_demultiplexed-seqs.qza \
    --o-visualization ${multi}_demultiplexed-seqs.qzv

  # Export the demultiplexed sequences
  qiime tools export \
    --input-path ${multi}_demultiplexed-seqs.qza \
    --output-path ${multi}_demultiplexed-seqs/

  # Return to the parent directory
  cd ..

done
```
  
Run script
```
chmod +x process_qiime.sh
./process_qiime.sh SRR045723 SRR047558 SRR057663 SRR058087 SRR058088 SRR058091 SRR058094 SRR058097 SRR058107 SRR058115
```
  
#### 3) Analysis in QIIME2 Env
import demultiplexed `fastq.gz` files to qiime2 artifact
- manifest file save in `/QIIME_preprocessing/`
```
conda activate qiime2-amplicon-2024.02
```
  

```
qiime tools import   \
--type 'SampleData[SequencesWithQuality]'   \
--input-path ./mapping_files/HMPV13_qiime2_manifest_total.txt   \
--output-path single-end-demux.qza   \
--input-format SingleEndFastqManifestPhred33V2 
```
  
Adapter trimming
```
qiime cutadapt trim-single \
  --i-demultiplexed-sequences single-end-demux.qza \
  --p-front ATTACCGCGGCTGCTGG  \
  --p-error-rate 0 \
  --p-discard-untrimmed \
  --o-trimmed-sequences single-end-trimmed.qza \
  --verbose
```
  
Denoising
```
qiime dada2 denoise-pyro \
  --i-demultiplexed-seqs single-end-trimmed.qza \
  --p-trunc-len 500 \#500 or 400
  --p-trim-left 20 \
  --output-dir dada2-out
```
  
```
# Assignment
qiime feature-classifier classify-sklearn \
  --i-classifier /data/Reference/16S/QIIME2/SILVA/138version/silva-138-99-full-length-nb-classifier.qza \
  --i-reads dada2-out/representative_sequences.qza \
  --o-classification dada2-out/taxonomy.qza
  
# Make phylogenetic tree 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences dada2-out/representative_sequences.qza \
   --output-dir dada2-out/tree

```
  
#### 4) QIIME2 to Phyloseq object
```
library(qiime2R)
library(phyloseq)
library(stringr)
library(dplyr)

ps<-qza_to_phyloseq(
    features="../input_V1V3_qiime/pyro-500/table.qza",
    tree="../input_V1V3_qiime/pyro-500/tree/rooted_tree.qza",
    taxonomy="../input_V1V3_qiime/pyro-500/taxonomy.qza",
    metadata = "../input_V1V3_qiime/1927_20230202-080822.txt"
    )
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6237 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 6237 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6237 tips and 6204 internal nodes ]
```

#### 5) Modify taxonomy 
- Reference: https://www.yanh.org/2021/01/01/microbiome-r/

Filtering Unassigned, Chloroplast, Mitochondria, Archaea
```
physeq <- subset_taxa(physeq, Kingdom %in% "Bacteria" ) 
physeq <- subset_taxa(physeq, Order %!in%  "Chloroplast" ) 
physeq <- subset_taxa(physeq, Family %!in%  "Mitochondria") 
```
    
Modify taxonomy format
- "NA", "_sp." to "_unclassified"
- "uncultured" to "Genus_uncultured"  
```
tax_clean <- function(TAX){ 
  # remove k_
  TAX.clean <- data.frame(row.names = row.names(TAX),
                          Kingdom = str_replace(TAX[,1], "d__",""),
                          Phylum = str_replace(TAX[,2], "p__",""),
                          Class = str_replace(TAX[,3], "c__",""),
                          Order = str_replace(TAX[,4], "o__",""),
                          Family = str_replace(TAX[,5], "f__",""),
                          Genus = str_replace(TAX[,6], "g__",""),
                          Species = str_replace(TAX[,7], "s__",""),
                          stringsAsFactors = FALSE)
  # replace NA to "" and other things
  TAX.clean[TAX.clean=="-"] <- ""
  TAX.clean[is.na(TAX.clean)] <- ""
  TAX.clean[TAX.clean == "NA"] <- ""
  
  # Spcies to Genus Species
  tax.clean <- TAX.clean
  
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,7] != ""){
      tax.clean$Species[i] <- paste(# tax.clean$Genus[i], 
        tax.clean$Species[i]# , sep = " "
        )
    } else if (tax.clean[i,2] == ""){
      kingdom <- paste(tax.clean[i,1], "unclassified",  sep = "_")
      tax.clean[i, 2:7] <- kingdom
    } else if (tax.clean[i,3] == ""){
      phylum <- paste(tax.clean[i,2], "unclassified", sep = "_")
      tax.clean[i, 3:7] <- phylum
    } else if (tax.clean[i,4] == ""){
      class <- paste(tax.clean[i,3], "unclassified", sep = "_")
      tax.clean[i, 4:7] <- class
    } else if (tax.clean[i,5] == ""){
      order <- paste(tax.clean[i,4], "unclassified", sep = "_")
      tax.clean[i, 5:7] <- order
    } else if (tax.clean[i,6] == ""){
      family <- paste(tax.clean[i,5], "unclassified",  sep = "_")
      tax.clean[i, 6:7] <- family
    } else if (tax.clean[i,7] == ""){
      tax.clean$Species[i] <- paste(tax.clean$Genus[i], "unclassified", sep = "_")
    }
  }

  tax.clean$Species <- gsub("*_sp.", "_unclassified", tax.clean$Species)
  
  tax.clean[tax.clean$Species %in% "uncultured_bacterium", "Species"] <- NA
  tax.clean[grepl("uncultured_", tax.clean$Species), "Species"] <- NA
  tax.clean[is.na(tax.clean)] <- ""
  
    tax.clean2 <- tax.clean
  
  for (i in 1:nrow(tax.clean2)){
    if (tax.clean2[i,7] != ""){
      tax.clean2$Species[i] <- paste(tax.clean2$Species[i])
    } else if (tax.clean2[i,2] == ""){
      kingdom <- paste(tax.clean2[i,1], "uncultured",  sep = "_")
      tax.clean2[i, 2:7] <- kingdom
    } else if (tax.clean2[i,3] == ""){
      phylum <- paste(tax.clean2[i,2], "uncultured", sep = "_")
      tax.clean2[i, 3:7] <- phylum
    } else if (tax.clean2[i,4] == ""){
      class <- paste(tax.clean2[i,3], "uncultured", sep = "_")
      tax.clean2[i, 4:7] <- class
    } else if (tax.clean2[i,5] == ""){
      order <- paste(tax.clean2[i,4], "uncultured", sep = "_")
      tax.clean2[i, 5:7] <- order
    } else if (tax.clean2[i,6] == ""){
      family <- paste(tax.clean2[i,5], "uncultured",  sep = "_")
      tax.clean2[i, 6:7] <- family
    } else if (tax.clean2[i,7] == ""){
      tax.clean2$Species[i] <- paste(tax.clean2$Genus[i], "uncultured", sep = "_")
    }
  }
    
      return(tax.clean2)
}



tax <- data.frame(tax_table(physeq))
tax.c <- tax_clean(tax)
tax_table(physeq) <- tax_table(as.matrix(tax.c))
```
