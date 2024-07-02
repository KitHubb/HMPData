This package contains data analyzed using the DADA2 pipeline, derived from HMP data originally generated through pyrosequencing and OTU-based methods.

  
## Install
```
library(devtools)
devtools::install_github("KitHubb/HMPData")
```
  
  
## Dataset
1. HMP13
- source: https://qiita.ucsd.edu/study/description/1927
- samples: 3,530
- object: Phyloseq
- Methods: QIIME2(2024.02), DADA2(1.32.0), SILVA Database(138)

    
2. HMP35
- source: https://qiita.ucsd.edu/study/description/1928
- samples: 6,346
- object: Phyloseq
- Methods: QIIME2, DADA2, SILVA Database
    
## Data
1. `V13p5`
- Using dada2-pyro
- Truncate the read length to 500 and remove 20 bp from the forward
2. `V13p4`
- Using dada2-pyro
- Truncate the read length to 450 and remove 20 bp from the forward
3. `V13s5`
- Using dada2-single
- Truncate the read length to 500 and remove 20 bp from the forward
4. `V13s4`
- Using dada2-single
- Truncate the read length to 450 and remove 20 bp from the forward

  
## Usage
```
library(phyloseq)
library(HMPData)

data(‘V13p5’) # dada2 pyrosequencing, trunc-length 500
V13p5
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 52512 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 52512 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 52512 tips and 51961 internal nodes ]

data(‘V13p4’) # dada2 pyrosequencing, trunc-length 450
V13p4
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 26053 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 26053 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 26053 tips and 25792 internal nodes ]


data(‘V13s5’) # dada2 single, trunc-length 500
V13s5
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 26773 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 26773 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 26773 tips and 26631 internal nodes ]


data(‘V13s4’) # dada2 single, trunc-length 450
V13s4
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 53991 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 53991 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 53991 tips and 53302 internal nodes ]


```

#### HMP35(not yes)


## Statistics
```


```


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
  --p-trunc-len 515 \
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

physeq<-qza_to_phyloseq(
    features="../input_V1V3_qiime/pyro-500/table.qza",
    tree="../input_V1V3_qiime/pyro-500/tree/rooted_tree.qza",
    taxonomy="../input_V1V3_qiime/pyro-500/taxonomy.qza",
    metadata = "../input_V1V3_qiime/1927_20230202-080822.txt"
    )
physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6237 taxa and 3530 samples ]
# sample_data() Sample Data:       [ 3530 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 6237 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6237 tips and 6204 internal nodes ]

```

Modify taxonomy format
- reference: https://www.yanh.org/2021/01/01/microbiome-r/
- `NA`, `_sp.` to `_unclassified`
- "uncultured" to "Genus_uncultured"  

#### 4) BLAST sequences
- Cut off: evalue< 1e-10, qcovus> 99,  pident> 98.75
- Exchange 32 unclassified ASVs as BLAST results

