정식 배포 XXX


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



## Comparison 

## Preprocessing
#### 1) Demultiplexing
```
conda activate qiime2-amplicon-2024.03

```

change `.fna` to `.fasta`
```
for file in *.fna; do # using root
    cp -- "$file" "${file%.fna}.fasta"
done 

```
Make function for demultiplexing
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
./process_qiime.sh SRR058107 SRR058107  ...
```

#### 2) Analysis in QIIME2 Env
import qiime2 artifect
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

#### 3) QIIME2 to Phyloseq object


#### 4) BLAST sequences

