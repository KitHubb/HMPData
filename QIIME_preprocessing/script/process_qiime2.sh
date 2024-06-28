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
  echo "Processing ${multi}..."

  # Create directory for each identifier
  echo "Creating directory for ${multi}..."
  mkdir -p ${multi}/

  # Copy fasta and qual files to the new directory
  echo "Copying fasta and qual files for ${multi}..."
  cp ./Rawdata/fastaqual/${multi}.fasta ${multi}/
  cp ./Rawdata/fastaqual/${multi}.qual ${multi}/

  # Change to the new directory
  cd ${multi}

  # Rename fasta and qual files
  echo "Renaming fasta and qual files for ${multi}..."
  mv ${multi}.fasta reads.fasta
  mv ${multi}.qual reads.qual

  # Import qiime artifact
  echo "Importing qiime artifact for ${multi}..."
  qiime tools import \
    --type MultiplexedSingleEndBarcodeInSequence \
    --input-format MultiplexedFastaQualDirFmt \
    --input-path ./ \
    --output-path ${multi}_seqs.qza

  # Demultiplex the reads
  echo "Demultiplexing reads for ${multi}..."
  qiime cutadapt demux-single \
    --i-seqs ${multi}_seqs.qza \
    --m-barcodes-file ../mapping_files/HMPV13_qiime2_mapping.txt \
    --m-barcodes-column barcode \
    --o-per-sample-sequences ${multi}_demultiplexed-seqs.qza \
    --o-untrimmed-sequences ${multi}_untrimmed.qza

  # Visualization
  echo "Creating visualization for ${multi}..."
  qiime demux summarize \
    --i-data ${multi}_demultiplexed-seqs.qza \
    --o-visualization ${multi}_demultiplexed-seqs.qzv

  # Export the demultiplexed sequences
  echo "Exporting demultiplexed sequences for ${multi}..."
  qiime tools export \
    --input-path ${multi}_demultiplexed-seqs.qza \
    --output-path ${multi}_demultiplexed-seqs/

  # Return to the parent directory
  cd ..

  echo "Finished processing ${multi}."
done

echo "All processes completed."

