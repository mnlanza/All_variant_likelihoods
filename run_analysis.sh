#!/bin/bash

# Exit on error
set -ex

# Ensure directories exist
mkdir -p output input/fasta input/logits

# Run the forward/reverse analysis
Rscript scripts/main_FR_83_all.R \
  --for_fasta input/fasta/codon83_variants.fasta \
  --rev_fasta input/fasta/reverse_complements_83.fasta \
  --for_raw_dir input/logits \
  --save \
  --output_dir output

# Optional: Check reverse complements
# Uncomment the following line to run the check
# python tests/check_reverse_complements.py \
#   input/fasta/codon83_variants.fasta \
#   input/fasta/reverse_complements_83.fasta

echo "✅ Analysis complete! Check output/ directory for results." 