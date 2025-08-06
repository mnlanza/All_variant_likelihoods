# All Variant Likelihoods

Analysis of codon variant likelihoods in the E. coli gyrase gene, comparing all 64 codon variants at codon 83 of the forward and reverse complement sequences.

## Overview

This project analyzes the likelihood of different codon variants at position 83 in the E. coli gyrase gene. It includes tools for:
- Generating all possible codon variants
- Computing likelihoods for both forward and reverse complement sequences
- Visualizing the results through various plots
- Validating reverse complement sequences

## Processing Pipeline

1. **Generate FASTA Files**:
   ```bash
   # First, generate forward strand variants
   python scripts/generate_cods_83.py  # Creates codon83_variants.fasta
   
   # Then, generate reverse complements
   python scripts/generate_r_cods_83.py  # Creates reverse_complements_83.fasta
   ```

2. **Run EVO Model**:
   - Input: Both FASTA files (`codon83_variants.fasta` and `reverse_complements_83.fasta`)
   - Output: Logits files in `input/logits/`
     - Forward: `input_83_*_logits.npy`
     - Reverse: `input_83_*_rc_logits.npy`

3. **Analyze Results**:
   ```bash
   # Run analysis script
   ./run_analysis.sh
   ```
   The script:
   - Loads FASTA sequences and logits
   - Calculates likelihoods
   - Generates comparison plots
   - Saves results to `output/`

## Directory Structure

```
.
├── input/                  # Input data directory
│   ├── fasta/             # FASTA sequence files
│   │   ├── codon83_variants.fasta        # Forward strand variants
│   │   └── reverse_complements_83.fasta  # Reverse complement variants
│   └── logits/            # Neural network logits files
│       ├── input_83_*.npy     # Forward strand logits
│       └── input_83_*_rc.npy  # Reverse complement logits
├── output/                # Generated files and results
│   ├── Norm_codon_ls.pdf     # Normalized codon likelihoods (forward)
│   ├── Norm_aa_ls.pdf        # Normalized amino acid likelihoods (forward)
│   ├── Rev_norm_codon_ls.pdf # Normalized codon likelihoods (reverse)
│   ├── Rev_norm_aa_ls.pdf    # Normalized amino acid likelihoods (reverse)
│   ├── FR_comp_ls.pdf        # Forward vs reverse likelihood comparison (normalized forward vs reversed strand)
│   └── FR_tot_lls.pdf        # Total log likelihoods comparison (forward vs reversed strand)
├── scripts/               # Helper scripts and utilities
│   ├── analyze_FR_83.R           # Main analysis script
│   ├── 83_all_fx.R              # Core analysis functions
│   ├── handle_EVO2_output.r      # Output handling utilities
│   └── evo2_analysis_functions.R # Analysis helper functions
├── tests/                # Validation and test scripts
│   └── check_reverse_complements.py  # checks that sequences in fasta is the reverse complement of another inputted fasta
└── requirements.txt      # Python dependencies
```

## Setup

1. Create and activate virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Unix/macOS
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Quick Start

Run the complete analysis using the provided shell script:
```bash
./run_analysis.sh
```

This will:
1. Create necessary directories if they don't exist
2. Run the forward/reverse analysis
3. Generate all plots in the output directory

## Manual Usage

### Main Analysis Script (scripts/analyze_FR_83.R)

Analyzes both forward and reverse complement sequences:

```bash
# Basic usage
Rscript scripts/analyze_FR_83.R \
  --for_fasta input/fasta/codon83_variants.fasta \
  --rev_fasta input/fasta/reverse_complements_83.fasta \
  --for_raw_dir input/logits \
  --save \
  --output_dir output

# All options
Rscript scripts/analyze_FR_83.R \
  --for_fasta <forward_fasta>        # Forward strand FASTA file
  --rev_fasta <reverse_fasta>        # Reverse complement FASTA file
  --for_raw_dir <forward_raw_dir>    # Directory with forward logits
  --rev_raw_dir <reverse_raw_dir>    # optional: default assumes for_raw_dir has rev logits
  --save                             # Save plots to PDF
  --output_dir <output_dir>          # Output directory (default: output)
  --highlight <positions>            # Comma-separated positions to highlight
```

## Utility Scripts

### scripts/generate_fasta.py

Generate FASTA files with codon variants.

```python
from scripts.generate_fasta import generate_all_83_variants

# Generate all possible codon variants at position 83
generate_all_83_variants(
    seq="ATGC...",                    # Reference sequence
    output_filename="variants.fasta",  # Output file
    codon_index=83,                   # Position to vary (1-based)
    same_dir_as_script=False          # If True, save next to script
)
```

### tests/check_reverse_complements.py

Validate reverse complement sequences:

```bash
python tests/check_reverse_complements.py \
  input/fasta/codon83_variants.fasta \
  input/fasta/reverse_complements_83.fasta
```

## File Requirements

1. FASTA Files (in `input/fasta/`):
   - Must contain sequence variants
   - Headers should follow format: `{position}_{codon}` (e.g., "83_AAA")
   - No description field required

2. Logits Files (in `input/logits/`):
   - Must be NumPy (.npy) files
   - Forward strand: `input_{position}_{codon}_logits.npy`
   - Reverse strand: `input_{position}_{codon}_rc_logits.npy`
