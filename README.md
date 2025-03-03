# viral_polyprotein_pipeline
Initial commit of viral polyprotein pipeline
# Viral Polyprotein Database Generator

A series of R scripts to create a comprehensive database of viral polyproteins and their mature peptide products from NCBI.

## Overview

This pipeline automates the following tasks:

1. Download all complete viral polyproteins from NCBI
2. Extract mature proteins from the nucleotide records
3. Track polyproteins that lack mature protein annotations
4. Use BLASTX to annotate unannotated polyproteins based on similarity to known mature proteins
5. Generate standardized FASTA files with custom formatted headers

## Requirements

- R (>= 4.0)
- Required R packages (will be automatically installed):
  - rentrez: For accessing NCBI databases
  - Biostrings: For biological sequence handling
  - stringr: For string manipulation
  - seqinr: For sequence handling
  - dplyr: For data manipulation
- NCBI BLAST+ (for the BLAST annotation step)

## Installation

1. Clone or download this repository:
```
git clone https://github.com/mihinduk/viral-polyprotein-database.git
cd viral-polyprotein-database
```

2. Make the scripts executable (optional, for Unix/Linux/Mac):
```
chmod +x *.R
```

## Usage

### Quick Start

To run the complete pipeline:

```
Rscript run_pipeline.R [max_records] [email] [api_key]
```

- `max_records`: (Optional) Limit to this many records (for testing)
- `email`: (Optional) Your email for NCBI E-utilities
- `api_key`: (Optional) Your NCBI API key for higher rate limits

If email/API key are not provided, you will be prompted to enter them.

### Step-by-Step Execution

You can also run each step of the pipeline separately:

1. First, set up your NCBI credentials:
```
Rscript 01_config.R
```

2. Download viral polyproteins:
```
Rscript 02_download_polyproteins.R [max_records]
```

3. Extract mature proteins:
```
Rscript 03_extract_mature_proteins.R [max_records]
```

4. Run BLAST annotation:
```
Rscript 04_blast_annotation.R [mature_file] [missing_file]
```

## Output Files

The pipeline generates several date-stamped files:

1. `YYYY_MM_DD_viral_polyprotein.fasta`: All viral polyproteins
2. `YYYY_MM_DD_viral_polyprotein_mature_proteins.fasta`: Mature proteins extracted from polyproteins
3. `YYYY_MM_DD_viral_polyprotein_missing_annotation.fasta`: Polyproteins without mature protein annotations
4. `YYYY_MM_DD_viral_polyprotein_missing_annotation.txt`: Metadata for unannotated polyproteins
5. `YYYY_MM_DD_blastx_results.txt`: BLAST results for unannotated polyproteins

## FASTA Header Format

All FASTA headers follow this standardized format:
```
>Accession|Viruses|phylum|class|order|family|genus|species|organism|protein_name
```

Missing taxonomy fields are replaced with placeholders (e.g., `Anelloviridae_no_genus`).

## Customization

You can modify the scripts to adjust:

- BLAST parameters (e-value, identity threshold, etc.)
- Search criteria for viral polyproteins
- FASTA header format

## Notes

- The pipeline uses the NCBI E-utilities API, which requires an email address and benefits from an API key
- BLAST annotation is resource-intensive for large datasets
- For very large datasets, consider running on a high-performance computing system

## License

This software is provided under the MIT License.

## Citation

If you use this pipeline in your research, please cite:

```
Your Name et al. (2025). Viral Polyprotein Database Generator: A pipeline for comprehensive viral protein database creation. Journal of Example Research, 1(1), 1-10.
```
