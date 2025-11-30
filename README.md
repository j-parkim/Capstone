# **GFF → Representative Protein FASTA Pipeline**

#### Author: Jeonghwa (J) Parkim
#### Version: 2025-11-30

## **Overview**
This is to provide a modular, lightweight, reproducible pipeline for:
1. **Inspecting** GFF/GTF annotation files.
2. (Optional) **Reformatting** inconsistent or non-standard GFF/GTF dialects
3. (Optioanl) **Cleaning** annotation files
    - remove unwanted genome regions and its associated info (e.g., chloroplast, plastID)
    - remove unwanted sources
    - keep only protein_coding genes with complete hierarchy (gene-mRNA-CDS)
4. **Extracting CDS + Protein FASTA** using `gffread`
5. **Selecting ONE representative isoform per gene**

The final output is:
<pre>representative_isoforms.faa
isoform_summary.tsv</pre>

## **Workflow**
Step 0: Inspect GFF<br>
    |<br>
    |-- ParseGFFinfo.py<br>
    |<br>
Step 0.5: (Optional) Reformat dialects<br>
    |
    |-- Reformat2GFF.py
    |
Step 1: (Optional) Cleanup
    |
    |-- Parse_whatUwant_gff.py
    |
Step 2: Generate FASTA
    |
    |-- Parse_ONE_Isoform.py (gffread)
    |
Step 3: Representative Isoform Selection
    |
    |-- Parse_ONE_Isoform.py
    V
Final representative FASTA + summary

#### **All run by `gff_to_rep_protein_fasta.py` master CLI pipeline**

## **Dependencies**
#### Required
- Python v3.13 or above
- `gffread`
- Python built-in packages:
    - argparse, os, subprecess, shutil, random, re, collections

## **Installation**
<pre>```
git clone https://github.com/j-parkim/Capstone.git
cd Capstone
```</pre>

Then run directly
`python gff_to_rep_protein_fasta.py -gff <annotation.gff3> -genome <genomefasta.fna> [options]`

## **Available Options**
Required:
  -gff                 Input GFF3 annotation file  
  -genome              Genome FASTA (for gffread)

Optional:
  --reformat           Standardize attributes (GFF3 format)
  --remove-genome      Remove genome regions (e.g., chloroplast mitochondrion)
  --remove-source      Remove sources (e.g., cmsearch, trnascan)
  --skip-cleanup       Do NOT clean → use original GFF

  --skip-gffread       Do NOT run gffread
  --protein-fasta      Provide protein FASTA manually when skipping gffread

  --criteria           longest | first | random
  --skip-isoform       Skip isoform selection (keep all proteins)

  --outdir             Output directory (default=results)
