# **GFF3 Cleaning, Standardization, and Protein Isoform Extraction**

#### Author: Jeonghwa (J) Parkim
#### Repository: https://github.com/j-parkim/Capstone
#### Version: 2025-11-30

## **Overview**
This repository contains a toolkit designed to solve real-world challenges when working with genome annotation files(GFF3) from sources such as RefSeq, GenBank, and Ensembl/Gencode. 

Real annotatino files are often **messy** - contatining inconsistent formats, invalid attributes, multi-source predictions, malformed hierarchies, and irrelevant organellar genome annotations. 

This toolkit is intended to provide a clean, modular, lightweight, reproducible pipeline for:
1. **Inspecting** GFF/GTF annotation structure.
2. (Optional) **Standardizing** inconsistent or non-standard GFF dialects
3. (Optioanl) **Cleaning** annotation files
    - removing unwanted genome regions and its associated info (e.g., chloroplast, plastID)
    - removing unwanted annotation sources (e.g., tRNAscan-SE, cmsearch)
    - keep only protein_coding genes with complete hierarchy (gene-mRNA-CDS)
4. **Extracting CDS + Protein FASTA** using `gffread`
5. **Selecting ONE representative isoform per gene**
6. Producing **summary** tables and cleaned FASTA outputs

The final output is:
<pre>representative_isoforms.faa
isoform_summary.tsv</pre>

This is the complete workflos used for my **NYU Capstone Project.**

## **Pipeline at a Glance**
 ┌──────────────────────────────┐
 │          Input GFF           │
 └───────────────┬──────────────┘
                 │
                 ▼
      (Step 0) Inspect Structure
      * ParseGffinfo.py
      - Sources
      - Feature Types
      - Genome
      - attribute summary for each source and featuretype
                 │
                 ▼
      (Step 0.5) (Optional) Standardize dialects
      * Reformat2GFF.py
      - Normalize separators, assigners, delimiter
      - Normalize keys (biotype, Name)
      - Fix missing IDs (if applicable)
                 │
                 ▼
      (Step 1) (Optional) Cleanup
      * Parse_whatUwant_gff.py
      - Remove unwanted genome regions and its seqIDs 
      - Remove data from unwanted sources
      - Keep only protein_coding gene → mRNA → CDS
                 │
                 ▼
      (2) Extract FASTA
      * Parse_ONE_Isoform.py
      - run `gffread`
      - output: CDS FASTA, Protein FASTA
                 │
                 ▼
      (3) Select Representative Isoform
      * Parse_ONE_Isoform.py
      Options:
        - longest  (default)
        - first
        - random
                 │
                 ▼
      ┌───────────────────────────────────┐
      │    Cleand GFFs                    |
      |    Protein FASTA                  |
      |    representative_isoforms.faa    │
      │    isoform_written_summary.tsv    │
      └───────────────────────────────────┘


## **Dependencies**
#### Required
- Python v3.13 or above
- `gffread`
- Python built-in packages:
    - argparse, os, subprecess, shutil, random, re, collections

## **Installation**
<pre>
git clone https://github.com/j-parkim/Capstone.git
cd Capstone
</pre>

Then, run directly<br>
```
python gff2faa.py \
    -gff annotationfile.gff3 \
    -genome genomefasta.fna \
    --reformat \
    --remove-genome mitochondrion,chloroplast \ # User input prompt
    --remove-source cmsearch \ # User input prompt
    --criteria longest \
    --outdir results
```

## **Available Options**
**Required**:
    -gff                 Input GFF3 annotation file  
    -genome              Genome FASTA (for gffread)

**Optional**:
    --reformat           Standardize attributes (GFF3 format)
    --remove-genome      Remove genome regions 
                        # option: enter in the commandline or user can input in prompt after summary is printed
    --remove-source      Remove sources 
                        # option: enter in the commandline or user can input in prompt after summary is printed
    --skip-cleanup       Do NOT clean → use original GFF

    --skip-gffread       Do NOT run gffread
    --protein-fasta      Provide protein FASTA manually when skipping gffread

    --criteria           longest | first | random
    --skip-isoform       Skip isoform selection (keep all proteins)

    --outdir             Output directory (default=results)

## Future Directions

#### 