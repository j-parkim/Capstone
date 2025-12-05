# Parse_ONE_Isoform.py
"""
Run gffread on a cleaned GFF (obtained from Parse_whatUwant_gff.py)
Then, parse FASTA and select ONE representative isoform per gene.

***ONE dependency*** : gffread

Representative logic supported: (11/30/2025)
    - longest
    - random

Outputs:
    - representative_isoforms.fasta
    - isoform_summary.tsv

"""
import os
import sys
import shutil
import subprocess
import random
import gzip
import tempfile
from collections import defaultdict

# ----------------------------------------------------------------------------
# 1. Check if gffread is installed

def check_gffread():
    exe = shutil.which("gffread")
    if exe is None:
        print("\n===[ERROR]=== gffread not found!")
        print("Please install gffread and ensure it is in your PATH.")
        print("Download: https://github.com/gpertea/gffread")
        print("via Conda: conda install bioconda::gffread")
        sys.exit(1) #exit out
    else:
        print(f"gffread found at: {exe}")

# ----------------------------------------------------------------------------
# _Helper function to unzip genome FASTA
def unzip_fasta(genome_fasta):
    """
    If genome FASTA provided is .gz, automatically decompress it
    gffread does not accept .gz
    """
    if genome_fasta.endswith(".gz"):
        print("[INFO] Decompressing gzipped FASTA......")

        tmpdir = tempfile.mkdtemp(prefix="genome_unzip_")
        out = os.path.join(tmpdir, os.path.basename(genome_fasta).replace(".gz",""))

        with gzip.open(genome_fasta,"rt") as infile, open(out,"w") as outfile:
            shutil.copyfileobj(infile,outfile)
        print(f"[INFO] Genome FASTA uncompressed to: {out}")
        return out
    return genome_fasta
# ----------------------------------------------------------------------------
# 2. Generate CDS and protein FASTA via gffread

def run_gffread(gff_file, genome_fasta, 
                out_cds = True, 
                out_protein=True, 
                out_prefix = "out"):
    """
    Run gffread to generate CDS FASTA(-x) and/or Protein FASTA(-y).

    Args:
        out_cds (bool): generate CDS FASTA (default:True)
        out_protein(bool): generate protein FASTA (default:True)
    
    Returns:
        (cds_out of None, protein_out or None)
    """
    check_gffread()

    if not out_cds and not out_protein:
        print("[ERROR]\n Neither out_cds nor out_protein is true.\n"
        "Nothing to generate. Aborting.")
        sys.exit(1)

    cds_out = f"{out_prefix}.cds.fasta" if out_cds else None
    protein_out = f"{out_prefix}.protein.fasta" if out_protein else None

    # Build gffread command

    genome = unzip_fasta(genome_fasta)
    cmd = ["gffread", gff_file, "-g", genome, "-F"]

    if out_cds:
        cmd += ["-x", cds_out]
    if out_protein:
        cmd += ["-y", protein_out]

    # Run it!
    print("[Processing.......]"," ".join(cmd))
    subprocess.run(cmd, check=True)

    if out_cds:
        print(f"CDS FASTA written: {cds_out}")
    if out_protein:
        print(f"Protein FASTA written: {out_protein}")

    return cds_out, protein_out

# ----------------------------------------------------------------------------
# 3. Read FASTA
def read_fasta(fn):
    """Returns dict: {header, sequence}"""
    seqs = {}
    header = None
    sequence = []

    with open(fn) as f:
        for l in f:
            l = l.rstrip()
            if not l:
                continue

            if l.startswith(">"):
                if header:
                    seqs[header] = "".join(sequence)
                header = l[1:].split()[0] # Remove ">" for dictionary keys
                sequence = []
            else:
                sequence.append(l.strip())

    if header:
        seqs[header] = "".join(sequence)

    return seqs

# ----------------------------------------------------------------------------
# 4. Parse cleaned GFF to extract hierarchy
def parse_cleaned_gff_hierarchy(gff_file,delimiter="\t",separator=";",assigner="="):
    """
    From cleaned protein_coding GFF (produced with Parse_whatUwant_gff.py)
    extract:
        - transcript_to_gene
        - transcript_to_protein_id
    
    To use in the final FASTA file.
    """
    transcript_to_gene = {}
    transcript_to_protein_id = defaultdict(set)

    def parse_attrs(attr_string):
        attrs = {}
        for kv in attr_string.split(separator):
            kv.strip()
            if not kv or assigner not in kv:
                continue
            key, val = kv.split(assigner,1)
            attrs[key.strip().lower()] = val.strip()
        return attrs
    
    with open(gff_file, "r") as gff:
        for l in gff:
            if not l.strip() or l.startswith("#"):
                continue

            fields = l.rstrip("\n").split(delimiter)
            if len(fields) < 9:
                continue

            ft = fields[2].lower()
            attrs = parse_attrs(fields[8])

            # mRNA/transceript lines: define transcript->gene
            if (ft.endswith("rna") or ft.endswith("transcript")):
                tid = attrs.get("id")
                gid = attrs.get("parent")
                if tid and gid:
                    transcript_to_gene[tid] = gid

            # CDS lines: define transcript->protein_id
            elif ft =="cds":
                tid = attrs.get("parent")
                pid = attrs.get("protein_id")
                if tid and pid:
                    transcript_to_protein_id[tid].add(pid)
    
    # Check Strict 1:1 transcript -> protein_id mapping
    transcript_to_protein = {}
    for tid, pidset in transcript_to_protein_id.items():
        if len(pidset) == 1:
            transcript_to_protein[tid] = next(iter(pidset))
        elif len(pidset) > 1:
            print(f"\n[ERROR] Multiple protein_id values found for {tid}.")
            print(f"Protein_ids: {pidset}")
            print("Please fix this issues before processing.")
            raise ValueError("Multiple protein_id entries for a single transcript.")

    return transcript_to_gene, transcript_to_protein

# ----------------------------------------------------------------------------
# 5. Select Representative isoform
def select_rep_isoform(cleaned_gff, protein_fasta, criteria = "longest",
                       out_fasta="representative_isoforms.faa", out_summary="isoform_summary.tsv"):
    """
    criteria: 
        - longest (Default)
        - random
        - first
    
    Requires:
        - cleaned_gff -> using Parse_whatUwant_gff.py
        - protein_Fast -> from run_gffread(..., out_protein=True)

    """

    # Step 1) Buld gene -> transcript -> protein_id
    transcript_to_gene, transcript_to_protein = parse_cleaned_gff_hierarchy(cleaned_gff)

    # Step 2) load protein sequences
    protein_seqs = read_fasta(protein_fasta)

    # Step 3) build gene -> transcripts
    gene_to_transcripts = defaultdict(list)
    for tid, gid in transcript_to_gene.items():
        if tid in protein_seqs: # only include transcripts with protein sequence
            gene_to_transcripts[gid].append(tid)
    
    # Step 4) select rep isoform per gene
    rep_seq = {}
    summary_rows = []

    for gid, tids in gene_to_transcripts.items():
        if criteria == "longest":
            best_tid = max(tids, key=lambda t: len(protein_seqs[t]))
        elif criteria == "random":
            best_tid = random.choice(tids)
        elif criteria == "first":
            best_tid = tids[0]
        else:
            raise ValueError("criteria supported is only: 'longest', 'random', or 'first'.")
        
        pid = transcript_to_protein.get(best_tid, "NA")
        seq = protein_seqs[best_tid]

        rep_seq[gid] = (best_tid,pid,seq)
        summary_rows.append([gid,best_tid,pid,len(seq)])

    # Step 5) write FASTA with only rep sequences
    with open(out_fasta, "w") as out:
        for gid, (tid, pid, seq) in rep_seq.items():
            out.write(f">{pid} gene={gid} transcript_id={tid} aa_len={len(seq)}\n")
            for i in range(0, len(seq), 60): # return every 60 sequences
                out.write(seq[i:i+60] + "\n")
    
    # Step 6) write summary table
    with open(out_summary, "w") as out:
        out.write("GeneID\tTranscriptID\tProteinID\tAA_length\n")
        for row in summary_rows:
            out.write("\t".join(map(str,row)) + "\n")

    print(f">>> Representative isoforms FASTA written: {out_fasta}")
    print(f">>> Summary table written: {out_summary}")

    return out_fasta, out_summary
