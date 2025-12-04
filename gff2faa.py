"""
Pipeline:
    0) Inspect GFF with *ParseGFFinfo.py*
        0.5) (Optional) Reformat dialects using *ReformatGFF.py*
    1) (Optional) *Parse_whatUwant_gff.py*: Clean up GFF (remove genome regions, sources, only keep protein_coding biotypes)
    2) *Parse_ONE_Isofom.py*: run gffread to produce CDS + protein FASTA 
    3) *Parse_ONE_Isofom.py*: Select one representative isoform per gene
                            (longest, random, first)

Author: Jeonghwa (J) Parkim
"""

import os
import argparse
import ParseGffinfo
import Reformat2GFF
import Parse_whatUwant_gff
import Parse_ONE_Isoform

# ------------------------------------------------------------------------
# Helper to store the final output files
# ------------------------------------------------------------------------
def ensure_outdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return os.path.abspath(path)

# ------------------------------------------------------------------------
# MAIN Command Line Interface
# ------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="GFF inspection "
        "+ GFF cleaning "
        "+ FASTA extraction "
        "+ Representative Isoform Selection",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-gff", required=True,
                        help="Input GFF3 annotation filepath")
    parser.add_argument("-genome", required=True,
                        help="Input Genome FASTA file for gffread")
    
    # ------ GFF format inspection / reformatting ------
    parser.add_argument("--reformat", action="store_true",
                        help="Run ReformatGff to standardize attributes")
    
    # ------ Cleanup -------
    parser.add_argument("--remove-genome", nargs="*", default=None,
                        help="Genome attribute(s) to remove (e.g., 'chloroplast', 'mitochondrion')")
    parser.add_argument("--remove-source",nargs="*",default=None,
                        help="Sources (column2) to remove (e.g., cmsearch)")
    parser.add_argument("--skip-cleanup", action="store_true",
                        help="Skip cleanup steps and use the original GFF")
    
    # ------ gffread ------
    parser.add_argument("--skip-gffread", action="store_true",
                        help="Skip gffread(must supply --protein-fasta manually)")
    parser.add_argument("--protein-fasta", default=None,
                        help="If --skip-gffread, MUST provide protein fasta file")
    
    # ------ Isoform selection -------
    parser.add_argument("--criteria", default="longest",
                        choices=["longest","first","random"],
                        help="Criteria for representative isoform selection")
    parser.add_argument("--skip-isoform", action="store_true",
                        help="Skip representative isoform selection, keep all sequences.")
    
    # ------ Output Directory ------
    parser.add_argument("--outdir", default="results",
                        help="Directory to store all output files.")
    
    args = parser.parse_args()
    outdir = ensure_outdir(args.outdir)

    # ------------------------------------------------------
    # 0) Inspect GFF with *ParseGFFinfo.py*
    # ------------------------------------------------------
    print("\n=== STEP 0: Inspecting GFF ===")
    inspector = ParseGffinfo.ParseGFFinfo(args.gff)
    inspector.source()
    inspector.featuretype()
    inspector.genome()
    inspector.show_attr_by_source_featuretype()

    # Add user prompt to enter any source/genome to remove
    # learned from the summary above
    remove_g = input("\nDo you want to remove any genome? \n (Enter comma-separated or 'none')").strip()
    remove_s = input("\nDo you want to remove data from any source? \n(Enter comma-separated or 'none)").strip()

    remove_genome_attr = [] if remove_g.lower() =="none" else [x.strip() for x in remove_g.split(",") if x.strip()]
    remove_sources = [] if remove_s.lower() == "none" else [x.strip() for x in remove_s.split(",") if x.strip()]
    # ------------------------------------------------------
    # 0.5) (Optional) Reformat dialects using *ReformatGFF.py*
    # ------------------------------------------------------
    gff_for_next = args.gff

    if args.reformat:
        print("\n=== STEP 0.5: Reformatting GFF ===")
        rf = Reformat2GFF.ReformatGFF(args.gff,
                                      remove_genome_attr=remove_genome_attr,
                                      remove_sources=remove_sources)
        rf.detect_format()
        out_reformed = os.path.join(outdir, "reformatted.gff3")
        rf.reformat_dialects("GFF3", outfile=out_reformed)
        gff_for_next = out_reformed
        print(f"[INFO] Reformatted GFF written to: {gff_for_next}")

    # ------------------------------------------------------
    # 1) (Optional) Clean up GFF 
    #    (remove genome regions, sources, only keep protein_coding biotypes)
    # ------------------------------------------------------
    if args.skip_cleanup:
        print("\n=== STEP 1: Not cleaning up (skipping) ===")
        cleaned_gff = gff_for_next
    else:
        print("\n=== STEP 1: Cleaning up in progress ===")
        cleaned_gff_temp = Parse_whatUwant_gff.full_cleanup(
            gff_for_next,
            genome_to_remove=args.remove_genome,
            source_to_remove=args.remove_source,
            protein_coding_only=True
        )
        # move output to outdir
        cleaned_gff = os.path.join(outdir, os.path.basename(cleaned_gff_temp))
        os.rename(cleaned_gff_temp, cleaned_gff)
    
    print(f"[INFO] Cleaned GFF: {cleaned_gff}")

    # ------------------------------------------------------
    # 2) run gffread to produce CDS + protein FASTA 
    # ------------------------------------------------------    
    if args.skip_gffread:
        print("\n=== STEP 2: gffread skipped ===")
        if not args.protein_fasta:
            raise ValueError("gffread SKIPPED, but --protein-fasta not provided. ")
        protein_fasta = args.protein_fasta
        print(f"[INFO] User provided protein FASTA {protein_fasta} will be used.")
    
    else:
        print("\n=== STEP 2: 'gffread' generating protein FASTA ===")
        cds_out, protein_out = Parse_ONE_Isoform.run_gffread(
            cleaned_gff,
            args.genome,
            out_prefix=os.path.join(outdir,"gffread_output")
        )
        protein_fasta = protein_out
        print(f"[INFO] Protein FASTA genearated: {protein_fasta}")

    # ------------------------------------------------------
    # 3) Select one representative isoform per gene
    #           (longest, random, first)
    # ------------------------------------------------------
    if args.skip_isoform:
        print(f"=== STEP 3: Representative Isoform not selected (SKIPPED) ===")
        print("[DONE]")
        return
    
    else:
        print(f"\n=== STEP 3: Generating FASTA with representative sequences only ===")

        out_faa = os.path.join(outdir,"representative_isoforms.faa")
        summary_tsv = os.path.join(outdir,"isoform_written_summary.tsv")

        Parse_ONE_Isoform.select_rep_isoform(
            cleaned_gff, protein_fasta, criteria=args.criteria,
            out_fasta=out_faa, out_summary=summary_tsv
        )

        print("\n=========== COMPLETE ============")
        print(f"Output Directory: {outdir}")
        print(f"- Cleaned GFF: {cleaned_gff}")
        print(f"- Protein FASTA: {protein_fasta}")
        print(f"- Representative FASTA: {out_faa}")
        print(f"- Selection Summary: {summary_tsv}")
        print("==============================================")

if __name__ == "__main__":
    main()
