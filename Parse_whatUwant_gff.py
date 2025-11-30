"""
Utility functions for filtering and cleaning GFF/GTF annotation files

Perform:
    - Remove regions by genome attributes (e..g, 'chloroplast', 'mitochondrion').
    - Remove entries by source (column2)
    - Retain only protein_coding biotype genes with complete hierarchy(gene-mRNA-CDS)
    - Combine both, producing cleaned GFFs

Written to work alongside ParseGFFinfo class

Example:
    import ParseGffinfo

    gff = ParseGFFinfo("thisisgff.gff3")
    parse_genome_info(gff.filepath)
    cleanup_by_genome_attr(gff.filepath,genome_to_exclude="chloroplast",outfile=True)
"""
from collections import defaultdict
from ParseGffinfo import ParseGFFinfo

# -------------------------------------------------------------------------
# # Clean up gff removing particular regions having given genome= attributes
# found using the "ParseGFFinfo.genome"
def cleanup_by_genome_attr(
        fn, 
        genome_to_exclude=None, 
        ft_4_genome_attr="region", 
        outfile=True, 
        delimiter = "\t", 
        separator = ";", 
        assigner="="):
    """
    Remove any regions whose genome attribute matches 'genome_to_exclude'
    (e.g. chroloplast)
    Remove the entire features belonged to that region based on column 1 (seqID)

    Args:
        fn (str): annotation(GFF/GTF) filepath
        genome_to_exclude (string or list): value(s) of genome attributes to remove
                            (found with 'ParseGFFinfo.genome', e.g., 'mitochondrion')
        ft_4_genome_attr (str):
            - default: region
            - add if any other feature types for genome annatation
        outfile (Bool or str):
            - True(Default): cleaned gff will be saved as original <fn>_cleaned.gff
            - False: print genome values excluded and their seq IDs, not write file
            - str: cleaned gff will be named as given str
        delimiter (str): column delimiter (default = tab)
        separator (str): attribute separator for each key-value pair (default = ';')
        assigner (str): Key-value assigned in attributes Default = '=')
    
    Returns:
        output file name.                  
    """
    if genome_to_exclude is None:
        genome_to_exclude = []
    if isinstance(genome_to_exclude, str):
        genome_to_exclude = [genome_to_exclude] # make it list if one string to work below

    exclude_set = {g.lower() for g in genome_to_exclude}
    excluded_seqids = set()

    # Identify seqIDs of regions to exclude
    with open(fn,'r') as gff:
        for l in gff:
            if l.startswith("#"):
                continue
            fields = l.rstrip().split(delimiter)
            if len(fields) <9:
                continue
            ft = fields[2].lower()

            # Only look at given feature type (like 'region')
            if ft != ft_4_genome_attr:
                continue

            # Parse attributes
            attrs = {
                kv.split(assigner,1)[0].strip().lower(): kv.split(assigner,1)[1].strip()
                for kv in fields[8].split(separator) if assigner in kv
            }
            if "genome" in attrs and attrs["genome"].lower() in exclude_set:
                excluded_seqids.add(fields[0])
    
    print(f"Found {len(excluded_seqids)} seq IDs were excluded with genome={exclude_set}")

    # Name output file of cleaned gff
    outfn = None
    if outfile:
        if isinstance(outfile, str) and outfile not in [True,False]:
            outfn = outfile
        else:
            outfn = fn + "_cleaned.gff"
    
    out = open(outfn, 'w') if outfn else None

    # Write or print filtered results
    with open(fn,'r') as gff:
        for l in gff:
            if l.startswith("#"):
                if out:
                    out.write(l)
                continue
            fields = l.rstrip().split(delimiter)
            if len(fields) <9:
                continue
            seqid = fields[0]
            if seqid in excluded_seqids:
                continue
            if out:                   
                out.write(l)
    
    if out:
        out.close()
        print(f">>>Cleaned file written: {outfn}")
    else:
        print("Output file not written as outfile=False.")

    print(f"SeqIDs for excluded regions: {excluded_seqids}")
    return outfn

# -------------------------------------------------------------------------
# Remove lines from particular source (Column2)
def cleanup_by_source(fn, sources_to_exclude, delimiter = "\t", outfile=True):
    """
    Remove lines where column2 (source) matches sources_to_exclude
    
    Args:
        fn(str) : input GFF file
        sources_to_exluce (list or str): unwanted data sources to remove
        delimiter (str): column delimiter (default = tab)
        outfile (Bool or str):
            - True(Default): cleaned gff will be saved as original <fn>_sourcecleaned.gff
            - False: print
            - str: cleaned gff will be named as given str
        """
    if isinstance(sources_to_exclude, str):
        sources_to_exclude = [sources_to_exclude]

    exclude_set = set(sources_to_exclude)
    removed = 0

    outfn = None
    if outfile:
        if isinstance(outfile, str) and outfile not in [True,False]:
            outfn = outfile
        else:
            outfn = fn + "_sourcecleaned.gff"
    
    out = open(outfn, 'w') if outfn else None

    # Write or print filtered results
    with open(fn,'r') as gff:
        for l in gff:
            if l.startswith("#"):
                if out:
                    out.write(l)
                continue
            fields = l.rstrip().split(delimiter)
            if len(fields) <2:
                continue
            if fields[1] in exclude_set:
                removed += 1
                continue
            if out:
                out.write(l)
    
    if out:
        out.close()
        print(f">>>Cleaned file written: {outfn}")
    else:
        print("Output file not written as outfile=False.")

    print(f"Sources cleaned: {removed} lines matching {exclude_set}")
    return outfn

# -------------------------------------------------------------------------
# Biotype = Protein_coding gene filter
# To confirm biotype=protein_coding is in the gff, 
# only in "gene" feature type or if anything else somehow,
# and sources that contain protein_coding biotype attributes 
class biotype_protein_coding(object):
    """
    Parse sources and feature types that has "protein_coding" biotype in attributes

    Then, keep the complete protein_coding hierarchy:
        gene (protein_coding) (parent)
            -> mRNA/transcript (children)
                -> CDS (grandchildren)
    """

    def __init__(self, filepath):
        self.filepath = filepath
            
    # Parse protein coding genes
    def parse_proteincoding_genes(self):
        """
        Parse protein coding genes only
        """
        genes = set()
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                fields = l.rstrip().split()
                if len(fields) < 9:
                    continue
                if "protein_coding" in fields[8].lower() and fields[2].lower() == "gene":
                    attrs = {
                        kv.split("=")[0].strip().lower(): kv.split("=",1)[1].strip() 
                        for kv in fields[8].split(";") 
                        if "=" in kv}
                    if "id" in attrs:
                        genes.add(attrs["id"])
                    else:
                        print("Check the file, there's no ID attributes.")
        return genes

    # parse childen for protein coding genes
    # gene (ID=gene-LOC1)
    #   - mRNA(ID=rna-XM_1;Parent=gene-LOC120112099)
    #       - exon(ID=exon-XM_1-1;Parent=rna-XM1)
    #       - CDS(ID=cds-XP_1;Parent=rna-XM1)
    def gene_to_children(self, outfile = False):
        """
        Parse children for protein coding genes

        Provide summary table at the end with counts of;
        - total protein coding genes
        - genes missing mRNA/transcript
        - genes missing CDS
        - genes with complete hierarchy

        Arg:
        outfile (bool or str):
            - False (default): print summary only
            - True : create a new GFF file(auto-named) 
                    only with complete hierarchy protein coding genes
                    (genes - mRNA - CDS)
        
        """
        protein_coding_genes = self.parse_proteincoding_genes()

        gene_to_transcripts = defaultdict(set)
        transcript_to_gene = {}
        transcript_to_cds = defaultdict(bool)
        gff_lines = [] # this is to parse when writing the filtered file

        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                fields = l.rstrip().split()
                if len(fields) < 9:
                    continue
                ft = fields[2].lower()
                attrs = {
                    kv.split("=")[0].strip().lower(): kv.split("=",1)[1].strip() 
                    for kv in fields[8].split(";") 
                    if "=" in kv
                    }
                gff_lines.append((fields, attrs, l))

                if ft in {"mrna", "transcript"} and "parent" in attrs:
                    parent_gene = attrs["parent"]
                    if parent_gene in protein_coding_genes:
                        if "id" in attrs:
                            transcript_id = attrs["id"]
                        else:
                            print(f"Check the file. ID missing. {l.strip()}")
                            continue
                        gene_to_transcripts[parent_gene].add(transcript_id)
                        transcript_to_gene[transcript_id] = parent_gene
                
                elif ft == "cds" and "parent" in attrs:
                    transcript_to_cds[attrs["parent"]] = True
        
        # Summary table
        total_genes = len(protein_coding_genes)
        genes_missing_transcripts = 0
        genes_missing_cds = 0
        genes_complete = set()
        
        for g in protein_coding_genes:
            transcripts = gene_to_transcripts.get(g,[])
            if not transcripts:
                genes_missing_transcripts += 1
                continue
            
            has_cds = any(transcript_to_cds[t] for t in transcripts)
            if not has_cds:
                genes_missing_cds += 1
            else:
                genes_complete.add(g)

        summary = {
            "Total protein coding genes": total_genes,
            "Genes missing transcripts": genes_missing_transcripts,
            "Genes missing CDS": genes_missing_cds,
            "Genes with complete children": len(genes_complete)
        }

        # Write filtered GFF 
        if outfile:
            if isinstance(outfile, str) and outfile.lower() not in {"true","false"}:
                outfn = outfile
            else:
                outfn = self.filepath + "_complete_protein_coding.gff"
        
            gene_count = 0
            transcript_count = 0
            cds_count = 0

            with open(outfn,'w') as out:
                # Grab the header lines from the original gff
                with open(self.filepath, 'r') as gff:
                    for l in gff:
                        if l.startswith("#"):
                            out.write(l)
                
                for fields, attrs, original_l in gff_lines:
                    ft = fields[2].lower()

                    if ft == "gene" and "id" in attrs and attrs["id"] in genes_complete:
                        out.write(original_l)
                        gene_count += 1
                    
                    elif ft in {"mrna", "transcript"} and "parent" in attrs:
                        if attrs["parent"] in genes_complete:
                            out.write(original_l)
                            transcript_count += 1
                    
                    elif ft == "cds" and "parent" in attrs:
                        parent_mrna = attrs["parent"]
                        if parent_mrna in transcript_to_gene and transcript_to_gene[parent_mrna] in genes_complete:
                            out.write(original_l)
                            cds_count += 1
            summary.update({
                "Genes written": gene_count,
                "mRNA/Transcripts written": transcript_count,
                "CDSs written": cds_count,
                "Output file": outfn
            })

        print("\n====Protein Coding Gene Summary====")
        for k,v in summary.items():
            print(f"{k:<35}: {v}")
        return outfn

# -------------------------------------------------------------------------
# Full clean up pipeline
def full_cleanup(
        gff_path,
        genome_to_remove = None,
        source_to_remove = None,
        protein_coding_only = True,
        suffix = "_cleaned"
    ):
    """
    Pipeline:
        1) Detect available genome values & sources via ParseGFFinfo
        2) Remove unwanted genome regions
        3) Remove unwanted sources
        4) Keep only protein coding genes with its children and grandchildren
    
    """
    print("\n=== Parsing file with ParseGFFinfo ===")
    parser = ParseGFFinfo(gff_path)
    genomes_found = parser.genome()["values"]
    sources_found = parser.source()

    print(f">Available genome values:{genomes_found}")
    print(f">Available sources : {sources_found}")

    # --- 1. Remove SeqIDs with unwanted genome region ---
    step1 = gff_path
    if genome_to_remove:
        step1 = cleanup_by_genome_attr(step1, genome_to_exclude=genome_to_remove,
                                       outfile=True)
    
    else:
        print("[FULL CLEANUP] SKIPPING genome cleanup.")

    
    # --- 2. Remove sources ---
    step2 = step1
    if source_to_remove:
        step2 = cleanup_by_source(step2, sources_to_exclude=source_to_remove,
                                  outfile=True)
    else:
        print("[FULL CLEANUP] SKIPPING source cleanup.")

    # --- 3. Protein coding gene only ---
    step3 = step2
    if protein_coding_only:
        pcg = biotype_protein_coding(step2)
        step3 = pcg.gene_to_children(outfile=True)

    print("\n=======Full Cleanup Complete========")
    print(f"Input file:               {gff_path}")
    print(f"After genome cleanup:     {step1}")
    print(f"After source cleanup:     {step2}")
    print(f"Protein coding output:    {step3}")
    print("=======================================")

    return step3