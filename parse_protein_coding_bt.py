from collections import defaultdict
       
# To confirm biotype=protein_coding is in the gff, 
# only in "gene" feature type or if anything else somehow,
# and sources that contain protein_coding biotype attributes 
class biotype_protein_coding(object):
    """
    Parse sources and feature types
    that has "protein_coding" biotype 
    in attributes

    Then, 
    for "protein_coding" biotype genes,
    parse their children(mRNA) and grandchildren(cds)
    """

    def __init__(self, filepath):
        self.filepath = filepath

    # Parse 1) all the sources that have protein coding biotype,
    # 2) "protein_coding" lines describing a gene
    # 3) non-gene "protein_coding" lines if any
    def source(self):
        source = set()
        protein_coding_lines = []
        non_gene_protein_coding_lines = []
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                if "protein_coding" in l:
                    fields = l.rstrip().split()
                    if fields[2].lower() == "gene":
                        source.add(fields[1])
                    else:
                        non_gene_protein_coding_lines.append(l)
                        print(f"Non-gene protein_coding line (feature={fields[2]}):\n{l.strip()}")
                    protein_coding_lines.append(l)
        
        if not protein_coding_lines:
            print("No protein_coding biotype in the file.")
        
        return source, protein_coding_lines, non_gene_protein_coding_lines
            
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
                        print("Check the file, there's no ID attributes")
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
                    if "=" in kv}
                gff_lines.append((fields, attrs, l))

                if ft in {"mrna", "transcript"} and "parent" in attrs:
                    parent_gene = attrs["parent"]
                    if parent_gene in protein_coding_genes and "id" in attrs:
                        transcript_id = attrs["id"]
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
                outfn = "complete_protein_coding.gff"
        
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
