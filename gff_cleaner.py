"""
Utility functions for filtering and cleaning GFF/GTF annotation files
based on genome attributes
(e..g, 'chloroplast', 'mitochondrion').

Written to work alongside ParseGFFinfo class

Example:
    import ParseGFFinfo

    gff = ParseGFFinfo("thisisgff.gff3")
    parse_genome_info(gff.filepath)
    cleanup_by_genome_attr(gff.filepath,genome_to_exclude="chloroplast",outfile=True)
"""

def parse_genome_info(fn, delimiter="\t", separator=";", assigner="="):
        """
        Parse genome information of region, or if any
        detecting if a line has 'genome=' in attributes.

        Args:
            fn (str) : GFF/GTF file path
            delimiter (str): column delimiter (default: tab'\\t')
            separator (str): attribute separator (default: ';')
            assigner (str) : key-value assigner for each attribute (default= '=')

        Return summary of 
        - number of lines that has 'genome'
        - Feature types(column 3) having genome in attributes
        - values of genome attribute
        """
        genome = set()
        genome_ft = set()
        genome_count= 0
        
        with open(fn, "r") as gff:
            for l in gff:
                if l.startswith("#"):
                     continue
                fields = l.rstip().split(delimiter)
                if len(fields) < 9:
                    continue
                attrs = {
                        kv.split(assigner,1)[0].strip().lower(): kv.split(assigner,1)[1].strip()
                        for kv in fields[8].split(separator) if assigner in kv
                        }
                if "genome" in attrs:
                    genome.add(attrs["genome"])
                    genome_ft.add(fields[2])
                    genome_count +=1

        print(f"# of lines containing genome : {genome_count}")
        print(f"Feature types having 'genome' in attributes: {genome_ft}")
        print(f"Values found for genome= : {genome}")

        return {
             "count": genome_count,
             "feature_types": genome_ft,
             "values": genome
        }

# To clean up gff exluding particular regions having given genome= attributes
# found using the "ParseGFFinfo.parse_genome_info"
def cleanup_by_genome_attr(
        fn, 
        genome_to_exclude=None, 
        ft_4_genome_attr="region", 
        outfile=True, 
        delimiter = "\t", 
        separator = ";", 
        assigner="="):
    """
    Exclude any region whose genome attribute matches given values
    (e.g. chroloplast)
    AND
    all features belongs to that region based on column 1 (seqID)

    Args:
        fn (str): annotation(GFF/GTF) filepath
        genome_to_exclude (string or list): 
                            add value(s) of genome attributes to exclude
                            (found with 'parse_genome_info')
                            e.g., mitochondrion
        ft_4_genome_attr (str):
            - default: region
            - add if any other feature types found parse_genome_info
        outfile (Bool or str):
            - True(Default): cleaned gff will be saved as original <fn>_cleaned.gff
            - False: print genome values excluded and their seq IDs, not write file
            - str: cleaned gff will be named as given str
        delimiter (str):
                    column delimiter (default = tab)
        separator (str):
                    attribute separator for each key-value pair
                    (default = ';')
        assigner (str):
                    Key-value assigned in attributes 
                    (Default = '=')
    
    Returns:
        set: excluded sequence IDs (regions removed).
                  
    """
    if genome_to_exclude is None:
        genome_to_exclude = []
    if isinstance(genome_to_exclude, str):
        genome_to_exclude = [genome_to_exclude] # make it list if one string to work below

    exclude_set = {e.lower() for e in genome_to_exclude}
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
            if ft != ft_4_genome_attr:
                continue
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
    return excluded_seqids

