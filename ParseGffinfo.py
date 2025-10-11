from collections import defaultdict

class ParseGFFinfo(object):
    """
    Parse unique values in following GFF columns:
    .source      : column 2
    .featuretype : column 3
    .attr        : attributes in column 9
                    for attributes, default arguments are
                    following standard RefSeq GFF3 format
                    - featuretype = "gene"
                    - separator = ";"
                    - assigner = "="   
    """
    
    def __init__(self, filepath):
        self.filepath = filepath
        
    
    # parse sources
    def source(self):
        source = set()
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                fields = l.rstrip().split()
                if len(fields) >=2:
                    source.add(fields[1])
                else:
                    print("Check the file. There is only one column.")
        return source
    
    # parse feature types
    def featuretype(self):
        features = set()
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                fields = l.rstrip().split()
                if len(fields) >= 3:
                    features.add(fields[2])
                else:
                    print("Line does not have 3 or more columns.")
        return features

    # dictionary for feature type(s) for each source
    def featuretype_by_source(self):
        ft_by_source = defaultdict(set)
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                fields = l.rstrip().split()
                if len(fields) >= 3:
                    ft_by_source[fields[1]].add(fields[2])
                else:
                    print("Line does not have 3 or more columns")
        return dict(ft_by_source)
    
    # parse attributes for assigned featuretype
    def attr(self, featuretype="gene", separator = ";", assigner = "="):
        keys = set() # to save all unique keys
        featuretype = featuretype.lower()
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue                
                fields = l.rstrip().split()
                if len(fields) < 9: 
                    continue # if there's no attribute column, skip
                if fields[2].lower() != featuretype: 
                    continue
                attributes = fields[8] # get the 9th column(attributes in general gff)
                for attr in attributes.split(separator):
                    attr = attr.strip()
                    if not attr: 
                        continue
                    if assigner in attr:
                        key, value = attr.split(assigner, 1)
                        keys.add(key.strip())
                    else:
                        keys.add(attr.strip())
               
        return keys
    
    # Parse genome information of region, or if any
    def parse_genome_info(self):
        genome = set()
        genome_ft = set()
        genome_count= 0
        with open(self.filepath, 'r') as gff:
            for l in gff:
                if l.startswith("#"):
                    continue
                fields = l.rstrip().split("\t")
                if len(fields) < 9:
                    continue
                attrs = {
                        kv.split("=")[0].strip().lower(): kv.split("=",1)[1].strip()
                        for kv in fields[8].split(";") if "=" in kv
                    }
                if "genome" in attrs:
                    genome.add(attrs["genome"])
                    genome_ft.add(fields[2])
                    genome_count +=1

        print(f"# of Lines contatining genome : {genome_count}")
        print(f"Feature types that have genome in attributes: {genome_ft}")
        print(f"Keys for genome= attributes: {genome}")

# To clean up gff exluding particular regions having given genome= attributes
# found using the "ParseGFFinfo.parse_genome_info"
def cleanup_by_genome_attr(fn, genome_to_exclude=None, ft_4_genome_attr="region", outfile=True):
    """
    Exclude any region whose genome attribute matches given values
    (e.g. chroloplast)
    AND
    all features belongs to that region
    based on column 1 (seqID)

    Args:
        -fn: gff filepath
        -genome_to_exclude: string or list 
                            add value(s) of genome attributes found .parse_genome_info
                            e.g., mitochondrion
        -ft_4_genome_attr: str
            - regions(default)
            - add if any other feature types found .parse_genome_info
        -outfile: Bool or str
            - True(Default): cleaned gff will be saved as original fn + _cleaned.gff
            - False: print genome values excluded and their seq IDs
            - str: cleaned gff will be named as given str
                  
    """
    if genome_to_exclude is None:
        genome_to_exclude = []
    exclude_set = {e.lower() for e in genome_to_exclude}

    excluded_seqids = set()

    # Identify region seqIDs to exclude
    with open(fn,'r') as gff:
        for l in gff:
            if l.startswith("#"):
                continue
            fields = l.rstrip().split("\t")
            if len(fields) <9:
                continue
            ft = fields[2].lower()
            if ft != ft_4_genome_attr:
                continue
            attrs = {
                kv.split("=")[0].strip().lower(): kv.split("=",1)[1].strip()
                for kv in fields[8].split(";") if "=" in kv
            }
            if "genome" in attrs and attrs["genome"].lower() in exclude_set:
                excluded_seqids.add(fields[0])
    
    print(f"Found {len(excluded_seqids)} seq IDs were excluded with genome={exclude_set}")

    outfn = None
    if outfile:
        if isinstance(outfile, str) and outfile not in [True,False]:
            outfn = outfile
        else:
            outfn = fn + "_cleaned.gff"
    
    out = open(outfn, 'w') if outfn else None

    with open(fn,'r') as gff:
        for l in gff:
            if l.startswith("#"):
                if out:
                    out.write(l)
                continue
            fields = l.rstrip().split("\t")
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

    return excluded_seqids



# Parse different feature types and attributes between files.

def find_diff_attributes(*inputs, outfile = False):
    """
    input argument format : 
    ("label for the gff1(e.g., rice)", filepath1) tuples

    outfile :
    False(Default): output printed to screen 
    True: A text file generated with output.

    Parse common and different feature types btw multiple GFFs

    Parse common and different attributes for each feature type btw multiple GFFs
    
    """
    # assign labels for the final output
    gffs = {label: ParseGFFinfo(path) for label, path in inputs}
    labels = list(gffs.keys())

    # if outfile is set, open file for writing.
    if outfile:
        filename = f"gff_comparison_{str(labels)}.txt"
        out = open(filename, "w")

        def write(*args, **kwargs):
            print(*args, **kwargs, file=out)
    else:
        def write(*args, **kwargs):
            print(*args, **kwargs)

    # Find distinctive feature type
    features = {label:g.featuretype() for label, g in gffs.items()}
    common_f = set.intersection(*features.values())

    write(f"Common feature types in all files: {sorted(common_f)}")
    write("!!!Unique features found!!!")
    for label, ft in features.items():
        unique_ft = ft - common_f
        if unique_ft:
            write(f"Only {label} has {sorted(unique_ft)}")
    write()

    # Find distinctive attributes for each feature type
    for ft in common_f:
        attrs = {label:g.attr(featuretype=ft) for label, g in gffs.items()}
        common_a = set.intersection(*attrs.values())
        all_a = set.union(*attrs.values())

        write(f"Feature type: {ft}")

        if all_a != common_a:
            write("!!!Different attributes found!!!")
            for label, attr in attrs.items():
                unique_a = attr - common_a
                if unique_a:
                    write(f"Attributes only in {label}: {sorted(unique_a)}")
        else:
            write(f"All files have the same attributes for {ft}")

        write(f"Number of common attribues for {ft}: {len(common_a)}")
        write(f"Common attributes for {ft}: {sorted(common_a)}\n")
        
    if outfile:
        out.close()
        print(f"Output written to {filename}")
