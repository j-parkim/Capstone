from collections import defaultdict

class ParseGFFinfo(object):
    """
    Parse unique values in following GFF columns:
    .source      : column 2
    .featuretype : column 3
    .attr        : column 9 (Attributes)

    Default arguments follow standard RefSeq GFF3:
        - delimiter: "\\t" (tab)
        - separator = ";"
        - assigner = "="   
    """
    
    def __init__(self, filepath, delimiter = "\t", separator = ";", assigner = "="):
        """
        Args:
        - filepath (str): annotation file path
        - delimiter (str): column delimiter 
                        (Default = tab ("\\t"))
        - separator (str): attribute separator for each key-value pair
                        (Default = ';')
        - assigner (str): Key-value assigned in attributes 
                        (Default = '=')
        """

        self.filepath = filepath
        self.delimiter = delimiter
        self.separator = separator
        self.assigner = assigner
        self.lines = []

        # Read file once and keep valid (non-comment) lines
        with open(filepath, 'r') as gff:
            for line in gff:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split(self.delimiter)
                if len(fields) >= 2:
                    self.lines.append(fields)
        if not self.lines:
            print("!!! Check the file. No valid lines found.")

    
    # parse sources
    def source(self):
        """Return unique values of column 2 (source)"""
        return {f[1] for f in self.lines if len(f)>= 2}
    
    # parse feature types
    def featuretype(self):
        """Return unique values of column 3 (feature type)"""
        return {f[2] for f in self.lines if len(f)>= 3}

    # dictionary for feature type(s) for each source
    def featuretype_by_source(self):
        """Return dictionary mapping each source to unique feature type(s)"""
        ft_by_source = defaultdict(set)
        for f in self.lines:
            if len(f) >= 3:
                ft_by_source[f[1]].add(f[2])
        return dict(ft_by_source)
    
    # parse set of unique attributes for a given featuretype
    def attr(self, featuretype="gene"):
        """
        Return set of unique attributes for a given feature type.
        Default featuretype = 'gene'
        """
        keys = set() # to save all unique keys
        featuretype = featuretype.lower()
        for f in self.lines:
            if len(f) < 9: 
                continue # if there's no attribute column, skip
            if f[2].lower() != featuretype: 
                continue
            attributes = f[8] # get the 9th column(attributes in general gff)
            for attr in attributes.split(self.separator):
                attr = attr.strip()
                if not attr: 
                    continue
                if self.assigner in attr:
                    key, value = attr.split(self.assigner, 1)
                    keys.add(key.strip())
                else:
                    keys.add(attr.strip())           
        return keys
    
    # Parse genome information of region, or if any
    def parse_genome_info(self):
        """
        Parse genome information of region, or if any
        detecting if a line has 'genome=' in attributes.
        Return summary of 
        - number of lines that has 'genome'
        - Feature types(column 3) having genome in attributes
        - values of genome attribute
        """
        genome = set()
        genome_ft = set()
        genome_count= 0
        
        for f in self.lines:
            if len(f) < 9:
                continue
            attrs = {
                    kv.split(self.assigner,1)[0].strip().lower(): kv.split(self.assigner,1)[1].strip()
                    for kv in f[8].split(self.separator) if self.assigner in kv
                    }
            if "genome" in attrs:
                genome.add(attrs["genome"])
                genome_ft.add(f[2])
                genome_count +=1

        print(f"# of lines containing genome : {genome_count}")
        print(f"Feature types having 'genome' in attributes: {genome_ft}")
        print(f"Values for genome= attributes: {genome}")

# To clean up gff exluding particular regions having given genome= attributes
# found using the "ParseGFFinfo.parse_genome_info"
def cleanup_by_genome_attr(fn, genome_to_exclude=None, ft_4_genome_attr="region", outfile=True, delimiter = "\t", separator = ";", assigner="="):
    """
    Exclude any region whose genome attribute matches given values
    (e.g. chroloplast)
    AND
    all features belongs to that region based on column 1 (seqID)

    Args:
        -fn: annotation(GFF/GTF) filepath
        -genome_to_exclude: string or list 
                            add value(s) of genome attributes found ParseGFFinfo.parse_genome_info
                            e.g., mitochondrion
        -ft_4_genome_attr: str
            - region(default)
            - add if any other feature types found ParseGFFinfo.parse_genome_info
        -outfile: Bool or str
            - True(Default): cleaned gff will be saved as original fn + _cleaned.gff
            - False: print genome values excluded and their seq IDs
            - str: cleaned gff will be named as given str
        -delimiter: str
                    column delimiter (default = tab)
        -separator: str
                    attribute separator for each key-value pair
                    (Default = ';')
        -assigner: str
                    Key-value assigned in attributes 
                    (Default = '=')
                  
    """
    if genome_to_exclude is None:
        genome_to_exclude = []
    if isinstance(genome_to_exclude, str):
        genome_to_exclude = [genome_to_exclude] # make it list if one string to work below

    exclude_set = {e.lower() for e in genome_to_exclude}
    excluded_seqids = set()

    # Identify region seqIDs to exclude
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
