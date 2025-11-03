import random
import re
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
            - delimiter (str): column delimiter (Default = tab ("\\t"))
            - separator (str): attribute separator for each key-value pair (Default = ';')
            - assigner (str): Key-value assigned in attributes (Default = '=')
        """

        self.filepath = filepath
        self.delimiter = delimiter
        self.separator = separator
        self.assigner = assigner
        self.format = "Unknown"
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

    def detect_attr_format(self, max_lines=100, apply=True):
        """
        This is to detect likely format that attribute string is using.
        
        Args:
            - max_lines (int) : number of lines to look at for format detection
            - apply (bool) : if True, update self.separator, self.assigner, self.format

        Returns:
            - dictionary: detected likely format of attributes
        """
        
        if not self.lines:
            print("No lines loaded. Please check the file or initialization")
            return
        
        sep_candidates = set()
        assign_candidates = set()
        any_quotes = False
        subformats = set()
        unknown_sep_count = 0
        unknown_asgn_count = 0
        bad_lines =[]
        sample_lines = []

        lines_to_check = random.sample(self.lines, min(max_lines, len(self.lines)))
        for f in lines_to_check:
            if len(f) < 9:
                continue
            attrs = f[8].strip()
            if not attrs:
                continue

            # Check separator
            found_sep = False
            if ";" in attrs:
                sep_candidates.add(";")
                found_sep = True

            if "," in attrs:
                # comma is tricky. 
                # First, check if comma is found inside quotes (GTFlike)
                if re.search(r'"[^"]*,[^"]*"', attrs): 
                    subformats.add("',' inside quoted value")
                        
                # Second, check if GFF3 standard format, then likely comma is subseparator
                elif ";" in sep_candidates:
                    if re.search(r';[^;]*,[^;]*;', attrs):
                        subformats.add("',' inside values with separator';'")
                        
                # Otherwise, assume comma as a separator, but adding comment in subformat
                else:
                    sep_candidates.add(",")
                    subformats.add("',' as main separator (PLEASE CONFIRM)")
                found_sep = True

            if "\t" in attrs:
                sep_candidates.add("\t")
                found_sep = True

            if not found_sep:
                unknown_sep_count += 1
                bad_lines.append(self.delimiter.join(f).strip()) # to present example lines with unknown separator
                    
            # Check assigner
            if "=" in attrs: # Standard GFF3
                assign_candidates.add("=")
            elif re.search(r'\w+\s+"[^"]+"', attrs): # standard GTF (e.g., gene_id "Gene1"
                assign_candidates.add(" ")
                any_quotes = True

            else: 
                unknown_asgn_count += 1
                bad_lines.append(self.delimiter.join(f).strip()) # to present example lines with unknown assigner
                    
            # Quotation mark 
            if '"' in attrs:
                any_quotes = True

        quoting = "present" if any_quotes else "absent"

        # Report unknowns 
        if unknown_sep_count > 0:
            print(f"!!! {unknown_sep_count} lines had no common separator(';', ',', or tab). Check the file.")
        if unknown_asgn_count > 0:
            print(f"!!! {unknown_asgn_count} lines had no common assigner('=' or space). Check the file.")

        # Now, update the foundings.
        if "=" in assign_candidates and " " not in assign_candidates:
            format_likely = "GFF3-like"
            assigner_likely = "="
        elif " " in assign_candidates:
            if any_quotes:
                format_likely = "GTF-like"
                assigner_likely = " "
        else:
            format_likely = "Unknown"
            assigner_likely = ", and/or ".join(sorted(assign_candidates)) if assign_candidates else "?"
            print(f"Check the file, Could not determine assigner - options: {self.assigner or 'none'}")

        if ";" in sep_candidates:
            separator_likely = ";"
        elif "," in sep_candidates:
            if any("inside" in sf for sf in subformats):
                separator_likely = ";"
                print("',' found inside attribute values - using ';' as main separator.")
            else:
                separator_likely = ","
                print("Verify if ',' is the main separator. No semicolon found, treating ',' as main separator.")
        elif "\t" in sep_candidates:
            separator_likely = "\t"
        else:
            separator_likely = ", and/or ".join(sorted(sep_candidates)) if sep_candidates else "?"
            print(f"!!!Check the file, Could not determine separator - options: {self.separator or 'none'}")
        
        if apply:
            self.format = format_likely
            self.assigner = assigner_likely
            self.separator = separator_likely


        # Print summary
        print("\n====== Attribute Format Detection Summary ======")
        print(f"Analyzed up to {max_lines} lines.")
        print(f"Format                : {format_likely}")
        print(f"Detected separator(s) : {', '.join(sep_candidates) or 'none'}")
        print(f"Detected assigner(s)  : {', '.join(assign_candidates) or 'none'}")
        print(f"Is there Quotes?      : {quoting}")
        print(f"Subformats observed   : {', '.join(subformats) or 'none'}")
        print(f"Lines with unknown format : {unknown_sep_count + unknown_asgn_count}")
        if bad_lines:
            print("Examples of unknown lines : ")
            for ex in bad_lines:
                print(f" {ex}")
        print("==================================================\n")
        
        result = {
        "Format": self.format,
        "Separators": sep_candidates,
        "Assigners": assign_candidates,
        "Quotes": any_quotes,
        "Subformats": subformats,
        "Unknown_count": unknown_sep_count + unknown_asgn_count,
        "Unknown_examples": bad_lines,
        }

        return result
    
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
