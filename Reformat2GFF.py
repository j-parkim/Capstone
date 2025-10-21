from collections import defaultdict
import re

class ReformatGFF(object):
    """
    To reformat inconsistent GFF/GTF files into standardized GFF3.
    
    Includes format detection:
    - Separator, assigner, quoting, and internal subformats.
    """

    def __init__(self, filepath, delimiter="\t"):
        """
        Args:
            - filepath (str): annotation file path
            - delimiter (str): column delimiter (default = tab)
        """
        self.filepath = filepath
        self.delimiter = delimiter
        self.separator = None # placeholder. TBD
        self.assigner = None # placeholder. TBD
        self.format = None # placeholder. TBD

    def example(self):
        """
        Present the first 10 rows of GFF files 
        as an example after comment lines
        """

        with open(self.filepath, 'r') as gff:
            ex_lines = []
            for l in gff:
                if l.startswith("#"):
                    continue
                ex_lines.append(l)
                if len(ex_lines) > 10:
                    break
            return ex_lines

    def detect_attr_format(self, max_lines=100, apply=True):
        """
        This is to detect likely format that attribute string is using.
        
        Args:
            - max_lines (int) : number of lines to look at for format detection
            - apply (bool) : if True, update self.separator, self.assigner, self.format

        Returns:
            - dictionary: detected likely format of attributes
        """
        sep_candidates = set()
        assign_candidates = set()
        any_quotes = False
        subformats = set()
        unknown_sep_count = 0
        unknown_asgn_count = 0
        bad_lines =[]

        with open(self.filepath, 'r') as gff:
                for i, l in enumerate(gff):
                    if i > max_lines:
                        break
                    if l.startswith("#"):
                        continue
                    fields = l.rstrip().split(self.delimiter)
                    if len(l.rstrip().split(self.delimiter)) < 9:
                        print("Check the file. There's no attribute column.")
                        continue # if there's no attribute column, skip
                    attrs = fields[8].strip()
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
                        bad_lines.append(l.strip()) # to present example lines with unknown separator
                    
                    # Check assigner
                    if "=" in attrs: # Standard GFF3
                        assign_candidates.add("=")
                    elif re.search(r'\w+\s+"[^"]+"', attrs): # standard GTF (e.g., gene_id "Gene1"
                        assign_candidates.add(" ")
                        any_quotes = True

                    else: 
                        unknown_asgn_count += 1
                        bad_lines.append(l.strip()) # to present example lines with unknown assigner
                    
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
            self.format = "GFF3-like"
            self.assigner = "="
        elif " " in assign_candidates:
            if any_quotes:
                self.format = "GTF-like"
                self.assigner = " "
        else:
            self.format = "Unknown"
            self.assigner = ", and/or ".join(sorted(assign_candidates)) if assign_candidates else "?"
            print(f"Check the file, Could not determine assigner - options: {self.assigner or 'none'}")

        if ";" in sep_candidates:
            self.separator = ";"
        elif "," in sep_candidates:
            if any("inside" in sf for sf in subformats):
                self.separator = ";"
                print("',' found inside attribute values - using ';' as main separator.")
            else:
                self.separator = ","
                print("Verify if ',' is the main separator. No semicolon found, treating ',' as main separator.")
        elif "\t" in sep_candidates:
            self.separator = "\t"
        else:
            self.separator = ", and/or ".join(sorted(sep_candidates)) if sep_candidates else "?"
            print(f"!!!Check the file, Could not determine separator - options: {self.separator or 'none'}")
        
        if not apply:
            format_likely = self.format
            assigner_likely = self.assigner
            separator_likely = self.separator
        else:
            format_likely = self.format
            assigner_likely = self.assigner
            separator_likely = self.separator

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
