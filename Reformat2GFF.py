from collections import defaultdict
class ReformatGFF(object):
    """
    To reformat inconsistent GFF/GTF files into standardized GFF3.
    using attribute detection from ParseGffinfo.detect_attr_format()
    """

    def __init__(self, filepath, delimiter="\t"):
        """
        Args:
            - filepath (str): annotation file path
            - delimiter (str): column delimiter (default = tab)
        """
        self.filepath = filepath
        self.parser = ParseGFFinfo(filepath)
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

    def detect_format(self, max_lines = 100):
        """Wrapper for ParseGFFinfo.detect_attr_format"""
        format_info = self.parser.detect_attr_format(max_lines=max_lines)

        seps = format_info.get("Separators", set())
        asgns= format_info.get("Assigners", set())

        # Decide separator
        if ";" in seps:
            self.separator = ";"
        elif "," in seps:
            self.separator = ","
        elif "\t" in seps:
            self.separator = "\t"
        else:
            self.separator = next(iter(seps),";")
            if seps:
                print(f"!!!!Multiple/Unknown separators deteccted: {seps}. Using '{self.separator}'")
        
        # Decide assigner
        if "=" in asgns:
            self.assigner = "="
        elif " " in asgns:
            self.assigner = " "
        else:
            self.assigner = next(iter(asgns), "=")
            if asgns:
                print(f"!!!!Multiple/Unknown assigners deteccted: {asgns}. Using '{self.assigner}'")

        
        self.format = format_info.get("Format","unknown")

        print("\n=== Format summary (ReformatGFF) ===")
        print(f"Format detected     : {self.format}")
        print(f"Separator chosen    : '{self.separator}'  (from {seps})")
        print(f"Assigner chosen     : '{self.assigner}'  (from {asgns})")
        print(f"Quoted attributes   : {format_info.get('Quotes', False)}")
        print(f"Subformats observed : {format_info.get('Subformats', 'none')}")
        print("====================================\n")

        return format_info
    
    # Internal helper function when needed
    def _ensure_valid_gff3_value(self, key, v):
        """
        Quote GFF3 attribute values only when necessary
            e.g., product=glucose-1-phosphate adenylyltransferase large subunit
                  -> product="glucose-1-phosphate adenylyltransferase large subunit"
        Skip quoting for knwon list-like attributes        
        """

        listlike_keys = {"dbxref", "ontology_term", "is_a", "derives_from", "belongs_to"}

        # quote text-like values
        if key.lower() not in listlike_keys:
            if any(symbol in v for symbol in [" ", ";", "=", ","]):
                if not (v.startswith('"') and v.endswith('"')): # check if already quoted
                    v = f'"{v}"'
        return v


    def reformat(self, tobe="GFF3", outfile = True):
        """
        Convert attributes to standardized GFF3 or GTF-like syntax.
        Uses detected separator/assigner if available.

        Args:
            tobe (str): choose the final file to be either "GFF3" or "GTF" format. (default: GFF3)
            outfile (Bool or str): generate file if true
                            True(Default) - save as filename_standardized.gff3
                            False - print to screen only
                            str - save using given filename
        """
        if self.separator is None or self.assigner is None: # if user didn't run detect_format in the class 
            print(">>>Detecting attribute format first...")
            self.detect_format()

        if tobe.lower() not in {"gff3","gtf"}:
            raise ValueError("tobe argument must be 'GFF3' or 'GTF'")
        
        sep = self.separator or ";"
        assigner = self.assigner or "="

        out = (f"{self.filepath}_standardized.{tobe.lower()}" 
               if outfile is True 
               else outfile if isinstance(outfile, str)
               else None)
        
        outfile_ = open(out, "w") if out else None

        with open(self.filepath, "r") as infile:
            for l in infile:
                if l.startswith("#"):
                    if outfile_:
                        outfile_.write(l)
                    else:
                        print(l.strip())
                    continue

                fields = l.rstrip().split(self.delimiter)
                if len(fields) < 9:
                    continue

                # Parse attribues using detected marks 
                # If line has no assigner returns, just print
                attrs = {}
                for kv in fields[8].split(sep):
                    kv = kv.strip()
                    if not kv:
                        continue
                    if assigner in kv:
                        k, v = kv.split(assigner,1)
                        attrs[k.strip()] = v.strip()
                    else:
                        attrs[kv.strip()] = ""

                # Convert to target syntax
                if tobe.lower() == "gtf":
                    new_attrs = "; ".join(f'{k} "{v}"' for k, v in attrs.items()) + ";"
                else: # gff3
                    valid_attrs = []
                    for k,v in attrs.items():
                        v = self._ensure_valid_gff3_value(k,v)
                        valid_attrs.append(f"{k}={v}")
                    new_attrs = ";".join(valid_attrs)
                
                if outfile_:
                    outfile_.write("\t".join(fields[:8] + [new_attrs]) + "\n")
                else:
                    print("\t".join(fields[:8] + [new_attrs]))
        
        if outfile_:
            outfile_.close()
            print(f"Reformatted file written to: {out}")
        else:
            print(">>> Reformatting completed. See below.")
