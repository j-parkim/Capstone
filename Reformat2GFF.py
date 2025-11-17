import gzip
import re
from collections import defaultdict
import ParseGffinfo

class ReformatGFF(object):
    """
    To reformat inconsistent GFF/GTF files into standardized GFF3.
    using attribute detection from ParseGffinfo.detect_attr_format()

    Performs:
        - dialect normalization(separator/assigner)
        - consolidation of biotype keys
        - controlled ID/Parent prefixing based on feature type
        - preservation of original information (not necessary modification made)

    """
    # -------------------------------------
    # -----Attribute normalization map-----
    # -------------------------------------

    gff3_attr_mapping = {
        # 1. "biotype" Key candidates
        "biotype" : "biotype",
        "gene_type": "biotype",
        "gene_biotype": "biotype",
        "transcript_biotype": "biotype",
        "rna_type": "biotype",

        # 2. MISC
        "gbkey": "Note"
    }

    # -------------------------------------------------------------------------------
    def __init__(self, filepath, delimiter="\t"):
        """
        Args:
            - filepath (str): annotation file path
            - delimiter (str): column delimiter (default = tab)
        """
        self.filepath = filepath
        self.parser = ParseGffinfo(filepath)
        self.delimiter = delimiter
        self.separator = None # placeholder. TBD
        self.assigner = None # placeholder. TBD
        self.format = None # placeholder. TBD

        self.exon_counters = defaultdict(int) # Tracks the last exon number for each transcript ID
        self.dynamic_exon_numbering = True # If already numbered, skip dynamic numbering

    # -------------------------------------------------------------------------------
    # Open File
    def _open_file(self):
        if self.filepath.endswith(".gz"):
            return gzip.open(self.filepath, "rt")
        return open(self.filepath, "r")

    # -------------------------------------------------------------------------------
    def example(self, n=10):
        """
        Present the first 10 rows (or as defined) of GFF files 
        as an example after comment lines
        """

        with self._open_file() as gff:
            ex_lines = []
            for l in gff:
                if l.startswith("#"):
                    continue
                ex_lines.append(l.rstrip("\n"))
                if len(ex_lines) > n:
                    break
            return ex_lines

    # -------------------------------------------------------------------------------
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
       
        # Decide assigner
        if "=" in asgns:
            self.assigner = "="
        elif " " in asgns:
            self.assigner = " "
        else:
            self.assigner = next(iter(asgns), "=")
       
        self.format = format_info.get("Format","unknown")

        print("\n=== Format summary (ReformatGFF) ===")
        print(f"Format detected     : {self.format}")
        print(f"Separator chosen    : '{self.separator}'  (from {seps})")
        print(f"Assigner chosen     : '{self.assigner}'  (from {asgns})")
        print(f"Quoted attributes   : {format_info.get('Quotes', False)}")
        print(f"Subformats observed : {format_info.get('Subformats', 'none')}")
        print("====================================\n")

        return format_info
    
    # -------------------------------------------------------------------------------
    # HELPER function
    def _select_best_id_for_ft(self, featuretype, attrs):
        """
        Returns the BEST ID value for feature type
        """

        ft = featuretype.lower()
        norm = {k.lower(): v.strip().strip('"') for k, v in attrs.items()}

        # Priority list in order
        priority = []
        # 1. Genes
        if ft == "gene":
            priority = ["id", "gene_id", "locus_tag", "name"]

        # 2. Transcripts (features contatin "rna" or "transcript" in their names)
        elif ft.endswith("rna") or ft.endswith("transcript"):
            priority = ["id", "transcript_id", "name", "gene_id"]

        # 3. Exons
        elif ft == "exon":
            priority = ["id", "exon_id", "transcript_id"]
                
        # 4. CDS
        elif ft == "cds":
            priority = ["id", "protein_id", "transcript_id"]

        # 5. Other features
        else:
            priority = ["id", "name"]
        
        for key in priority:
            if key in norm and norm[key]:
                return norm[key]
            
        return None

    # -------------------------------------------------------------------------------
    def _normalize_attr_names(self, attrs, featuretype):
        """
        - Applies correct ID prefix based on feature type
        - Applies corret Parent prefix based on hierarchy
        - Preserve identifiers
        - Preserve other attributes as-is
        """
        normalized = {}
        ft = featuretype.lower()
        is_rna = (ft.endswith("rna")) or (ft.endswith("transcript"))

        #----Define prefix based on featuretype---
        id_prefix = (
            "gene-" if (ft == "gene") else
            "rna-" if is_rna else
            "exon-" if (ft == "exon") else
            "cds-" if (ft == "cds") else
            f"{ft}-" # if not in map, just use feature type.lower
        )

        # Prefixes to check for to prevent double-prefixing(rna-rna-XM.dflskdjf)
        known_prefixes = {"gene-", "rna-", "exon-", "cds-","transcript-"}

        #-----Mapping, Parent extraction -----
        parent_val = None

        # Keys that will be used to select the best ID.
        id_to_keep = {"id", "gene_id", "transcript_id", "exon_id", "protein_id", "locus_tag"}

        for k, v in attrs.items():
            k_low = k.lower()
            v_clean = v.strip('"') # Remove quotations from values
            
            # 0. Extract Parent
            if k_low == "parent":
                parent_val = v_clean
                continue

            # 1. Check if the key is in the normalization map
            mapped = self.gff3_attr_mapping.get(k_low)

            if mapped:
                # Handle mapped keys
                if mapped in normalized:
                    # Keep the longest/most descriptive for biotype
                    if mapped == "biotype":
                        if len(v_clean) > len(normalized[mapped]):
                            normalized[mapped] = v_clean
                else:
                    normalized[mapped] = v_clean
            
            if k_low in id_to_keep:
                # preserve as is.
                if k not in normalized:
                    normalized[k] = v_clean
            
            # Preserve other info as is
            elif mapped is None and k_low not in id_to_keep:
                if k not in normalized:
                    normalized[k] = v_clean
        
        # 1. Normalize ID
        id_val = self._select_best_id_for_ft(featuretype, attrs)

        if id_val:
            already_prefixed = any(id_val.startswith(p) for p in known_prefixes)

            # ----EXON numbering check-----
            if (ft == "exon") and self.dynamic_exon_numbering:

                if not already_prefixed:
                    #Check if the selected ID already has a suffix number
                    has_number = re.search(r'[.-]\d+$', id_val)

                    if not has_number:
                        # If missing, use parent Id to track the count
                        parent_id_base = attrs.get("transcript_id")

                        if parent_id_base:
                            # increment the counter for this specific transcript
                            parent_id_base = parent_id_base.strip('"')
                            self.exon_counters[parent_id_base] += 1
                            new_suffix = f"-{self.exon_counters[parent_id_base]}"

                            # Use the Parent/Transcript ID as the base for the new exon ID
                            id_val = parent_id_base + new_suffix

                            # We just reconstructed the ID without a prefix, so we'll add it
                            already_prefixed = False

            normalized["ID"] = id_val if already_prefixed else id_prefix + id_val

        # 2. Normalize Parent
        if parent_val:
            parent_prefix = "rna-" # as parent attributes exist for transcripts and cds
            if is_rna:
                parent_prefix = "gene-"
            
            already_prefixed = any(parent_val.startswith(p) for p in known_prefixes)

            normalized["Parent"] = parent_val if already_prefixed else parent_prefix + parent_val

        return normalized
    # -------------------------------------------------------------------------------
    def reformat_dialects(self, tobe="GFF3", outfile = True):
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

        with self._open_file() as infile:
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
                
                featuretype = fields[2]

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

                # Apply standardization and prefixing
                if tobe.lower() == "gff3":
                    normed_attrs = self._normalize_attr_names(attrs,featuretype)
                else: # GTF
                    normed_attrs = {k: v.strip('"') for k,v in attrs.items()}

                # Convert to target syntax
                if tobe.lower() == "gtf":
                    new_attrs = "; ".join(f'{k} "{v}"' for k, v in normed_attrs.items()) + ";"
                else: # gff3
                    new_attrs = ";".join(f"{k}={v}" if v!= "" else k for k,v in normed_attrs.items())
                
                if outfile_:
                    outfile_.write("\t".join(fields[:8] + [new_attrs]) + "\n")
                else:
                    print("\t".join(fields[:8] + [new_attrs]))
        
        if outfile_:
            outfile_.close()
            print(f"Reformatted file written to: {out}")
        else:
            print(">>> Reformatting completed. See below.")
