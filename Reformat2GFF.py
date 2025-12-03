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
        - Check if ID/Parent attributes are there
        - WARN if essential structural attributes missing
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
    def __init__(self, filepath, delimiter="\t", 
                 remove_genome_attr=None, remove_sources=None):
        """
        Args:
            - filepath (str): annotation file path
            - delimiter (str): column delimiter (default = tab)
        """
        self.filepath = filepath
        self.parser = ParseGffinfo.ParseGFFinfo(filepath)
        self.delimiter = delimiter
        self.separator = None # placeholder. TBD
        self.assigner = None # placeholder. TBD
        self.format = None # placeholder. TBD
        self.remove_genome_attr = set(remove_genome_attr or [])
        self.remove_sources = set(remove_sources or [])

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
    def _normalize_attr_names(self, attrs, featuretype, summary, l_idx):
        """
        - Normalize some keys (such as biotype)
        - Detect missing ID and Parent where expected
        - Add ID ONLY when an alternative valid ID exits defined by _select_best_id_for_ft
        """
        normalized = {}
        ft = featuretype.lower()
        is_rna = (ft.endswith("rna") or ft.endswith("transcript"))

        # ------ Normalize mapped keys (biotype, Note, etc) ------
        for k, v in attrs.items():
            k_low = k.lower()
            v_clean = v.strip('"') # Remove quotations from values

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
            else:
                normalized[k] = v_clean
        
        # ------ ID handling (WARNING AND COPYING) ------
        has_id = any(k.lower() == "id" for k in normalized)

        if not has_id:
                candidate = self._select_best_id_for_ft(featuretype, attrs)
                if candidate:
                    print(f"[----INFO---- line {l_idx}] adding ID={candidate} for {featuretype} as missing")
                    normalized["ID"]=candidate
                    summary["ID_added_cases"] += 1
                else:
                    print(f"[WARNING!!! line {l_idx}] Missing ID and no alternative ID candidate for {featuretype}. Check the file.")
                    summary["ID_MISSING_cases"] += 1

        # ----- Checking Parent attribute where expected (WARNING ONLY) ------
        has_parent = any(k.lower() == "parent" for k in attrs)

        if is_rna and not has_parent:
            summary["RNA/transcript_missing_Parentkey"] += 1
            print(f"[WARNING!!! LINE {l_idx}] Missing Parent key \n {summary["CURRENT_LINE_RAW"]}")
        elif ft == "exon" and not has_parent:
            summary["EXON_missing_Parentkey"] += 1
            print(f"[WARNING!!! LINE {l_idx}] Missing Parent key \n {summary["CURRENT_LINE_RAW"]}")            
        elif ft == "cds" and not has_parent:
            summary["CDS_missing_Parentkey"] += 1
            print(f"[WARNING!!! LINE {l_idx}] Missing Parent key \n {summary["CURRENT_LINE_RAW"]}")

        return normalized
    # -------------------------------------------------------------------------------
    def reformat_dialects(self, tobe="GFF3", outfile = True):
        """
        - Convert attributes to standardized GFF3 or GTF-like syntax.
        - Uses detected separator/assigner if available.

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

        summary = {
            "Total_lines": 0,
            "ID_added_cases": 0,
            "ID_MISSING_cases": 0,
            "RNA/transcript_missing_Parentkey": 0,
            "EXON_missing_Parentkey": 0,
            "CDS_missing_Parentkey": 0,
            "strand_fixed_?_to_.": 0
        }

        out = (f"{self.filepath}_standardized.{tobe.lower()}" 
               if outfile is True 
               else outfile if isinstance(outfile, str)
               else None)
        
        outfile_ = open(out, "w") if out else None

        with self._open_file() as infile:
            for l_idx, l in enumerate(infile, start=1):

                if l.startswith("#"):
                    if outfile_:
                        outfile_.write(l)
                    else:
                        print(l.strip())
                    continue

                summary["Total_lines"] += 1
                summary["CURRENT_LINE_RAW"] = l.rstrip("\n")

                fields = l.rstrip().split(self.delimiter)
                if len(fields) < 9:
                    continue
                
                featuretype = fields[2]

                # ------ Skip if source and/or Genome match ------
                source = fields[1]
                if source in self.remove_sources:
                    continue

                if "genome" in attrs and attrs["genome"] in self.remove_genome_attr:
                    continue

                # ------ Fix strand marked as '?'  -> '.' -------
                # As it can throw an error in downstream processes
                if len(fields) > 6 and fields[6] == "?":
                    fields[6] = "."
                    summary["strand_fixed_?_to_."] +=1

                # ------ Parse attributes using detected symbols ------
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

                # ------ Normalize the keys ------
                if tobe.lower() == "gff3":
                    normed_attrs = self._normalize_attr_names(attrs, featuretype, summary, l_idx)
                else: # GTF
                    normed_attrs = {k: v.strip('"') for k,v in attrs.items()}

                # ------ Convert to target syntax ------
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

        # ------ Summary ------
        print("\n========Normalization Summary========")
        for k,v in summary.items():
            if k =="CURRENT_LINE_RAW":
                continue
            print(f"{k:35} : {v}")
        print("========================================")

        return summary
