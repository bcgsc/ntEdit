#!/usr/bin/env python3
'''
Given the result of LOJ between VCF files, output valid ntEdit VCF file
'''
import argparse
import shlex
import subprocess
import gzip

class Vcf:
    "Represents a VCF entry"
    def __init__(self, chr, position, idenfitier, ref, alt, qual, filter, info, format=None, integration=None):
        self.chr = chr
        self.position = int(position)
        self.idenfitier = idenfitier
        self.ref = ref
        self.alt = alt.split(",")
        self.qual = qual
        self.filter = filter
        self.info_alts = self.parse_info(info)
        self.format = format
        self.integration = integration

    def remove_info(self):
        "Clear out the info and position attributes"
        self.info_alts = {}
        self.position = -1

    def get_info(self, alt_base):
        "Return the info for the given alt_base"
        info_str = ""
        if alt_base in self.alt and alt_base in self.info_alts:
            info_str += ";".join(self.info_alts[alt_base])
        return info_str

    def strip_string_from_info(self, info):
        "Strip ^NA from the info string" 
        return info.replace("^NA", "")

    def add_info(self, alt_dict):
        "Add INFO from the dictionary"
        for alt_base in alt_dict:
            if alt_base in self.info_alts:
                continue
            self.info_alts[alt_base] = alt_dict[alt_base]
            self.alt.append(alt_base)

    def parse_info(self, info):
        "Parse INFO string to data structures for: alt alleles, both"
        info_formatted = self.strip_string_from_info(info) if "^" in info else info
        alt_alleles = {}
        for item in info_formatted.split(";"):
            if "=" in item:
                key, value = item.split("=")
                split_val = value.split(",")
                if len(split_val) > 1 and len(split_val) == len(self.alt):
                    for i, v in enumerate(split_val):
                        if self.alt[i] not in alt_alleles:
                            alt_alleles[self.alt[i]] = []
                        alt_alleles[self.alt[i]].append(f"{key}={v}")
                else:
                    for alt_base in self.alt:
                        if alt_base not in alt_alleles:
                            alt_alleles[alt_base] = []
                        alt_alleles[alt_base].append(f"{key}={value}")

        return alt_alleles

    def get_info_dict(self):
        "Return the INFO in a dictionary form"
        return {alt_base: {info.split("=")[0]: info.split("=")[1] for info in info_list} \
                for alt_base, info_list in self.info_alts.items()}


    def __str__(self):
        strs = []
        for alt_base in self.alt:
            if self.format is None or self.integration is None:
                strs.append(f"{self.chr}\t{self.position}\t{self.idenfitier}\t{self.ref}\t{alt_base}\t{self.qual}"\
                        f"\t{self.filter}\t{self.get_info(alt_base)}")
            else:
                strs.append(f"{self.chr}\t{self.position}\t{self.idenfitier}\t{self.ref}\t{alt_base}\t{self.qual}"\
                            f"\t{self.filter}\t{self.get_info(alt_base)}\t{self.format}\t{self.integration}")
        return "\n".join(strs)

def write_header(vcffile, outfile, info_only=False):
    "Write header from input VCF file to standard out"
    vcf_open = open(vcffile, 'r', encoding="utf8") if not vcffile.endswith(".gz") else gzip.open(vcffile, 'rt')
    for line in vcf_open:
        line = line.strip()
        if info_only and line.startswith("##INFO"):
            outfile.write(f"{line}\n")
        elif not info_only and line.startswith("##"):
            outfile.write(f"{line}\n")
        elif not info_only and line.startswith("#"):
            vcf_open.close()
            return line
        elif line.startswith("#"):
            continue
        else:
            vcf_open.close()
            return
    vcf_open.close()

def get_all_vcf_info(vcf1, vcf2, alt_base):
    "Concatenate INFO strings from both VCF files"
    return f"{vcf1.get_info(alt_base)};{vcf2.get_info(alt_base)}".strip(";")

def are_compatible(vcf1, vcf2):
    "Ensure that at least one alt in vcf1 matches vcf2, references match, positions match"
    if vcf1.ref != vcf2.ref:
        return False
    if vcf1.position != vcf2.position:
        return False
    for alt_base in vcf1.alt:
        if alt_base in vcf2.alt:
            return True
    return False

def print_vcf_line(ntedit_vcf, l_vcf, outfile):
    "Print the info into a VCF"
    if l_vcf.position == -1:
        outfile.write(f"{ntedit_vcf}\n")
    else:
        if not are_compatible(ntedit_vcf, l_vcf):
            outfile.write(f"{ntedit_vcf}\n")
            return
        if len(ntedit_vcf.alt) > 1:
            out_strs = []
            for alt_base in ntedit_vcf.alt:
                out_strs.append(f"{ntedit_vcf.chr}\t{ntedit_vcf.position}\t{ntedit_vcf.idenfitier}\t{ntedit_vcf.ref}\t"
                        f"{alt_base}\t{ntedit_vcf.qual}\t{ntedit_vcf.filter}\t"
                        f"{get_all_vcf_info(ntedit_vcf, l_vcf, alt_base)}\t{ntedit_vcf.format}\t{ntedit_vcf.integration}")
            out_str = "\n".join(out_strs)
            outfile.write(f"{out_str}\n")
        else:
            alt_base = ntedit_vcf.alt[0]
            out_str = (f"{ntedit_vcf.chr}\t{ntedit_vcf.position}\t{ntedit_vcf.idenfitier}"
                    f"\t{ntedit_vcf.ref}\t{alt_base}\t{ntedit_vcf.qual}\t{ntedit_vcf.filter}\t"
                    f"{get_all_vcf_info(ntedit_vcf, l_vcf, alt_base)}\t{ntedit_vcf.format}\t"
                    f"{ntedit_vcf.integration}")
            outfile.write(f"{out_str}\n")

def parse_bedtools_loj(infile, outfile):
    "Parse the LOJ from bedtools to ntEdit-formatted VCF"
    ntedit_vcf = None
    l_vcf = None

    with open(infile, 'r', encoding="utf8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            ntedit_vcf_new = Vcf(*line[:10])
            l_vcf_new = Vcf(*line[10:18])
            if ntedit_vcf is not None and ntedit_vcf.position == ntedit_vcf_new.position and ntedit_vcf.chr == ntedit_vcf_new.chr:
                # This is the same position as before, tally extra INFO from l_vcf. Know that this cannot be due to ntEdit VCF.
                if l_vcf.position == -1 and are_compatible(ntedit_vcf, l_vcf_new):
                    # Currently only storing incompatible entry, OK to overwrite
                    l_vcf = l_vcf_new
                elif are_compatible(ntedit_vcf, l_vcf_new):
                    l_vcf.add_info(l_vcf_new.info_alts)
                    l_vcf.position = l_vcf_new.position
                else:
                    continue

            else:
                if ntedit_vcf is not None:
                    # This is not the same position, print previous and update
                    print_vcf_line(ntedit_vcf, l_vcf, outfile)
                ntedit_vcf = ntedit_vcf_new
                l_vcf = l_vcf_new
                if l_vcf_new.position != -1 and not are_compatible(ntedit_vcf_new, l_vcf_new):
                    l_vcf.remove_info()

        print_vcf_line(ntedit_vcf, l_vcf, outfile)

def fold_attributes(vcf_list):
    "Fold the attributes"
    attribute_lists = []
    alts = []
    for vcf in vcf_list:
        assert len(vcf.alt) == 1
        alts.append(vcf.alt[0])
        attribute_dict = vcf.get_info_dict()
        list_attributes = [key for _, base_dict in attribute_dict.items() for key in base_dict]
        attribute_lists.extend(list_attributes)
    attribute_set = set(attribute_lists)

    new_dict = {} # Key: [] to track lists of INFO
    for attribute in attribute_set:
        new_dict[attribute] = []
        for vcf in vcf_list:
            vcf_attributes = vcf.get_info_dict()[vcf.alt[0]]
            if attribute in vcf_attributes:
                new_dict[attribute].append(vcf_attributes[attribute])
            else:
                new_dict[attribute].append("NA")

    new_info_line = []
    visited_attributes = set()
    for attribute in attribute_lists:
        if attribute not in visited_attributes:
            new_info_line.append(f"{attribute}={','.join(new_dict[attribute])}")
            visited_attributes.add(attribute)
    rep_vcf = vcf_list[0]
    out_str = f"{rep_vcf.chr}\t{rep_vcf.position}\t{rep_vcf.idenfitier}\t{rep_vcf.ref}\t{','.join(alts)}\t"\
            f"{rep_vcf.qual}\t{rep_vcf.filter}\t{';'.join(new_info_line)}\t{rep_vcf.format}\t{rep_vcf.integration}"
    return out_str


def print_folded_vcf_line(vcf_list, outfile):
    "Print the lines to the VCF, folding variants and INFO if applicable"
    if len(vcf_list) == 1:
        outfile.write(f"{vcf_list[0]}\n")
    else:
        new_info_line = fold_attributes(vcf_list)
        outfile.write(f"{new_info_line}\n")

def refold_variants(in_vcf, out_prefix):
    "Re-interleave or refold the variants, so a given position is on one line"
    prev_vcf = None
    vcfs = []
    with open(in_vcf, 'r', encoding="utf8") as fin:
        with open(f"{out_prefix}.vcf", 'w', encoding="utf8") as fout:
            for line in fin:
                line = line.strip()
                if line.startswith("#"):
                    fout.write(f"{line}\n")
                else:
                    vcf_line = Vcf(*line.split("\t"))
                    if prev_vcf is not None and \
                        (prev_vcf.position != vcf_line.position or prev_vcf.chr != vcf_line.chr):
                        # Different, ready to print
                        print_folded_vcf_line(vcfs, fout)
                        vcfs = []
                    vcfs.append(vcf_line)
                    prev_vcf = vcf_line
            print_folded_vcf_line(vcfs, fout)
            
def remove_tmp_vcf(filename):
    "Remove this temporary file"
    cmd = f"rm {filename}"
    ret_code = subprocess.run(shlex.split(cmd), check=True)
    assert ret_code.returncode == 0


def main():
    """
    Given the result of LOJ between VCF files, output valid ntEdit VCF file
    """
    parser = argparse.ArgumentParser(
        description="Given the result of LOJ between VCF files, output valid ntEdit VCF file")
    parser.add_argument("-b", "--bedtools", help="Bedtools intersect file, or - to read from standard in",
                        required=True, type=str)
    parser.add_argument("--vcf", help="ntEdit VCF file (for header)", required=True, type=str)
    parser.add_argument("--vcf_l", help="-l VCF file (for header)", required=True, type=str)
    parser.add_argument("-p", "--prefix", help="Prefix for out files", required=False,
                        default="ntedit_snv_vcf", type=str)

    args = parser.parse_args()

    args.bedtools = "/dev/stdin" if args.bedtools == "-" else args.bedtools

    with open(f"{args.prefix}.tmp.vcf", 'w', encoding="utf8") as fout:
        header = write_header(args.vcf, fout, info_only=False)
        write_header(args.vcf_l, fout, info_only=True)
        fout.write(f"{header}\n")

        parse_bedtools_loj(args.bedtools, fout)

    refold_variants(f"{args.prefix}.tmp.vcf", args.prefix)

    remove_tmp_vcf(f"{args.prefix}.tmp.vcf")

if __name__ == "__main__":
    main()
