#!/usr/bin/env python

# Copyright (C) 2019  Shengwei Hou
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import re
import sys
import Bio
import argparse
from pybedtools import BedTool


class Fasta(object):

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq.upper()

    def get_gc_content(self):
        ret = 0
        if len(self.seq) == 0:
            return ret
        for c in self.seq:
            if c == "G" or c == "C":
                ret += 1
        return float(ret)/len(self.seq)

    def __str__(self):
        return ">"+self.header+"\n"+self.seq+"\n"


def parse_fasta(fasta_file):
    """
    :param fasta_file: input fasta file
    :return:           yield Fasta record as a generator
    """
    header = ""
    seq = []
    with open(fasta_file, "r") as ih:
        for line in ih:
            if line.startswith(">"):
                if header:
                    yield Fasta(header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield Fasta(header, "".join(seq))


def substract_genes_from_chroms(chrom_bed, gene_bed):

    chroms = BedTool(chrom_bed)
    genes = BedTool(gene_bed)
    chrom_sub_genes = chroms.subtract(genes)

    return chrom_sub_genes 


def get_intergenic_sequences(contig_file, chrom_bed, gene_bed, output_igr_fa):

    # subtract gene bed and write sequences 
    chrom_sub_genes = substract_genes_from_chroms(chrom_bed, gene_bed)
    chrom_sub_genes.sequence(contig_file)
    chrom_sub_genes.save_seqs(output_igr_fa)


def merge_intergenic_sequences(igr_fa, output_merged_fa, fill_N_length=100):
    with open(output_merged_fa, "w") as oh:
        igr_records = parse_fasta(igr_fa)
        first_igr = next(igr_records)
        header = first_igr.header.split(':')[0]
        seq = [first_igr.seq]
        for rec in igr_records:
            _header = rec.header.split(':')[0]
            if _header == header:
                seq.append('N'*fill_N_length)
                seq.append(rec.seq)
            else:
                oh.write(">" + header +"\n")
                oh.write("".join(seq)+"\n")
                header = rec.header.split(':')[0]
                seq = [rec.seq]


def generate_bed_files(prodigal_gene_file, output_chrom_bed, output_gene_bed):
    p = re.compile(r'seqlen=(\d+);seqhdr=\"(\S+)\"$')
    with open(prodigal_gene_file, "r") as ih, open(output_chrom_bed, "w") as chrom_bed, open(output_gene_bed, "w") as gene_bed:
        for line in ih:
            if line.startswith("# Sequence Data: "):
                match = p.findall(line)
                seqlen = match[0][0]
                seqhdr = match[0][1]
                chrom_bed.write("{chrom}\t{start}\t{end}\n".format(chrom=seqhdr, start=1, end=seqlen))
            elif line.startswith("# Run Data:") or line.startswith("Beg") or not line.strip():
                continue
            else:
                line = line.strip().split()
                gene_start = line[0]
                gene_end = line[1]
                gene_bed.write("{chrom}\t{start}\t{end}\n".format(chrom=seqhdr, start=gene_start, end=gene_end))


def main():

    # main parser
    parser = argparse.ArgumentParser(description="get intergenic regions of input contigs")
    parser.add_argument("input_contig_file", help="input contig file in fasta format")
    parser.add_argument("input_prodigal_gene_file", help="input all potential gene file predicted by prodigal")
    parser.add_argument("--fill_N_length", type=int, default=100, help="The number of N's to merge two intergenic regions [100]")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--output_dir", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    if not args.prefix:
        basename = os.path.basename(args.input_contig_file)
        args.prefix = os.path.splitext(basename)[0]
       
    output_chrom_bed = os.path.join(args.output_dir, args.prefix+"_chrom.bed")
    output_gene_bed = os.path.join(args.output_dir, args.prefix+"_gene.bed")
    output_igr_fa = os.path.join(args.output_dir, args.prefix+"_igr.fa")
    output_merged_fa = os.path.join(args.output_dir, args.prefix+"_igr_merged.fa")
    if os.path.exists(output_gene_bed):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # generate bed files
    generate_bed_files(args.input_prodigal_gene_file, output_chrom_bed, output_gene_bed)

    # get intergenic region in fasta format 
    get_intergenic_sequences(args.input_contig_file, output_chrom_bed, output_gene_bed, output_igr_fa)

    # merge sequences in each contig 
    merge_intergenic_sequences(output_igr_fa, output_merged_fa, args.fill_N_length)


if __name__ == "__main__":
    main()

