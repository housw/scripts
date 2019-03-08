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
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from MulticoreTSNE import MulticoreTSNE as TSNE
from freqgen import k_mer_frequencies, codon_frequencies, genetic_codes


def get_codon_frequency_per_contig(contig2seqs):

    contig2frequencies = {} # {contig1:{'ATG':001, ...}, contig2:{'ATG':002, ... }}

    for contig, seqs in contig2seqs.items():
        frequency = codon_frequencies(seq=seqs, mode="absolute", genetic_code=11)
        contig2frequencies[contig] = frequency

    return contig2frequencies


def get_genes_per_contig(input_prodigal_nucl_file):
    """ return a dict which contains list of genes for each contig
    """

    contig2seqs = {} # {contig1:[seq1, seq2], contig2:[seq1, seq2]}

    for gene in SeqIO.parse(input_prodigal_nucl_file, "fasta"):
        header = gene.name
        contig = "_".join(header.split("_")[0:-1])
        seq = gene.seq
        # don't use genes if there are 'N's 
        #if "N" in seq:
        #    print("[WARNING]: {header} contains 'N' in the seq, ignored!".format(header=header))
        #    continue
        if len(seq)  % 3 != 0:
            print("[WARNING]: The length of {header} is not divisible by 3, ignored!".format(header=header))
            continue
        if contig not in contig2seqs:
            contig2seqs[contig] = [seq]
        else:
            contig2seqs[contig].append(seq)

    return contig2seqs 


def get_codon_frequencies_for_contigs(input_prodigal_nucl_file, output_file, genetic_code=11):

    # get seqs per contig
    contig2genes = get_genes_per_contig(input_prodigal_nucl_file)

    # get frequencies
    contig2frequencies = get_codon_frequency_per_contig(contig2genes)

    # write to output
    with open(output_file, "w") as oh:
        oh.write("Contig_ID"+"\t"+"\t".join(genetic_codes[genetic_code])+"\n")
        for contig, frequency_dict in contig2frequencies.items():
            line = [contig]
            for codon in genetic_codes[genetic_code]:
                frequency = frequency_dict.get(codon, 0)
                line.append(str(frequency))
            oh.write("\t".join(line)+"\n")


def compute_codon_tSNE_coordinates(input_codon_freq, output_file, threads=20):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """

    # run t-SNE 
    df = pd.read_csv(input_codon_freq, sep="\t", index_col=0, header=0)
    arr = np.array(df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=['CodonFreq_tSNE_X', 'CodonFreq_tSNE_Y'])
    tSNE_df.to_csv(output_file, sep="\t", header=True, index=True)


def main():

    # main parser
    parser = argparse.ArgumentParser(description="get codon frequencies of contigs")
    parser.add_argument("input_prodigal_ncul_file", help="input prodigal nucleotide file for predicted genes in fasta format")
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
        basename = os.path.basename(args.input_prodigal_ncul_file)
        args.prefix = os.path.splitext(basename)[0]
    CodonFreq_out_file = os.path.join(args.output_dir, args.prefix+"_CodonFreq.tsv")
    CodonFreq_tSNE_file = os.path.join(args.output_dir, args.prefix+"_CodonFreq_tSNE.tsv")
    if os.path.exists(CodonFreq_out_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # convert
    get_codon_frequencies_for_contigs(args.input_prodigal_ncul_file, CodonFreq_out_file)

    # t-SNE
    compute_codon_tSNE_coordinates(CodonFreq_out_file, CodonFreq_tSNE_file, threads=20)


if __name__ == "__main__":
    main()
