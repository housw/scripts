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
import pymer
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
from MulticoreTSNE import MulticoreTSNE as TSNE


def get_kmer_count_per_contig(contig, ksize=5):

    contig_name = contig.name
    contig_seq = str(contig.seq)

    # initialize ExactKmerCounter
    kc = pymer.ExactKmerCounter(ksize)
    kc.consume(contig_seq)
    kmer_dict = kc.to_dict(sparse=False)

    return (contig_name, kmer_dict)


def get_kmer_counts_for_contigs(input_contig_file, output_file, k=5, cpus=10):

    # do kmer count with multiple cores using mp
    contig_iter = SeqIO.parse(input_contig_file, "fasta")
    pool = mp.Pool(processes=cpus)
    result_iter = pool.imap(get_kmer_count_per_contig, contig_iter)

    # write 
    first_contig_name, first_kmer_dict = next(result_iter)
    header = sorted(first_kmer_dict.keys())
    first_count = [str(first_kmer_dict[kmer]) for kmer in header]

    with open(output_file, 'w') as oh:
        oh.write("Contig_ID" +"\t"+ "\t".join(header)+"\n")
        oh.write(first_contig_name +"\t"+ "\t".join(first_count)+"\n")
        # ask for forgiveness than permission
        while True:
            try:
                contig_name, kmer_dict = next(result_iter)
                count = [str(kmer_dict[kmer]) for kmer in header]
                oh.write(contig_name +"\t"+ "\t".join(count)+"\n")
            except IndexError:
                break
            except StopIteration:
                break

def get_kmer_frequencies_for_contigs(input_kmer_count_file, output_file):
    """ divide each kmer count by the sum of kmer counts in each row
    """

    kmercount = pd.read_csv('test_KmerCount.tsv', sep='\t', header=0, index_col=0)
    kmerfreq = kmercount.div(kmercount.sum(axis=1), axis=0)
    kmerfreq.to_csv(output_file, sep="\t", header=True, index=True)


def compute_kmer_tSNE_coordinates(input_kmer_freq, output_file, threads=20):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """

    # run t-SNE 
    df = pd.read_csv(input_kmer_freq, sep="\t", index_col=0, header=0)
    arr = np.array(df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=['KmerFreq_tSNE_X', 'KmerFreq_tSNE_Y'])
    tSNE_df.to_csv(output_file, sep="\t", header=True, index=True)


def main():

    # main parser
    parser = argparse.ArgumentParser(description="get kmer frequencies of contigs")
    parser.add_argument("input_contig_file", help="input contig file in fasta format")
    parser.add_argument("-k", "--kmer_size", type=int, default=5, help="kmer size [5]")
    parser.add_argument('-c', "--cpus", type=int, default=10, help="number of cpus to count kmer [10]")
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

    KmerCount_out_file = os.path.join(args.output_dir, args.prefix+"_KmerCount.tsv")
    KmerFreq_out_file = os.path.join(args.output_dir, args.prefix+"_KmerFreq.tsv")
    KmerFreq_tSNE_file = os.path.join(args.output_dir, args.prefix+"_KmerFreq_tSNE.tsv")
    if os.path.exists(KmerCount_out_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # count exact kmers
    get_kmer_counts_for_contigs(args.input_contig_file, KmerCount_out_file, k=args.kmer_size, cpus=args.cpus)

    # calculate kmer frequencies
    get_kmer_frequencies_for_contigs(KmerCount_out_file, KmerFreq_out_file)

    # t-SNE
    compute_kmer_tSNE_coordinates(KmerFreq_out_file, KmerFreq_tSNE_file, threads=args.cpus)


if __name__ == "__main__":
    main()
