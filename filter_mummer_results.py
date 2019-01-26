#!/usr/bin/env python

# Copyright (C) 2016  Shengwei Hou: housw2010'at'gmail'dot'com
#
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
from Bio import SeqIO


rev_dict = {'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C', 
            'N': 'N'}


class MUMMER(object):

    def __init__(self, header, hits=[]):
        self.header = header
        self.hits = hits

    def filter_hits_by_length(self, length=20):
        ret = []
        for hit in self.hits:
            if int(hit.split()[-1]) >= length:
                ret.append(hit)
        return ret

    def __str__(self):
        return ">"+self.header+"\n"+"\n".join(self.hits)+"\n"


def parse_mummer_results(mummer_output):
    """
    :param mummer_output: input mummer output file
    :return:           yield MUMMER record as a generator
    """
    header = ""
    hits = []
    with open(mummer_output, "r") as ih:
        for line in ih:
            if line.startswith(">"):
                if header:
                    yield MUMMER(header, hits)
                header = line[1:].strip()
                hits = []
            else:
                hits.append(line.strip())
        yield MUMMER(header, hits)


def filter_mummer_by_length(input_mummer_result, length_cutoff, out_file):
    mummer_records = parse_mummer_results(input_mummer_result)
    with open(out_file, "w") as oh:
        for rec in mummer_records:
            hits = rec.filter_hits_by_length(length=length_cutoff)
            if hits:
                oh.write(">"+rec.header+"\n")
                for hit in hits:
                    oh.write(hit+"\n")

def main():

    # main parser
    parser = argparse.ArgumentParser(description="parse and filter mummer results")
    parser.add_argument('input_mummer_result', help='input mummer 4 column result file')
    parser.add_argument('-l', "--length_cutoff", type=int, default=20,
                        help="minimum length cutoff [default=20]")
    parser.add_argument('-q', "--query_list", nargs='+',
                        help="only write hits in this name list")
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
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if not args.prefix:
        basename = os.path.basename(args.input_mummer_result)
        args.prefix = os.path.splitext(basename)[0]
    out_file = os.path.join(args.output_dir, args.prefix+"_filtered.txt")
    if os.path.exists(out_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # filter mummer result 
    filter_mummer_by_length(args.input_mummer_result, args.length_cutoff, out_file)

if __name__ == "__main__":
    main()
