#!/usr/bin/env python
# Copyright (C) 2016  Shengwei Hou
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


class CmscanTbl(object):

    def __init__(self, ncRNA_name, ncRNA_accession, ncRNA_clanname, ncRNA_trunc, ncRNA_pass, mdl_from, mdl_to, query_name, query_from, query_to, query_strand, score, evalue, *args, **kwargs):

        self.ncRNA_name = ncRNA_name
        self.ncRNA_accession = ncRNA_accession
        self.ncRNA_clanname = ncRNA_clanname
        self.ncRNA_trunc = ncRNA_trunc
        self.ncRNA_pass = ncRNA_pass
        self.mdl_from = int(mdl_from)
        self.mdl_to = int(mdl_to) 
        self.query_name = query_name
        self.query_from = int(query_from)
        self.query_to = int(query_to)
        self.query_strand = query_strand
        self.score = float(score)
        self.evalue = float(evalue)

    def get_ncRNA_seq(self, seq_record):
        if self.query_strand == "+":
            ret_seq = seq_record[self.query_from-1:self.query_to]
        else:
            ret_seq = seq_record[self.query_to-1:self.query_from].reverse_complement()
        
        return ret_seq

    def __str__(self):
        return self.ncRNA_name +"\t"+ self.ncRNA_accession +"\t"+ \
               self.ncRNA_clanname +"\t"+ self.ncRNA_trunc +"\t"+ \
               self.ncRNA_pass +"\t"+ str(self.mdl_from) +"\t"+ \
               str(self.mdl_to) +"\t"+ self.query_name +"\t"+ \
               str(self.query_from) +"\t"+ str(self.query_to) +"\t"+ \
               self.query_strand +"\t"+ str(self.score) +"\t"+ str(self.evalue)


def parse_cmscan_tblout(cmscan_tblout_file):
    """
    :param cmscan_tblout_file: input cmscan domtblout file
    :return:           yield CmscanTbl record as a generator
    """

    # handle gzip compressed file
    if cmscan_tblout_file.endswith(".gz"):
        ih = gzip.open(cmscan_tblout_file, "r")
    else:
        ih = open(cmscan_tblout_file, "r")

    # parse records
    with open(cmscan_tblout_file, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            # interested 0-based column index 
            ncRNA_name = line[1]
            ncRNA_accession = line[2]
            ncRNA_clanname = line[5]
            ncRNA_trunc = line[12]
            ncRNA_pass = line[13]
            mdl_from = line[7]
            mdl_to = line[8]
            query_name = line[3]
            query_from = line[9]
            query_to = line[10]
            query_strand = line[11]
            score = line[16]
            evalue  = line[17]
            yield CmscanTbl(ncRNA_name, ncRNA_accession, ncRNA_clanname, ncRNA_trunc, ncRNA_pass, 
                            mdl_from, mdl_to, query_name, query_from, query_to, query_strand, score, evalue)


def get_fasta_from_cmscan_tblout_result(cmscan_domtblout, query_fasta,
                                    out_file, evalue_cutoff=1e-10, 
                                    name_list=[], accession_list=[]):
    """given a query fasta, and cmscan_domtblout, find the hits for each
       model that meet the evalue cutoff. 
    """

    # parse input query fasta 
    header2seq_record = {} # {header:seq_record}
    for fasta_rec in SeqIO.parse(query_fasta, "fasta"):
        header = fasta_rec.name
        seq = fasta_rec.seq
        assert header not in header2seq_record, "Duplicated fasta record header: {h}".format(h=header)
        header2seq_record[header] = seq

    with open(out_file, "w") as oh:
        for rec in parse_cmscan_tblout(cmscan_domtblout):
            query_name = rec.query_name
            query_seq = header2seq_record.get(query_name, None)
            if not query_seq:
                raise Exception("contig {q} was not found in query fasta file {f}".format(q=query_name, f=query_fasta))
            ncRNA_seq = rec.get_ncRNA_seq(query_seq)
            print(ncRNA_seq)
            if (not name_list) and (not accession_list):
                if rec.evalue <= evalue_cutoff:
                    oh.write(">" + rec.ncRNA_name +"\n")
                    oh.write(str(ncRNA_seq)+"\n")
            else:
                name_list = [] if name_list is None else name_list
                accession_list = [] if accession_list is None else accession_list
                if (rec.ncRNA_name in name_list) or (rec.ncRNA_accession in accession_list):
                    if rec.evalue <= evalue_cutoff:
                        oh.write(">" + rec.ncRNA_name +"\n")
                        oh.write(str(ncRNA_seq)+"\n")


def main():

    # main parser
    parser = argparse.ArgumentParser(description="get fasta from hmmseach result")
    parser.add_argument('input_cmscan_tblout', help='input cmscan domtblout file')
    parser.add_argument('input_query_fasta', help="query fasta file,used in cmscan")
    parser.add_argument('-e', "--evalue_cutoff", type=float, default=1e-10,
                        help="maximum evalue cutoff [default=1e-10]")
    parser.add_argument('-n', "--name_list", nargs='+',
                        help="only extract ncRNAs in this name list")
    parser.add_argument('-a', "--accession_list", nargs='+', 
                        help="only extract ncRNAs in this accession list")
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
        basename = os.path.basename(args.input_cmscan_tblout)
        args.prefix = os.path.splitext(basename)[0]
    out_file = os.path.join(args.output_dir, args.prefix+"_ncRNAs.fa")
    if os.path.exists(out_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # get fasta from cmscan result
    get_fasta_from_cmscan_tblout_result(args.input_cmscan_tblout, args.input_query_fasta,
                                    out_file, args.evalue_cutoff, args.name_list, args.accession_list)


if __name__ == "__main__":
    main()
