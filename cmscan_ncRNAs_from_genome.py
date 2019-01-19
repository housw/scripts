#!/usr/bin/python


import os
import sys
import subprocess
import argparse
from Bio import SeqIO


def do_cmscan(input_fna, rfam_cm, rfam_clanin, output_dir, prefix):
    """ this function will be used to scan ncRNAs based on input_cm database, 
        from input_fna sequences, the infernal v1.1.2 needs to be installed 
    """

    # output 
    output_prefix = os.path.join(output_dir, prefix)

    # compute database size
    # For the purposes of Infernal, the total database size is the number of nucleotides that will be searched, in units of megabases (Mb, millions of nucleotides). So, it is the total number of nucleotides in all sequences that make up the genome, multiplied by two (because both strands will be searched), and divided by 1,000,000 (to convert to millions of nucleotides).
    database_size = 0
    for contig in SeqIO.parse(input_fna, "fasta"):
        contig_size = len(contig.seq)
        database_size += contig_size 
    database_size = (2 * database_size) / 1000000.0

    # run cmscan
    try:
        cmscan_proc = subprocess.Popen(['cmscan',
                                        "-Z", str(database_size),
                                        "--cut_ga", 
                                        "--rfam", 
                                        "--nohmmonly",
                                        "--tblout", output_prefix+".tblout",
                                        "--fmt", "2", 
                                        "--clanin", rfam_clanin,
                                        "--cpu", "20",
                                        "-o",  output_prefix+".cmscan",
                                        rfam_cm,
                                        input_fna 
                                       ])
        cmscan_proc.wait()
    except Exception as e:
        print "Can't run cmscan due to the following error:\n\t{err}".format(err=e)


def main():

    # main parser
    parser = argparse.ArgumentParser(description="scan ncRNAs against rfam using cmscan for input genome in fasta format")
    parser.add_argument("input_fna", help="input genome in fasta format")
    parser.add_argument("--path_to_rfam_cm", default="/home/shengwei/db/Rfam/current/Rfam.cm", help="absolute path to rfam.cm")
    parser.add_argument("--path_to_rfam_clanin", default="/home/shengwei/db/Rfam/current/Rfam.clanin", help="absolute path to rfam.clanin")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--output_dir", help="output directory")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    args = parser.parse_args()

    # prefix handeling
    if not args.prefix:
        basename = os.path.basename(args.input_fna)
        args.prefix = os.path.splitext(basename)[0]
    
    # output_dir handeling
    if not args.output_dir:
        args.output_dir = os.path.dirname(os.path.abspath(args.input_fna))

    # do cmscan
    do_cmscan(args.input_fna, args.path_to_rfam_cm, args.path_to_rfam_clanin, args.output_dir, args.prefix)    


if __name__ == "__main__":
    main()
