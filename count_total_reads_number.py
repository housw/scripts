#!/usr/bin/env python

# Copyright (C) 2018  Shengwei Hou
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


from __future__ import print_function
import os
import sys
import argparse
import subprocess


class ReadCounter(object):

    _grep_map = {"none":"grep",
                "gz":"zgrep",
                "zip":"zgrep",
                "bz2":"bzgrep"
                }

    def __init__(self, input_file, format, compress_type):
        self.input_file = input_file
        self.format = format
        self.compress_type = compress_type

    def count_read_number(self):
        pass


class FastaReadCounter(ReadCounter):

    def count_read_number(self):
        grep_prog = ReadCounter._grep_map.get(self.compress_type)
        cmd = "{grep_prog} -cP '^>' {input_file}".format(grep_prog=grep_prog, input_file=self.input_file)
        return subprocess.check_output(cmd, shell=True)


class FastqReadCounter(ReadCounter):

    def count_read_number(self):
        grep_prog = ReadCounter._grep_map.get(self.compress_type)
        cmd = "{grep_prog} -cP '^@' {input_file}".format(grep_prog=grep_prog, input_file=self.input_file)
        return subprocess.check_output(cmd, shell=True)


class FastqcReadCounter(ReadCounter):

    def count_read_number(self):
        if self.compress_type == "zip":
            fastqc_zip = self.input_file
            fastqc_folder=fastqc_zip.rstrip(".zip")
            cmd = "unzip -c {fastqc_zip} {fastqc_folder}/fastqc_data.txt | grep '^Total Sequences'".format(
                fastqc_zip=fastqc_zip, fastqc_folder=fastqc_folder)
        else:
            fastqc_folder = self.input_file
            cmd = "cat {fastqc_folder}/fastqc_data.txt | grep '^Total Sequences'".format(fastqc_folder=fastqc_folder)
        line = subprocess.check_output(cmd, shell=True)
        count = line.strip().split()[-1]
        return count



class CounterDispatcher(ReadCounter):

    counter_map = {"fasta":FastaReadCounter,
                   "fastq":FastqReadCounter,
                   "fastqc":FastqcReadCounter
                   }

    def __init__(self, input_file, format, compress_type):
        try:
            super().__init__(input_file, format, compress_type)
        except Exception as e:
            print("WARNING: Python3 is not supported by your interpreter: {err_msg}, using Python2 instead".format(err_msg=e))
            super(CounterDispatcher, self).__init__(input_file, format, compress_type)
        self.counter = self._get_read_counter()
        self.read_count = 0

    def _get_read_counter(self):
        return CounterDispatcher.counter_map.get(self.format, None)

    def count_read_number(self):
        counter = self.counter(self.input_file, self.format, self.compress_type)
        self.read_count = counter.count_read_number()

    def __str__(self):
        _basename = os.path.basename(self.input_file)
        _filestem = os.path.splitext(_basename)[0]
        return "{filestem} : {read_count:d}".format(filestem=_filestem, read_count=int(self.read_count))


def main():

    # main parser
    parser = argparse.ArgumentParser(description="count total read number in given input file")
    parser.add_argument("input_file", help="input file")
    parser.add_argument("-t", "--format", choices=["infer", "fasta", "fastq", "fastqc"],
                        default="infer", help="input file format, \
                             can be infer/fasta/fastq/fastqc, default is infer.")
    parser.add_argument("-c", "--compress_type", choices=["none", "gz", "zip", "bz2"],
                        default="infer", help="input file compress type, \
                                 can be infer/none/gz/zip/bz2, default is infer, \
                                 it will guess gz/zip/bz2 from the suffix")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()
    suffix = args.input_file.split(".")[-1]

    # compress_type handeling
    compress_type = args.compress_type
    compress_map = { "gz": "gz", "gzip":"gz",
                     "bz2": "bz2", "bzip2": "bz2",
                     "zip": "zip"}
    if compress_type == "infer":
        if suffix in compress_map:
            compress_type = compress_map.get(suffix)
        else:
            compress_type = "none"

    # input file format handeling
    format_map = { "fa": "fasta", "fas": "fasta", "fasta": "fasta", "fna": "fasta", "faa": "fasta",
                   "fq": "fastq", "fastq": "fastq",
                   "fastqc": "fastqc"}
    format = args.format
    if format == "infer":
        if "fastqc" in args.input_file:
            suffix = "fastqc"
        elif suffix in compress_map:
            if len(args.input_file.split(".")) > 2:
                suffix = args.input_file.split(".")[-2]
            else:
                sys.stderr.write("\nERROR: {suffix} file format is not supported!\n".format(suffix=suffix))
                sys.exit(1)

        if suffix in format_map:
            format = format_map.get(suffix)
        else:
            sys.stderr.write("\nERROR: {suffix} file format is not supported!\n".format(suffix=suffix))
            sys.exit(1)

    # input and output handeling
    if args.prefix:
        prefix = args.prefix
        out_file = os.path.join(args.out_folder, prefix + "_read_number.txt")
    else:
        basename = os.path.basename(args.input_file)
        prefix = os.path.splitext(basename)[0]
        out_file = os.path.join(args.out_folder, prefix + "_read_number.txt")

    if os.path.exists(out_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    with open(out_file, "w") as oh:
        counter = CounterDispatcher(args.input_file, format, compress_type)
        counter.count_read_number()
        oh.write(str(counter) + '\n')




if __name__ == "__main__":
    main()
