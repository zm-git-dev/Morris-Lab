# trim the 3' end of reads.
#
# Some more complicated read trimmers will trim partial adapter from
# 3' end or will consider read quality in deciding how much to trim.
# This routine is much simpler and just trims a fixed amount from the
# 3' end of each read.
#
# usage: python trim_reads.py -t 16 <infile.fastq> <out.fastq>
#


from __future__ import with_statement
import sys
import os
from optparse import OptionParser

from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator



def main(in_file, out_file, trim=0):
    trim = int(trim)

    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):

                trim_seq = seq[:len(seq)-trim]
                trim_qual = qual[:len(qual)-trim]
                out_handle.write("@%s\n%s\n+\n%s\n" % (title, trim_seq,
                                                       trim_qual))


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--trim", dest="trim", default=0)
    options, args = parser.parse_args()
    if len(args) == 0:
        sys.exit(run_tests(sys.argv))
    else:
        kwd = dict(trim = options.trim)
        main(*args, **kwd)
