#!/bin/env python
import sys
import pysam
import argparse
def get_args():
    '''Define arguments.'''
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--bam',
                        help="The bam file",
                        required=True)
    parser.add_argument('-o', '--out',
                        help="The outfile",
                        default='outfile.bam')
    args = parser.parse_args()
    return args

# set the tag names - take a look at SAM spec to pick an appropriate one
args = get_args()
infile = pysam.AlignmentFile(args.bam, "rb")
out = pysam.AlignmentFile(args.out, "wb", template=infile)
for read in infile.fetch():
      read.set_tag('RX', read.qname.split(":")[-1])
      out.write(read)

infile.close()
out.close()
