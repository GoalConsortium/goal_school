#!/usr/bin/env python3

'''Convert tagAlign files for further processing.'''

import os
import argparse
import shutil
import subprocess
import shlex
import logging
from multiprocessing import cpu_count
from python_utils import utils

EPILOG = '''
For more details:
        %(prog)s --help
'''

# SETTINGS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)


def get_args():
    '''Define arguments.'''

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--bam',
                        help="The bam file to convert.",
                        required=True)

    parser.add_argument('-p', '--paired',
                        help="True/False if paired-end or single end.",
                        default=False,
                        action='store_true')

    args = parser.parse_args()
    return args

# Functions


def check_tools():
    '''Checks for required componenets on user system'''

    logger.info('Checking for required libraries and components on this system')

    bedtools_path = shutil.which("bedtools")
    if bedtools_path:
        logger.info('Found bedtools: %s', bedtools_path)
    else:
        logger.error('Missing bedtools')
        raise Exception('Missing bedtools')

    samtools_path = shutil.which("samtools")
    if samtools_path:
        logger.info('Found samtools: %s', samtools_path)
    else:
        logger.error('Missing samtools')
        raise Exception('Missing samtools')


def convert_mapped(bam, tag_filename):
    '''Use bedtools to convert to tagAlign.'''

    out, err = utils.run_pipe([
        "bamToBed -i %s" % (bam),
        r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'""",
        "gzip -nc"], outfile=tag_filename)


def convert_mapped_pe(bam, bam_basename):
    '''Use bedtools to convert to tagAlign PE data.'''

    bedpe_filename = bam_basename + ".bedpe.gz"

    # Name sort bam to make BEDPE
    nmsrt_bam_filename = bam_basename + ".nmsrt.bam"
    samtools_sort_command = \
        "samtools sort -n -@%d -o %s %s" \
        % (cpu_count(), nmsrt_bam_filename, bam)

    logger.info(samtools_sort_command)
    subprocess.check_output(shlex.split(samtools_sort_command))

    out, err = utils.run_pipe([
        "bamToBed -bedpe -mate1 -i %s" % (nmsrt_bam_filename),
        "gzip -nc"], outfile=bedpe_filename)


def main():
    args = get_args()
    paired = args.paired
    bam = args.bam

    # Create a file handler
    handler = logging.FileHandler('convert_reads.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Convert PE or SE to tagAlign
    bam_basename = os.path.basename(
        utils.strip_extensions(bam, ['.bam']))

    tag_filename = bam_basename + '.tagAlign.gz'
    convert_mapped(bam, tag_filename)

    if paired:  # paired-end data
        convert_mapped_pe(bam, bam_basename)
    else:
        bedse_filename = bam_basename + ".bedse.gz"
        shutil.copy(tag_filename, bedse_filename)


if __name__ == '__main__':
    main()
