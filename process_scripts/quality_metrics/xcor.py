#!/usr/bin/env python3

'''Compute cross-correlation analysis.'''

import os
import argparse
import shutil
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

    parser.add_argument('-t', '--tag',
                        help="The tagAlign file to qc on.",
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

    r_path = shutil.which("R")
    if r_path:
        logger.info('Found R: %s', r_path)
    else:
        logger.error('Missing R')
        raise Exception('Missing R')


def xcor(tag, paired):
    '''Use spp to calculate cross-correlation stats.'''

    tag_basename = os.path.basename(utils.strip_extensions(tag, ['.gz']))

    uncompressed_tag_filename = tag_basename
    out, err = utils.run_pipe([
        'gzip -dc %s > ' % (tag)], outfile=uncompressed_tag_filename)

    # Subsample tagAlign file
    numReads = 15000000

    subsampled_tag_filename = \
        tag_basename + ".%d.tagAlign.gz" % (numReads/1000000)
    steps = [
        'grep -v "chrM" %s' % (uncompressed_tag_filename),
        'shuf -n %d --random-source=%s' % (numReads, uncompressed_tag_filename)
    ]
    if paired:
        steps.extend([r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'"""])

    steps.extend(['gzip -nc'])

    out, err = utils.run_pipe(steps, outfile=subsampled_tag_filename)

    # Calculate Cross-correlation QC scores
    cc_scores_filename = subsampled_tag_filename + ".cc.qc"
    cc_plot_filename = subsampled_tag_filename + ".cc.plot.pdf"

    # CC_SCORE FILE format
    # Filename <tab>
    # numReads <tab>
    # estFragLen <tab>
    # corr_estFragLen <tab>
    # PhantomPeak <tab>
    # corr_phantomPeak <tab>
    # argmin_corr <tab>
    # min_corr <tab>
    # phantomPeakCoef <tab>
    # relPhantomPeakCoef <tab>
    # QualityTag

    run_spp_command = shutil.which("run_spp.R")
    out, err = utils.run_pipe([
        "Rscript %s -c=%s -p=%d -filtchr=chrM -savp=%s -out=%s" %
        (run_spp_command, subsampled_tag_filename, cpu_count(),
         cc_plot_filename, cc_scores_filename)
    ])

    return cc_scores_filename


def main():
    args = get_args()
    paired = args.paired
    tag = args.tag

    # Create a file handler
    handler = logging.FileHandler('xcor.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Calculate Cross-correlation
    xcor_filename = xcor(tag, paired)


if __name__ == '__main__':
    main()
