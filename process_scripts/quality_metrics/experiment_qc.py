#!/usr/bin/env python3

'''Experiment correlation and enrichment of reads over genome-wide bins.'''


import argparse
import logging
import subprocess
import shutil
from multiprocessing import cpu_count
import pandas as pd


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

    parser.add_argument('-d', '--design',
                        help="The design file to run QC (tsv format).",
                        required=True)

    args = parser.parse_args()
    return args


def check_tools():
    '''Checks for required componenets on user system.'''

    logger.info('Checking for required libraries and components on this system')

    deeptools_path = shutil.which("deeptools")
    if deeptools_path:
        logger.info('Found deeptools: %s', deeptools_path)
    else:
        logger.error('Missing deeptools')
        raise Exception('Missing deeptools')


def generate_read_summary(design):
    '''Generate summary of data based on read counts.'''

    bam_files = ' '.join(design['bam_reads'])
    labels = ' '.join(design['sample_id'])
    mbs_filename = 'sample_mbs.npz'

    mbs_command = (
        "multiBamSummary bins -p %d --bamfiles %s --labels %s -out %s"
        % (cpu_count(), bam_files, labels, mbs_filename)
        )

    logger.info("Running deeptools with %s", mbs_command)

    read_summary = subprocess.Popen(mbs_command, shell=True)
    out, err = read_summary.communicate()

    return mbs_filename


def check_correlation(mbs):
    '''Check Spearman pairwise correlation of samples based on read counts.'''

    spearman_filename = 'heatmap_SpearmanCorr.png'
    spearman_params = "--corMethod spearman --skipZero" + \
                    " --plotTitle \"Spearman Correlation of Read Counts\"" + \
                    " --whatToPlot heatmap --colorMap RdYlBu --plotNumbers"

    spearman_command = (
        "plotCorrelation -in %s %s -o %s"
        % (mbs, spearman_params, spearman_filename)
    )

    logger.info("Running deeptools with %s", spearman_command)

    spearman_correlation = subprocess.Popen(spearman_command, shell=True)
    out, err = spearman_correlation.communicate()


def check_coverage(design):
    '''Asses the sequencing depth of samples.'''

    bam_files = ' '.join(design['bam_reads'])
    labels = ' '.join(design['sample_id'])
    coverage_filename = 'coverage.png'
    coverage_params = "-n 1000000 --plotTitle \"Sample Coverage\"" + \
                    " --ignoreDuplicates --minMappingQuality 10"

    coverage_command = (
        "plotCoverage -b %s --labels %s %s --plotFile %s"
        % (bam_files, labels, coverage_params, coverage_filename)
        )

    logger.info("Running deeptools with %s", coverage_command)

    coverage_summary = subprocess.Popen(coverage_command, shell=True)
    out, err = coverage_summary.communicate()


def update_controls(design):
    '''Update design file to append controls list.'''

    logger.info("Running control file update.")

    file_dict = design[['sample_id', 'bam_reads']] \
                .set_index('sample_id').T.to_dict()

    design['control_reads'] = design['control_id'] \
                                .apply(lambda x: file_dict[x]['bam_reads'])

    logger.info("Removing rows that are there own control.")

    design = design[design['control_id'] != design['sample_id']]

    return design


def check_enrichment(sample_id, control_id, sample_reads, control_reads):
    '''Asses the enrichment per sample.'''

    fingerprint_filename = sample_id + '_fingerprint.png'

    fingerprint_command = (
        "plotFingerprint -b %s %s --labels %s %s --plotFile %s"
        % (sample_reads, control_reads, sample_id, control_id, fingerprint_filename)
        )

    logger.info("Running deeptools with %s", fingerprint_command)

    fingerprint_summary = subprocess.Popen(fingerprint_command, shell=True)
    out, err = fingerprint_summary.communicate()


def main():
    args = get_args()
    design = args.design

    # Create a file handler
    handler = logging.FileHandler('experiment_qc.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Read files
    design_df = pd.read_csv(design, sep='\t')

    # Run correlation
    mbs_filename = generate_read_summary(design_df)
    check_correlation(mbs_filename)

    # Run coverage
    check_coverage(design_df)

    # Run enrichment
    new_design_df = update_controls(design_df)
    for index, row in new_design_df.iterrows():
        check_enrichment(
            row['sample_id'],
            row['control_id'],
            row['bam_reads'],
            row['control_reads'])


if __name__ == '__main__':
    main()
