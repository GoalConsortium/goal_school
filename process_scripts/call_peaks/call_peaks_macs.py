#!/usr/bin/env python3

'''Generate peaks from data.'''

import os
import argparse
import shutil
import logging
from multiprocessing import cpu_count
from python_utils import utils
from quality_metrics.xcor import xcor as calculate_xcor

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
                        help="The tagAlign file to perform peak calling on.",
                        required=True)

    parser.add_argument('-x', '--xcor',
                        help="The cross-correlation file (if already calculated).",
                        required=True)

    parser.add_argument('-c', '--con',
                        help="The control tagAling file used for peak calling.",
                        required=True)

    parser.add_argument('-s', '--sample',
                        help="The sample id to name the peak files.",
                        required=True)

    parser.add_argument('-g', '--genome',
                        help="The genome size of reference genome.",
                        required=True)

    parser.add_argument('-z', '--size',
                        help="The file with chromosome sizes of reference genome.",
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

    macs_path = shutil.which("macs2")
    if r_path:
        logger.info('Found MACS2: %s', macs_path)
    else:
        logger.error('Missing MACS2')
        raise Exception('Missing MACS2')

    bg_bw_path = shutil.which("bedGraphToBigWig")
    if bg_bw_path:
        logger.info('Found bedGraphToBigWig: %s', bg_bw_path)
    else:
        logger.error('Missing bedGraphToBigWig')
        raise Exception('Missing bedGraphToBigWig')

    bedtools_path = shutil.which("bedtools")
    if bedtools_path:
        logger.info('Found bedtools: %s', bedtools_path)
    else:
        logger.error('Missing bedtools')
        raise Exception('Missing bedtools')


def call_peaks_macs(experiment, xcor, control, prefix, genome_size, chrom_sizes):

    # Extract the fragment length estimate from column 3 of the
    # cross-correlation scores file
    with open(xcor, 'r') as xcor_fh:
        firstline = xcor_fh.readline()
        frag_lengths = firstline.split()[2]  # third column
        fragment_length = frag_lengths.split(',')[0]  # grab first value
        logger.info("Fraglen %s" % (fragment_length))


    # Generate narrow peaks and preliminary signal tracks

    command = 'macs2 callpeak ' + \
              '-t %s -c %s ' % (experiment, control) + \
              '-f BED -n %s ' % (prefix) + \
              '-g %s -p 1e-2 --nomodel --shift 0 --extsize %s --keep-dup all -B --SPMR' % (genome_size, fragment_length)

    logger.info(command)
    returncode = utils.block_on(command)
    logger.info("MACS2 exited with returncode %d" % (returncode))
    assert returncode == 0, "MACS2 non-zero return"

    # MACS2 sometimes calls features off the end of chromosomes.
    # Remove coordinates outside chromosome sizes

    narrowpeak_fn = '%s_peaks.narrowPeak' % (prefix)
    clipped_narrowpeak_fn = 'clipped-%s' % (narrowpeak_fn)


    steps = ['slopBed -i %s -g %s -b 0' % (narrowpeak_fn, chrom_sizes),
             'bedClip stdin %s %s' % (chrom_sizes, clipped_narrowpeak_fn)]

    out, err = utils.run_pipe(steps)

    # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format
    # (score must be <1000)
    rescaled_narrowpeak_fn = utils.rescale_scores(
        clipped_narrowpeak_fn, scores_col=5)

    # Sort by Col8 in descending order and replace long peak names in Column 4
    # with Peak_<peakRank>
    steps = [
        'sort -k 8gr,8gr %s' % (rescaled_narrowpeak_fn),
        r"""awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'"""
    ]

    out, err = utils.run_pipe(steps, '%s' % (narrowpeak_fn))

    # For Fold enrichment signal tracks

    # This file is a tab delimited file with 2 columns Col1 (chromosome name),
    # Col2 (chromosome size in bp).

    command = 'macs2 bdgcmp ' + \
          '-t %s_treat_pileup.bdg ' % (prefix) + \
          '-c %s_control_lambda.bdg ' % (prefix) + \
          '-o %s_FE.bdg ' % (prefix) + \
          '-m FE'

    logger.info(command)
    returncode = utils.block_on(command)
    logger.info("MACS2 exited with returncode %d" % (returncode))
    assert returncode == 0, "MACS2 non-zero return"

    # Remove coordinates outside chromosome sizes (MACS2 bug)
    fc_bedgraph_fn = '%s.fc.signal.bedgraph' % (prefix)
    fc_bedgraph_sorted_fn = 'sorted-%s' % (fc_bedgraph_fn)
    fc_signal_fn = "%s.fc_signal.bw" % (prefix)
    steps = ['slopBed -i %s_FE.bdg -g %s -b 0' % (prefix, chrom_sizes),
             'bedClip stdin %s %s' % (chrom_sizes, fc_bedgraph_fn)]

    out, err = utils.run_pipe(steps)

    # Sort file
    out, err = utils.run_pipe([
        'bedSort %s %s' % (fc_bedgraph_fn, fc_bedgraph_sorted_fn)])

    # Convert bedgraph to bigwig
    command = 'bedGraphToBigWig ' + \
          '%s ' % (fc_bedgraph_sorted_fn) + \
          '%s ' % (chrom_sizes) + \
          '%s' % (fc_signal_fn)

    logger.info(command)
    returncode = utils.block_on(command)
    logger.info("bedGraphToBigWig exited with returncode %d" % (returncode))
    assert returncode == 0, "bedGraphToBigWig non-zero return"

    # For -log10(p-value) signal tracks

    # Compute sval =
    # min(no. of reads in ChIP, no. of reads in control) / 1,000,000
    out, err = utils.run_pipe(['gzip -dc %s' % (experiment), 'wc -l'])
    chip_reads = out.strip()
    out, err = utils.run_pipe(['gzip -dc %s' % (control), 'wc -l'])
    control_reads = out.strip()
    sval = str(min(float(chip_reads), float(control_reads)) / 1000000)

    logger.info("chip_reads = %s, control_reads = %s, sval = %s" %
                (chip_reads, control_reads, sval))

    command = 'macs2 bdgcmp ' + \
          '-t %s_treat_pileup.bdg ' % (prefix) + \
          '-c %s_control_lambda.bdg ' % (prefix) + \
          '-o %s_ppois.bdg ' % (prefix) + \
          '-m ppois -S %s' % (sval)

    logger.info(command)
    returncode = utils.block_on(command)
    assert returncode == 0, "MACS2 non-zero return"

    # Remove coordinates outside chromosome sizes (MACS2 bug)
    pvalue_bedgraph_fn = '%s.pval.signal.bedgraph' % (prefix)
    pvalue_bedgraph_sorted_fn = 'sort-%s' % (pvalue_bedgraph_fn)
    pvalue_signal_fn = "%s.pvalue_signal.bw" % (prefix)
    steps = ['slopBed -i %s_ppois.bdg -g %s -b 0' % (prefix, chrom_sizes),
             'bedClip stdin %s %s' % (chrom_sizes, pvalue_bedgraph_fn)]

    out, err = utils.run_pipe(steps)

    # Sort file
    out, err = utils.run_pipe([
        'bedSort %s %s' % (fc_bedgraph_fn, pvalue_bedgraph_sorted_fn)])

    # Convert bedgraph to bigwig
    command = 'bedGraphToBigWig ' + \
          '%s ' % (pvalue_bedgraph_sorted_fn) + \
          '%s ' % (chrom_sizes) + \
          '%s' % (pvalue_signal_fn)

    logger.info(command)
    returncode = utils.block_on(command)
    logger.info("bedGraphToBigWig exited with returncode %d" % (returncode))
    assert returncode == 0, "bedGraphToBigWig non-zero return"

    # Remove temporary files
    os.remove(clipped_narrowpeak_fn)
    os.remove(rescaled_narrowpeak_fn)


def main():
    args = get_args()
    tag = args.tag
    xcor = args.xcor
    con = args.con
    sample = args.sample
    genome_size = args.genome
    chrom_size = args.size
    paired = args.paired

    # Create a file handler
    handler = logging.FileHandler('call_peaks.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Calculate Cross-correlation if not already calcualted
    if xcor == 'Calculate':
        xcor_file = calculate_xcor(tag, paired)
    else:
        xcor_file = xcor

    # Call Peaks using MACS2
    call_peaks_macs(tag, xcor_file, con, sample, genome_size, chrom_size)


if __name__ == '__main__':
    main()
