#!/usr/bin/env python3

'''Generate pooled and pseudoreplicate from data.'''

import argparse
import logging
import pandas as pd
import numpy as np
import os
import sys
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

    parser.add_argument('-d', '--design',
                        help="The design file to make experiemnts (tsv format).",
                        required=True)

    parser.add_argument('-p', '--paired',
                        help="True/False if paired-end or single end.",
                        default=False,
                        action='store_true')

    parser.add_argument('-c', '--cutoff',
                        help="Cutoff ratio used for choosing controls.",
                        default=1.2)

    args = parser.parse_args()
    return args


def check_replicates(design):
    '''Check the number of replicates for the experiment.'''

    no_rep = design.shape[0]

    return no_rep


def check_controls(design):
    '''Check the number of controls for the experiment.'''

    no_controls = len(design.control_tag_align.unique())

    return no_controls


def get_read_count_ratio(design):
    '''Get the ratio of read counts for unique controls.'''

    controls = design.control_tag_align.unique()
    control_dict = {}
    for con in controls:
        no_reads = utils.count_lines(con)
        control_dict[con] = no_reads

    control_matrix = {c: control_dict for c in controls}
    control_df = pd.DataFrame.from_dict(control_matrix)

    control_ratio = control_df.divide(list(control_dict.values()), axis=0)
    return control_ratio


def pool(tag_files, outfile, paired):
    '''Pool files together.'''

    if paired:
        file_extension = '.bedpe.gz'
    else:
        file_extension = '.bedse.gz'

    pooled_filename = outfile + file_extension

    # Merge files
    out, err = utils.run_pipe([
        'gzip -dc %s' % (' '.join(tag_files)),
        'gzip -cn'], outfile=pooled_filename)

    return pooled_filename


def self_psuedoreplication(tag_file, prefix, paired):
    '''Make 2 self-psuedoreplicates.'''

    # Get total number of reads
    no_lines = utils.count_lines(tag_file)

    # Number of lines to split into
    lines_per_rep = (no_lines+1)/2

    # Make an array of number of psuedoreplicatesfile names
    pseudoreplicate_dict = {r: prefix + '.pr' + str(r) + '.bedse.tagAlign.gz'
                            for r in [0, 1]}

    # Shuffle and split file into equal parts
    # by using the input to seed shuf we ensure multiple runs with the same
    # input will produce the same output
    # Produces two files named splits_prefix0n, n=1,2

    splits_prefix = 'temp_split'

    out, err = utils.run_pipe([
        'gzip -dc %s' % (tag_file),
        'shuf --random-source=%s' % (tag_file),
        'split -d -l %d - %s' % (lines_per_rep, splits_prefix)])

    # Convert read pairs to reads into standard tagAlign file

    for i, index in enumerate([0, 1]):
        string_index = '0' + str(index)
        steps = ['cat %s' % (splits_prefix + string_index)]
        if paired:
            steps.extend([r"""awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}'"""])
        steps.extend(['gzip -cn'])
        out, err = utils.run_pipe(steps, outfile=pseudoreplicate_dict[i])
        os.remove(splits_prefix + string_index)

    return pseudoreplicate_dict


def main():
    args = get_args()
    paired = args.paired
    design = args.design
    cutoff_ratio = args.cutoff

    # Create a file handler
    handler = logging.FileHandler('experiment_generation.log')
    logger.addHandler(handler)

    # Read files as dataframes
    design_df = pd.read_csv(design, sep='\t')

    # Get current directory to build paths
    cwd = os.getcwd()

    # Check Number of replicates and replicates
    no_reps = check_replicates(design_df)
    no_unique_controls = check_controls(design_df)

    if no_reps == 1:
        logger.info("No other replicate specified "
                    "so processing as an unreplicated experiment.")
        replicated = False

    else:
        logger.info("Multiple replicates specified "
                    "so processing as a replicated experiment.")
        replicated = True

    if no_unique_controls == 1 and replicated:
        logger.info("Only a single control was specified "
                    "so using same control for replicates, pool and psuedoreplicates.")
        single_control = True
    else:
        logger.info("Will merge only unique controls for pooled.")
        single_control = False

    # Pool the controls for checking
    if not single_control:
        control_df = get_read_count_ratio(design_df)
        control_files = design_df.control_tag_align.unique()
        pool_control = pool(control_files, "pool_control", paired)
    else:
        pool_control = design_df.control_tag_align.unique()[0]

    # Psuedoreplicate and update design accordingly
    if not replicated:

        # Duplicate rows and update for pool and psuedoreplicates
        experiment_id = design_df.at[0, 'experiment_id']
        replicate = design_df.at[0, 'replicate']
        design_new_df = design_df.loc[np.repeat(design_df.index, 4)].reset_index()
        design_new_df['replicate'] = design_new_df['replicate'].astype(str)
        design_new_df.at[1, 'sample_id'] = experiment_id + '_pr1'
        design_new_df.at[1, 'replicate'] = '1_pr'
        design_new_df.at[1, 'xcor'] = 'Calculate'
        design_new_df.at[2, 'sample_id'] = experiment_id + '_pr2'
        design_new_df.at[2, 'replicate'] = '2_pr'
        design_new_df.at[2, 'xcor'] = 'Calculate'
        design_new_df.at[3, 'sample_id'] = experiment_id + '_pooled'
        design_new_df.at[3, 'replicate'] = 'pooled'
        design_new_df.at[3, 'xcor'] = 'Calculate'

        # Make 2 self psuedoreplicates
        self_pseudoreplicates_dict = {}
        for rep, tag_file in zip(design_df['replicate'], design_df['tag_align']):
            replicate_prefix = experiment_id + '_' + str(rep)
            self_pseudoreplicates_dict = \
                self_psuedoreplication(tag_file, replicate_prefix, paired)

        # Update design to include new self pseudo replicates
        for rep, pseudorep_file in self_pseudoreplicates_dict.items():
            path_to_file = cwd + '/' + pseudorep_file
            replicate = rep + 1
            design_new_df.loc[replicate, 'tag_align'] = path_to_file

        # Drop index column
        design_new_df.drop(labels='index', axis=1, inplace=True)

    else:
        # Make pool of replicates
        replicate_files = design_df.tag_align.unique()
        experiment_id = design_df.at[0, 'experiment_id']
        pool_experiment = pool(replicate_files, experiment_id + "_pooled", paired)

        # Make self psuedoreplicates equivalent to number of replicates
        pseudoreplicates_dict = {}
        for rep, tag_file in zip(design_df['replicate'], design_df['tag_align']):
            replicate_prefix = experiment_id + '_' + str(rep)
            pr_dict = self_psuedoreplication(tag_file, replicate_prefix, paired)
            pseudoreplicates_dict[rep] = pr_dict

        # Merge self psuedoreplication for each true replicate
        pseudoreplicates_df = pd.DataFrame.from_dict(pseudoreplicates_dict)
        pool_pseudoreplicates_dict = {}
        for index, row in pseudoreplicates_df.iterrows():
            replicate_id = index + 1
            pr_filename = experiment_id + ".pr" + str(replicate_id) + '.bedse.gz'
            pool_replicate = pool(row, pr_filename, False)
            pool_pseudoreplicates_dict[replicate_id] = pool_replicate

        design_new_df = design_df
        # Check controls against cutoff_ratio
        # if so replace with pool_control
        # unless single control was used

        if not single_control:
            path_to_pool_control = cwd + '/' + pool_control
            if control_df.values.max() > 1.2:
                logger.info("Number of reads in controls differ by " +
                            " > factor of %f. Using pooled controls." % (cutoff_ratio))
                design_new_df['control_tag_align'] = path_to_pool_control
            else:
                for index, row in design_new_df.iterrows():
                    exp_no_reads = utils.count_lines(row['tag_align'])
                    con_no_reads = utils.count_lines(row['control_tag_align'])
                    if con_no_reads < exp_no_reads:
                        logger.info("Fewer reads in control than experiment." +
                                    "Using pooled controls for replicate %s."
                                    % row['replicate'])
                        design_new_df.loc[index, 'control_tag_align'] = \
                                                            path_to_pool_control
        else:
            path_to_pool_control = pool_control

        # Add in pseudo replicates
        tmp_metadata = design_new_df.loc[0].copy()
        tmp_metadata['control_tag_align'] = path_to_pool_control
        for rep, pseudorep_file in pool_pseudoreplicates_dict.items():
            tmp_metadata['sample_id'] = experiment_id + '_pr' + str(rep)
            tmp_metadata['replicate'] = str(rep) + '_pr'
            tmp_metadata['xcor'] = 'Calculate'
            path_to_file = cwd + '/' + pseudorep_file
            tmp_metadata['tag_align'] = path_to_file
            design_new_df = design_new_df.append(tmp_metadata)

        # Add in pool experiment
        tmp_metadata['sample_id'] = experiment_id + '_pooled'
        tmp_metadata['replicate'] = 'pooled'
        tmp_metadata['xcor'] = 'Calculate'
        path_to_file = cwd + '/' + pool_experiment
        tmp_metadata['tag_align'] = path_to_file
        design_new_df = design_new_df.append(tmp_metadata)

    # Write out new dataframe
    design_new_df.to_csv(experiment_id + '_ppr.tsv',
                         header=True, sep='\t', index=False)


if __name__ == '__main__':
    main()
