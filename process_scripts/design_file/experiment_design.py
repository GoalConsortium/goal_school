#!/usr/bin/env python3

'''Generate experiment design files for downstream processing.'''

import argparse
import logging
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
                        help="The design file to make experiemnts (tsv format).",
                        required=True)

    args = parser.parse_args()
    return args


def update_controls(design):
    '''Update design file to append controls list.'''

    logger.info("Running control file update.")

    file_dict = design[['sample_id', 'tag_align']] \
                .set_index('sample_id').T.to_dict()

    design['control_tag_align'] = design['control_id'] \
                                .apply(lambda x: file_dict[x]['tag_align'])

    logger.info("Removing rows that are there own control.")

    design = design[design['control_id'] != design['sample_id']]

    return design


def make_experiment_design(design):
    '''Make design file by grouping for each experiment'''

    logger.info("Running experiment design generation.")

    for experiment, df_experiment in design.groupby('experiment_id'):
        experiment_file = experiment + '.tsv'
        df_experiment.to_csv(experiment_file, header=True, sep='\t', index=False)


def main():
    args = get_args()
    design = args.design

    # Create a file handler
    handler = logging.FileHandler('experiment_generation.log')
    logger.addHandler(handler)

    # Read files as dataframes
    design_df = pd.read_csv(design, sep='\t')

    # Update design file for check_controls
    new_design_df = update_controls(design_df)

    # write out experiment design files
    make_experiment_design(new_design_df)


if __name__ == '__main__':
    main()
