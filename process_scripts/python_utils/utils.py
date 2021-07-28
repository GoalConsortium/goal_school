#!/usr/bin/env python3

'''General utilities.'''


import shlex
import logging
import subprocess
import sys
import os


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = True


def run_pipe(steps, outfile=None):
    # TODO:  capture stderr
    from subprocess import Popen, PIPE
    p = None
    p_next = None
    first_step_n = 1
    last_step_n = len(steps)
    for n, step in enumerate(steps, start=first_step_n):
        logger.debug("step %d: %s" % (n, step))
        if n == first_step_n:
            if n == last_step_n and outfile:  # one-step pipeline with outfile
                with open(outfile, 'w') as fh:
                    print("one step shlex: %s to file: %s" % (shlex.split(step), outfile))
                    p = Popen(shlex.split(step), stdout=fh)
                break
            print("first step shlex to stdout: %s" % (shlex.split(step)))

            p = Popen(shlex.split(step), stdout=PIPE)
        elif n == last_step_n and outfile:  # only treat the last step specially if you're sending stdout to a file
            with open(outfile, 'w') as fh:
                print("last step shlex: %s to file: %s" % (shlex.split(step), outfile))
                p_last = Popen(shlex.split(step), stdin=p.stdout, stdout=fh)
                p.stdout.close()
                p = p_last
        else:  # handles intermediate steps and, in the case of a pipe to stdout, the last step
            print("intermediate step %d shlex to stdout: %s" % (n, shlex.split(step)))
            p_next = Popen(shlex.split(step), stdin=p.stdout, stdout=PIPE)
            p.stdout.close()
            p = p_next
    out, err = p.communicate()
    return out, err


def block_on(command):
    process = subprocess.Popen(shlex.split(command), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    for line in iter(process.stdout.readline, b''):
        sys.stdout.write(line.decode('utf-8'))
    process.communicate()
    return process.returncode


def strip_extensions(filename, extensions):
    '''Strips extensions to get basename of file.'''

    basename = filename
    for extension in extensions:
        basename = basename.rpartition(extension)[0] or basename

    return basename


def count_lines(filename):
    import mimetypes
    compressed_mimetypes = [
        "compress",
        "bzip2",
        "gzip"
        ]
    mime_type = mimetypes.guess_type(filename)[1]
    if mime_type in compressed_mimetypes:
        catcommand = 'gzip -dc'
    else:
        catcommand = 'cat'
    out, err = run_pipe([
        '%s %s' % (catcommand, filename),
        'wc -l'
        ])
    return int(out)


def rescale_scores(filename, scores_col, new_min=10, new_max=1000):
    sorted_fn = 'sorted-%s' % (filename)
    rescaled_fn = 'rescaled-%s' % (filename)

    out, err = run_pipe([
        'sort -k %dgr,%dgr %s' % (scores_col, scores_col, filename),
        r"""awk 'BEGIN{FS="\t";OFS="\t"}{if (NF != 0) print $0}'"""],
        sorted_fn)

    out, err = run_pipe([
        'head -n 1 %s' % (sorted_fn),
        'cut -f %s' % (scores_col)])
    max_score = float(out.strip())
    logger.info("rescale_scores: max_score = %s" % (max_score))

    out, err = run_pipe([
        'tail -n 1 %s' % (sorted_fn),
        'cut -f %s' % (scores_col)])
    min_score = float(out.strip())
    logger.info("rescale_scores: min_score = %s" % (min_score))

    a = min_score
    b = max_score
    x = new_min
    y = new_max

    if min_score == max_score:  # give all peaks new_min
        rescale_formula = "x"
    else:  # n is the unscaled score from scores_col
        rescale_formula = "((n-a)*(y-x)/(b-a))+x"

    out, err = run_pipe(
        [
            'cat %s' % (sorted_fn),
            r"""awk 'BEGIN{OFS="\t"}{n=$%d;a=%d;b=%d;x=%d;y=%d}"""
            % (scores_col, a, b, x, y) +
            r"""{$%d=int(%s) ; print $0}'"""
            % (scores_col, rescale_formula)
        ],
        rescaled_fn)

    os.remove(sorted_fn)

    return rescaled_fn
