from __future__ import print_function
import glob
import os
import sys
from random import shuffle
import pandas as pd
import numpy as np


def listFiles(path_to_fast5s):
    files = [x for x in glob.glob(path_to_fast5s) if os.stat(x).st_size != 0]
    shuffle(files)
    return files


def fivemer_histogram(kmer, path_to_fast5s, threshold, limit, outpath):
    def parse_table(assignments):
        try:
            data = pd.read_table(assignments,
                                 usecols=(0, 1, 2, 3),
                                 dtype={"5mer": np.str,
                                        "mean": np.float64,
                                        "prob": np.float64,
                                        "strand": np.str},
                                 header=None,
                                 names=['5mer', 'mean', 'prob', 'strand'])
            return data
        except Exception as e:
            print("Error opening file {filename}".format(filename=assignments), file=sys.stderr)
            return None

    def collect_means_with_threshold():
        files = listFiles(path_to_fast5s=path_to_fast5s)
        assert len(files) != 0, "Didn't find any files here {}".format(path_to_fast5s)
        total_template = 0
        total_complement = 0
        outfile = open(outpath + "{kmer}.histogram".format(kmer=kmer), 'a')
        for f in files:
            data = parse_table(f)
            if data is None:
                continue
            if total_template < limit:
                d = data.ix[(data['prob'] >= threshold) &
                            (data['5mer'] == kmer) &
                            (data['strand'] == 't')]
                d = d.drop('5mer', 1)
                d.to_csv(outfile, sep='\t', header=False, index=False)
                total_template += d.shape[0]
            if total_complement < limit:
                d = data.ix[(data['prob'] >= threshold) &
                            (data['5mer'] == kmer) &
                            (data['strand'] == 'c')]
                d = d.drop('5mer', 1)
                d.to_csv(outfile, sep='\t', header=False, index=False)
                total_complement += d.shape[0]
            if total_template >= limit and total_complement >= limit:
                break
        outfile.close()
        return total_template, total_complement
    t, c = collect_means_with_threshold()
    print("Got {t} template and {c} complement assignments for {kmer}".format(t=t, c=c, kmer=kmer))
    return
