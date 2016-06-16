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
        total = 0
        outfile = open(outpath + "{kmer}.histogram".format(kmer=kmer), 'a')
        for f in files:
            data = parse_table(f)
            if data is None:
                continue
            d = data.ix[(data['prob'] >= threshold) &
                        (data['5mer'] == kmer)]
            d = d.drop('5mer', 1)
            d.to_csv(outfile, sep='\t', header=False, index=False)
            total += d.shape[0]
            if total >= limit:
                break
        outfile.close()
        return total
    t = collect_means_with_threshold()
    print("Got {n} assignments for {kmer}".format(n=t, kmer=kmer))
    return
