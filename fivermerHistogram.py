from __future__ import print_function
import h5py
import glob
import os
import sys
from random import shuffle


def fivemer_histogram(kmer, path_to_fast5s, threshold, outpath):
    def listFiles():
        files = [x for x in glob.glob(path_to_fast5s) if os.stat(x).st_size != 0]
        shuffle(files)
        return files

    def open_fast5(fast5):
        try:
            f = h5py.File(fast5, 'r')
            return f
        except Exception as e:
            print("Error opening file {filename}".format(filename=fast5), file=sys.stderr)
            return None

    def collect_means_with_threshold(address, limit=500):
        files = listFiles()
        assert len(files) != 0, "Didn't find any files here {}".format(path_to_fast5s)
        means = []
        for f in files:
            minion_read = open_fast5(f)
            if minion_read is None:
                continue
            means += [x[0] for x in minion_read[address] if x[4] == kmer and x[7] >= threshold]
            minion_read.close()
            if len(means) >= limit:
                break
        return means

    def write_means(means, filename):
        if len(means) == 0:
            print("Didn't find any means for {}".format(kmer))
            return
        with open(outpath + filename, 'w') as f:
            for u in means:
                f.write("{u}\t".format(u=u))
            f.write("\n")
            f.close()

    template_address = "Analyses/Basecall_1D_000/BaseCalled_template/Events"
    complement_address = "Analyses/Basecall_1D_000/BaseCalled_complement/Events"

    template_means = collect_means_with_threshold(template_address)
    complement_means = collect_means_with_threshold(complement_address)

    write_means(template_means, "{kmer}.template.tsv".format(kmer=kmer))
    write_means(complement_means, "{kmer}.complement.tsv".format(kmer=kmer))
    return
