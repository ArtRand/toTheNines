from __future__ import print_function
import h5py
import glob
import os
import sys


def outputTables(fast5, outpath):
    def writeLine(strand):
        outfile.write(outline.format(kmer=_[4], eventMean=_[0], prob=_[7], strand=strand))

    try:
        f = h5py.File(fast5, 'r')
    except Exception as e:
        print("Error opening file {filename}".format(filename=fast5), file=sys.stderr)
        return None

    filename = fast5.split("/")[-1]

    template_address = "Analyses/Basecall_1D_000/BaseCalled_template/Events"

    complement_address = "Analyses/Basecall_1D_000/BaseCalled_complement/Events"

    outline = "{kmer}\t{eventMean}\t{prob}\t{strand}\n"
    outfile = open(outpath + filename + ".assignments", 'w')

    if template_address in f:
        for _ in f[template_address]:
            writeLine(strand='t')
    if complement_address in f:
        for _ in f[complement_address]:
            writeLine(strand='c')
    f.close()
