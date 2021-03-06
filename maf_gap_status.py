import argparse
import warnings
import copy
from collections import OrderedDict
import logging
import ipdb
import sys
import matplotlib
matplotlib.use('Agg')
import pylab as plt
matplotlib.style.use('ggplot')
import maf_iterate
import numpy as np
import seaborn


def main():
    # maf = sys.argv[1]
    # anchor = sys.argv[2]
    anchor = 'hg19'
    files = ['/bar/twlab-shared/Epigenome_Evolution/genomes/maf_100way/maf/chr' + str(x) + '.4species.strict2.maf' for x in range(1, 22)]
    files.append('/bar/twlab-shared/Epigenome_Evolution/genomes/maf_100way/maf/chrX.4species.strict2.maf')
    files.append('/bar/twlab-shared/Epigenome_Evolution/genomes/maf_100way/maf/chrY.4species.strict2.maf')
    length = []
    ins = []
    dele = []
    for maf in files:
        for block in maf_iterate.maf_iterator(maf):
            block_len = block['req'][anchor]['length']
            block_ins = 0
            block_del = 0
            for species in block['req']:
                # ipdb.set_trace()
                species_dict = block['req'][species]
                if species_dict['aln'] == 1:
                    if 'leftCount' in species_dict and species_dict['leftCount'] > block_ins:
                        block_ins = species_dict['leftCount']
                elif species_dict['aln'] == 0:
                    if species_dict['deletion'] > block_del:
                        block_del = species_dict['deletion']
            length.append(block_len)
            if block_ins:
                ins.append(block_ins)
            if block_del:
                dele.append(block_del)

    ipdb.set_trace()
    seaborn.set()
    lengthplot = seaborn.kdeplot(np.array(length))
    lengthplot.set(xlim=(0, 1000))
    fig = lengthplot.get_figure()
    fig.savefig("all.4species_1klength.png")
    plt.clf()
    insplot = seaborn.kdeplot(np.array(ins))
    insplot.set(xlim=(0, 1000))
    fig = insplot.get_figure()
    fig.savefig("all.4species_1kinsertion.png")
    plt.clf()
    delplot = seaborn.kdeplot(np.array(dele))
    delplot.set(xlim=(0, 1000))
    fig = delplot.get_figure()
    fig.savefig("all.4species_1kdeletion.png")
    plt.close(fig)


if __name__ == '__main__':
    main()
