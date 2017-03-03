import argparse
import warnings
import copy
from collections import OrderedDict
import logging
import ipdb
import sys
import maf_iterate
import numpy as np
import seaborn

def main():
    # maf = sys.argv[1]
    # anchor = sys.argv[2]
    anchor = 'hg19'
    files = ['chr' + str(x) + '.4species.strict.maf' for x in range(1, 22)]
    files.append('chrX.4species.strict.maf')
    files.append('chrY.4species.strict.maf')
    length = []
    ins = []
    dele = []
    for maf in files:
        print(maf)
        for block in maf_iterate.maf_iterator(maf[22]):
            length.append(block['req'][anchor][length])
            block_ins = 0
            block_del = 0
            for species in block['req']:
                if species['aln'] == 1:
                    if 'leftCount' in species and species['leftCount'] > block_ins:
                        block_ins = species['leftCount']
                elif species['aln'] == 0:
                    if species['deletion'] > block_del:
                        block_del = species['deletion']
            if block_ins:
                ins.append(block_ins)
            if block_del:
                dele.append(block_del)
    seaborn.set_style('whitegrid')
    lengthplot = seaborn.kdeplot(np.array(length))
    lengthplot.savefig("length.png")
    insplot = seaborn.kdeplot(np.array(ins))
    insplot.savefig("insertion.png")
    delplot = seaborn.kdeplot(np.array(dele))
    delplot.savefig("deletion.png")
if __name__ == '__main__':
    main()
