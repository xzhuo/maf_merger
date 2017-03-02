import argparse
import warnings
import copy
from collections import OrderedDict
import logging
import ipdb


def _get_args():
    parser = argparse.ArgumentParser(description='simple arguments')
    parser.add_argument(
        '--inmaf',
        '-i',
        action="store",
        dest="maf",
        help='The input maf file.',
    )
    parser.add_argument(
        '--outmaf',
        '-o',
        action="store",
        dest="out",
        help='The onput maf file.',
    )
    parser.add_argument(
        '--assemblies',
        '-a',
        action="store",
        dest="assemblies",
        help="the genomes included in the output, comma seperated. The reference genome must be the first in the list.",
    )
    parser.add_argument(
        '--log',
        '-l',
        action="store_true",
        dest="log",
        default=False,
        help='print a log file',
    )
    parser.add_argument(
        '--directory',
        '-d',
        action="store",
        dest="dir",
        help='The directory with all the maf files,',
    )
    parser.add_argument(
        '--format',
        '-f',
        action="store",
        dest="format",
        help='The output format, could be maf or bed.',
    )
    parser.add_argument(
        '--gap',
        '-g',
        type=int,
        action="store",
        dest="gap",
        help='the threshold to separate small and big indels. blocks separated by small indels are merged, separated by big indels are kept separated.',
    )
    return parser.parse_args()


def main():
    args = _get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    anchor = genomes[0]
    if args.log:
        logging.basicConfig(filename=args.out + '.log', filemode='w', level=logging.INFO)
    holding_blocks = []
    for block in maf_iterator(args.maf):
        if _is_complete(block):
            curr_block = _clean_block(block, genomes)
            if len(holding_blocks) == 0:
                holding_blocks.append(curr_block)
            elif:
                # compare holding_block[-1] with curr_block
                # merge holding_block[-1] and curr_block into holding_block[-1] if mergable
                # if gap in any holding_block, evaluate if the gap is too long to hold. if yes, break.


def maf_iterator(in_maf):
    curr_block = {}
    with open(in_maf, 'r') as Fh:
        for num, line in enumerate(Fh, 1):
            if line.startswith('#'):
                # Out.write(line + "\n")
                continue
            linelist = line.split()
            if len(linelist) < 2:
                continue
            lead = linelist.pop(0)
            if lead == 'a':  # new alignmentblock

                if curr_block:
                    yield curr_block

                curr_block = {}
                for item in linelist:
                    try:
                        (key, value) = item.split("=")
                        curr_block['anno'] = {}
                        curr_block['anno'][key] = float(value)
                    except ValueError:
                        print("Found abnormal a line in %d!" % (num))
            elif lead == 's':
                (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                (start, length, strand, chrlenth, seq) = linelist
                if 'req' not in curr_block:
                    print("assign anchor here")
                    anchor = assembly  # assign anchor as the first assembly in block
                    curr_block['req'] = OrderedDict()
                curr_block['req'][assembly] = {}
                curr_block['req'][assembly]['aln'] = 1  # it is an alignment seq in the block
                curr_block['req'][assembly]['chrom'] = chrom
                curr_block['req'][assembly]['start'] = int(start)
                curr_block['req'][assembly]['length'] = int(length)
                curr_block['req'][assembly]['strand'] = strand
                curr_block['req'][assembly]['chrlenth'] = int(chrlenth)
                curr_block['req'][assembly]['seq'] = seq
            elif lead == 'i':
                (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                (leftStatus, leftCount, rightStatus, rightCount) = linelist
                curr_block['req'][assembly]['chrom'] = chrom
                curr_block['req'][assembly]['leftStatus'] = leftStatus
                curr_block['req'][assembly]['leftCount'] = int(leftCount)
                curr_block['req'][assembly]['rightStatus'] = rightStatus
                curr_block['req'][assembly]['rightCount'] = int(rightCount)
            elif lead == 'e':
                (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                (start, length, strand, chrlenth, gapStatus) = linelist
                curr_block['req'][assembly] = {}
                curr_block['req'][assembly]['aln'] = 0  # it is not an alignment seq in the block
                curr_block['req'][assembly]['chrom'] = chrom
                curr_block['req'][assembly]['start'] = int(start)
                curr_block['req'][assembly]['length'] = int(length)
                curr_block['req'][assembly]['strand'] = strand
                curr_block['req'][assembly]['chrlenth'] = int(chrlenth)
                curr_block['req'][assembly]['gapStatus'] = gapStatus
                curr_block['req'][assembly]['deletion'] = curr_block['req'][anchor]['length']
            elif lead == 'q':
                (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                quality = linelist[0]
                curr_block['req'][assembly]['quality'] = quality

        else:
            yield curr_block


def _is_complete(curr_block, genomes):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block['req']:
            complete = 0
            logging.info("%s not found" % assembly)
            break
    return complete


def _is_aln(curr_block, genomes):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block['req']:
            complete = 0
            logging.info("%s not found" % assembly)
            break
        elif curr_block['req'][assembly]['aln'] == 0:
            complete = 0
            logging.info("%s is a gap" % assembly)
    return complete


def _clean_block(block, genomes):
    clean_block = {}
    clean_block['anno'] = block['anno']
    for assembly in block['req']:
        if assembly in genomes:
            clean_block['req'][assembly] = block['req'][assembly]
    return clean_block


def _can_merge(last_assembly, curr_assembly):
    merged_assembly = {}
    if (curr_assembly['chrom'] == last_assembly['chrom'] and curr_assembly['strand'] == curr_assembly['strand']):
        if last_assembly['aln'] == 0 and curr_assembly['aln'] == 0:
            if (curr_assembly['start'] == last_assembly['start'] and
                    curr_assembly['length'] == last_assembly['length'] and
                    curr_assembly['strand'] == last_assembly['strand'] and
                    curr_assembly['gapStatus'] == last_assembly['gapStatus']):
                mergability = 1
                merged_assembly = merge_eachblocks(last_assembly, curr_assembly, 'gap')
            else:
                mergability = 0
        elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 1:
            if (last_assembly['rightStatus'] == 'C' and
                    last_assembly['rightCount'] == 0 and
                    curr_assembly['leftStatus'] == 'C' and
                    curr_assembly['leftCount'] == 0):
                if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
                    mergability = 1
                    merged_assembly = merge_eachblocks(last_assembly, curr_assembly, 'aln')
                else:
                    mergability = 0
            else:
                mergability = 0
        else:
            mergability = 0
    else:
        mergability = 0
    return mergability, merged_assembly


def merge_eachblocks(last_assembly, curr_assembly, kind):
    merged = {}
    for key in ('chrom', 'start', 'strand', 'chrlenth', 'aln'):
        merged[key] = last_assembly[key]
    if kind == 'anchor':
        merged['length'] = last_assembly['length'] + curr_assembly['length']
        merged['seq'] = last_assembly['seq'] + curr_assembly['seq']
        if 'quality' in last_assembly and 'quality' in curr_assembly:
            merged['quality'] = last_assembly['quality'] + curr_assembly['quality']

    if kind == 'aln':
        merged['length'] = last_assembly['length'] + last_assembly['rightCount'] + curr_assembly['length']
        merged['seq'] = last_assembly['seq'] + 'N' * last_assembly['rightCount'] + curr_assembly['seq']
        if 'quality' in last_assembly and 'quality' in curr_assembly:
            merged['quality'] = last_assembly['quality'] + '0' * last_assembly['rightCount'] + curr_assembly['quality']
        for key in ('leftStatus', 'leftCount'):
            merged[key] = last_assembly[key]
        for key in ('rightStatus', 'rightCount'):
            merged[key] = curr_assembly[key]

    if kind == 'mix':
        merged['length'] = last_assembly['length'] + curr_assembly['length']
        if curr_assembly['gapStatus'] == 'C':
            merged['seq'] = last_assembly['seq'] + '-' * curr_assembly['length']
        if curr_assembly['gapStatus'] == 'I' or curr_assembly['gapStatus'] == 'M':
            merged['seq'] = last_assembly['seq'] + 'N' * curr_assembly['length']
        if 'quality' in last_assembly and 'quality' in curr_assembly:
            merged['quality'] = last_assembly['quality'] + curr_assembly['quality']

        for key in ('leftStatus', 'leftCount'):
            merged[key] = last_assembly[key]
        merged['rightCount'] = last_assembly['rightCount'] - curr_assembly['length']
        if last_assembly['rightStatus'] == 'N' or last_assembly['rightStatus'] == 'n':
            merged['rightStatus'] = last_assembly['rightStatus']
        elif merged['rightCount'] > 0:
            merged['rightStatus'] = 'I'
        elif merged['rightCount'] == 0:
            merged['rightStatus'] = 'C'
        else:
            ipdb.set_trace()
            logging.error('length error during merge blocks!')
    if kind == 'mix2':
        merged['aln'] = curr_assembly['aln']
        if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
            merged['length'] = last_assembly['length'] + curr_assembly['length']
            if last_assembly['gapStatus'] == 'C':
                merged['seq'] = curr_assembly['seq']
            if last_assembly['gapStatus'] == 'I' or last_assembly['gapStatus'] == 'M':
                merged['seq'] = 'N' * last_assembly['length'] + curr_assembly['seq']
            if 'quality' in last_assembly and 'quality' in curr_assembly:
                merged['quality'] = last_assembly['quality'] + curr_assembly['quality']
        else:
            print('merge blocks failed!!')
        for key in ('rightStatus', 'rightCount'):
            merged[key] = curr_assembly[key]
        merged['leftCount'] = last_assembly['length'] - curr_assembly['leftCount']
        if curr_assembly['leftStatus'] == 'N' or curr_assembly['leftStatus'] == 'n':
            merged['leftStatus'] = curr_assembly['leftStatus']
        elif merged['leftCount'] > 0:
            merged['leftStatus'] = 'I'
        elif merged['leftCount'] == 0:
            merged['leftStatus'] = 'C'
        else:
            ipdb.set_trace()
            logging.error('length error during merge blocks!')
    if kind == 'gap':
        for key in ('gapStatus', 'length', 'quality'):
            if key in last_assembly and key in curr_assembly:
                merged[key] = curr_assembly[key]
        merged['deletion'] = last_assembly['deletion'] + curr_assembly['deletion']

    return merged


if __name__ == '__main__':
    main()
