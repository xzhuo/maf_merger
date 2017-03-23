import argparse
import warnings
import copy
from collections import OrderedDict
import logging
import ipdb


def get_args():
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
    parser.add_argument(  # TODO, incoporate this.
        '--level',
        '-v',
        action="store",
        dest="lvl",
        default=False,
        help='numeric logging level. 0 for NOTSET, 10 for DEBUG, 20 for INFO, 30 for WARNING, 40 for ERROR, 50 for CRITICAL.',
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
    args = get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    # anchor = genomes[0]
    if args.log:
        logging.basicConfig(filename=args.out + '.log', filemode='w', level=logging.INFO)
    holding_blocks = []
    for block in maf_iterator(args.maf):
        if _is_complete(block, genomes):
            curr_block = clean_block(block, genomes)
            if len(holding_blocks) == 0:
                holding_blocks.append(curr_block)
            else:
                # compare holding_block[-1] with curr_block
                # merge holding_block[-1] and curr_block into holding_block[-1] if mergable
                # if gap in any holding_block, evaluate if the gap is too long to hold. if yes, break.
                last_block = holding_blocks[-1]
                status = _compare_blocks(last_block, curr_block)  # return break, totally mergable or hold


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
                    # print("assign anchor here")
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


def clean_block(block, genomes):
    clean_block = {}
    clean_block['anno'] = block['anno']
    clean_block['req'] = OrderedDict()
    for assembly in block['req']:
        if assembly in genomes:
            clean_block['req'][assembly] = block['req'][assembly]

    if len(clean_block['req']) == 0:
        clean_block = {}
    return clean_block


# def _compare_blocks(last_block, curr_block):
#     if len(last_block['req']) == len(curr_block['req']):
#         for assembly in last_block['req']:
#             if assembly in curr_block['req']:
#                 slkdj
#             else:
#                 status = "break"
#     else:
#         status = "break"


# def _can_merge(last_assembly, curr_assembly):
#     merged_assembly = {}
#     if (curr_assembly['chrom'] == last_assembly['chrom'] and curr_assembly['strand'] == curr_assembly['strand']):
#         if last_assembly['aln'] == 0 and curr_assembly['aln'] == 0:
#             if (curr_assembly['start'] == last_assembly['start'] and
#                     curr_assembly['length'] == last_assembly['length'] and
#                     curr_assembly['strand'] == last_assembly['strand'] and
#                     curr_assembly['gapStatus'] == last_assembly['gapStatus']):
#                 mergability = 1
#                 merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'gap')
#             else:
#                 mergability = 0
#         elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 1:
#             if (last_assembly['rightStatus'] == 'C' and
#                     last_assembly['rightCount'] == 0 and
#                     curr_assembly['leftStatus'] == 'C' and
#                     curr_assembly['leftCount'] == 0):
#                 if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
#                     mergability = 1
#                     merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'aln')
#                 else:
#                     mergability = 0
#             else:
#                 mergability = 0
#         elif curr_assembly['aln'] == 0 and last_assembly['aln'] == 1:  # combine small indel gap block with aln block.
#             if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
#                 if curr_assembly['length'] < indel_length:
#                     status = 'diff10'
#                 else:
#                     status = 'b'
#                     logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#             elif last_assembly['start'] + last_assembly['length'] == curr_assembly['start'] + curr_assembly['length']:
#                 if curr_assembly['length'] < indel_length:
#                     status = 'diff10'
#                 else:
#                     status = 'b'
#                     logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#             else:
#                 status = 'b'
#                 ipdb.set_trace()
#                 logging.error("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#         elif curr_assembly['aln'] == 1 and last_assembly['aln'] == 0:  # combine small indel gap block with aln block.
#             if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
#                 if last_assembly['length'] < indel_length:
#                     status = 'diff01'
#                 else:
#                     status = 'b'
#                     logging.info("starting gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#             elif last_assembly['start'] + last_assembly['length'] == curr_assembly['start'] + curr_assembly['length']:
#                 logging.warning("do I need mix2 here?")
#                 if curr_assembly['length'] < indel_length:
#                     status = 'diff01'
#                 else:
#                     status = 'b'
#                     logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#             else:
#                 status = 'b'
#                 # ipdb.set_trace()
#                 logging.info("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#         else:
#             status = 'b'  # b for break
#             # ipdb.set_trace()
#             logging.info("breaking blocks! %s, %s, %d, %d and %d" % (assembly, curr_assembly['chrom'], curr_assembly['start'], last_assembly['aln'], curr_assembly['aln']))





#     else:
#         mergability = 0
#     return mergability, merged_assembly


def strict_can_merge(last_assembly, curr_assembly):
    merged_assembly = {}
    if (curr_assembly['chrom'] == last_assembly['chrom'] and curr_assembly['strand'] == curr_assembly['strand']):
        if last_assembly['aln'] == 0 and curr_assembly['aln'] == 0:
            if (curr_assembly['start'] == last_assembly['start'] and
                    curr_assembly['length'] == last_assembly['length'] and
                    curr_assembly['strand'] == last_assembly['strand'] and
                    curr_assembly['gapStatus'] == last_assembly['gapStatus']):
                mergability = 1
                merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'gap')
            else:
                mergability = 0
        elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 1:
            if 'rightStatus' in last_assembly:
                if (last_assembly['rightStatus'] == 'C' and
                        last_assembly['rightCount'] == 0 and
                        curr_assembly['leftStatus'] == 'C' and
                        curr_assembly['leftCount'] == 0):
                    if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
                        mergability = 1
                        merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'aln')
                    else:
                        mergability = 0
                else:
                    mergability = 0
            else:  # the anchor assembly
                if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
                    mergability = 1
                    merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'anchor')
                else:
                    mergability = 0

        else:
            mergability = 0
    else:
        mergability = 0
    return mergability, merged_assembly


def block_can_merge(last_assembly, curr_assembly, indel_length):
    merged_assembly = {}
    if (curr_assembly['chrom'] == last_assembly['chrom'] and curr_assembly['strand'] == curr_assembly['strand']):
        if last_assembly['aln'] == 0 and curr_assembly['aln'] == 0:
            if (curr_assembly['start'] == last_assembly['start'] and
                    curr_assembly['length'] == last_assembly['length'] and
                    curr_assembly['strand'] == last_assembly['strand'] and
                    curr_assembly['gapStatus'] == last_assembly['gapStatus']):
                mergability = 1
                indel_dict = {'insertion': curr_assembly['length'], 'deletion': last_assembly['deletion'] + curr_assembly['deletion']}
                merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'gap')
            else:
                mergability = 0
                indel_dict = {'insertion': curr_assembly['length'], 'deletion': curr_assembly['deletion']}
                logging.warning("Two different continuous e blocks?")
        elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 1:
            if 'rightStatus' in last_assembly:
                indel_dict = {'insertion': curr_assembly['rightCount'], 'deletion': 0}
                if (last_assembly['rightStatus'] == 'C' and
                        last_assembly['rightCount'] == 0 and
                        curr_assembly['leftStatus'] == 'C' and
                        curr_assembly['leftCount'] == 0):
                    if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
                        mergability = 1
                        merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'aln')
                    else:
                        mergability = 0
                        logging.warning("space between 2 continuous s blocks?")
                elif (last_assembly['rightStatus'] == 'I' and
                        last_assembly['rightCount'] < indel_length and
                        curr_assembly['leftStatus'] == 'I' and
                        curr_assembly['leftCount'] < indel_length and
                        last_assembly['start'] + last_assembly['length'] + last_assembly['rightCount'] == curr_assembly['start']):
                    mergability = 1
                    merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'aln')
                elif (last_assembly['rightStatus'] == 'M' and
                        last_assembly['rightCount'] < indel_length and
                        curr_assembly['leftStatus'] == 'M' and
                        curr_assembly['leftCount'] < indel_length and
                        last_assembly['start'] + last_assembly['length'] + last_assembly['rightCount'] == curr_assembly['start']):
                    mergability = 1
                    merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'aln')

                else:
                    mergability = 0
                    logging.warning("Big insertion in s blocks?")
            else:  # the anchor assembly
                indel_dict = {'insertion': 0, 'deletion': 0}
                if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
                    mergability = 1
                    merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'anchor')
                else:
                    mergability = 0
        elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 0:
            indel_dict = {'insertion': curr_assembly['length'], 'deletion': curr_assembly['deletion']}
            if curr_assembly['length'] < indel_length:
                if (last_assembly['start'] + last_assembly['length'] == curr_assembly['start']):
                    mergability = 2
                else:
                    mergability = 0
            else:
                mergability = 0
        elif last_assembly['aln'] == 0 and curr_assembly['aln'] == 1:
            # the hard work ahead here:
            indel_dict = {'insertion': curr_assembly['rightCount'], 'deletion': 0}
            if last_assembly['length'] < indel_length:
                if (last_assembly['start'] + last_assembly['length'] == curr_assembly['start']):
                    mergability = 3
                else:
                    mergability = 0
            else:
                mergability = 0

        else:
            mergability = 0
    else:
        mergability = 0
        if curr_assembly['aln'] == 1:
            if 'rightStatus' in last_assembly:
                indel_dict = {'insertion': curr_assembly['rightCount'], 'deletion': 0}
            else:
                indel_dict = {'insertion': 0, 'deletion': 0}
        elif curr_assembly['aln'] == 0:
            indel_dict = {'insertion': curr_assembly['length'], 'deletion': curr_assembly['deletion']}
    # return mergability
    return mergability, indel_dict


def merge_assemblies(last_assembly, curr_assembly, kind):
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
        if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
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
                logging.error('length error during merge blocks!')
        elif last_assembly['start'] + last_assembly['length'] == curr_assembly['start'] + curr_assembly['length']:
            merged = last_assembly
        else:
            logging.error('merge blocks failed in mix!')

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
        elif last_assembly['start'] == curr_assembly['start']:
            merged = curr_assembly
        else:
            logging.error('merge blocks failed in mix2!')
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


def print_blocks(blocks, indel_length, Out):  # combine blocks to units, then print them
    if len(blocks) > 0:
        units = []
        for curr_block in blocks:
            if len(units) > 0:
                last_unit = units[-1]
                if last_unit['stat'] == curr_block['stat']:
                    # merge
                    units[-1] = combine_blocks_to_units(last_unit, curr_block)
                else:
                    units.append(curr_block)
                    if len(units) > 2:
                        # merge 3 elements in units to 2 if len < 50, or pop and print out the first one
                        if unit_length(units[-3]) > 50 and unit_length(units[-2]) > 50:
                            print_block(units[-3], Out)
                            units.pop(-3)
                        else:
                            units = merge_units(units[-3], units[-2], units[-1], Out)



                    # if (unit_length(units[-2])) < indel_length:
                    #     ipdb.set_trace()
                    #     # merge_units(units[-3], units[-2], units[-1])
                    # # else:
                    # #     print_pop_unit(units[-3])
                    # print_block(units[-2], Out)
                    # units.pop(-2)
            else:
                units.append(curr_block)
        for unit in units:
            print_block(unit, Out)


def combine_blocks_to_units(unit, block):
    merged = {}
    merged['anno'] = {}
    merged['req'] = {}
    merged['stat'] = unit['stat']
    merged['anno']['score'] = unit['anno']['score'] + block['anno']['score']
    for assembly in unit['req']:
        merged['req'][assembly] = merge_all_assemblies(unit['req'][assembly], block['req'][assembly])
    return merged


def merge_all_assemblies(last_assembly, curr_assembly):
    if last_assembly['aln'] == 1 and curr_assembly['aln'] == 1:
        if 'leftStatus' in last_assembly:
            merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'aln')
        else:
            merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'anchor')
    elif last_assembly['aln'] == 0 and curr_assembly['aln'] == 0:
        merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'gap')
    elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 0:
        merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'mix')
    elif last_assembly['aln'] == 0 and curr_assembly['aln'] == 1:
        merged_assembly = merge_assemblies(last_assembly, curr_assembly, 'mix2')
    return merged_assembly


def unit_length(unit):
        try:
            return max(unit['req'][assembly]['length'] for assembly in unit['req'])
        except:
            ipdb.set_trace()


def merge_units(first_unit, mid_unit, last_unit, Out):  # TODO here returns a list with 2 units
    merged1 = first_unit
    merged2 = last_unit

    for assembly in mid_unit['stat']:
        if first_unit['stat'][assembly]['num'] == mid_unit['stat'][assembly]['num']:
            if first_unit['stat'][assembly]['status'] != 'N':
                merged1['req'][assembly] = merge_all_assemblies(first_unit['req'][assembly], mid_unit['req'][assembly])
        elif last_unit['stat'][assembly]['num'] == mid_unit['stat'][assembly]['num']:
            if last_unit['stat'][assembly]['status'] != 'N':
                merged2['req'][assembly] = merge_all_assemblies(mid_unit['req'][assembly], last_unit['req'][assembly])
        else:
            merged1 = mid_unit
            merged2 = last_unit
            print_block(first_unit, Out)
            logging.warning("hmmmmm")
            break

    return [merged1, merged2]


# def print_unit(unit):


def print_block(block, Out):
    refList = []
    alnList = []
    gapList = []
    if block:
        Out.write("a\tscore=%.6f\n" % (block['anno']['score']))
        for key in block['req']:
            if 'seq' in block['req'][key] and 'leftStatus' not in block['req'][key]:
                refList.append(key)
            if 'seq' in block['req'][key] and 'leftStatus' in block['req'][key]:
                alnList.append(key)
            if 'gapStatus' in block['req'][key]:
                gapList.append(key)
        # if len(refList) > 1:
        #     print("missing i line?")
        for key in refList:
            # Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block['req'][key]['chrom'], block['req'][key]['start'], block['req'][key]['length'], block['req'][key]['strand']))
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block['req'][key]['chrom'], block['req'][key]['start'], block['req'][key]['length'], block['req'][key]['strand'], block['req'][key]['chrlenth'], block['req'][key]['seq']))
        for key in alnList:
            # Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block['req'][key]['chrom'], block['req'][key]['start'], block['req'][key]['length'], block['req'][key]['strand']))
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block['req'][key]['chrom'], block['req'][key]['start'], block['req'][key]['length'], block['req'][key]['strand'], block['req'][key]['chrlenth'], block['req'][key]['seq']))
            Out.write("i\t%s.%s\t%s\t%d\t%s\t%d\n" % (key, block['req'][key]['chrom'], block['req'][key]['leftStatus'], block['req'][key]['leftCount'], block['req'][key]['rightStatus'], block['req'][key]['rightCount']))
        for key in gapList:
            # Out.write("%s.%s\t%d\t%d\t%s\t-\n" % (key, block['req'][key]['chrom'], block['req'][key]['start'], block['req'][key]['length'], block['req'][key]['strand']))
            Out.write("e\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block['req'][key]['chrom'], block['req'][key]['start'], block['req'][key]['length'], block['req'][key]['strand'], block['req'][key]['chrlenth'], block['req'][key]['gapStatus']))


        # for key in block:
        #     if key == 'score':
        #         continue
        #     Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))


        Out.write("\n")


if __name__ == '__main__':
    main()
