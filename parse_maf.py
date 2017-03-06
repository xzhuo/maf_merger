import sys
import argparse
import warnings
import copy
from collections import OrderedDict
import ipdb
import logging


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
    curr_block = OrderedDict()
    last_block = OrderedDict()
    with open(args.maf, 'r') as Fh:
        with open(args.out, 'w') as Out:
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
                        last_block = do_strick_merge_stuff(last_block, curr_block, genomes, Out)

                    curr_block = OrderedDict()
                    for item in linelist:
                        try:
                            (key, value) = item.split("=")
                            curr_block[key] = float(value)
                        except ValueError:
                            print("Found abnormal a line in %d!" % (num))
                elif lead == 's':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (start, length, strand, chrlenth, seq) = linelist
                    curr_block[assembly] = {}
                    curr_block[assembly]['aln'] = 1  # it is an alignment seq in the block
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['start'] = int(start)
                    curr_block[assembly]['length'] = int(length)
                    curr_block[assembly]['strand'] = strand
                    curr_block[assembly]['chrlenth'] = int(chrlenth)
                    curr_block[assembly]['seq'] = seq
                elif lead == 'i':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (leftStatus, leftCount, rightStatus, rightCount) = linelist
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['leftStatus'] = leftStatus
                    curr_block[assembly]['leftCount'] = int(leftCount)
                    curr_block[assembly]['rightStatus'] = rightStatus
                    curr_block[assembly]['rightCount'] = int(rightCount)
                elif lead == 'e':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (start, length, strand, chrlenth, gapStatus) = linelist
                    curr_block[assembly] = {}
                    curr_block[assembly]['aln'] = 0  # it is not an alignment seq in the block
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['start'] = int(start)
                    curr_block[assembly]['length'] = int(length)
                    curr_block[assembly]['strand'] = strand
                    curr_block[assembly]['chrlenth'] = int(chrlenth)
                    curr_block[assembly]['gapStatus'] = gapStatus
                    curr_block[assembly]['deletion'] = curr_block[anchor]['length']
                elif lead == 'q':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    quality = linelist[0]
                    curr_block[assembly]['quality'] = quality

            else:
                last_block = do_strick_merge_stuff(last_block, curr_block, genomes, Out)
                print_block(last_block, Out)


def do_strick_merge_stuff(last_block, curr_block, genomes, Out):  # merge alignment blocks that are continued without any gap. C, 0 or same gap.
    if _is_complete(curr_block, genomes):
        # ipdb.set_trace()
        if _is_complete(last_block, genomes):
            anchor = genomes[0]
            merged_block = OrderedDict()
            merge = 1
            for assembly in genomes:
                if assembly == anchor:
                    merged_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'anchor')
                else:
                    mergability, merged_assembly = _can_merge(last_block[assembly], curr_block[assembly])
                    if mergability == 1:
                        merged_block[assembly] = merged_assembly
                    else:
                        merge = 0
                        break
            if merge == 0:
                print_block(last_block, Out)
                last_block = copy.deepcopy(curr_block)
            else:
                merged_block['score'] = last_block['score'] + curr_block['score']
                last_block = copy.deepcopy(merged_block)
        else:
            last_block = copy.deepcopy(curr_block)
    return last_block


    # get gap in between 2 blocks for each assembly
    # get gap in the curr_block for each assembly
    # add up gaps from last block to in between gap to curr_block gap
    # if gaps < the threshold,


def _is_complete(curr_block, genomes):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block:
            complete = 0
            logging.info("%s not found" % assembly)
            break
    return complete


def _is_aln(curr_block, genomes):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block:
            complete = 0
            logging.info("%s not found" % assembly)
            break
        elif curr_block[assembly]['aln'] == 0:
            complete = 0
            logging.info("%s is a gap" % assembly)
    return complete


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





# def _compare_assembly(last_assembly, curr_assembly, indel_length):
#     merged_assembly = OrderedDict()
#     if (curr_assembly['chrom'] == last_assembly['chrom'] and curr_assembly['strand'] == curr_assembly['strand']):
#         if last_assembly['aln'] == 0 and curr_assembly['aln'] == 0:
#             if (curr_assembly['start'] == last_assembly['start'] and
#                     curr_assembly['length'] == last_assembly['length'] and
#                     curr_assembly['strand'] == last_assembly['strand'] and
#                     curr_assembly['gapStatus'] == last_assembly['gapStatus']):
#                 status = 'same'
#                 merged_assembly = last_assembly
#             else:
#                 status = 'break'
#                 ipdb.set_trace()
#                 logging.error("I don't expect this %s, %d" % (curr_assembly['chrom'], curr_assembly['start']))
#         elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 1:
#             if curr_assembly['leftStatus'] == 'C':
#                 if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
#                     status = 'same'
#                 else:
#                     logging.warning("C warning! Possible SV in %s, %d, %s, %d!" % (last_assembly['chrom'], last_assembly['start'], curr_assembly['chrom'], curr_assembly['start']))
#                     status = 'break'
#             elif curr_assembly['leftStatus'] == 'I' or curr_assembly['leftStatus'] == 'M':
#                 if last_assembly['start'] + last_assembly['length'] + last_assembly['rightCount'] == curr_assembly['start']:
#                     status = 'same'
#                 else:
#                     logging.warning("%s warning! Possible SV in %s, %d, %s, %d!" % (curr_assembly['leftStatus'], last_assembly['chrom'], last_assembly['start'], curr_assembly['chrom'], curr_assembly['start']))
#                     status = 'break'
#             else:
#                 status = 'break'
#                 ipdb.set_trace()
#                 logging.info("Different Status? break blocks %s, %d, %s and %s" % (curr_assembly['chrom'], curr_assembly['start'], last_assembly['rightStatus'], curr_assembly['leftStatus']))

#         elif last_assembly['aln'] == 1 and curr_assembly['aln'] == 0:
#             status = 'break'
#         elif last_assembly['aln'] == 0 and curr_assembly['aln'] == 1:
#             if last_assembly['start'] + last_assembly['length'] == curr_assembly['start']:
#                 if curr_assembly['length'] < indel_length:
#                     status[assembly] = 'diff10'
#                 else:
#                     status[assembly] = 'b'
#                     logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#             elif last_assembly['start'] + last_assembly['length'] == curr_assembly['start'] + curr_assembly['length']:
#                 if curr_assembly['length'] < indel_length:
#                     status[assembly] = 'diff10'
#                 else:
#                     status[assembly] = 'b'
#                     logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))
#             else:
#                 status[assembly] = 'b'
#                 ipdb.set_trace()
#                 logging.error("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_assembly['chrom'], last_assembly['start'], last_assembly['length'], curr_assembly['start']))

#         else:
#             logging.error("I don't I am going to get here...what is 'aln' of last_assembly and curr_assembly?")
#     else:
#         status = 'break'
#     return status, merged_assembly


def compare_blocks(last_block, curr_block, genomes, indel_length):  # return "merge_now", "hold" or "break"
    anchor = genomes[0]
    status = {}
    merge_block = OrderedDict()
    if last_block:
        for assembly in curr_block:
            if assembly in last_block:
                if (curr_block[assembly]['chrom'] == last_block[assembly]['chrom'] and curr_block[anchor]['strand'] == last_block[anchor]['strand']):
                    # the basic continous requirement
                    if assembly == anchor:
                        if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                            merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'anchor')
                            merge_block[assembly]['insertion'] = 0
                            merge_block[assembly]['deletion'] = 0
                            merge_block[assembly]['total_insertion'] = last_block[assembly]['total_insertion']
                            merge_block[assembly]['total_deletion'] = last_block[assembly]['total_deletion']
                        # else:
                            # discontinue anchor?
                    elif curr_block[assembly]['aln'] and last_block[assembly]['aln']:
                        if curr_block[assembly]['leftStatus'] == 'C':
                            if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                                merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                                merge_block[assembly]['insertion'] = 0
                                merge_block[assembly]['deletion'] = 0
                                merge_block[assembly]['total_insertion'] = last_block[assembly]['total_insertion']
                                merge_block[assembly]['total_deletion'] = last_block[assembly]['total_deletion']
                            else:
                                logging.warning("C warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                        elif curr_block[assembly]['leftStatus'] == 'I' or curr_block[assembly]['leftStatus'] == 'M':
                            if curr_block[assembly]['leftCount'] < indel_length:
                                if last_block[assembly]['start'] + last_block[assembly]['length'] + last_block[assembly]['rightCount'] == curr_block[assembly]['start']:
                                    merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                                    merge_block[assembly]['insertion'] = 0
                                    merge_block[assembly]['deletion'] = 0
                                    merge_block[assembly]['total_insertion'] = last_block[assembly]['total_insertion'] + curr_block[assembly]['leftCount']
                                    merge_block[assembly]['total_deletion'] = last_block[assembly]['total_deletion']
                                else:
                                    logging.warning("%s warning! Possible SV in %s? %s, %d, %s, %d!" % (curr_block[assembly]['leftStatus'], assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                                    status[assembly] = 'b'
                            # else:
                                # a big insertion here
                        # elif curr_block[assembly]['leftStatus'] == 'N' or curr_block[assembly]['leftStatus'] == 'n':
                            # a new sequence
                            # ipdb.set_trace()
                        else:
                            ipdb.set_trace()
                    elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 0:
                        if curr_block[assembly]['length'] < indel_length:
                            if (curr_block[assembly]['start'] == last_block[assembly]['start']
                                    and curr_block[assembly]['length'] == last_block[assembly]['length']
                                    and curr_block[assembly]['strand'] == last_block[assembly]['strand']
                                    and curr_block[assembly]['gapStatus'] == last_block[assembly]['gapStatus']):
                                merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'gap')
                                merge_block[assembly]['insertion'] = curr_block[assembly]['length']
                            else:
                                status[assembly] = 'b'
                                ipdb.set_trace()
                                logging.error("I don't expect this %s, %d" % (curr_block[anchor]['chrom'], curr_block[anchor]['start']))
                        # else:
                            # a big insertion here

                # else:
                    # Do whatever necessary to last_block
            # else:
                # assembly is new in curr_block
    # else:
        #  whatever I should do if last_block is empty



    if _is_complete(last_block) and _is_complete(curr_block):
        for assembly in genomes:
            if (curr_block[assembly]['chrom'] == last_block[assembly]['chrom'] and curr_block[anchor]['strand'] == last_block[anchor]['strand']):
                if assembly == anchor:
                    if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                        status[assembly] = 'same'
                    else:
                        status[assembly] = 'b'
                elif curr_block[assembly]['aln'] == 1 and last_block[assembly]['aln'] == 1:
                    if curr_block[assembly]['leftStatus'] == 'C':
                        if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                            status[assembly] = 'same'
                        else:
                            logging.warning("C warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                            status[assembly] = 'b'
                    elif curr_block[assembly]['leftStatus'] == 'I' or curr_block[assembly]['leftStatus'] == 'M':
                        if last_block[assembly]['start'] + last_block[assembly]['length'] + last_block[assembly]['rightCount'] == curr_block[assembly]['start']:
                            status[assembly] = 'same'
                        else:
                            logging.warning("%s warning! Possible SV in %s? %s, %d, %s, %d!" % (curr_block[assembly]['leftStatus'], assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                            status[assembly] = 'b'
                    else:
                        status[assembly] = 'b'
                        ipdb.set_trace()
                        logging.info("Different Status? break blocks %s, %s, %d, %s and %s" % (assembly, curr_block[assembly]['chrom'], curr_block[assembly]['start'], last_block[assembly]['rightStatus'], curr_block[assembly]['leftStatus']))
                elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 0:
                    if (curr_block[assembly]['start'] == last_block[assembly]['start']
                            and curr_block[assembly]['length'] == last_block[assembly]['length']
                            and curr_block[assembly]['strand'] == last_block[assembly]['strand']
                            and curr_block[assembly]['gapStatus'] == last_block[assembly]['gapStatus']):
                        status[assembly] = 'same'
                    else:
                        status[assembly] = 'b'
                        ipdb.set_trace()
                        logging.error("I don't expect this %s, %d" % (curr_block[anchor]['chrom'], curr_block[anchor]['start']))
                elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 1:  # combine small indel gap block with aln block.
                    if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                        if curr_block[assembly]['length'] < indel_length:
                            status[assembly] = 'diff10'
                        else:
                            status[assembly] = 'b'
                            logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                    elif last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start'] + curr_block[assembly]['length']:
                        if curr_block[assembly]['length'] < indel_length:
                            status[assembly] = 'diff10'
                        else:
                            status[assembly] = 'b'
                            logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                    else:
                        status[assembly] = 'b'
                        ipdb.set_trace()
                        logging.error("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                elif curr_block[assembly]['aln'] == 1 and last_block[assembly]['aln'] == 0:  # combine small indel gap block with aln block.
                    if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                        if last_block[assembly]['length'] < indel_length:
                            status[assembly] = 'diff01'
                        else:
                            status[assembly] = 'b'
                            logging.info("starting gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                    elif last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start'] + curr_block[assembly]['length']:
                        logging.warning("do I need mix2 here?")
                        if curr_block[assembly]['length'] < indel_length:
                            status[assembly] = 'diff01'
                        else:
                            status[assembly] = 'b'
                            logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                    else:
                        status[assembly] = 'b'
                        # ipdb.set_trace()
                        logging.info("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                else:
                    status[assembly] = 'b'  # b for break
                    # ipdb.set_trace()
                    logging.info("breaking blocks! %s, %s, %d, %d and %d" % (assembly, curr_block[assembly]['chrom'], curr_block[assembly]['start'], last_block[assembly]['aln'], curr_block[assembly]['aln']))

            else:
                blockstatus = "break"
    else:
        blockstatus = "break"
    return blockstatus


def compare_merge_blocks(last_block, curr_block, genomes, indel_length):
    blockstatus = "merge"
    status = {}
    merge_block = OrderedDict()
    anchor = genomes[0]
    if last_block:
        if (curr_block[anchor]['chrom'] == last_block[anchor]['chrom']
                and curr_block[anchor]['start'] == last_block[anchor]['start'] + last_block[anchor]['length']
                and curr_block[anchor]['strand'] == last_block[anchor]['strand']):
            status[anchor] = 'same'  # continue
            merge_block[anchor] = merge_eachblocks(last_block[anchor], curr_block[anchor], 'anchor')

            for assembly in genomes[1:]:
                if assembly not in curr_block or assembly not in last_block:
                    ipdb.set_trace()
                    logging.error("hmmmm")
                if curr_block[assembly]['chrom'] == last_block[assembly]['chrom']:
                    if curr_block[assembly]['aln'] == 1 and last_block[assembly]['aln'] == 1:  # combine two aln blocks
                        if curr_block[assembly]['leftStatus'] == 'C':
                            if curr_block[assembly]['leftCount'] > 0:
                                ipdb.set_trace()
                                logging.error("Error! C not associated with 0 at %s, %d!" % (curr_block[anchor]['chrom'], curr_block[anchor]['start']))
                            elif last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                                status[assembly] = 'same'
                                merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                            else:
                                logging.warning("C warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                                status[assembly] = 'b'
                        elif curr_block[assembly]['leftStatus'] == 'I':
                            if (last_block[assembly]['start'] + last_block[assembly]['length'] + last_block[assembly]['rightCount'] == curr_block[assembly]['start']):
                                if curr_block[assembly]['leftCount'] > indel_length:
                                    status[assembly] = 'b'
                                else:
                                    status[assembly] = 'same'
                                    merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                            else:
                                # ipdb.set_trace()
                                logging.warning("I warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                                status[assembly] = 'b'
                        elif curr_block[assembly]['leftStatus'] == 'M':
                            if (last_block[assembly]['start'] + last_block[assembly]['length'] + last_block[assembly]['rightCount'] == curr_block[assembly]['start']):
                                if curr_block[assembly]['leftCount'] > indel_length:
                                    status[assembly] = 'b'
                                else:
                                    status[assembly] = 'same'
                                    merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                            else:
                                # ipdb.set_trace()
                                logging.warning("M warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                                status[assembly] = 'b'
                        else:
                            status[assembly] = 'b'
                            # ipdb.set_trace()
                            logging.info("Different Status? break blocks %s, %s, %d, %s and %s" % (assembly, curr_block[assembly]['chrom'], curr_block[assembly]['start'], last_block[assembly]['rightStatus'], curr_block[assembly]['leftStatus']))
                    elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 0:  # combine two gap blocks.
                        if (curr_block[assembly]['start'] == last_block[assembly]['start']
                                and curr_block[assembly]['length'] == last_block[assembly]['length']
                                and curr_block[assembly]['strand'] == last_block[assembly]['strand']
                                and curr_block[assembly]['gapStatus'] == last_block[assembly]['gapStatus']):
                            status[assembly] = 'same'
                            merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'gap')
                        else:
                            status[assembly] = 'b'
                            ipdb.set_trace()
                            logging.error("I don't expect this %s, %d" % (curr_block[anchor]['chrom'], curr_block[anchor]['start']))
                    elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 1:  # combine small indel gap block with aln block.
                        if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                            if curr_block[assembly]['length'] < indel_length:
                                status[assembly] = 'diff10'
                                # print("merge gap and aln! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                                merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'mix')
                            else:
                                status[assembly] = 'b'
                                logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                        elif last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start'] + curr_block[assembly]['length']:
                            if curr_block[assembly]['length'] < indel_length:
                                status[assembly] = 'diff10'
                                merge_block[assembly] = last_block[assembly]
                            else:
                                status[assembly] = 'b'
                                logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                        else:
                            status[assembly] = 'b'
                            ipdb.set_trace()
                            logging.error("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                    elif curr_block[assembly]['aln'] == 1 and last_block[assembly]['aln'] == 0:  # combine small indel gap block with aln block.
                        if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                            if last_block[assembly]['length'] < indel_length:
                                status[assembly] = 'diff01'
                                # print("merge gap and aln! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                                merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'mix2')
                            else:
                                status[assembly] = 'b'
                                logging.info("starting gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                        elif last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start'] + curr_block[assembly]['length']:
                            logging.warning("do I need mix2 here?")
                            if curr_block[assembly]['length'] < indel_length:
                                status[assembly] = 'diff01'
                                # merge_block[assembly] = last_block[assembly]
                            else:
                                status[assembly] = 'b'
                                logging.info("gap too big to merge! %s, %s. %d and %d to %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                        else:
                            status[assembly] = 'b'
                            # ipdb.set_trace()
                            logging.info("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                    else:
                        status[assembly] = 'b'  # b for break
                        # ipdb.set_trace()
                        logging.info("breaking blocks! %s, %s, %d, %d and %d" % (assembly, curr_block[assembly]['chrom'], curr_block[assembly]['start'], last_block[assembly]['aln'], curr_block[assembly]['aln']))
                else:
                    status[assembly] = 'b'
                    logging.info("different chromosomes! %s and %s" % (last_block[assembly]['chrom'], curr_block[assembly]['chrom']))

            merge_block['score'] = last_block['score'] + curr_block['score']
        else:
            status[anchor] = 'b'  # break
    else:  # return "merge" and curr_block if last_block is empty
        merge_block = curr_block
    if status:  # skip empty merge
        blocklist = []
        maxgap = 0
        last_maxgap = 0
        curr_maxgap = 0
        for assembly in genomes:
            if status[assembly] == 'b':
                blockstatus = 'break'
                break
            elif status[assembly] not in blocklist:
                blocklist.append(status[assembly])
            if curr_block[assembly]['aln']:
                curr_maxgap = max(curr_maxgap, curr_block[assembly]['length'])
            if last_block[assembly]['aln']:
                last_maxgap = max(last_maxgap, last_block[assembly]['length'])
        if not blockstatus == 'break':
            if "diff01" in blocklist:
                maxgap = max(maxgap, last_maxgap)
            if "diff10" in blocklist:
                maxgap = max(maxgap, curr_maxgap)
            # if len(blocklist) > 1 and maxgap > indel_length:
            if len(blocklist) > 1:
                if maxgap <= indel_length:
                    blockstatus = 'hold'
                else:
                    blockstatus = 'break'
    if not _is_complete(merge_block, genomes):
        if not blockstatus == 'break':
            ipdb.set_trace()
    return blockstatus, merge_block


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


def stat_block(curr_block, genomes):
    anchor = genomes[0]
    if anchor in curr_block and "length" in curr_block[anchor]:
        return curr_block[anchor]['length']
    else:
        logging.error("%s not found?" % anchor)


# def __merge_blocks(last_block, curr_block, genomes, Out):
#     if _is_complete(curr_block, genomes):  # if the curr_block contains all the species in alignment
#         if _is_continue(last_block, curr_block, genomes):  # all species in curr_block are continue, merge curr_block into last_block and return it.
#             last_block['score'] += curr_block['score']
#             for assembly in genomes:
#                 for key in ('length', 'seq', 'quality'):
#                     if key in last_block[assembly]:
#                         last_block[assembly][key] += curr_block[assembly][key]
#                 for key in ('rightStatus', "rightCount"):
#                     if key in last_block[assembly]:
#                         last_block[assembly][key] = curr_block[assembly][key]
#             return(last_block)
#         else:  # not all speceis are continue in the curr_block, print last_block, and return curr_block as next last_block.
#             print_block(last_block, Out)
#             return(curr_block)

#     else:  # skip curr_block and print last_block because some species in missing in the curr_block.
#         print_block(last_block, Out)
#         curr_block.clear()
#         return (curr_block)


def print_block(block, Out):
    refList = []
    alnList = []
    gapList = []
    if bool(block):
        Out.write("a\tscore=%.6f\n" % (block['score']))
        # Out.write("a\tscore= %.6f\n" % (6.66))
        for key in block:
            if key == 'score':
                continue
            if 'seq' in block[key] and 'leftStatus' not in block[key]:
                refList.append(key)
            if 'seq' in block[key] and 'leftStatus' in block[key]:
                alnList.append(key)
            if 'gapStatus' in block[key]:
                gapList.append(key)
        # if len(refList) > 1:
        #     print("missing i line?")
        for key in refList:
            # Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
        for key in alnList:
            # Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))
            Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
            Out.write("i\t%s.%s\t%s\t%d\t%s\t%d\n" % (key, block[key]['chrom'], block[key]['leftStatus'], block[key]['leftCount'], block[key]['rightStatus'], block[key]['rightCount']))
        for key in gapList:
            # Out.write("%s.%s\t%d\t%d\t%s\t-\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))
            Out.write("e\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['gapStatus']))


        # for key in block:
        #     if key == 'score':
        #         continue
        #     Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))


        Out.write("\n")

if __name__ == '__main__':
    main()
