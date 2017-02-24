import sys
import argparse
import warnings
import copy
from collections import OrderedDict
import pprint
import ipdb

def main():
    args = _get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    curr_block = OrderedDict()
    last_block = OrderedDict()
    with open(args.maf, 'r') as Fh:
        with open(args.out, 'w') as Out:
            for num, line in enumerate(Fh, 1):
                if line.startswith('#'):
                    Out.write(line + "\n")
                    continue
                linelist = line.split()
                if len(linelist) < 2:
                    continue
                lead = linelist.pop(0)
                if lead == 'a':  # new alignmentblock
                    if _is_complete(curr_block, genomes) and _is_complete(last_block, genomes):
                        Two_blocks, merged_block = compare_merge_blocks(last_block, curr_block, genomes, Out, args.gap)
                    else:
                        Two_blocks = 'break'
                    if Two_blocks == 'new':
                        last_block = copy.deepcopy(curr_block)
                    elif Two_blocks == "aln":
                        last_block = copy.deepcopy(merged_block)
                    elif Two_blocks == "gap":
                        last_block = copy.deepcopy(merged_block)
                    else:
                        if _is_complete(last_block, genomes):
                            print_block(last_block, Out)
                        if _is_complete(curr_block, genomes):
                            last_block = copy.deepcopy(curr_block)
                    curr_block = OrderedDict()
                    for item in linelist:
                        try:
                            (key, value) = item.split("=")
                            curr_block[key] = float(value)
                        except ValueError:
                            print("Found abnormal a line in %d!" % (num))
                if lead == 's':
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
                if lead == 'i':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    if assembly not in genomes:
                        continue
                    (leftStatus, leftCount, rightStatus, rightCount) = linelist
                    curr_block[assembly]['chrom'] = chrom
                    curr_block[assembly]['leftStatus'] = leftStatus
                    curr_block[assembly]['leftCount'] = int(leftCount)
                    curr_block[assembly]['rightStatus'] = rightStatus
                    curr_block[assembly]['rightCount'] = int(rightCount)
                if lead == 'e':
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
                if lead == 'q':
                    (assembly, chrom) = linelist.pop(0).split(".", maxsplit=1)
                    quality = linelist[0]
                    curr_block[assembly]['quality'] = quality

            else:
                Two_blocks, merged_block = compare_merge_blocks(last_block, curr_block, genomes, Out, args.gap)
                if Two_blocks == "aln":
                    last_block = copy.deepcopy(merged_block)
                    print_block(last_block, Out)
                elif Two_blocks == "gap":
                    last_block = copy.deepcopy(merged_block)
                    print_block(last_block, Out)
                else:
                    print_block(last_block, Out)
                    print_block(curr_block, Out)


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
        '--reference',
        '-r',
        action="store",
        dest="ref",
        help='the reference genome used to anchor the alignment,',
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


def _is_complete(curr_block, genomes):
    complete = 1  # if curr_block contains all the species in genomes
    for assembly in genomes:
        if assembly not in curr_block:
            complete = 0
            break
    return complete


def compare_merge_blocks(last_block, curr_block, genomes, Out, indel_length):
    blockstatus = "aln"
    status = {}
    merge_block = OrderedDict()
    anchor = genomes[0]
    if (curr_block[anchor]['chrom'] == last_block[anchor]['chrom']
            and curr_block[anchor]['start'] == last_block[anchor]['start'] + last_block[anchor]['length']
            and curr_block[anchor]['strand'] == last_block[anchor]['strand']):
        status[anchor] = 'C'  # continue
        merge_block[anchor] = merge_eachblocks(last_block[anchor], curr_block[anchor], 'anchor')

        for assembly in genomes[1:]:
            if assembly not in curr_block or assembly not in last_block:
                ipdb.set_trace()
                print("hmmmm")
            if curr_block[assembly]['aln'] == 1 and last_block[assembly]['aln'] == 1:  # combine two aln blocks
                if curr_block[assembly]['leftStatus'] == 'C' and last_block[assembly]['rightStatus'] == 'C':
                    if curr_block[assembly]['leftCount'] > 0 or last_block[assembly]['rightCount'] > 0:
                        ipdb.set_trace()
                        print("Error! C not associated with 0 at %s, %d!" % (curr_block[anchor]['chrom'], curr_block[anchor]['start']))
                    elif last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                        status[assembly] = 'C'
                        merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                    else:
                        ipdb.set_trace()
                        print("something wrong with C aln lines")
                elif curr_block[assembly]['leftStatus'] == 'I' and last_block[assembly]['rightStatus'] == 'I':
                    if (curr_block[assembly]['leftCount'] == last_block[assembly]['rightCount']
                            and last_block[assembly]['start'] + last_block[assembly]['length'] + last_block[assembly]['rightCount'] == curr_block[assembly]['start']):
                        if curr_block[assembly]['leftCount'] > indel_length:
                            status[assembly] = 'I'
                        else:
                            status[assembly] = 'C'
                            merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                    else:
                        # import ipdb; ipdb.set_trace()
                        print("I warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                        status[assembly] = 'b'
                elif curr_block[assembly]['leftStatus'] == 'M' and last_block[assembly]['rightStatus'] == 'M':
                    if (curr_block[assembly]['leftCount'] == last_block[assembly]['rightCount']
                            and last_block[assembly]['start'] + last_block[assembly]['length'] + last_block[assembly]['rightCount'] == curr_block[assembly]['start']):
                        if curr_block[assembly]['leftCount'] > indel_length:
                            status[assembly] = 'M'
                        else:
                            status[assembly] = 'C'
                            merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'aln')
                    else:
                        # import ipdb; ipdb.set_trace()
                        print("M warning! Possible SV in %s? %s, %d, %s, %d!" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], curr_block[assembly]['chrom'], curr_block[assembly]['start']))
                        status[assembly] = 'b'
                else:
                    status[assembly] = 'e'
            elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 0:  # combine two gap blocks.
                if (curr_block[assembly]['chrom'] == last_block[assembly]['chrom']
                        and curr_block[assembly]['start'] == last_block[assembly]['start']
                        and curr_block[assembly]['length'] == last_block[assembly]['length']
                        and curr_block[assembly]['strand'] == last_block[assembly]['strand']
                        and curr_block[assembly]['gapStatus'] == last_block[assembly]['gapStatus']):
                    status[assembly] = curr_block[assembly]['gapStatus']
                    merge_block[assembly] = merge_eachblocks(last_block[assembly], curr_block[assembly], 'gap')
                else:
                    status[assembly] = 'b'
                    ipdb.set_trace()
                    # pprint.PrettyPrinter(indent=4).pprint(curr_block)
                    print("I don't expect this %s, %d" % (curr_block[anchor]['chrom'], curr_block[anchor]['start']))
            elif curr_block[assembly]['aln'] == 0 and last_block[assembly]['aln'] == 1:  # combine small indel gap block with aln block.
                if curr_block[assembly]['chrom'] == last_block[assembly]['chrom']:
                    if last_block[assembly]['start'] + last_block[assembly]['length'] == curr_block[assembly]['start']:
                        if curr_block[anchor]['length'] < indel_length or curr_block[assembly]['length'] < indel_length:
                            status[assembly] = 'C'
                            merge_eachblocks(last_block[assembly], curr_block[assembly], 'mix')
                        else:
                            status[assembly] = 'b'
                    else:
                        status[assembly] = 'b'
                        print("length discrepency! %s, %s. %d and %d is not %d" % (assembly, last_block[assembly]['chrom'], last_block[assembly]['start'], last_block[assembly]['length'], curr_block[assembly]['start']))
                else:
                    status[assembly] = 'b'
                    print("different chrom in %s, %s, %s!" % (assembly, last_block[assembly]['chrom'], curr_block[assembly]['chrom']))
            else:
                status[assembly] = 'b'  # b for break
        merge_block['score'] = last_block['score'] + curr_block['score']
    else:
        status[anchor] = 'b'  # break
    for assembly in genomes:
        try:
            if status[assembly] == 'b':
                blockstatus = 'break'
                break
            elif curr_block[assembly]['aln'] == 0:
                blockstatus = 'gap'
        except KeyError:
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
        merged['length'] = last_assembly['length'] + last_assembly['rightCount'] + curr_assembly['length']
        if curr_assembly['gapStatus'] == 'C':
            merged['seq'] = last_assembly['seq'] + 'N' * last_assembly['rightCount'] + '-' * curr_assembly['length']
        if curr_assembly['gapStatus'] == 'I' or curr_assembly['gapStatus'] == 'M':
            merged['seq'] = last_assembly['seq'] + 'N' * last_assembly['rightCount'] + 'N' * curr_assembly['length']
        if 'quality' in last_assembly and 'quality' in curr_assembly:
            merged['quality'] = last_assembly['quality'] + '0' * last_assembly['rightCount'] + curr_assembly['quality']
        for key in ('leftStatus', 'leftCount'):
            merged[key] = last_assembly[key]
            merged['rightCount'] = last_assembly['rightCount'] - curr_assembly['length']
            if merged['rightCount'] > 0:
                merged['rightStatus'] = 'I'
            elif merged['rightCount'] == 0:
                merged['rightStatus'] = 'C'
            else:
                ipdb.set_trace()
                print('length error during merge blocks!')
    if kind == 'gap':
        for key in ('gapStatus', 'length', 'quality'):
            if key in last_assembly and key in curr_assembly:
                merged[key] = curr_assembly[key]

    return merged


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
        Out.write("a\tscore: %.6f\n" % (block['score']))
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
        if len(refList) > 1: print("missing i line?")
        for key in refList:
            Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))
            #Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
        for key in alnList:
            Out.write("%s.%s\t%d\t%d\t%s\t=\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))
            #Out.write("s\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['seq']))
            #Out.write("i\t%s.%s\t%s\t%d\t%s\t%d\n" % (key, block[key]['chrom'], block[key]['leftStatus'], block[key]['leftCount'], block[key]['rightStatus'], block[key]['rightCount']))
        for key in gapList:
            Out.write("%s.%s\t%d\t%d\t%s\t-\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand']))
            # Out.write("e\t%s.%s\t%d\t%d\t%s\t%d\t%s\n" % (key, block[key]['chrom'], block[key]['start'], block[key]['length'], block[key]['strand'], block[key]['chrlenth'], block[key]['gapStatus']))
        Out.write("\n")

if __name__ == '__main__':
    main()
