import maf_iterate
import copy
from collections import OrderedDict
import ipdb
import logging


def main():
    args = maf_iterate.get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    indel = args.gap
    # level = args.lvl
    # logging.basicConfig(level=logging.WARNING)
    if args.log:
        logging.basicConfig(filename=args.out + '.log', filemode='w', level=logging.INFO)
    holding_blocks = []

    for block in maf_iterate.maf_iterator(args.maf):
        curr_block = maf_iterate.clean_block(block, genomes)
        # if maf_iterate._is_complete(block, genomes):
        if 'req' in curr_block:
            curr_block['stat'] = {}
            if len(holding_blocks) > 0:
                last_block = holding_blocks[-1]
                for assembly in genomes:
                    curr_block['stat'][assembly] = {}
                    if assembly in curr_block['req']:
                        if assembly in last_block['req']:
                            mergability, indel_dict = maf_iterate.block_can_merge(last_block['req'][assembly], curr_block['req'][assembly], indel)
                            curr_block['req'][assembly]["deletion"] = indel_dict['deletion']
                            curr_block['req'][assembly]["insertion"] = indel_dict['insertion']
                            if mergability == 1:
                                curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                                curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num']
                                if 'tempstatus' in last_block['stat'][assembly]:
                                    curr_block['stat'][assembly]['tempstatus'] = last_block['stat'][assembly]['tempstatus']
                                    curr_block['stat'][assembly]['tempnum'] = last_block['stat'][assembly]['tempnum']

                            elif mergability == 2:
                                curr_block['stat'][assembly]['tempstatus'] = 1
                                curr_block['stat'][assembly]['tempnum'] = last_block['stat'][assembly]['num']
                                curr_block['stat'][assembly]['status'] = 0
                                curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] + 1
                            elif mergability == 3:
                                # ipdb.set_trace()
                                curr_block['stat'][assembly]['status'] = 1
                                curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] - 1
                                i = -1
                                while ('tempstatus' in holding_blocks[i]['stat'][assembly] and holding_blocks[i]['stat'][assembly]['tempstatus'] == 1):
                                    holding_blocks[i]['stat'][assembly]['status'] = holding_blocks[i]['stat'][assembly].pop('tempstatus')
                                    holding_blocks[i]['stat'][assembly]['num'] = holding_blocks[i]['stat'][assembly].pop('tempnum')
                                    i -= 1
                            else:
                                curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                                curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] + 1
                        else:
                            curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                            curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] + 1
                    elif assembly not in last_block['req']:
                        curr_block['stat'][assembly]['status'] = 'N'
                        curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num']
                    else:
                        curr_block['stat'][assembly]['status'] = 'N'
                        curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] + 1

                holding_blocks.append(curr_block)

            else:
                for assembly in genomes:
                    curr_block['stat'][assembly] = {}
                    if assembly in curr_block['req']:
                        curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                        curr_block['stat'][assembly]['num'] = 1
                    else:
                        curr_block['stat'][assembly]['status'] = "N"
                        curr_block['stat'][assembly]['num'] = 1
                holding_blocks.append(curr_block)

    print("ok, everything in holding_blocks, now print them out:")
    with open(args.out, 'w') as Out:
        maf_iterate.print_blocks(holding_blocks, indel, Out)


if __name__ == '__main__':
    main()
