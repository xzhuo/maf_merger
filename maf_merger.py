import maf_iterate
import copy
from collections import OrderedDict


def main():
    args = maf_iterate.get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    indel = args.gap
    holding_blocks = []
    with open(args.out, 'w') as Out:
        for block in maf_iterate.maf_iterator(args.maf):
            curr_block = maf_iterate.clean_block(block, genomes)
            if 'req' in curr_block:
                if len(holding_blocks) > 0:
                    last_block = holding_blocks[-1]
                    for assembly in genomes:
                        if assembly in curr_block['req']:
                            if assembly in last_block['req']:
                                mergability, indel_dict = maf_iterate.block_can_merge(last_block['req'][assembly], curr_block['req'][assembly], indel)
                                curr_block['req'][assembly]["deletion"] = indel_dict['deletion']
                                curr_block['req'][assembly]["insertion"] = indel_dict['insertion']
                                if mergability == 1:
                                    curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                                    curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num']

                                elif mergability == 2:
                                    curr_block['stat'][assembly]['tempstatus'] = 1
                                    curr_block['stat'][assembly]['tempnum'] = last_block['stat'][assembly]['num']
                                    curr_block['stat'][assembly]['status'] = 0
                                    curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] + 1
                                elif mergability == 3:
                                    curr_block['stat'][assembly]['status'] = 1
                                    curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num'] - 1
                                    i = -1
                                    while ('tempstatus' in holding_blocks[i]['stat'][assembly] and holding_blocks[i]['stat'][assembly]['tempstatus'] == 0):
                                        holding_blocks[i]['stat'][assembly]['status'] = 1
                                        holding_blocks[i]['stat'][assembly]['num'] -= 1
                                        holding_blocks[i]['stat'][assembly].pop('tempstatus')
                                        holding_blocks[i]['stat'][assembly].pop('tempnum')
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
                        if assembly in curr_block['req']:
                            curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                            curr_block['stat'][assembly]['num'] = 1
                        else:
                            curr_block['stat'][assembly]['status'] = "N"
                            curr_block['stat'][assembly]['num'] = 1
                    holding_blocks.append(curr_block)
        maf_iterate.print_blocks(holding_blocks, Out)


if __name__ == '__main__':
    main()
