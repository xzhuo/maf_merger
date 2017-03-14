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
                                mergability, merged_assembly = maf_iterate.block_can_merge(last_block['req'][assembly], curr_block['req'][assembly], indel)
                                if mergability == 1:
                                    curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                                    curr_block['stat'][assembly]['num'] = last_block['stat'][assembly]['num']
                                elif mergability == 2:
                                    sdf
                                elif mergability == 3:
                                    ssf
                                else:
                                    break
                                    for assembly in genomes:
                                        if assembly in curr_block['req']:
                                            curr_block['stat'][assembly]['status'] = curr_block['req'][assembly]['aln']
                                            curr_block['stat'][assembly]['num'] = 1
                                        else:
                                            curr_block['stat'][assembly]['status'] = "N"
                                            curr_block['stat'][assembly]['num'] = 1
                                    maf_iterate.print_blocks(holding_blocks, Out)
                                    holding_blocks = []
                                    holding_blocks.append(curr_block)
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
                    maf_iterate.print_blocks(holding_blocks, Out)
                    holding_blocks = []
                    holding_blocks.append(curr_block)
        
        maf_iterate.print_blocks(holding_blocks, Out)


if __name__ == '__main__':
    main()
