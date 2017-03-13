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
                    if ('req' in last_block and
                            len(curr_block['req']) == len(last_block['req'])):
                        merged_block = {}
                        merged_block['req'] = OrderedDict()
                        merged_block['anno'] = {}
                        merge = 1
                        for assembly in curr_block['req']:
                            if assembly in last_block['req']:
                                holding_assemblies = [x['req'][assembly] for x in holding_blocks]
                                mergability, merged_assembly = maf_iterate.indel_can_merge(holding_assemblies, curr_block['req'][assembly], indel)
                                if mergability == 1:
                                    merged_block['req'][assembly] = merged_assembly
                                elif mergability == 2:
                                    merge = 2
                                else:
                                    merge = 0
                                    break
                            else:
                                merge = 0
                                break
                    else:
                        merge = 0

                    if merge == 0:
                        # print all the blocks in holding block, and initiate holding_block, append curr_block to it
                        for block in holding_blocks:
                            maf_iterate.print_block(block, Out)

                        holding_blocks = []
                        last_block = copy.deepcopy(curr_block)
                    elif merge == 2:
                        append_block = copy.deepcopy(curr_block)
                        holding_blocks.append(append_block)
                    else:
                        merged_block['anno']['score'] = last_block['anno']['score'] + curr_block['anno']['score']
                        last_block = copy.deepcopy(merged_block)
                else:
                    holding_blocks.append(curr_block)
        maf_iterate.print_block(last_block, Out)

if __name__ == '__main__':
    main()
