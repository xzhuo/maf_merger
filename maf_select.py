import maf_iterate
import copy
from collections import OrderedDict


def main():
    args = maf_iterate.get_args()
    genomes = args.assemblies.split(",")  # The list to store all included genomes
    last_block = {}
    with open(args.out, 'w') as Out:
        for block in maf_iterate.maf_iterator(args.maf):
            curr_block = maf_iterate.clean_block(block, genomes)
            if 'req' in curr_block:
                if 'req' in last_block:
                    if len(curr_block['req']) == len(last_block['req']):
                        merged_block = {}
                        merged_block['req'] = OrderedDict()
                        merged_block['anno'] = {}
                        merge = 1
                        for assembly in curr_block['req']:
                            if assembly in last_block['req']:
                                mergability, merged_assembly = maf_iterate.strict_can_merge(last_block['req'][assembly], curr_block['req'][assembly])
                                if mergability == 1:
                                        merged_block['req'][assembly] = merged_assembly
                                else:
                                    merge = 0
                                    break
                            else:
                                merge = 0
                                break
                    else:
                        merge = 0
                else:
                    merge = 0
            else:
                merge = 0

            if merge == 0:
                maf_iterate.print_block(last_block, Out)
                last_block = copy.deepcopy(curr_block)
            else:
                merged_block['anno']['score'] = last_block['anno']['score'] + curr_block['anno']['score']
                last_block = copy.deepcopy(merged_block)
        maf_iterate.print_block(last_block, Out)

if __name__ == '__main__':
    main()
