import maf_iterate


args = maf_iterate.get_args()
genomes = args.assemblies.split(",")  # The list to store all included genomes
with open(args.out, 'w') as Out:
    for block in maf_iterate.maf_iterator(args.maf):
        curr_block = maf_iterate.clean_block(block, genomes)
        if maf_iterate._is_complete(curr_block, genomes):
            maf_iterate.print_block(curr_block, Out)
