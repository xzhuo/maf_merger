import maf_iterate
from intervaltree import Interval, IntervalTree
import logging


# start with finding indels between hg19 and mm10.

args = maf_iterate.get_args()
anchor = 'hg19'
source = 'mm10'
flank5k = 'chr1_hg19_mm10_5k.maf'
exact = 'chr1_hg19_mm10_0.maf'
flanking = 5000
curr_insertion = {}
curr_deletion = {}
insertiontree = IntervalTree()
deletiontree = IntervalTree()
for block in maf_iterate.maf_iterator(args.maf):
    if source in block['req']:
        # only + strain in anchor genome
        start = block['req'][anchor]['start']
        end = block['req'][anchor]['start'] + block['req'][anchor]['length']

        if block['req'][source]['aln'] == 1 and block['req'][source]['leftCount'] > 0:
            if 'end' in curr_insertion or 'start' not in curr_insertion:
                logging.warning("another half of insertion incorrect!")
                curr_insertion = {}
            else:
                curr_insertion['end'] = end
                curr_insertion['part2'] = block
                insertiontree[curr_insertion['start']:curr_insertion['end']] = curr_insertion
                curr_insertion = {}
        elif block['req'][source]['aln'] == 1 and block['req'][source]['rightCount'] > 0:
            if curr_insertion == {}:
                curr_insertion['start'] = start
                curr_insertion['part1'] = block
            else:
                logging.warning("insertion part already exist?")
                curr_insertion = {}

        elif block['req'][source]['aln'] == 0:
            if 'part' in curr_deletion and block['req'][source] == curr_deletion['part']['req'][source]:
                deletiontree.removei(curr_deletion['start'], curr_deletion['end'], curr_deletion)
                curr_deletion['end'] = end
                deletiontree[curr_deletion['start']:curr_deletion['end']] = curr_deletion
            else:
                curr_deletion['start'] = start
                curr_deletion['end'] = end
                curr_deletion['part'] = block
                deletiontree[curr_deletion['start']:curr_deletion['end']] = curr_deletion

with open(flank5k, 'w') as flanking_outfile:
    with open(exact, 'w') as exact_outfile:
        for Iobj in sorted(insertiontree):
            Dset = sorted(deletiontree.search(Iobj.begin, Iobj.end))
            closeset = sorted(deletiontree.search(Iobj.begin - flanking, Iobj.end + flanking))
            for closeobj in closeset:
                if closeobj in Dset:
                    maf_iterate.print_block(Iobj.data['part1'], flanking_outfile)
                    maf_iterate.print_block(Iobj.data['part2'], flanking_outfile)
                    maf_iterate.print_block(closeobj.data['part'], flanking_outfile)
                else:
                    maf_iterate.print_block(Iobj.data['part1'], exact_outfile)
                    maf_iterate.print_block(Iobj.data['part2'], exact_outfile)
                    maf_iterate.print_block(closeobj.data['part'], exact_outfile)
