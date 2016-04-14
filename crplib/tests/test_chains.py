

import re as re
import numpy as np
import unittest as unittest
import itertools as itt

from crplib.auxiliary.text_parsers import get_chain_iterator, _read_chain_header
from crplib.commands.convert_chain import build_index_structures


class TestChainHandling(unittest.TestCase):
    def setUp(self):
        self.chainfile = TESTCHAINS.split('\n')
        self.target_sizes = {'chr16': 1000, 'chr21': 1000, 'chr14': 1000, 'chrX': 155270560}
        self.query_sizes = {'chr8': 2500, 'chr10': 1000, 'chr14': 4050, 'chr2': 10000, 'chrX_random': 1785075}
        self.qcheck = re.compile('(chr)?[0-9]+(\s|$)')
        self.splits = {'chr16': [581, 582, 756, 765, 862, 863, 948, 949],
                       'chr21': [826, 831, 844, 845, 875, 878, 884, 886, 914, 916],
                       'chr14': [845, 848, 890, 892, 893, 895, 912, 915, 218, 219, 259, 260, 397, 398, 498, 508]}

    def tearDown(self):
        pass

    def test_target_missing(self):

        chainit = get_chain_iterator(self.chainfile, tselect='chr1')
        with self.assertRaises(AssertionError):
            for block in chainit:
                pass
        return

    def test_parse_chainheader(self):

        for line in self.chainfile:
            if not line or not line.startswith('chain'):
                continue
            chain = _read_chain_header(line)
            tchrom = chain[0]
            tsize = chain[1]
            qchrom = chain[5]
            qsize = chain[6]
            self.assertIn(tchrom, self.target_sizes)
            self.assertIn(qchrom, self.query_sizes)
            self.assertEqual(tsize, self.target_sizes[tchrom])
            self.assertEqual(qsize, self.query_sizes[qchrom])
        return

    def test_make_target_mask(self):

        for chrom in self.target_sizes.keys():
            if chrom == 'chrX':
                # this is checked in test_target_missing
                continue
            csize = self.target_sizes[chrom]
            chainit = get_chain_iterator(self.chainfile, tselect=chrom, qcheck=self.qcheck)
            mask, splits, select = build_index_structures(chainit, csize)
            self.assertListEqual(sorted(splits), sorted(self.splits[chrom]), 'Splits mismatch for chrom {}'.format(chrom))
            self.assertEqual(0, select[0], 'Select not 0 at beginning')
            self.assertEqual(1, select[1], 'Select not 1 at position 1')
            self.assertEqual(0, select[-1], 'Select not 0 at last position')
            indices = np.arange(csize, dtype=np.int32)
            check_indices = []
            select_blocks = 0
            for i in range(0, len(splits), 2):
                check_indices.extend(list(range(splits[i], splits[i+1], 1)))
                select_blocks += 1
            select_indices = np.compress(~mask, indices)
            self.assertListEqual(sorted(select_indices), sorted(check_indices), 'Wrong index set masked regions')
            self.assertEqual(select_blocks, np.sum(select), 'Mismatch for number of selected blocks')
        return

# chains extracted from hg19_to_mm9.rbest.chain.gz, chromosome sizes adapted
# to be modest on memory consumption for the test

TESTCHAINS = "chain 12 chr16 1000 + 581 949 chr8 2500 - 1820 2341 612448\n\
1   174 305\n\
4   0   1\n\
4   0   1\n\
1   97  95\n\
1   85  107\n\
1\n\
\n\
chain 0 chr21 1000 + 826 916 chr10 1000 + 555 644 50921\n\
5   13  11\n\
1   30  33\n\
3   6   6\n\
2   28  26\n\
2\n\
\n\
\n\
chain 8 chr14 1000 + 845 915 chr14 4050 + 3937 4001 256715\n\
3   42  41\n\
2   1   0\n\
2   17  13\n\
3\n\
\n\
chain 44 chr14 1000 + 218 508 chr2 10000 + 8864 9142 188346\n\
1       40      40\n\
1       137     124\n\
1       100     101\n\
10\n\
\n\
\n\
chain 9 chrX 155270560 + 114559207 114559208 chrX_random 1785075 - 1331258 1331259 326536\n\
1"