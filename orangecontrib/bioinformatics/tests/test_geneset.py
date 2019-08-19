import os
import unittest
from tempfile import mkstemp

from orangecontrib.bioinformatics.geneset import GeneSet, GeneSets, GeneSetException, filename, filename_parse


class TestGeneSets(unittest.TestCase):
    test_gs_id = 'test_gs'
    test_name = 'test_name'
    test_organism = '9606'
    test_hierarchy = ('GO', 'biological_process')
    test_file = 'GO-biological_process-9606.gmt'
    test_genes = [123, 321, 111, 222, 333]  # ncbi-like ids

    def test_file_name(self):
        file_name = filename(self.test_hierarchy, self.test_organism)
        self.assertEqual(file_name, self.test_file)

        hierarchy, org = filename_parse(file_name)
        self.assertEqual(hierarchy, self.test_hierarchy)
        self.assertEqual(org, self.test_organism)

    def test_gmt_file_format(self):

        gs = GeneSet(
            gs_id=self.test_gs_id,
            name=self.test_name,
            genes=self.test_genes,
            hierarchy=self.test_hierarchy,
            organism=self.test_organism,
            link='',
        )

        fd, file_name = mkstemp()

        # write to file
        write_sets = GeneSets([gs])
        write_sets.to_gmt_file_format(file_name)

        with open(file_name, 'r') as temp_f:
            line = temp_f.readline()
            columns = line.strip().split('\t')
            self.assertGreater(len(columns), 0)

        # read from file
        read_sets = GeneSets.from_gmt_file_format(file_name)
        self.assertIsNotNone(read_sets)
        self.assertGreater(len(read_sets), 0)
        self.assertEqual(read_sets.common_hierarchy(), self.test_hierarchy)
        self.assertEqual(read_sets.common_org(), self.test_organism)

        # clean-up
        os.close(fd)
        os.remove(file_name)

    def test_gene_set(self):

        gs1 = GeneSet(
            gs_id=self.test_gs_id,
            name=self.test_name,
            genes=self.test_genes,
            hierarchy=self.test_hierarchy,
            organism=self.test_organism,
            link='',
        )

        gs2 = GeneSet(gs_id='test2', name='test_name2')

        self.assertEqual(gs1.gmt_description(), 'test_gs,GO-biological_process,9606,test_name,_,_,_')
        self.assertFalse(gs1 == gs2)
        self.assertNotEqual(gs1, gs2)
        self.assertTrue(gs1 == gs1)

    def test_gene_sets(self):
        gs1 = GeneSet(
            gs_id=self.test_gs_id,
            name=self.test_name,
            genes=self.test_genes,
            hierarchy=self.test_hierarchy,
            organism=self.test_organism,
            link='',
        )

        gs2 = GeneSet(gs_id='test2', name='test_name2', hierarchy=('Test', 'test'), organism='3702')
        gs3 = GeneSet(gs_id='test3', name='test_name3', hierarchy=('Test', 'test'), organism='3702')

        sets = GeneSets([gs1, gs2, gs3])
        self.assertIsNotNone(sets)

        self.assertRaises(GeneSetException, sets.common_org)
        self.assertRaises(GeneSetException, sets.common_hierarchy)

        self.assertGreater(len(sets.hierarchies()), 1)

        split_by_hierarchy = sets.split_by_hierarchy()
        self.assertLess(len(split_by_hierarchy), len(sets))


if __name__ == '__main__':
    unittest.main()
