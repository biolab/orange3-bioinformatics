import unittest

import six

from orangecontrib.bioinformatics.kegg import pathway, databases


class TestGenome(unittest.TestCase):
    def test_genome(self):
        genome = databases.Genome()
        entry_keys = list(genome.keys())

        for key in entry_keys[:3] + entry_keys[-3:]:
            self.assertTrue(key in genome)
            self.assertTrue(key in genome)
            entry = genome[key]
            self.assertEqual(entry.entry_key, key)
            self.assertIsInstance(entry, genome.ENTRY_TYPE)
            self.assertIsInstance(entry.name, six.string_types)
            self.assertIsInstance(entry.taxid, six.string_types)

        self.assertTrue(genome.search("homo sapiens")[0] == "hsa")
        entry = genome['hsa']
        self.assertEqual(entry.taxid, "9606")


class TestGenes(unittest.TestCase):
    def _tester(self, org):
        genes = databases.Genes(org)
        keys = list(genes.keys())
        keys = keys[:3] + keys[-3:]
        all_entries = []
        for gene in keys:
            self.assertTrue(gene in genes)
            entry = genes[gene]
            self.assertEqual(entry.entry_key, genes.get(gene).entry_key, "__getitem__ and get return different result")

            self.assertTrue(gene.endswith(entry.entry_key))

            self.assertIsInstance(entry, genes.ENTRY_TYPE)
            self.assertIsInstance(entry.aliases(), list)

            self.assertTrue(all(isinstance(a, six.string_types) for a in entry.aliases()))
            all_entries.append(entry)

        self.assertSequenceEqual(
            [(e.name, e.entry_key) for e in all_entries],
            [(e.name, e.entry_key) for e in genes.batch_get(keys)],
            "batch_get returns different result",
        )

    def test_hsa(self):
        self._tester("hsa")

    def test_sce(self):
        self._tester("sce")

    def test_ddi(self):
        self._tester("ddi")


class TestPathways(unittest.TestCase):
    def _tester(self, path_id):
        pathways = databases.Pathway()

        entry = pathways[path_id]
        self.assertTrue(path_id.endswith(entry.entry_key))
        self.assertIsInstance(entry, pathways.ENTRY_TYPE)

        genes = entry.gene or []

        path = pathway.Pathway(path_id)
        self.assertEqual(sorted(genes), sorted(path.genes()))

    def test(self):
        self._tester("path:hsa05130")


class TestUtils(unittest.TestCase):
    def test_batch_iter(self):
        iter = range(25)
        expected = [list(range(10)), list(range(10, 20)), list(range(20, 25))]
        for exp, batch in zip(expected, databases.batch_iter(iter, 10)):
            self.assertEqual(exp, batch)


if __name__ == '__main__':
    unittest.main()
