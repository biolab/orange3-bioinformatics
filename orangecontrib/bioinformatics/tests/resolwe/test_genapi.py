import unittest

import Orange

from orangecontrib.bioinformatics.resolwe import genapi


class TestGenapi(unittest.TestCase):
    def test_login(self):
        from orangecontrib.bioinformatics import resolwe

        email = 'anonymous@genialis.com'
        password = 'anonymous'
        url = genapi.DEFAULT_URL

        self.assertTrue(resolwe.connect(email, password, url, 'genesis'))

    def test_objects(self):
        from orangecontrib.bioinformatics import resolwe
        from orangecontrib.bioinformatics.resolwe.utils import etc_to_table

        email = 'anonymous@genialis.com'
        password = 'anonymous'
        url = genapi.DEFAULT_URL

        gen = resolwe.connect(email, password, url, 'genesis')
        etc_objects = gen.fetch_etc_objects()

        self.assertEqual(type(etc_objects), list)  # test if return type is list
        self.assertTrue(etc_objects)  # test if list is not empty
        self.assertEqual(etc_objects[0].type, 'data:etc:')  # test if it contains correct objects

        #  Test experiment D. purpureu

        for obj in etc_objects:
            if obj.name == 'D. purpureum':
                test_experiment = obj

        self.assertTrue(test_experiment)
        self.assertEqual(test_experiment.id, 626)

        json, table_name = gen.download_etc_data(test_experiment.id)

        self.assertEqual(type(json), dict)
        self.assertEqual(len(json["etc"].keys()), 2)
        self.assertEqual(len(json["etc"]["genes"].keys()), 12410)
        self.assertEqual(len(json["etc"]["genes"]["DPU_G0071544"]), 7)
        self.assertEqual(json["etc"]["timePoints"], [0, 4, 8, 12, 16, 20, 24])

        self.assertEqual(type(etc_to_table(json)), Orange.data.table.Table)


if __name__ == '__main__':
    unittest.main()
