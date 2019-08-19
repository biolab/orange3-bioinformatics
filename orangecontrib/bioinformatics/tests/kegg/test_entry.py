import unittest

from six import StringIO

from orangecontrib.bioinformatics.kegg.entry import DBEntry, parser, entry_decorate

TEST_ENTRY = """\
ENTRY       test_id    something else
NAME        test
DESCRIPTION This is a test's description.
            It spans
            multiple lines
  SUB       This is a description's sub
            section
///
"""


@entry_decorate
class Entry(DBEntry):
    pass


class TestEntry(unittest.TestCase):
    def test_entry(self):
        """
        Test basic DBEntry class.
        """
        entry = Entry(TEST_ENTRY)
        self.assertEqual(entry.entry_key, "test_id")
        self.assertEqual(entry.ENTRY.TITLE, "ENTRY")
        self.assertEqual(str(entry), TEST_ENTRY[:-4])


class TestParser(unittest.TestCase):
    def test_parser(self):
        parse = parser.DBGETEntryParser()
        stream = StringIO(TEST_ENTRY)

        expected = [
            (parse.ENTRY_START, None, None),
            (parse.SECTION_START, "ENTRY", "test_id    something else\n"),
            (parse.SECTION_END, "ENTRY", None),
            (parse.SECTION_START, "NAME", "test\n"),
            (parse.SECTION_END, "NAME", None),
            (parse.SECTION_START, "DESCRIPTION", "This is a test's description.\n"),
            (parse.TEXT, None, "It spans\n"),
            (parse.TEXT, None, "multiple lines\n"),
            (parse.SUBSECTION_START, "SUB", "This is a description's sub\n"),
            (parse.TEXT, None, "section\n"),
            (parse.SUBSECTION_END, "SUB", None),
            (parse.SECTION_END, "DESCRIPTION", None),
            (parse.ENTRY_END, None, None),
        ]
        self.assertSequenceEqual(list(parse.parse(stream)), expected)


if __name__ == '__main__':
    unittest.main()
