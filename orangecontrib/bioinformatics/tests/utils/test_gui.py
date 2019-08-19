import unittest

from AnyQt.QtCore import QStringListModel

from orangecontrib.bioinformatics.widgets.utils.gui import TokenListCompleter


class TestCompleter(unittest.TestCase):
    def test_token_list_completer(self):
        from Orange.widgets.utils.itemmodels import PyListModel

        completer = TokenListCompleter()
        completer.setTokenList(["foo", "bar", "baz"])
        completer.setCompletionPrefix("foo b")

        def completions(completer):
            current = completer.currentRow()
            items = []
            for i in range(completer.completionCount()):
                completer.setCurrentRow(i)
                items.append(completer.currentCompletion())
            completer.setCurrentRow(current)
            return items

        self.assertSequenceEqual(completions(completer), ["foo bar", "foo baz"])
        completer.setModel(None)
        self.assertSequenceEqual(completer.tokenList(), [])
        self.assertSequenceEqual(completions(completer), [])

        completer.setModel(QStringListModel(["a", "ab", "b"]))
        self.assertSequenceEqual(completer.tokenList(), ["a", "ab", "b"])
        completer.setCompletionPrefix("a a")
        self.assertSequenceEqual(completions(completer), ["a a", "a ab"])

        completer.setModel(PyListModel(["a", "aa", "ab", "ad"]))
        self.assertSequenceEqual(completer.tokenList(), ["a", "aa", "ab", "ad"])

        completer.model()[-1] = "z"
        self.assertSequenceEqual(completer.tokenList(), ["a", "aa", "ab", "z"])

        completer.model()[-2:] = ["ax", "az"]
        self.assertSequenceEqual(completer.tokenList(), ["a", "aa", "ax", "az"])

        completer.setSeparator(",")
        completer.setCompletionPrefix("a, a")
        self.assertSequenceEqual(completions(completer), ["a, a", "a, aa", "a, ax", "a, az"])
