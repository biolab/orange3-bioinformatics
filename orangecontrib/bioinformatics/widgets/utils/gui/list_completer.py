""" List Completer """
import re

from AnyQt.QtCore import Qt, QObject, QStringListModel
from AnyQt.QtWidgets import QCompleter


class TokenListCompleter(QCompleter):
    """
    A completer allowing completion of multiple words (tokens)

    Example
    -------
    completer = TokenListCompleter()
    completer.setTokenList(["foo", "bar", "baz"])
    completer.setCompletionPrefix("foo b")
    for i in range(completer.completionCount()):
         completer.setCurrentRow(i)
         print(completer.currentCompletion())

    """

    def __init__(self, *args, **kwargs):
        separator = kwargs.pop("separator", " ")
        super().__init__(*args, **kwargs)
        self.__tokenList = []
        self.__completerModel = None
        self.__separator = separator
        # The current 'known' completion prefix (tracked in splitPath)
        self.__currentKnownPrefix = ""
        self.setModelSorting(QCompleter.CaseSensitivelySortedModel)
        super().setModel(QStringListModel(self))

    def setTokenList(self, tokenList):
        """
        Set the token word completion list

        Equivalent to
        `completer.setModel(QStringListModel(tokenList, completer))`

        Note that this is the preferred method for setting the completion
        model.

        Parameters
        ----------
        tokenList : List[str]
        """
        tokenList = list(sorted(set(tokenList)))
        self.setModel(QStringListModel(tokenList, self))

    def tokenList(self):
        """
        Return the current token word completion list.

        Return
        ------
        tokens : List[str]
        """
        return list(self.__tokenList)

    def setModel(self, model):
        """
        Reimplemented.

        Parameters
        ----------
        model : QAbstractItemModel
        """
        if model is self.__completerModel:
            return

        if self.__completerModel is not None:
            self.__completerModel.dataChanged.disconnect(self.__initDynamicModel)
            self.__completerModel.rowsInserted.disconnect(self.__initDynamicModel)
            self.__completerModel.rowsRemoved.disconnect(self.__initDynamicModel)

            if QObject.parent(self.__completerModel) is self:
                self.__completerModel.deleteLater()
            self.__completerModel = None

        self.__completerModel = model

        if self.__completerModel is not None:
            self.__completerModel.dataChanged.connect(self.__initDynamicModel)
            self.__completerModel.rowsInserted.connect(self.__initDynamicModel)
            self.__completerModel.rowsRemoved.connect(self.__initDynamicModel)

        self.__initDynamicModel()

    def model(self):
        return self.__completerModel

    def setSeparator(self, separator):
        if self.__separator != separator:
            self.__separator = separator

    def separator(self):
        return self.__separator

    # QCompleter::setCompleterPrefix is not a virtual method, and is called
    # directly by the QLineEdit to set the prefix. There does not appear to
    # be any other direct way that a QCompleter subclass could be notified
    # of the completion prefix change.
    # However in QCompleter::setCompleterPrefix
    # the prefix is immediately split by splitPath which is virtual.
    # So we use this to track changes in the current prefix.
    def splitPath(self, prefix):
        """reimplemented."""
        items = super().splitPath(prefix)
        if self.__currentKnownPrefix != self.completionPrefix():
            self.__currentKnownPrefix = self.completionPrefix()
            self.__updateEffectiveCompleterModel()
        return items

    def __updateEffectiveCompleterModel(self):
        """
        Prefix all items in the effective completer model with the current
        prefix to enable the completion of multiple keywords.
        """
        prefix = self.completionPrefix()
        model = super().model()
        assert isinstance(model, QStringListModel)
        separator = self.__separator
        if not prefix.rstrip().endswith(separator) and separator in prefix:
            # We are in the middle of a second (or later) word completion
            prefix, rest = prefix.rsplit(separator, 1)
            # preserve any whitespace after the seperator
            match = re.match(r"^(?P<whitespace>\s*)(?P<rest>.*)", rest)
            if match:
                whitespace, rest = match.groups()
            else:
                whitespace = ""
            items = [prefix + separator + whitespace + item for item in self.__tokenList]
        else:
            # Either completing the the first word or immediately after
            # the separator. In both cases the original tokens list is
            # the active completion list
            items = self.__tokenList
        itemset = set(items)
        current = model.stringList()

        if itemset != set(current):
            model.setStringList(list(sorted(itemset)))

    def __initDynamicModel(self):
        """
        [Re]Initialize the effective completion model from the client
        supplied model
        """
        model = self.__completerModel
        if model is not None:
            if isinstance(model, QStringListModel):
                tokens = model.stringList()
            else:
                tokens = [model.data(model.index(row, 0), Qt.DisplayRole) for row in range(model.rowCount())]
                tokens = [str(token) for token in filter(None, tokens)]
        else:
            tokens = []

        tokenList = list(sorted(set(tokens)))
        if tokenList != self.__tokenList:
            self.__tokenList = tokenList
            self.__updateEffectiveCompleterModel()
