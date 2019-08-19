"""
Wrapper classes for db entry fields to support pythonic
interface.

"""


class DBField(object):
    """
    Base DBGET entry field
    """

    __SLOTS__ = ["text"]

    def __init__(self, text):
        self.text = text

    def _convert(self):
        """
        Convert the contents into python representation using builtin types.
        """
        return self.text.rstrip("\n")


class DBSimpleField(DBField):
    """
    Simple field (with no subsections).
    """

    __SLOTS__ = ["text"]
    # TITLE must be set in subclasses or object instances
    TITLE = None

    def __str__(self):
        return self.format()

    def format(self, section_indent=12, subsection_indent=0):
        fmt = (" " * subsection_indent) + "%-" + str(section_indent - subsection_indent) + "s%s"
        text = self._indent(self.text, section_indent)
        text = fmt % (self.TITLE, text)
        return text

    def _indent(self, text, section_indent=12):
        indent_str = "\n" + " " * section_indent
        nl_count = text.count("\n")
        return text.replace("\n", indent_str, nl_count - 1)


class DBEntryField(DBSimpleField):
    """
    ENTRY field (all entries start with this field)
    """

    __SLOTS__ = ["text"]
    TITLE = "ENTRY"


class DBNameField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "NAME"


class DBDefinitionField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "DEFINITION"


class DBFieldWithSubsections(DBSimpleField):
    """
    A field with subsections (for instance REFERENCE in genome)
    """

    __SLOTS__ = ["text", "subsections"]
    TITLE = None
    SUBSECTIONS = None

    def __init__(self, text, subsections=None):
        self.text = text
        self.subsections = subsections or []

    def format(self, section_indent=12, subsection_indent=2):
        text = DBSimpleField.format(self, section_indent, subsection_indent=0)
        subsections = [sub.format(section_indent, subsection_indent) for sub in self.subsections]
        return "".join([text] + subsections)

    def _convert(self):
        my = DBSimpleField._convert(self)
        subs = [(s.TITLE.lower(), s._convert()) for s in self.subsections]
        return (my, subs)


class DBTaxonomyField(DBFieldWithSubsections):
    __SLOTS__ = ["text", "subsections"]
    TITLE = "TAXONOMY"
    SUBSECTIONS = ["LINEAGE"]

    @property
    def taxid(self):
        return DBSimpleField._convert(self).split(":")[1]


class DBDataSourceField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "DATA_SOURCE"


class DBReference(DBFieldWithSubsections):
    __SLOTS__ = ["text", "subsections"]
    TITLE = "REFERENCE"
    SUBSECTIONS = ["AUTHORS", "TITLE", "JOURNAL"]

    @property
    def authors(self):
        return self.subsections[0]

    @property
    def title(self):
        return self.subsections[1]

    @property
    def journal(self):
        return self.subsections[2]


class DBDBLinks(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "DBLINKS"

    @property
    def links(self):
        return [tuple(s.split(": ", 1)) for s in self.text.splitlines()]

    def _convert(self):
        # Some dblinks can span multiple lines but are always 'indented'
        links = DBSimpleField._convert(self).replace("\n ", "").splitlines()
        links = [tuple(link.split(": ", 1)) for link in links]
        links = [(key, [v for v in values.split(" ") if v]) for key, values in links]
        return dict(links)


class DBPathway(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "PATHWAY"

    @property
    def pathways(self):
        return self._convert()

    def _convert(self):
        text = DBSimpleField._convert(self)
        return [line.split(" ", 1)[0] for line in text.splitlines()]


class DBAASeq(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "AASEQ"

    @property
    def sequence(self):
        return self.split("\n", 1)[1].replace("\n", "")

    @property
    def sequence_lenght(self):
        return int(self.text.split("\n", 1)[0])

    def _convert(self):
        text = DBSimpleField._convert(self)
        count, seq = text.split("\n", 1)
        return seq.replace("\n", "")


class DBNTSeq(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "NTSEQ"

    @property
    def sequence(self):
        return self.split("\n", 1)[1].replace("\n", "")

    @property
    def sequence_lenght(self):
        return int(self.text.split("\n", 1)[0])

    def _convert(self):
        text = DBSimpleField._convert(self)
        count, seq = text.split("\n", 1)
        return seq.replace("\n", "")


class DBPathwayMapField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "PATHWAY_MAP"

    def kgml_url(self):
        return "http://www.genome.jp/kegg-bin/download?entry={0}&format=kgml".format(self.pathway_id)

    @property
    def pathway_id(self):
        return self.text.split(" ", 1)[0]


class DBGeneField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "GENE"

    def _convert(self):
        text = DBSimpleField._convert(self)
        lines = text.splitlines()
        return [line.split(" ", 1)[0] for line in lines]


class DBEnzymeField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "ENZYME"

    def _convert(self):
        text = DBSimpleField._convert(self)
        lines = text.splitlines()
        return lines


class DBCompoundField(DBSimpleField):
    __SLOTS__ = ["text"]
    TITLE = "COMPOUND"

    def _convert(self):
        text = DBSimpleField._convert(self)
        lines = text.splitlines()
        return [line.split(" ", 1)[0] for line in lines]
