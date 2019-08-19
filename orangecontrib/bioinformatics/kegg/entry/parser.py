"""
A parser for DBGET database entries

"""
from __future__ import print_function

from six import StringIO


class DBGETEntryParser(object):
    r"""
    A DBGET entry parser (inspired by ``xml.dom.pulldom``).

    Example
    -------
    >>> stream = StringIO(
    ...     "ENTRY foo\n"
    ...     "NAME  foo's name\n"
    ...     "  BAR A subsection of 'NAME'\n"
    ... )
    >>> parser = DBGETEntryParser()
    >>> for event, title, contents_part in parser.parse(stream):
    ...    print(parser.EVENTS[event], title, repr(contents_part))
    ...
    ENTRY_START None None
    SECTION_START ENTRY 'foo\n'
    SECTION_END ENTRY None
    SECTION_START NAME "foo's name\n"
    SUBSECTION_START BAR "A subsection of 'NAME'\n"
    SUBSECTION_END BAR None
    SECTION_END NAME None
    ENTRY_END None None
    """
    #: Entry start events
    ENTRY_START = 0

    #: Entry end event
    ENTRY_END = 1

    #: Section start event
    SECTION_START = 2

    #: Section end event
    SECTION_END = 3

    #: Subsection start event
    SUBSECTION_START = 4

    #: Subsection end event
    SUBSECTION_END = 5

    #: Text element event
    TEXT = 6

    EVENTS = ["ENTRY_START", "ENTRY_END", 'SECTION_START', 'SECTION_END', 'SUBSECTION_START', 'SUBSECTION_END', 'TEXT']

    def __init__(self):
        pass

    def parse(self, stream):
        entry_offset = None
        section_title = None
        subsection_title = None
        textline_start = None

        for line in stream:
            startswith = line.startswith
            # TODO: Reorder by frequency (for faster fallthrough)
            if startswith("ENTRY"):
                # Start parsing new entry
                yield (self.ENTRY_START, None, None)
                title, rest = self._partition_section_title(line)
                entry_offset = len(line) - len(rest)
                textline_start = " " * entry_offset
                yield (self.SECTION_START, title, rest)
                yield (self.SECTION_END, title, None)

            elif startswith("///"):
                # End entry
                if subsection_title is not None:
                    # End current subsection if any
                    yield (self.SUBSECTION_END, subsection_title, None)
                    subsection_title = None

                if section_title is not None:
                    # End current section if any
                    yield (self.SECTION_END, section_title, None)
                    section_title = None

                yield (self.ENTRY_END, None, None)
                entry_offset = None
                textline_start = None

            elif not startswith(" "):
                # Start new section
                if subsection_title is not None:
                    # End current subsection if any
                    yield (self.SUBSECTION_END, subsection_title, None)
                    subsection_title = None

                if section_title is not None:
                    # End current section if any
                    yield (self.SECTION_END, section_title, None)
                    section_title = None

                title, rest = self._partition_section_title(line)
                section_title = title
                yield (self.SECTION_START, section_title, rest)

            elif startswith(textline_start):
                # A line of text
                # TODO: pass the current subsection/section title
                yield (self.TEXT, None, line[entry_offset:])

            elif startswith(" "):
                # Start a new subsection
                if subsection_title is not None:
                    # End current subsection
                    yield (self.SUBSECTION_END, subsection_title, None)
                title, rest = self._partition_subsection_title(line)
                subsection_title = title
                yield (self.SUBSECTION_START, subsection_title, rest)

        # Close any remaining sections/entries
        if subsection_title is not None:
            yield (self.SUBSECTION_END, subsection_title, None)
        if section_title is not None:
            yield (self.SECTION_END, section_title, None)
        if entry_offset is not None:
            yield (self.ENTRY_END, None, None)

    def parse_string(self, string):
        return self.parse(StringIO(string))

    def _partition_section_title(self, line):
        """
        Split the section title from the rest of the line
        """
        try:
            title, rest = line.split(" ", 1)
        except ValueError:
            # no contents, only section title
            title = line.rstrip()
            rest = ""
        rest = rest.lstrip(" ")
        return title, rest

    def _partition_subsection_title(self, line):
        """
        Split the subsection title from the rest of the line
        """
        line = line.lstrip(" ")
        return self._partition_section_title(line)
