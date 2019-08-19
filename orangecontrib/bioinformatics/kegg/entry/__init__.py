"""
DBGET entry
"""
from __future__ import absolute_import

import warnings
from collections import defaultdict

from . import fields
from .parser import DBGETEntryParser

__all__ = ["parser", "fields"]

# TODO: Remove the use of entry_decorate decorator
# for constructing a DBEntry subclass, make fields
# properties with __get__ method, and explicit assignment
# and meaningful docstrings


def entry_decorate(cls):
    """
    Decorate the DBEntry subclass with properties for accessing
    the fields through the 'DBField._convert' interface

    """
    reserved_names_map = {"class": "class_", "def": "def_"}

    def construct_one(name):
        def get(self):
            field = getattr(self, name, None)
            if field is not None:
                return field._convert()
            else:
                return None

        return property(get, doc=name)

    def construct_multiple(name):
        def get(self):
            field = getattr(self, name, None)
            if field is not None:
                return [f._convert() for f in field]
            else:
                return None

        return property(get, doc=name)

    for name, field in cls.FIELDS:
        name_lower = name.lower()
        if not hasattr(cls, name_lower):
            if name in cls.MULTIPLE_FIELDS:
                prop = construct_multiple(name)
            else:
                prop = construct_one(name)
            setattr(cls, reserved_names_map.get(name_lower, name_lower), prop)

    return cls


class DBEntry(object):
    """
    A DBGET entry object.
    """

    FIELDS = [("ENTRY", fields.DBEntryField)]
    MULTIPLE_FIELDS = []

    def __init__(self, text=None):
        self._sections = {}
        self.fields = []
        if text is not None:
            self.parse(text)

    @property
    def entry_key(self):
        """
        Primary entry key used for identifying the entry.
        """
        return self.entry.split(" ", 1)[0]

    def parse(self, text):
        """
        Parse `text` string containing a formated DBGET entry.
        """
        parser = DBGETEntryParser()
        gen = parser.parse_string(text)
        field_constructors = dict(self.FIELDS)

        current = None
        current_subfield = None
        entry_fields = []
        for (event, title, text) in gen:
            if event == DBGETEntryParser.SECTION_START:
                if title in field_constructors:
                    ftype = field_constructors[title]
                else:
                    ftype = fields.DBSimpleField
                current = ftype(text)
                if current.TITLE is None:
                    current.TITLE = title
            elif event == DBGETEntryParser.SECTION_END:
                entry_fields.append(current)
                current = None
            elif event == DBGETEntryParser.SUBSECTION_START:
                current_subfield = fields.DBSimpleField(text)
                current_subfield.TITLE = title
                if not isinstance(current, fields.DBFieldWithSubsections):
                    # Upgrade simple fields to FieldWithSubsection
                    new = fields.DBFieldWithSubsections(current.text)
                    new.TITLE = current.TITLE
                    current = new

            elif event == DBGETEntryParser.SUBSECTION_END:
                current.subsections.append(current_subfield)
                current_subfield = None
            elif event == DBGETEntryParser.TEXT:
                if current_subfield is not None:
                    current_subfield.text += text
                elif current is not None:
                    current.text += text
            elif event == DBGETEntryParser.ENTRY_END:
                break

        self.fields = entry_fields
        self._consolidate()

    def _consolidate(self):
        """
        Update mapping to field entries.
        """
        registered_fields = dict(self.FIELDS)
        multiple_fields = set(self.MULTIPLE_FIELDS)

        for field in self.fields:
            title = field.TITLE
            if title not in registered_fields:
                import warnings

                warnings.warn("Nonregisterd field %r in %r" % (title, type(self)))

            if title in multiple_fields:
                if not hasattr(self, title):
                    setattr(self, title, [])
                getattr(self, title).append(field)
            else:
                setattr(self, title, field)

    def __str__(self):
        return self.format()

    def format(self, section_indent=12):
        """
        Return a DBGET formated string representation.
        """
        return "".join(f.format(section_indent) for f in self.fields)

    def get(self, key, default=None):
        raise NotImplementedError

        f = getattr(self, key, None)
        if f is not None:
            if key in self.MULTIPLE_FIELDS:
                return [f.text for f in f]
            else:
                return f.text
        else:
            return None
