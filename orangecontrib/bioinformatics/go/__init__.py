"""  Gene Ontology module """
import os
import re
import sys
import tarfile
import warnings
from collections import namedtuple, defaultdict

import six

from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.utils import statistics, serverfiles, progress_bar_milestones
from orangecontrib.bioinformatics.go.config import DOMAIN, FILENAME_ONTOLOGY, FILENAME_ANNOTATION

intern = sys.intern
default_database_path = os.path.join(serverfiles.localpath(), DOMAIN)

_CVS_REVISION_RE = re.compile(r"^(rev)?(\d+\.\d+)+$")

evidence_types = {
    # Experimental
    'EXP': 'Inferred from Experiment',
    'IDA': 'Inferred from Direct Assay',
    'IPI': 'Inferred from Physical Interaction',  # [with <database:protein_name>]',
    'IMP': 'Inferred from Mutant Phenotype',
    'IGI': 'Inferred from Genetic Interaction',  # [with <database:gene_symbol[allele_symbol]>]',
    'IEP': 'Inferred from Expression Pattern',
    # Computational Analysis Evidence Codes
    'ISS': 'Inferred from Sequence Similarity',  # [with <database:sequence_id>] ',
    'ISA': 'Inferred from Sequence Alignment',
    'ISO': 'Inferred from Sequence Orthology',
    'ISM': 'Inferred from Sequence Model',
    'IGC': 'Inferred from Genomic Context',
    'RCA': 'Inferred from Reviewed Computational Analysis',
    # Author Statement Evidence Codes
    'TAS': 'Traceable author statement',
    'NAS': 'Non-traceable author statement',
    # Curatorial Statement Evidence Codes
    'IC': 'Inferred by curator',
    'ND': 'No biological data available',
    # Computationally-assigned Evidence Codes
    'IEA': 'Inferred from electronic annotation',  # [to <database:id>]',
    # Obsolete Evidence Codes
    'NR': 'Not Recorded(Obsolete)',
}

evidence_dict = defaultdict(int, [(e, 2 ** i) for i, e in enumerate(evidence_types.keys())])

evidence_types_ordered = [
    'EXP',
    'IDA',
    'IPI',
    'IMP',
    'IGI',
    'IEP',
    # Computational Analysis Evidence Codes
    'ISS',
    'ISA',
    'ISO',
    'ISM',
    'IGC',
    'RCA',
    # Author Statement Evidence Codes
    'TAS',
    'NAS',
    # Curatorial Statement Evidence Codes
    'IC',
    'ND',
    # Computationally-assigned Evidence Codes
    'IEA',
    # Obsolete Evidence Codes
    'NR',
]

multiplicity_set = {
    "alt_id",
    "is_a",
    "subset",
    "synonym",
    "related_synonym",
    "exact_synonym",
    "broad_synonym",
    "narrow_synonym",
    "xref_analog",
    "xref_unknown",
    "relationship",
}

multiple_tag_set = multiplicity_set

builtin_obo_objects = [
    """
[Typedef]
id: is_a
name: is_a
range: OBO:TERM_OR_TYPE
domain: OBO:TERM_OR_TYPE
definition: The basic subclassing relationship [OBO:defs]""",
    """[Typedef]
id: disjoint_from
name: disjoint_from
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that two classes are disjoint [OBO:defs]""",
    """[Typedef]
id: instance_of
name: instance_of
range: OBO:TERM
domain: OBO:INSTANCE
definition: Indicates the type of an instance [OBO:defs]""",
    """[Typedef]
id: inverse_of
name: inverse_of
range: OBO:TYPE
domain: OBO:TYPE
definition: Indicates that one relationship type is the inverse of another [OBO:defs]""",
    """[Typedef]
id: union_of
name: union_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the union of several others [OBO:defs]""",
    """[Typedef]
id: intersection_of
name: intersection_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the intersection of several others [OBO:defs]""",
]


class OBOObject:
    """ Represents a generic OBO object (e.g. Term, Typedef, Instance, ...)
    """

    _INTERN_TAGS = ["id", "name", "namespace", "alt_id", "is_a"]

    def __init__(self, stanza=None, ontology=None):
        self.ontology = ontology
        self._lines = []
        self.values = {}
        self.related = set()
        self.related_to = set()
        if stanza:
            self.parse_stanza(stanza)

    def parse_stanza(self, stanza):
        intern_tags = set(self._INTERN_TAGS)
        for line in stanza.splitlines():
            if ":" not in line:
                continue
            tag, rest = line.split(":", 1)
            value, modifiers, comment = "", "", ""
            if "!" in rest:
                rest, comment = rest.split("!")
            if "{" in rest:
                value, modifiers = rest.split("{", 1)
                modifiers = modifiers.strip("}")
            else:
                value = rest
            tag = intern(tag)
            value = value.strip()
            comment = comment.strip()
            if tag in intern_tags:
                value, comment = intern(value), intern(comment)
            self._lines.append((tag, value, modifiers, comment))
            if tag in multiple_tag_set:
                self.values.setdefault(tag, []).append(value)
            else:
                self.values[tag] = value
        self.related = set(self.related_objects())
        self.__dict__.update(self.values)
        if "def" in self.__dict__:
            self.__dict__["def_"] = self.__dict__["def"]

    def related_objects(self):
        """Return a list of tuple pairs where the first element is relationship
        type_id and the second id of object to whom the relationship applies to.

        """
        # TODO: add other defined Typedef ids
        type_ids = [intern("is_a")]
        result = [(type_id, id) for type_id in type_ids for id in self.values.get(type_id, [])]
        result = result + [tuple(map(intern, r.split(None, 1))) for r in self.values.get("relationship", [])]
        return result

    def __repr__(self):
        """ Return a string representation of the object in OBO format
        """
        _repr = "[%s]\n" % type(self).__name__
        for tag, value, modifiers, comment in self._lines:
            _repr = _repr + tag + ": " + value
            if modifiers:
                _repr = _repr + "{ " + modifiers + " }"
            if comment:
                _repr = _repr + " ! " + comment
            _repr = _repr + "\n"
        return _repr

    def __str__(self):
        """ Return the OBO object id entry
        """
        return "%s: %s" % (self.id, self.name)

    def __iter__(self):
        """ Iterates over sub terms
        """
        for type_id, id in self.related_to:
            yield (type_id, self.ontology[id])


class Term(OBOObject):
    pass


class Typedef(OBOObject):
    pass


class Instance(OBOObject):
    pass


class Ontology:
    """
    :class:`Ontology` is the class representing a gene ontology.

    :param str filename:
        A filename of an .obo formated file.
    :param progress_callback:
        Optional `float -> None` function.


    Example
    --------
        >>> # Load the current ontology (downloading it if necessary)
        >>> ontology = Ontology()

        >>> term_ids = list(ontology)
        >>> term = ontology[term_ids[0]]

    """

    version = 1

    def __init__(self, filename=None, progress_callback=None):
        self.terms = {}
        self.typedefs = {}
        self.instances = {}
        self.slims_subset = set()
        self.alias_mapper = {}
        self.reverse_alias_mapper = defaultdict(set)
        self.header = ""

        if filename is not None:
            self.parse_file(filename, progress_callback)
        else:
            filename = serverfiles.localpath_download(DOMAIN, FILENAME_ONTOLOGY)
            self.parse_file(filename, progress_callback)

    @classmethod
    def load(cls, progress_callback=None):
        """
        A class method that tries to load the ontology file from
        default_database_path. It looks for a filename starting with
        'gene_ontology'. If not found it will download it.

        """
        filename = serverfiles.localpath_download(DOMAIN, FILENAME_ONTOLOGY)

        return cls(filename, progress_callback=progress_callback)

    Load = load

    def parse_file(self, file, progress_callback=None):
        """
        Parse the file. file can be a filename string or an open file like
        object. The optional progressCallback will be called with a single
        argument to report on the progress.
        """
        if isinstance(file, str):
            if os.path.isfile(file) and tarfile.is_tarfile(file):
                f = tarfile.open(file).extractfile("gene_ontology_edit.obo")
            elif os.path.isfile(file):
                f = open(file)
            elif os.path.isdir(file):
                f = open(os.path.join(file, "gene_ontology_edit.obo"))
            else:
                raise ValueError("Cannot open %r for parsing" % file)
        else:
            f = file

        data = [line.decode() if not isinstance(line, str) else line for line in f.readlines()]
        data = "".join([line for line in data if not line.startswith("!")])
        self.header = data[: data.index("[Term]")]
        c = re.compile(r"\[.+?\].*?\n\n", re.DOTALL)
        data = c.findall(data)

        milestones = progress_bar_milestones(len(data), 90)
        for i, block in enumerate(builtin_obo_objects + data):
            if block.startswith("[Term]"):
                term = Term(block, self)
                self.terms[term.id] = term
            elif block.startswith("[Typedef]"):
                typedef = Typedef(block, self)
                self.typedefs[typedef.id] = typedef
            elif block.startswith("[Instance]"):
                instance = Instance(block, self)
                self.instances[instance.id] = instance
            if progress_callback and i in milestones:
                progress_callback(90.0 * i / len(data))

        self.alias_mapper = {}
        self.reverse_alias_mapper = defaultdict(set)
        milestones = progress_bar_milestones(len(self.terms), 10)
        for i, (id, term) in enumerate(six.iteritems(self.terms)):
            for type_id, parent in term.related:
                self.terms[parent].related_to.add((type_id, id))
            try:
                self.alias_mapper.update([(alt_id, id) for alt_id in term.alt_id])
                self.reverse_alias_mapper[id].update(term.alt_id)
            except AttributeError:
                pass
            if progress_callback and i in milestones:
                progress_callback(90.0 + 10.0 * i / len(self.terms))

    def defined_slims_subsets(self):
        """
        Return a list of defined subsets in the ontology.

        :rtype: :class:`list` of :class:`str`

        """
        return [line.split()[1] for line in self.header.splitlines() if line.startswith("subsetdef:")]

    def named_slims_subset(self, subset):
        """
        Return all term IDs in a named `subset`.

        :param str subset: A string naming a subset in the ontology.
        :rtype: :class:`list` of :class:`str`

        .. seealso:: :func:`defined_slims_subsets`

        """
        return [id for id, term in self.terms.items() if subset in getattr(term, "subset", set())]

    def set_slims_subset(self, subset):
        """
        Set the `slims_subset` term subset to `subset`.

        :param set subset: A subset of GO term IDs.

        `subset` may also be a string, in which case the call is equivalent
        to ``ont.set_slims_subsets(ont.named_slims_subset(subset))``

        """
        if isinstance(subset, str):
            self.slims_subset = set(self.named_slims_subset(subset))
        else:
            self.slims_subset = set(subset)

    def slims_for_term(self, term):
        """
        Return a list of slim term IDs for `term`.

        This is a list of `most specific` slim terms to which `term` belongs.

        :param str term: Term ID.

        """
        queue = {term}
        visited = set()
        slims = set()
        while queue:
            term = queue.pop()
            visited.add(term)
            if term in self.slims_subset:
                slims.add(term)
            else:
                queue.update({tid for _, tid in self[term].related} - visited)
        return slims

    def extract_super_graph(self, terms):
        """
        Return all super terms of `terms` up to the most general one.

        :param list terms: A list of term IDs.

        """
        terms = [terms] if isinstance(terms, str) else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update({tid for _, tid in self[term].related} - visited)
        return visited

    def extract_sub_graph(self, terms):
        """
        Return all sub terms of `terms`.

        :param list terms: A list of term IDs.

        """
        terms = [terms] if type(terms) == str else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update({tid for _, tid in self[term].related_to} - visited)
        return visited

    def term_depth(self, term, cache_={}):
        """
        Return the minimum depth of a `term`.

        (length of the shortest path to this term from the top level term).

        """
        if term not in cache_:
            cache_[term] = min([self.term_depth(parent) + 1 for _, parent in self[term].related] or [1])
        return cache_[term]

    def __getitem__(self, termid):
        """
        Return a :class:`Term` object with `termid`.

        :param str term: An id of a 'Term' in the ontology.
        :rtype: :class:`Term`

        """
        if termid in self.terms:
            return self.terms[termid]
        elif termid in self.alias_mapper:
            return self.terms[self.alias_mapper[termid]]
        else:
            raise KeyError(termid)

    def __iter__(self):
        """
        Iterate over all term ids in ontology.
        """
        return iter(self.terms)

    def __len__(self):
        """
        Return number of terms in ontology.
        """
        return len(self.terms)

    def __contains__(self, termid):
        """
        Return `True` if a term with `termid` is present in the ontology.
        """
        return termid in self.terms or termid in self.alias_mapper


annotation_fields = ["tax_id", "gene_id", "go_id", "evidence", "qualifier", "go_term", "pubMed", "aspect"]


_AnnotationRecordBase = namedtuple("AnnotationRecord", annotation_fields)


class AnnotationRecord(_AnnotationRecordBase):
    """ An annotation record mapping a gene to a term.

    See ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/README for description
    if individual fields under <gene2go> section.

    """

    def __new__(cls, *args):
        if len(args) == 1 and isinstance(args[0], str):
            args = map(intern, args[0].split("\t"))
        return super(AnnotationRecord, cls).__new__(cls, *args)

    @classmethod
    def from_string(cls, string):
        """ Create an instance from a line in a annotations file format from serverfiles.
        """
        return AnnotationRecord._make(map(intern, string.strip().split('\t')))


class Annotations:
    """ :class:`Annotations` object holds the annotations.

    :param str organism:
        an organism specifier (e.g. ``'9606'``). Annotations for that organism will be loaded.

    :param ontology: :class:`Ontology` object for annotations
    :type ontology: :class:`Ontology`

    """

    def __init__(self, organism, ontology=None, progress_callback=None, filename=None):
        #: A dictionary mapping a gene (gene_id) to a set of all annotations of that gene.
        self.gene_annotations = defaultdict(list)

        #: A dictionary mapping a GO term id to a set of annotations that are directly annotated to that term
        self.term_anotations = defaultdict(list)

        self.all_annotations = defaultdict(list)

        self._gene_names = None
        self._gene_names_dict = None

        #: A list of all :class:`AnnotationRecords` instances.
        self.annotations = []
        self.header = ''
        self.taxid = organism

        self._ontology = ontology

        if filename is None:
            try:
                filename = serverfiles.localpath_download(
                    DOMAIN, FILENAME_ANNOTATION.format(organism), progress_callback=progress_callback
                )
            except FileNotFoundError:
                raise taxonomy.UnknownSpeciesIdentifier(organism)

        self._parse_file(filename)

    @property
    def ontology(self):
        return self._ontology

    @ontology.setter
    def ontology(self, ontology):
        """ Set the ontology to use in the annotations mapping.
        """
        self.all_annotations = defaultdict(list)
        self._ontology = ontology

    def _ensure_ontology(self):
        if self.ontology is None:
            self.ontology = Ontology()

    def _parse_file(self, file_path):

        with open(file_path, 'r') as anno_file:
            self.header = anno_file.readline()

            for line in anno_file.readlines():
                self.add_annotation(AnnotationRecord.from_string(line))

    def add_annotation(self, a):
        """ Add a single :class:`AnotationRecord` instance to this object.
        """
        if not isinstance(a, AnnotationRecord):
            a = AnnotationRecord(a)
        if not a.gene_id or not a.go_id or a.qualifier == 'NOT':
            return

        self.gene_annotations[a.gene_id].append(a)
        self.term_anotations[a.go_id].append(a)

        self.annotations.append(a)
        self.all_annotations = defaultdict(list)

    def get_genes_with_known_annotation(self, genes):
        """ Return only genes with known annotation

        :param genes: List of genes

        """
        return {gene for gene in genes if self.gene_annotations[gene]}

    def _collect_annotations(self, go_id, visited):
        """ Recursive function collects and caches all annotations for id
        """
        if go_id not in self.all_annotations and go_id not in visited:
            if go_id in self.ontology.reverse_alias_mapper:
                annotations = [
                    self.term_anotations.get(alt_id, []) for alt_id in self.ontology.reverse_alias_mapper[go_id]
                ] + [self.term_anotations[go_id]]
            else:
                annotations = [self.term_anotations[go_id]]  # annotations for this term alone
            visited.add(go_id)

            for type_id, child in self.ontology[go_id].related_to:
                aa = self._collect_annotations(child, visited)
                if type(aa) == set:
                    annotations.append(aa)  # if it was already reduced in get_all_annotations
                else:
                    annotations.extend(aa)
            self.all_annotations[go_id] = annotations
        return self.all_annotations[go_id]

    def get_annotations_by_go_id(self, go_id):
        """ Return a set of all annotations (instances of :obj:`AnnotationRecord`)
        for GO term `id` and all it's subterms.

        :param :class:`str` go_id: GO term id

        """
        self._ensure_ontology()
        id = self.ontology.alias_mapper.get(go_id, go_id)
        if id not in self.all_annotations or type(self.all_annotations[id]) == list:
            annot_set = set()
            for annots in self._collect_annotations(id, set()):
                annot_set.update(annots)
            self.all_annotations[id] = annot_set
        return self.all_annotations[id]

    def get_genes_by_go_term(self, go_id, evidence_codes=None):
        """ Return a list of genes annotated by specified `evidence_codes`
        to GO term 'id' and all it's subterms."

        :param str go_id: GO term id

        :param list-of-strings evidence_codes:
               List of evidence codes to consider when matching annotations to terms.

        """
        evidence_codes = set(evidence_codes or evidence_dict.keys())
        annotations = self.get_annotations_by_go_id(go_id)
        return list({int(ann.gene_id) for ann in annotations if ann.evidence in evidence_codes})

    def genes(self):
        return {ann.gene_id for ann in self.annotations}

    def get_enriched_terms(
        self,
        genes,
        reference=None,
        evidence_codes=None,
        slims_only=False,
        aspect=None,
        prob=statistics.Binomial(),
        use_fdr=True,
        progress_callback=None,
    ):
        """
        Return a dictionary of enriched terms, with tuples of
        (list_of_genes, p_value, reference_count) for items and term
        ids as keys. P-Values are FDR adjusted if use_fdr is True (default).

        :param genes: List of genes
        :param reference: List of genes (if None all genes included in the annotations will be used).
        :param evidence_codes:  List of evidence codes to consider.
        :param slims_only: If `True` return only slim terms.
        :param aspect: Which aspects to use. Use all by default;
                       one of Process (biological process),
                       Function (molecular function) or Component (cellular component)
        :param prob:
        :param use_fdr:
        :param progress_callback:
        """

        all_genes = set(genes)

        if aspect is None:
            aspects_set = {'Process', 'Component', 'Function'}
        elif isinstance(aspect, str):
            aspects_set = {aspect}
        else:
            aspects_set = aspect

        if reference is None:
            reference = self.genes()

        evidence_codes = set(evidence_codes or evidence_dict.keys())
        annotations = [
            ann
            for gene in genes
            for ann in self.gene_annotations[gene]
            if ann.evidence in evidence_codes and ann.aspect in aspects_set
        ]

        ref_annotations = {
            ann
            for gene in reference
            for ann in self.gene_annotations[gene]
            if ann.evidence in evidence_codes and ann.aspect in aspects_set
        }

        annotations_dict = defaultdict(set)
        for ann in annotations:
            annotations_dict[ann.go_id].add(ann)

        self._ensure_ontology()

        if slims_only and not self.ontology.slims_subset:
            warnings.warn("Unspecified slims subset in the ontology! " "Using 'goslim_generic' subset", UserWarning)
            self.ontology.set_slims_subset('goslim_generic')

        terms = annotations_dict.keys()
        filtered_terms = [term for term in terms if term in self.ontology]

        if len(terms) != len(filtered_terms):
            term_diff = set(terms) - set(filtered_terms)
            warnings.warn(
                "%s terms in the annotations were not found in the " "ontology." % ",".join(map(repr, term_diff)),
                UserWarning,
            )

        terms = self.ontology.extract_super_graph(filtered_terms)
        res = {}

        milestones = progress_bar_milestones(len(terms), 100)

        for i, term in enumerate(terms):
            if slims_only and term not in self.ontology.slims_subset:
                continue
            all_annotations = self.get_annotations_by_go_id(term).intersection(ref_annotations)
            all_annotated_genes = {ann.gene_id for ann in all_annotations}
            mapped_genes = all_genes.intersection(all_annotated_genes)

            if len(reference) > len(all_annotated_genes):
                mapped_reference_genes = reference.intersection(all_annotated_genes)
            else:
                mapped_reference_genes = all_annotated_genes.intersection(reference)

            res[term] = (
                [gene for gene in mapped_genes],
                prob.p_value(len(mapped_genes), len(reference), len(mapped_reference_genes), len(genes)),
                len(mapped_reference_genes),
            )

            if progress_callback and i in milestones:
                progress_callback(100.0 * i / len(terms))

        if use_fdr:
            res = sorted(res.items(), key=lambda x: x[1][1])
            res = {
                id: (genes, p, ref)
                for (id, (genes, _, ref)), p in zip(res, statistics.FDR([p for _, (_, p, _) in res]))
            }
        return res

    def get_annotated_terms(self, genes, direct_annotation_only=False, evidence_codes=None, progress_callback=None):
        """ Return all terms that are annotated by genes with evidence_codes.
        """

        genes = [genes] if type(genes) == str else genes
        genes = {gene for gene in genes}

        evidence_codes = set(evidence_codes or evidence_dict.keys())
        annotations = [ann for gene in genes for ann in self.gene_annotations[gene] if ann.evidence in evidence_codes]

        dd = defaultdict(set)
        for ann in annotations:
            dd[ann.go_id].add(ann.gene_id)

        if not direct_annotation_only:
            self._ensure_ontology()
            terms = dd.keys()
            filtered_terms = [term for term in terms if term in self.ontology]
            if len(terms) != len(filtered_terms):
                term_diff = set(terms) - set(filtered_terms)
                warnings.warn(
                    "%s terms in the annotations were not found in the " "ontology." % ",".join(map(repr, term_diff)),
                    UserWarning,
                )

            terms = self.ontology.extract_super_graph(filtered_terms)
            for i, term in enumerate(terms):
                term_annotations = self.get_annotations_by_go_id(term).intersection(annotations)
                dd[term].update([ann.gene_id for ann in term_annotations])
        return dict(dd)

    def __add__(self, iterable):
        """ Return a new Annotations object with combined annotations
        """
        return Annotations([a for a in self] + [a for a in iterable], ontology=self.ontology)

    def __iadd__(self, iterable):
        """ Add annotations to this instance
        """
        self.extend(iterable)
        return self

    def __contains__(self, item):
        return item in self.annotations

    def __iter__(self):
        """ Iterate over all AnnotationRecord objects in annotations
        """
        return iter(self.annotations)

    def __len__(self):
        """ Return the number of annotations
        """
        return len(self.annotations)

    def __getitem__(self, index):
        """ Return the i-th annotation record
        """
        return self.annotations[index]

    def __getslice__(self, *args):
        return self.annotations.__getslice__(*args)

    def add(self, line):
        """ Add one annotation
        """
        self.add_annotation(line)

    def append(self, line):
        """ Add one annotation
        """
        self.add_annotation(line)

    def extend(self, lines):
        """ Add multiple annotations
        """
        for line in lines:
            self.add_annotation(line)


def filter_by_p_value(terms, p_value=0.01):
    """ Filters the terms by the p-value. Assumes terms is a dict with
    the same structure as returned from get_enriched_terms.
    """
    return dict(filter(lambda x: x[1][1] <= p_value, terms.items()))


if __name__ == "__main__":
    pass
