""" Orange data table and other data manipulation utilities """
import numbers
from types import SimpleNamespace
from typing import Tuple, Sequence
from functools import wraps

from Orange.data import Table, Variable
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.utils.messages import UnboundMsg

__missing_annotation = UnboundMsg('Missing annotation on gene IDs and organism in the input data.')
__missing_gene_id = UnboundMsg('Missing gene ID information. Make sure that Table is properly annotated.')
__missing_tax_id = UnboundMsg('Missing organism information. Make sure that Table is properly annotated.')
__unable_to_locate_genes = UnboundMsg('Unable to locate genes. Make sure that Table is properly annotated.')


class TableAnnotation(SimpleNamespace):
    """ Data Table hints """

    # Organism in data table
    tax_id: str = 'taxonomy_id'

    # This indicates position of genes in data table
    gene_as_attr_name: str = 'gene_as_attribute_name'

    # This indicates a column name (if genes are in rows)
    gene_id_column: str = 'gene_id_column'

    # This indicates attribute name (if genes are in columns)
    gene_id_attribute: str = 'gene_id_attribute'


def check_table_annotation(f):
    """Wrapper for widget's input method that checks if the data on the input is correctly annotated.

    A widget in bioinformatics add-on expects that every Table has additional
    information stored as table attributes:
       - taxonomy_id = 'taxonomy id for given organism'
       - gene_as_attribute_name = 'location of gene names (rows/columns)'
       - gene_id_attribute/gene_id_column = 'attribute/column name'

    """

    @wraps(f)
    def wrapper(widget, data: Table, *args, **kwargs):
        widget.Error.add_message('missing_annotation', __missing_annotation)
        widget.Error.add_message('missing_gene_id', __missing_gene_id)
        widget.Error.add_message('missing_tax_id', __missing_tax_id)
        widget.Error.add_message('unable_to_locate_genes', __unable_to_locate_genes)

        widget.Error.missing_annotation.clear()
        widget.Error.missing_gene_id.clear()
        widget.Error.missing_tax_id.clear()
        widget.Error.unable_to_locate_genes.clear()

        if data is not None and isinstance(data, Table):
            attributes: dict = data.attributes

            tax_id: str = TableAnnotation.tax_id
            gene_id_column: str = TableAnnotation.gene_id_column
            gene_id_attribute: str = TableAnnotation.gene_id_attribute
            gene_as_attr_name: str = TableAnnotation.gene_as_attr_name

            if not attributes:
                widget.Error.missing_annotation()
                data = None

            elif tax_id not in attributes:
                widget.Error.missing_tax_id()
                data = None

            elif gene_as_attr_name not in attributes:
                widget.Error.unable_to_locate_genes()
                data = None

            elif gene_as_attr_name in attributes:
                if (attributes[gene_as_attr_name] and gene_id_attribute not in attributes) or (
                    not attributes[gene_as_attr_name] and gene_id_column not in attributes
                ):

                    widget.Error.unable_to_locate_genes()
                    data = None

        return f(widget, data, *args, **kwargs)

    return wrapper


# TODO: remove this and replace with TableAnnotation namespace

# species
TAX_ID = 'taxonomy_id'

# Will be set to True if gene names are represented as attribute names.
# If gene names are in rows, we set this value to False. (user must select proper column index)
GENE_AS_ATTRIBUTE_NAME = 'gene_as_attribute_name'

# Name of the column where rows are gene ids
GENE_ID_COLUMN = 'gene_id_column'

# Name of the variable attribute that holds gene id
GENE_ID_ATTRIBUTE = 'gene_id_attribute'

# Error strings
ERROR_ON_MISSING_ANNOTATION = 'Missing annotation on gene IDs and organism in the input data.'
ERROR_ON_MISSING_GENE_ID = 'Missing gene ID information in the input data'
ERROR_ON_MISSING_TAX_ID = 'Missing organism information in the input data'
