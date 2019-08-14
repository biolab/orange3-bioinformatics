""" Orange data table and other data manipulation utilities """
import numbers
from types import SimpleNamespace
from typing import Tuple, Sequence

import Orange.data


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

col_spec = Sequence[Tuple[Orange.data.Variable, Sequence[numbers.Real]]]
metas_spec = Sequence[Tuple[Orange.data.Variable, Sequence[str]]]
