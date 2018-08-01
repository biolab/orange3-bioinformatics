""" Orange data table and other data manipulation utilities """
import numbers
import numpy
import Orange.data

from operator import itemgetter
from typing import Sequence, Tuple

from Orange.widgets.widget import OWWidget, Msg

# Data hints variables

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


def append_columns(data, attributes=(), class_vars=(), metas=()):
    # type: (Orange.data.Table, col_spec, col_spec, metas_spec) -> Orange.data.Table
    """ Append a set of columns to a data table.


    Parameters
    ----------
    data : Orange.data.Table
        Primary table.
    attributes : Sequence[Tuple[Orange.data.Variable], Sequence[float]]
        A Sequence of variable and column data tuples to append to the
        `data`.
    class_vars : Sequence[Tuple[Orange.data.Variable], Sequence[float]]
        A Sequence of variable and column data tuples to append to the
        `data` in the
    metas : Sequence[Tuple[Orange.data.Variable], Sequence[str]]
        A Sequence of variable and column data tuples to append to the
        `data`

    Returns
    -------
    data : Orange.data.Table
        A copy of the original `data` input extended with all columns from
        `attributes`, `class_vars`, `metas` parameters

    Note
    ----
    All variables in the original and new columns should be distinct.
    """
    domain = data.domain
    new_attributes = tuple(map(itemgetter(0), attributes))
    new_class_vars = tuple(map(itemgetter(0), class_vars))
    new_metas = tuple(map(itemgetter(0), metas))

    new_domain = Orange.data.Domain(
        domain.attributes + new_attributes,
        domain.class_vars + new_class_vars,
        domain.metas + new_metas
    )

    def ascolumn(array, n):
        # type: (Sequence[float], int) -> numpy.ndarray
        array = numpy.asarray(array)
        if array.ndim < 2:
            array = array.reshape((n, 1))
        return array
    N = len(data)

    attr_cols = [ascolumn(col, N) for _, col in attributes]
    class_cols = [ascolumn(col, N) for _, col in class_vars]
    meta_cols = [ascolumn(col, N) for _, col in metas]

    new_data = data.from_table(new_domain, data)

    for i, (var, col) in enumerate(zip(new_attributes, attr_cols),
                                   start=len(domain.attributes)):
        assert new_data.domain.attributes[i] is var
        new_data.X[:, i] = col.ravel()

    for i, (var, col) in enumerate(zip(new_class_vars, class_cols),
                                   start=len(domain.class_vars)):
        assert new_data.domain.class_vars[i] is var
        new_data._Y[:, i] = col.ravel()

    for i, (var, col) in enumerate(zip(new_metas, meta_cols),
                                   start=len(domain.metas)):
        assert new_data.domain.metas[i] is var
        new_data.metas[:, i] = col.ravel()

    return new_data
