""" Utils """
import io
import os
import gzip
import json

from Orange.data import Table, Domain, TimeVariable, StringVariable, ContinuousVariable

from orangecontrib.bioinformatics.utils import local_cache

#  Support cache with requests_cache module
LOCAL_CACHE_DIR = os.path.join(local_cache, 'resolwe')
try:
    os.makedirs(LOCAL_CACHE_DIR)
except OSError:
    pass

CACHE_BACKEND = 'sqlite'
GENAPI_CACHE = os.path.join(LOCAL_CACHE_DIR, 'GenAPI_cache')
RESOLWEAPI_CACHE = os.path.join(LOCAL_CACHE_DIR, 'ResolweAPI_cache')


def transpose_table(table):
    """Transpose the rows and columns of the table.

    Args:
        table: Data in :obj:`Orange.data.Table`

    Returns:
         Transposed :obj:`Orange.data.Table`. (Genes as columns)
    """

    # TODO: remove this and use Orange.data.Table.transpose

    attrs = table.domain.attributes
    attr = [ContinuousVariable.make(ex['Gene'].value) for ex in table]
    #  Set metas
    new_metas = [
        StringVariable.make(name) if name != 'Time' else TimeVariable.make(name)
        for name in sorted(table.domain.variables[0].attributes.keys())
    ]
    domain = Domain(attr, metas=new_metas)
    meta_values = [[exp.attributes[var.name] for var in domain.metas] for exp in attrs]

    return Table(domain, table.X.transpose(), metas=meta_values)


def etc_to_table(etc_json, time_var=False):
    """Converts data from Json to :obj:`Orange.data.table`

    Args:
        etc_json (dict): Data in json like format from genesis
        time_var (bool): Create column of time points. Default is set to False.

    Returns:
        :obj:`Orange.data.Table`
    """

    variables = []
    time_point = 1
    for time in etc_json['etc']['timePoints']:
        var = ContinuousVariable('TP ' + str(time_point))
        var.attributes['Time'] = str(time)
        variables.append(var)
        time_point += 1

    meta_attr = StringVariable.make('Gene')
    domain = Domain(variables, metas=[meta_attr])

    table = []
    for row in etc_json['etc']['genes']:
        gene_expression = list(etc_json['etc']['genes'][row])
        gene_expression.append(row)
        table.append(gene_expression)

    orange_table = Table.from_list(domain, table)

    if time_var:
        orange_table = transpose_table(orange_table)

    return orange_table


def response_to_json(response):
    if not response.ok:
        response.raise_for_status()

    response_gzipped = io.BytesIO(response.content)
    response_content = io.TextIOWrapper(gzip.GzipFile(fileobj=response_gzipped), encoding='utf-8')

    try:
        json_object = json.load(response_content)
    except ValueError as e:
        raise ValueError('Downloaded data is not a valid JSON') from e

    return json_object
