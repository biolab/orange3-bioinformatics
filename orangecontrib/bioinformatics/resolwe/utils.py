""" Utils """
import io
import os
import gzip
import json

from Orange.data import Table, Domain, StringVariable, ContinuousVariable

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


def etc_to_table(etc_json, transpose=False):
    """Converts data from Json to :obj:`Orange.data.table`

    Args:
        etc_json (dict): Data in json like format from genesis
        transpose (bool): Transpose table so that genes are in columns.
                          Default is set to False.

    Returns:
        :obj:`Orange.data.Table`
    """

    variables = []
    time_point = 1
    for time in etc_json['etc']['timePoints']:
        var = ContinuousVariable(str(time))
        var.attributes['Time point'] = str(time_point)
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
    if transpose:
        orange_table = Table.transpose(
            orange_table,
            feature_names_column=meta_attr.name,
            meta_attr_name='Time',
            remove_redundant_inst=True,
        )

        old_domain = orange_table.domain
        new_domain = Domain(
            old_domain.attributes,
            old_domain.class_vars,
            [ContinuousVariable(var.name) for var in old_domain.metas],
        )

        time_as_floats = orange_table.get_column('Time').astype(int).reshape(-1, 1)
        orange_table = orange_table.transform(new_domain)
        orange_table[:, new_domain['Time']] = time_as_floats
    return orange_table


def response_to_json(response):
    if not response.ok:
        response.raise_for_status()

    response_gzipped = io.BytesIO(response.content)
    response_content = io.TextIOWrapper(
        gzip.GzipFile(fileobj=response_gzipped), encoding='utf-8'
    )

    try:
        json_object = json.load(response_content)
    except ValueError as e:
        raise ValueError('Downloaded data is not a valid JSON') from e

    return json_object
