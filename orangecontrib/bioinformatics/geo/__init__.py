"""
NCBI GEO module


GEO datasets are hosted on datasets.biolab.si and are not part of
a bioinformatics knowledge base (we can't use bioinformatics.utils.serverfiles module).

"""
import os
import json

import numpy as np
from serverfiles import LocalFiles, ServerFiles

from Orange.data import Table, Domain, StringVariable, DiscreteVariable
from Orange.data import filter as table_filter
from Orange.misc.environ import data_dir

from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation

domain = 'geo'
_local_cache_path = os.path.join(data_dir(), domain)
_all_info_file = os.path.join(_local_cache_path, '__INFO__')
_server_url = 'http://download.biolab.si/datasets/geo/'
pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed/{}'

server_files = ServerFiles(server=_server_url)
local_files = LocalFiles(_local_cache_path, serverfiles=server_files)


def is_cached(gds_id):
    return os.path.exists(os.path.join(_local_cache_path, gds_id + '.tab'))


def info_cache(f):
    """Store content of __INFO__ file locally."""

    def wrapper():
        if not os.path.isfile(_all_info_file) or os.path.getsize(_all_info_file) == 0:
            with open(_all_info_file, 'w') as fp:
                json.dump(list(server_files.allinfo().items()), fp)

        return f()

    return wrapper


@info_cache
def dataset_all_info():
    """
    Return a list of all GEO datasets from databases.biolab.si.
    """

    with open(_all_info_file, 'r') as fp:
        return json.load(fp)


def dataset_download(gds_id, samples=None, transpose=False, callback=None):
    file_name = '{}.tab'.format(gds_id)
    file_path = local_files.localpath_download(file_name, extract=True, callback=callback)

    table = Table(file_path)
    title = table.name
    gds_info = local_files.info(file_name)
    table_annotations = {TableAnnotation.tax_id: gds_info['taxid']}

    if callback:
        callback()

    if samples is not None:
        filters = [table_filter.FilterStringList(sample, sample_types) for sample, sample_types in samples.items()]
        table = table_filter.Values(filters)(table)

        column_values = []
        for meta_var in samples.keys():
            column_values.append(table.get_column_view(table.domain[meta_var])[0])

        class_values = list(map('|'.join, zip(*column_values)))

        _class_values = list(set(class_values))
        map_class_values = {value: key for (key, value) in enumerate(_class_values)}
        class_var = DiscreteVariable(name='class', values=_class_values)
        _domain = Domain(table.domain.attributes, table.domain.class_vars + (class_var,), table.domain.metas)

        table = table.transform(_domain)
        col, _ = table.get_column_view(class_var)
        col[:] = [map_class_values[class_val] for class_val in class_values]

    if transpose:
        table = Table.transpose(table, feature_names_column='sample_id', meta_attr_name='genes')

        # When transposing a table, variable.attributes get picked up as numerical values instead of strings.
        # We need to convert from Continuous to StringVariable
        _genes = [
            [str(int(gene)) if not np.isnan(gene) else '?']
            for gene in table.get_column_view('Entrez ID')[0].astype(np.float64)
        ]
        new_var = StringVariable('Entrez ID')
        metas = [var for var in table.domain.metas if var.name != 'Entrez ID'] + [new_var]
        new_domain = Domain(table.domain.attributes, table.domain.class_vars, metas)
        table = table.transform(new_domain)
        table[:, new_var] = _genes

        # table name is lost after transpose
        table.name = title

        table_annotations[TableAnnotation.gene_as_attr_name] = not gds_info[TableAnnotation.gene_as_attr_name]
        table_annotations[TableAnnotation.gene_id_column] = gds_info[TableAnnotation.gene_id_attribute]
    else:
        table_annotations[TableAnnotation.gene_as_attr_name] = gds_info[TableAnnotation.gene_as_attr_name]
        table_annotations[TableAnnotation.gene_id_attribute] = gds_info[TableAnnotation.gene_id_attribute]

    if callback:
        callback()

    table.attributes = table_annotations
    return table
