import json
from typing import Dict, List, Tuple, Optional
from collections import Counter, defaultdict

import numpy as np
from resdk.resources.data import Sample
from resdk.resources.relation import Relation

from Orange.data import Table, Domain, Variable, TimeVariable, DiscreteVariable
from Orange.data.io_util import guess_data_type


def flatten_descriptor(descriptor: Dict[str, dict]) -> Dict[str, dict]:
    """
    Flattens sample descriptor one level deep. Adjust this
    method if there is a need to support descriptors with
    multiple nested choices.

    For example:
        {
          'general': { ... },
          'experiment': {'chip_seq': { ... }, ...}
        }
    Is split into:
        {
          'general.field_name1': { ... },
          'experiment.field_name2': { ... },
          'chip_seq.field_name3': { ... },
          ...
        }

    """
    flat_descriptor = {}
    for section, values in descriptor.items():
        for name, field_value in values.items():
            if isinstance(field_value, dict):
                # its a subsection
                flat_descriptor.update({f'{name}.{key}': val for key, val in field_value.items()})
            else:
                flat_descriptor.update({f'{section}.{name}': field_value})
    return flat_descriptor


def flatten_descriptor_schema(schema: List[dict]) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """
    Flattens descriptor schema one level deep. Adjust this
    method if there is a need to support schemas which have
    multiple layers of nested subgroups.

    For example:
        [
          {'name': 'general', 'group': [ ... ]},
          {'name': 'experiment', 'group: [ {'name': 'chip_seq', 'group': [ ... ]}, ... ]}
        ]
    Is split into:
        {
          'general.field_name1': { ... },
          'experiment.field_name2': { ... },
          'chip_seq.field_name3': { ... },
          ...
        }

    """

    flat_schema = {}
    section_name_to_label = {}
    for section in schema:
        section_name_to_label.update({section['name']: section['label']})

        for field in section['group']:
            if 'group' in field:
                # field is a subgroup
                flat_schema.update({f"{field['name']}.{sub_field['name']}": sub_field for sub_field in field['group']})
                section_name_to_label.update({sub_field['name']: sub_field['label'] for sub_field in field['group']})
            else:
                flat_schema.update({f"{section['name']}.{field['name']}": field})

    return flat_schema, section_name_to_label


def handle_field_type(field: dict, field_values: np.array) -> Optional[Tuple[Variable, list]]:
    field_name = field['label']
    field_type = field['type']

    # Multiple choices fields are always discrete
    if 'choices' in field:
        discrete_values = sorted((choice['value'] for choice in field['choices']))
        var = DiscreteVariable(field_name, values=(str(val) for val in discrete_values))
        return var, [var.to_val(value) for value in field_values]
    else:
        discrete_values, _field_values, col_type = guess_data_type(field_values)
        if discrete_values is not None:
            # indexes = np.nonzero(_field_values[:, None] == discrete_values)[1]
            var = DiscreteVariable(field_name, values=(str(val) for val in discrete_values))
            return var, [var.to_val(str(value)) for value in _field_values]

        var = col_type(field_name)
        # Orange does not handle Time variables correctly
        if field_type in ('basic:date:', 'basic:datetime:'):
            _var = TimeVariable(field_name)
            _field_values = [_var.parse(time) for time in field_values]
            var = _var

        return var, _field_values


def descriptors(samples: List[Sample]) -> Dict[Variable, List[List[str]]]:
    def _drop_duplicate_descriptor_schemas(schemas: List[Dict[str, dict]]) -> List[Dict[str, dict]]:
        """
        Return a list where descriptor schemas are ensured to be unique.
        """
        return [json.loads(_schema) for _schema in {json.dumps(_schema) for _schema in schemas}]

    flatten_descriptor_schemas = [flatten_descriptor_schema(sample.descriptor_schema.schema) for sample in samples]
    schema, *tail = _drop_duplicate_descriptor_schemas([schema for schema, _ in flatten_descriptor_schemas])

    if len(tail):
        raise ValueError('Descriptor schema is not unique for selected samples.')

    section_name_to_label = [name_to_label for _, name_to_label in flatten_descriptor_schemas][0]
    meta_key_to_field = {key: field for key, field in schema.items()}
    meta_key_to_field.update(
        {
            'other.sample_name': {'label': 'Sample Name', 'type': 'basic:string:'},
            'other.sample_slug': {'label': 'Sample Slug', 'type': 'basic:string:'},
            'other.sample_id': {'label': 'Sample ID', 'type': 'basic:string:'},
        }
    )

    metadata = defaultdict(list)
    for sample in samples:
        descriptor = flatten_descriptor(sample.descriptor)

        for field_name in schema.keys():
            field = schema[field_name]
            descriptor_value = descriptor.get(field_name, '')

            if field['type'] == 'basic:boolean:' and descriptor_value is not None:
                descriptor_value = str(descriptor_value)

            metadata[field_name].append(descriptor_value)

        metadata['other.sample_name'].append(sample.name)
        metadata['other.sample_slug'].append(sample.slug)
        metadata['other.sample_id'].append(sample.id)

    # Filter fields with empty values
    metadata = {k: v for k, v in metadata.items() if not all(x == '' for x in v)}

    # handle duplicates
    split_section_field = [label.split(sep='.') for label in metadata.keys()]
    field_names = [field_name for _, field_name in split_section_field]
    sections = [section for section, _ in split_section_field]
    duplicates = {item for item, count in Counter(field_names).items() if count > 1}
    for key, section, field_name in zip(metadata.keys(), sections, field_names):
        if field_name in duplicates:
            meta_key_to_field[key]['label'] = f"{section_name_to_label[section]} {meta_key_to_field[key]['label']}"

    return dict(handle_field_type(meta_key_to_field[key], np.array(values)) for key, values in metadata.items())


def relations(samples: List[Sample], _relations: List[Relation]) -> Dict[Variable, List[List[str]]]:
    sample_ids = [sample.id for sample in samples]
    metadata = defaultdict(list)

    for relation in _relations:
        discrete_values = sorted({p['label'] for p in relation.partitions})
        sample_id_to_index = {p['entity']: discrete_values.index(p['label']) for p in relation.partitions}

        for samp_id in sample_ids:
            metadata[DiscreteVariable(relation.category, values=discrete_values)].append([sample_id_to_index[samp_id]])

    return metadata


def get_target_column(source: Table, incoming: Table) -> Variable:
    """
    Return matching column from both tables if there is any.
    """

    target_column_names = ('Sample ID', 'Sample slug', 'Sample name', 'Sample')

    column_candidates = [
        column
        for column in set(incoming.domain.metas) & set(source.domain.metas)
        if column.name in target_column_names
    ]

    if not column_candidates:
        raise ValueError(f'No matching columns detected. Use any of {target_column_names}')

    # Use the first column candidate. It does not really matter which one is used.
    return column_candidates[0]


def merge_clinical_data(source: Table, incoming: Table) -> Table:
    target_column = get_target_column(source, incoming)
    source_target_column, _ = source.get_column_view(source.domain.index(target_column))
    incoming_target_column, _ = incoming.get_column_view(incoming.domain.index(target_column))

    # We can not be sure that clinical data is in the same order (sample wise) as source data.
    # Map incoming row indexes to row indexes in source data
    # Example: X = ['B', 'C', 'A']
    #          Y = ['C', 'A', 'B']
    #          mask = [2, 0, 1] -> where should an element from array Y be to match array X
    source_index_by_row = {index: str(row_value) for index, row_value in enumerate(source_target_column)}
    incoming_row_by_index = {str(row_value): index for index, row_value in enumerate(incoming_target_column)}
    source_row_order = [
        incoming_row_by_index[source_index_by_row[row_idx]] for row_idx in range(len(source_target_column))
    ]

    # Columns from incoming table to include in source table
    incoming_columns = tuple(column for column in incoming.domain.metas if column not in source.domain.metas)

    # create new domain and transform data table
    domain = Domain(source.domain.variables, source.domain.class_vars, source.domain.metas + incoming_columns)
    # table = source.copy()
    table = source.transform(domain)

    for col in incoming_columns:
        col_view, _ = incoming.get_column_view(incoming.domain.index(col))
        table[:, col] = col_view[source_row_order].reshape(col_view.shape[0], 1)

    return table
