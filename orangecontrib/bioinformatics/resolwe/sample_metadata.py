import json
from typing import Dict, List, Optional
from collections import defaultdict

from resdk.resources.data import Sample
from resdk.resources.relation import Relation

from Orange.data import Table, Domain, Variable, TimeVariable, StringVariable, DiscreteVariable, ContinuousVariable


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
          'field_name1': { ... },
          'field_name2': { ... },
          'field_name3': { ... },
          ...
        }

    """
    flat_descriptor = {}
    for section, values in descriptor.items():
        for name, value in values.items():
            if isinstance(value, dict):
                # its a subsection
                flat_descriptor.update(value)
            else:
                flat_descriptor.update({name: value})
    return flat_descriptor


def flatten_descriptor_schema(schema: List[dict]) -> Dict[str, dict]:
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
          'field_name1': { ... },
          'field_name2': { ... },
          'field_name3': { ... },
          ...
        }

    """

    def map_choices_with_value(_field: dict) -> dict:
        _field = _field.copy()
        if 'choices' in _field:
            _field['choices'] = {
                choice['value']: (index, choice['label']) for index, choice in enumerate(_field['choices'])
            }
        return _field

    flat_schema = {}
    for section in schema:
        for field in section['group']:
            if 'group' in field:
                # field is a subgroup
                flat_schema.update({field['name']: map_choices_with_value(field) for field in field['group']})
            else:
                flat_schema.update({field['name']: map_choices_with_value(field)})

    return flat_schema


def field_type_to_variable(field: dict) -> Optional[Variable]:
    field_type = field['type']

    # Check for multiple choices first.
    if 'choices' in field:
        # sort by index
        discrete_values = (label for _, label in sorted(field['choices'].values(), key=lambda x: x[0]))
        return DiscreteVariable(field['label'], values=discrete_values)

    if field_type == 'basic:string:':
        return StringVariable(field['label'])

    if field_type in ('basic:decimal:', 'basic:integer:'):
        return ContinuousVariable(field['label'])

    if field_type == 'basic:boolean:':
        return DiscreteVariable(field['label'], values=['False', 'True'])

    if field_type in ('basic:date:', 'basic:datetime:'):
        return TimeVariable(field['label'])


def descriptors(samples: List[Sample]) -> Dict[Variable, List[List[str]]]:
    def _drop_duplicate_descriptor_schemas(schemas: List[Dict[str, dict]]) -> List[Dict[str, dict]]:
        """
        Return a list where descriptor schemas are ensured to be unique.
        """
        return [json.loads(_schema) for _schema in {json.dumps(_schema, sort_keys=True) for _schema in schemas}]

    descriptor_schemas = [flatten_descriptor_schema(sample.descriptor_schema.schema) for sample in samples]
    schema, *tail = _drop_duplicate_descriptor_schemas(descriptor_schemas)

    if len(tail):
        raise ValueError('Descriptor schema is not unique for selected samples.')

    label_to_variable = {field['label']: field_type_to_variable(field) for _, field in schema.items()}

    metadata = defaultdict(list)
    for sample in samples:
        descriptor = flatten_descriptor(sample.descriptor)

        for field_name in schema.keys() & descriptor.keys():
            field = schema[field_name]
            field_choice = descriptor[field_name]

            # if multiple choices or boolean type store discrete value index.
            if 'choices' in field:
                meta_value = field['choices'][field_choice][0]
            elif field['type'] == 'basic:boolean:':
                meta_value = int(field_choice)
            elif field['type'] in ('basic:date:', 'basic:datetime:'):
                time_var = TimeVariable('placeholder')
                meta_value = time_var.parse(field_choice)
            else:
                meta_value = field_choice
            metadata[field['label']].append([meta_value])

        metadata['Sample name'].append([sample.name])
        metadata['Sample slug'].append([sample.slug])
        metadata['Sample ID'].append([sample.id])

    return {label_to_variable.get(key, StringVariable(key)): value for key, value in metadata.items()}


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
