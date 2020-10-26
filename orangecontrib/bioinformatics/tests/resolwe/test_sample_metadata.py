import unittest
from unittest.mock import MagicMock

from resdk.resources.sample import Sample

from Orange.data import TimeVariable, StringVariable, DiscreteVariable, ContinuousVariable

from orangecontrib.bioinformatics.resolwe.sample_metadata import (
    descriptors,
    flatten_descriptor,
    field_type_to_variable,
    flatten_descriptor_schema,
)

SCHEMA = [
    {
        'name': 'section_1',
        'group': [
            {'name': 'field_1', 'type': 'basic:string:', 'label': 'Field 1'},
            {
                'name': 'field_2',
                'type': 'basic:string:',
                'label': 'Field 2',
                'choices': [
                    {'label': 'Choice 1', 'value': 'choice_1'},
                    {'label': 'Choice 2', 'value': 'choice_2'},
                    {'label': 'Choice 3', 'value': 'choice_3'},
                    {'label': 'Choice 4', 'value': 'choice_4'},
                ],
            },
        ],
    },
    {
        'name': 'section_2',
        'group': [
            {
                'name': 'hidden_section_3',
                'group': [
                    {'name': 'field_3', 'type': 'basic:decimal:', 'label': 'Field 3'},
                    {'name': 'field_4', 'type': 'basic:integer:', 'label': 'Field 4'},
                ],
            },
            {'name': 'field_5', 'type': 'basic:boolean:', 'label': 'Field 5'},
            {'name': 'field_6', 'type': 'basic:date:', 'label': 'Field 6'},
            {'name': 'field_7', 'type': 'basic:datetime:', 'label': 'Field 7'},
        ],
    },
]

DESCRIPTOR_SCHEMA = {'id': 1, 'slug': 'test-schema', 'name': 'Test schema', 'version': '1.0.0', 'schema': SCHEMA}


class TestSampleMetadata(unittest.TestCase):
    def test_flatten_descriptor_schema(self):
        sample = Sample(id=1, descriptor_schema=DESCRIPTOR_SCHEMA, resolwe=MagicMock())
        schema = flatten_descriptor_schema(sample.descriptor_schema.schema)

        expected_keys = ['field_1', 'field_2', 'field_3', 'field_4', 'field_5', 'field_6', 'field_7']
        self.assertTrue(len(schema.keys()) == 7)
        self.assertTrue(all(key in expected_keys for key in schema.keys()))

        # test if fields with multiple choices are properly handled
        discrete_field = schema['field_2']
        self.assertTrue(len(discrete_field.keys()) == 4)
        self.assertEqual(tuple(index for index, _ in discrete_field['choices'].values()), (0, 1, 2, 3))
        self.assertEqual(
            tuple(val for _, val in discrete_field['choices'].values()),
            ('Choice 1', 'Choice 2', 'Choice 3', 'Choice 4'),
        )

    def test_flatten_descriptor(self):
        descriptor = {
            'section_1': {'field_1': 'field_1_value'},
            'section_2': {'hidden_section_3': {'field_3': 'field_3_value'}, 'field_5': 'field_5_value'},
        }
        sample = Sample(id=1, descriptor=descriptor, resolwe=MagicMock())
        schema = flatten_descriptor(sample.descriptor)

        expected_keys = ['field_1', 'field_3', 'field_5']
        self.assertTrue(len(schema.keys()) == 3)
        self.assertTrue(all(key in expected_keys for key in schema.keys()))

    def test_field_type_to_variable(self):
        sample = Sample(id=1, descriptor_schema=DESCRIPTOR_SCHEMA, resolwe=MagicMock())
        schema = flatten_descriptor_schema(sample.descriptor_schema.schema)

        field = schema['field_1']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, StringVariable)
        self.assertEqual(var.name, field['label'])

        field = schema['field_2']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, DiscreteVariable)
        self.assertEqual(var.name, field['label'])
        self.assertEqual(var.values, tuple(val for _, val in field['choices'].values()))

        field = schema['field_3']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, ContinuousVariable)
        self.assertEqual(var.name, field['label'])

        field = schema['field_4']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, ContinuousVariable)
        self.assertEqual(var.name, field['label'])

        field = schema['field_5']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, DiscreteVariable)
        self.assertEqual(var.name, field['label'])
        self.assertEqual(var.values, ('False', 'True'))

        field = schema['field_6']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, TimeVariable)
        self.assertEqual(var.name, field['label'])

        field = schema['field_7']
        var = field_type_to_variable(field)
        self.assertIsInstance(var, TimeVariable)
        self.assertEqual(var.name, field['label'])

    def test_descriptors(self):
        samples = []
        descriptor = {
            'section_1': {'field_1': 'field_1_value', 'field_2': 'choice_1'},
            'section_2': {
                'hidden_section_3': {'field_3': 123.456, 'field_4': 123},
                'field_5': True,
                'field_6': '2020-11-01',
                'field_7': '2020-11-02 12:00:00',
            },
        }
        samples.append(
            Sample(name='sample1', descriptor=descriptor, descriptor_schema=DESCRIPTOR_SCHEMA, resolwe=MagicMock())
        )

        descriptor = {
            'section_1': {'field_1': 'field_1_value', 'field_2': 'choice_2'},
            'section_2': {
                'hidden_section_3': {'field_3': 123.456, 'field_4': 123},
                'field_5': False,
                'field_6': '2020-11-01',
                'field_7': '2020-11-02 13:00:00',
            },
        }
        samples.append(
            Sample(name='sample2', descriptor=descriptor, descriptor_schema=DESCRIPTOR_SCHEMA, resolwe=MagicMock())
        )

        descriptor = {
            'section_1': {'field_1': 'field_1_value', 'field_2': 'choice_3'},
            'section_2': {
                'hidden_section_3': {'field_3': 123.456, 'field_4': 123},
                'field_5': False,
                'field_6': '2020-11-01',
                'field_7': '2020-11-02 14:00:00',
            },
        }
        samples.append(
            Sample(name='sample3', descriptor=descriptor, descriptor_schema=DESCRIPTOR_SCHEMA, resolwe=MagicMock())
        )

        metadata = descriptors(samples)

        self.assertTrue(len(metadata.keys()) == 10)
        for var in metadata.keys():
            self.assertIsInstance(var, (ContinuousVariable, DiscreteVariable, StringVariable, TimeVariable))

        self.assertTrue(StringVariable('Sample name') in metadata)
        self.assertTrue(StringVariable('Sample slug') in metadata)
        self.assertTrue(StringVariable('Sample ID') in metadata)

        # ensure that discrete variables match sample descriptor
        self.assertEqual(metadata.get(DiscreteVariable('Field 5')), [[1], [0], [0]])
        self.assertEqual(metadata.get(DiscreteVariable('Field 2')), [[0], [1], [2]])
