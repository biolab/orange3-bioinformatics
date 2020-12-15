import unittest
from unittest.mock import MagicMock

from resdk.resources.sample import Sample

from Orange.data import TimeVariable, StringVariable, DiscreteVariable, ContinuousVariable

from orangecontrib.bioinformatics.resolwe.sample_metadata import (
    descriptors,
    flatten_descriptor,
    flatten_descriptor_schema,
)

SCHEMA = [
    {
        'name': 'section_1',
        'label': 'Section 1',
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
        'label': 'Section 2',
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
            {'name': 'field_1', 'type': 'basic:string:', 'label': 'Field 1'},
        ],
    },
]

DESCRIPTOR_SCHEMA = {'id': 1, 'slug': 'test-schema', 'name': 'Test schema', 'version': '1.0.0', 'schema': SCHEMA}


class TestSampleMetadata(unittest.TestCase):
    def test_flatten_descriptor_schema(self):
        sample = Sample(id=1, descriptor_schema=DESCRIPTOR_SCHEMA, resolwe=MagicMock())
        schema, _ = flatten_descriptor_schema(sample.descriptor_schema.schema)

        expected_keys = [
            'section_1.field_1',
            'section_1.field_2',
            'hidden_section_3.field_3',
            'hidden_section_3.field_4',
            'section_2.field_5',
            'section_2.field_6',
            'section_2.field_7',
            'section_2.field_1',
        ]
        self.assertTrue(len(schema.keys()) == 8)
        self.assertTrue(all(key in expected_keys for key in schema.keys()))

        # test if fields with multiple choices are properly handled
        discrete_field = schema['section_1.field_2']
        self.assertTrue(len(discrete_field.keys()) == 4)
        self.assertEqual(
            ('choice_1', 'choice_2', 'choice_3', 'choice_4'),
            tuple(choice['value'] for choice in discrete_field['choices']),
        )

    def test_flatten_descriptor(self):
        descriptor = {
            'section_1': {'field_1': 'field_1_value'},
            'section_2': {'hidden_section_3': {'field_3': 'field_3_value'}, 'field_5': 'field_5_value'},
        }
        sample = Sample(id=1, descriptor=descriptor, resolwe=MagicMock())
        schema = flatten_descriptor(sample.descriptor)

        expected_keys = ['section_1.field_1', 'hidden_section_3.field_3', 'section_2.field_5']
        self.assertTrue(len(schema.keys()) == 3)
        self.assertTrue(all(key in expected_keys for key in schema.keys()))

    def test_descriptors(self):
        samples = []
        descriptor = {
            'section_1': {'field_1': 'field_1_value1', 'field_2': 'choice_1'},
            'section_2': {
                'hidden_section_3': {'field_3': 123.456, 'field_4': 123},
                'field_5': True,
                'field_6': '2020-11-01',
                'field_7': '2020-11-02 12:00:00',
                'field_1': 'test_dup_field1',
            },
        }
        samples.append(
            Sample(
                name='sample1',
                slug='sample2',
                id=1,
                descriptor=descriptor,
                descriptor_schema=DESCRIPTOR_SCHEMA,
                resolwe=MagicMock(),
            )
        )

        descriptor = {
            'section_1': {'field_1': 'field_1_value2', 'field_2': 'choice_2'},
            'section_2': {
                'hidden_section_3': {'field_3': 123.456, 'field_4': 123},
                'field_5': False,
                'field_6': '2020-11-01',
                'field_7': '2020-11-02 13:00:00',
                'field_1': 'test_dup_field2',
            },
        }
        samples.append(
            Sample(
                name='sample2',
                slug='sample2',
                id=2,
                descriptor=descriptor,
                descriptor_schema=DESCRIPTOR_SCHEMA,
                resolwe=MagicMock(),
            )
        )

        descriptor = {
            'section_1': {'field_1': 'field_1_value3', 'field_2': 'choice_3'},
            'section_2': {
                'hidden_section_3': {'field_3': 123.456, 'field_4': 123},
                'field_5': False,
                'field_6': '2020-11-01',
                'field_7': '2020-11-02 14:00:00',
                'field_1': 'test_dup_field3',
            },
        }
        samples.append(
            Sample(
                name='sample3',
                slug='sample3',
                id=3,
                descriptor=descriptor,
                descriptor_schema=DESCRIPTOR_SCHEMA,
                resolwe=MagicMock(),
            )
        )

        metadata = descriptors(samples)
        self.assertTrue(len(metadata.keys()) == 11)
        for var in metadata.keys():
            self.assertIsInstance(var, (ContinuousVariable, DiscreteVariable, StringVariable, TimeVariable))

        # test for prefix in duplicate fields
        variable_names = {var.name for var in metadata.keys()}
        self.assertTrue('Field 1' not in variable_names)
        self.assertTrue('Section 1 Field 1' in variable_names)
        self.assertTrue('Section 2 Field 1' in variable_names)

        # ensure that discrete variables match sample descriptor
        self.assertEqual(list(metadata.get(DiscreteVariable('Field 5'))), [1, 0, 0])
        self.assertEqual(list(metadata.get(DiscreteVariable('Field 2'))), [0, 1, 2])
