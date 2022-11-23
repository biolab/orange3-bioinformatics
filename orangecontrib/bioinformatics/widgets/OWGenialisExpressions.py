import asyncio
from enum import IntEnum
from typing import Any, Dict, List, Tuple, Optional, NamedTuple
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
import resdk.tables
from resdk.resources.data import Data
from resdk.resources.utils import iterate_schema

from AnyQt.QtCore import Qt, QSize, pyqtSignal
from AnyQt.QtWidgets import QTableView, QHeaderView, QAbstractItemView

from Orange.data import Table, Domain, ContinuousVariable
from Orange.util import wrap_callback
from Orange.widgets import gui, widget, settings
from Orange.widgets.widget import Msg, Output
from Orange.data.pandas_compat import vars_from_df
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.itemmodels import PyTableModel

from orangecontrib.bioinformatics import resolwe
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher
from orangecontrib.bioinformatics.preprocess import (
    ZScore,
    LogarithmicScale,
    QuantileTransform,
    QuantileNormalization,
)
from orangecontrib.bioinformatics.ncbi.taxonomy import species_name_to_taxid
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation
from orangecontrib.bioinformatics.widgets.components.resolwe import (
    SignIn,
    SortBy,
    ItemsPerPage,
    PaginationComponent,
    FilterByDateModified,
    QuantileTransformDist,
    NormalizationComponent,
    CollapsibleFilterComponent,
    get_credential_manager,
)


class TableHeader(IntEnum):
    id = 0
    slug = 1
    name = 2
    samples = 3
    species = 4
    created = 5
    modified = 6
    contributor = 7
    description = 8
    tags = 9

    @staticmethod
    def labels():
        return [
            'Id',
            'Slug',
            'Name',
            'Samples',
            'Species',
            'Created',
            'Modified',
            'Contributor',
            'Description',
            'Tags',
        ]


class GenialisExpressionsModel(PyTableModel):
    def __init__(self, parent_widget, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent = parent_widget

    def flags(self, index):
        """
        Disable the row selection by clicking on the first column.
        """
        return (
            Qt.ItemIsEnabled
            if index.column() == 0
            else Qt.ItemIsEnabled | Qt.ItemIsSelectable
        )

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return

        row, column = self.mapToSourceRows(index.row()), index.column()

        role_value = self._roleData.get(row, {}).get(column, {}).get(role)
        if role_value is not None:
            return role_value

        try:
            value = self[row][column]
            tag = self[row][TableHeader.tags]
            slug = self[row][TableHeader.slug]
            url = f'{self.parent.res.url}/{tag}/search/collection/{slug}'
        except IndexError:
            return

        if role == Qt.DisplayRole:
            if isinstance(value, datetime):
                return str(value.date().strftime('%m/%d/%Y'))
            elif column == TableHeader.description:
                if value:
                    return f'{str(value[:100])} ...'
            return str(value)

        elif role == Qt.ToolTipRole:
            if isinstance(value, datetime):
                return str(value.strftime('%b %d, %Y %H:%M'))
            if tag and column == TableHeader.id:
                return url
            return str(value)

        # elif role == Qt.TextAlignmentRole and isinstance(value, Number):
        #     return Qt.AlignRight | Qt.AlignVCenter

        elif role == gui.LinkRole:
            if tag:
                return url

    def set_data(
        self, collections: List[Dict[str, str]], col_to_species: Dict[str, str]
    ):
        def model_row(collection: Dict[str, str]) -> List[Any]:
            first_name = collection.get('contributor', {}).get('first_name', '')
            last_name = collection.get('contributor', {}).get('last_name', '')
            user_name = collection.get('contributor', {}).get('username', '')
            contributor = f'{first_name} {last_name}'.strip()

            tags = collection.get('tags', [])
            tags = [tag for tag in tags if 'community' in tag]
            if len(tags) == 1:
                _, tag = tags[0].split(':')
            else:
                tag = ''

            return [
                collection['id'],
                collection.get('slug', ''),
                collection.get('name', ''),
                collection.get('entity_count', ''),
                col_to_species.get(collection['id'], ''),
                datetime.strptime(
                    collection.get('created', ''), '%Y-%m-%dT%H:%M:%S.%f%z'
                ),
                datetime.strptime(
                    collection.get('modified', ''), '%Y-%m-%dT%H:%M:%S.%f%z'
                ),
                contributor if contributor else user_name,
                collection.get('description', ''),
                tag,
            ]

        self.wrap([model_row(result) for result in collections])


class Expression(NamedTuple):
    type: str
    name: str


class Process(NamedTuple):
    slug: str
    name: str


class DataOutputOptions(NamedTuple):
    expression_type: Tuple[Expression]
    expression_sources: Tuple[str]
    process: Tuple[Process]


def available_data_output_options(data_objects: List[Data]) -> DataOutputOptions:
    """
    Traverse the data objects in the selected collection and store the
    information regarding available expression types, process sources/slugs.

    The method returns a named tuple (`DataOutputOptions`) which used for
    creating radio buttons in the control area.
    """
    expression_types = sorted({data.output['exp_type'] for data in data_objects})
    expression_types = (Expression('rc', 'Read Counts'),) + tuple(
        Expression(exp_type, exp_type) for exp_type in expression_types
    )

    process_slugs = sorted(
        {(data.process.slug, data.process.name) for data in data_objects}
    )
    process_slugs = tuple(
        Process(proc_slug, proc_name) for proc_slug, proc_name in process_slugs
    )

    expression_sources = tuple({data.output['source'] for data in data_objects})

    return DataOutputOptions(
        expression_type=expression_types,
        expression_sources=expression_sources,
        process=process_slugs,
    )


class OWGenialisExpressions(widget.OWWidget, ConcurrentWidgetMixin):
    name = 'Genialis Expressions'
    priority = 30
    want_main_area = True
    want_control_area = True
    icon = '../widgets/icons/OWGenialisExpressions.svg'

    pagination_availability = pyqtSignal(bool, bool)

    norm_component = settings.SettingProvider(NormalizationComponent)
    pagination_component = settings.SettingProvider(PaginationComponent)
    filter_component = settings.SettingProvider(CollapsibleFilterComponent)

    exp_type: int
    exp_type = settings.Setting(1, schema_only=True)

    proc_slug: int
    proc_slug = settings.Setting(0, schema_only=True)

    exp_source: int
    exp_source = settings.Setting(0, schema_only=True)

    append_qc_data: bool
    append_qc_data = settings.Setting(False, schema_only=True)

    auto_commit: bool
    auto_commit = settings.Setting(False, schema_only=True)

    class Outputs:
        table = Output('Expressions', Table)

    class Warning(widget.OWWidget.Warning):
        no_expressions = Msg('Expression data objects not found.')
        no_data_objects = Msg('No expression data matches the selected options.')
        unexpected_feature_type = Msg(
            'Can not import expression data, unexpected feature type "{}".'
        )
        multiple_feature_type = Msg(
            'Can not import expression data, multiple feature types found.'
        )

    def __init__(self):
        super().__init__()
        ConcurrentWidgetMixin.__init__(self)

        self._res: Optional[resolwe.resapi.ResolweAPI] = None

        # Store collection ID from currently selected row
        self.selected_collection_id: Optional[str] = None
        # Store data output options
        self.data_output_options: Optional[DataOutputOptions] = None
        # Cache output data table
        self.data_table: Optional[Table] = None
        # Cache clinical metadata
        self.clinical_metadata: Optional[Table] = None

        self.exp_type_combo = gui.comboBox(
            self.controlArea,
            self,
            'exp_type',
            label='Expression Type',
            callback=self.on_output_option_changed,
        )
        self.proc_slug_combo = gui.comboBox(
            self.controlArea,
            self,
            'proc_slug',
            label='Process Name',
            callback=self.on_output_option_changed,
        )
        self.exp_source_combo = gui.comboBox(
            self.controlArea,
            self,
            'exp_source',
            label='Expression source',
            callback=self.on_output_option_changed,
        )

        self.norm_component = NormalizationComponent(self, self.controlArea)
        self.norm_component.options_changed.connect(self.on_normalization_changed)

        box = gui.widgetBox(self.controlArea, 'Sample QC')
        gui.checkBox(
            box,
            self,
            'append_qc_data',
            'Append QC data',
            callback=self.on_output_option_changed,
        )

        gui.rubber(self.controlArea)
        box = gui.widgetBox(self.controlArea, 'Sign in')
        self.user_info = gui.label(box, self, '')
        self.server_info = gui.label(box, self, '')

        box = gui.widgetBox(box, orientation=Qt.Horizontal)
        self.sign_in_btn = gui.button(
            box, self, 'Sign In', callback=self.sign_in, autoDefault=False
        )
        self.sign_out_btn = gui.button(
            box, self, 'Sign Out', callback=self.sign_out, autoDefault=False
        )

        self.commit_button = gui.auto_commit(
            self.controlArea, self, 'auto_commit', '&Commit', box=False
        )
        self.commit_button.button.setAutoDefault(False)

        # Main area
        self.table_view = QTableView()
        self.table_view.setAlternatingRowColors(True)
        self.table_view.viewport().setMouseTracking(True)
        self.table_view.setShowGrid(False)
        self.table_view.verticalHeader().hide()
        self.table_view.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents
        )
        self.table_view.horizontalHeader().setStretchLastSection(True)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.SingleSelection)

        self.model = GenialisExpressionsModel(self)
        self.model.setHorizontalHeaderLabels(TableHeader.labels())
        self.table_view.setModel(self.model)
        self.table_view.selectionModel().selectionChanged.connect(
            self.on_selection_changed
        )

        self.filter_component = CollapsibleFilterComponent(self, self.mainArea)
        self.filter_component.options_changed.connect(self.on_filter_changed)
        self.mainArea.layout().addWidget(self.table_view)
        self.pagination_component = PaginationComponent(self, self.mainArea)
        self.pagination_component.options_changed.connect(self.update_collections_view)

        self.sign_in(silent=True)

    @property
    def res(self):
        return self._res

    @res.setter
    def res(self, value: resolwe.resapi.ResolweAPI):
        if isinstance(value, resolwe.resapi.ResolweAPI):
            self._res = value
            self.update_user_status()
            self.update_collections_view()
            self.__invalidate()
            self.Outputs.table.send(None)

    def __invalidate(self):
        self.data_table = None
        self.selected_collection_id = None
        self.clinical_metadata = None

        self.data_output_options = None
        self.exp_type_combo.clear()
        self.proc_slug_combo.clear()
        self.exp_source_combo.clear()

        self.Outputs.table.send(None)
        self.Warning.no_expressions.clear()
        self.Warning.multiple_feature_type.clear()
        self.Warning.unexpected_feature_type.clear()
        self.Warning.no_data_objects.clear()

    def update_user_status(self):
        user = self.res.get_currently_logged_user()

        if user:
            user_info = (
                f"{user[0].get('first_name', '')} "
                f"{user[0].get('last_name', '')}".strip()
            )
            user_info = (
                f"User: {user_info if user_info else user[0].get('username', '')}"
            )
            self.sign_in_btn.setEnabled(False)
            self.sign_out_btn.setEnabled(True)
        else:
            user_info = 'User: Anonymous'
            self.sign_in_btn.setEnabled(True)
            self.sign_out_btn.setEnabled(False)

        self.user_info.setText(user_info)
        self.server_info.setText(f'Server: {self.res.url[8:]}')

    def sign_in(self, silent=False):
        dialog = SignIn(self, server_type=resolwe.RESOLWE_PLATFORM)

        if silent:
            dialog.sign_in()
            if dialog.resolwe_instance is not None:
                self.res = dialog.resolwe_instance
            else:
                self.res = resolwe.connect(
                    url=resolwe.resapi.DEFAULT_URL, server_type=resolwe.RESOLWE_PLATFORM
                )

        if not silent and dialog.exec():
            self.res = dialog.resolwe_instance

    def sign_out(self):
        # Use public credentials when user signs out
        self.res = resolwe.connect(
            url=resolwe.resapi.DEFAULT_URL, server_type=resolwe.RESOLWE_PLATFORM
        )
        # Remove username and password
        cm = get_credential_manager(resolwe.RESOLWE_PLATFORM)
        if cm.username:
            del cm.username
        if cm.password:
            del cm.password

    def on_filter_changed(self):
        self.pagination_component.reset_pagination()
        self.update_collections_view()

    def get_query_parameters(self) -> Dict[str, str]:
        params = {
            'limit': ItemsPerPage.values()[self.pagination_component.items_per_page],
            'offset': self.pagination_component.offset,
            'ordering': SortBy.values()[self.filter_component.sort_by],
        }

        if self.filter_component.filter_by_full_text:
            params.update({'text': self.filter_component.filter_by_full_text})

        if self.filter_component.filter_by_name:
            params.update({'name__icontains': self.filter_component.filter_by_name})

        if self.filter_component.filter_by_contrib:
            params.update({'contributor_name': self.filter_component.filter_by_contrib})

        if self.filter_component.filter_by_owner:
            params.update({'owners_name': self.filter_component.filter_by_owner})

        last_modified = FilterByDateModified.values()[
            self.filter_component.filter_by_modified
        ]
        if last_modified:
            params.update({'modified__gte': last_modified.isoformat()})

        return params

    def get_collections(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        # Get response from the server
        collections = self.res.get_collections(**self.get_query_parameters())
        # Loop trough collections and store ids
        collection_ids = [
            collection['id'] for collection in collections.get('results', [])
        ]
        # Get species by collection ids
        collection_to_species = self.res.get_species(collection_ids)

        return collections, collection_to_species

    def update_collections_view(self):
        collections, collection_to_species = self.get_collections()

        # Pass the results to data model
        self.model.set_data(collections.get('results', []), collection_to_species)
        self.table_view.setItemDelegateForColumn(
            TableHeader.id, gui.LinkStyledItemDelegate(self.table_view)
        )
        self.table_view.setColumnHidden(TableHeader.slug, True)
        self.table_view.setColumnHidden(TableHeader.tags, True)

        # Check pagination parameters and emit pagination_availability signal
        next_page = True if collections.get('next') else False
        previous_page = True if collections.get('previous') else False
        self.pagination_availability.emit(next_page, previous_page)

    def normalize(self, table: Table) -> Optional[Table]:
        if not table:
            return

        if self.norm_component.quantile_norm:
            table = QuantileNormalization()(table)

        if self.norm_component.log_norm:
            table = LogarithmicScale()(table)

        if self.norm_component.z_score_norm:
            table = ZScore(axis=self.norm_component.z_score_axis)(table)

        if self.norm_component.quantile_transform:
            axis = self.norm_component.quantile_transform_axis
            quantiles = table.X.shape[int(not axis)]
            distribution = QuantileTransformDist.values()[
                self.norm_component.quantile_transform_dist
            ]
            table = QuantileTransform(
                axis=axis, n_quantiles=quantiles, output_distribution=distribution
            )(table)

        return table

    @gui.deferred
    def commit(self):
        self.Warning.no_data_objects.clear()
        self.cancel()
        self.start(self.runner)

    def on_output_option_changed(self):
        self.data_table = None
        self.commit.deferred()

    def on_clinical_data_changed(self):
        self.clinical_metadata = self.fetch_clinical_metadata()
        self.commit.deferred()

    def on_normalization_changed(self):
        self.commit.deferred()

    def on_selection_changed(self):
        self.__invalidate()

        collection_id: str = self.get_selected_row_data(TableHeader.id)
        if not collection_id:
            return

        self.selected_collection_id = collection_id
        data_objects = self.res.get_expression_data_objects(collection_id)
        self.data_output_options = available_data_output_options(data_objects)

        self.exp_type_combo.addItems(
            exp_name for _, exp_name in self.data_output_options.expression_type
        )
        if self.exp_type >= len(self.data_output_options.expression_type):
            self.exp_type = 0
        self.exp_type_combo.setCurrentIndex(self.exp_type)

        self.proc_slug_combo.addItems(
            proc_name for _, proc_name in self.data_output_options.process
        )
        if self.proc_slug >= len(self.data_output_options.process):
            self.proc_slug = 0
        self.proc_slug_combo.setCurrentIndex(self.proc_slug)

        self.exp_source_combo.addItems(self.data_output_options.expression_sources)
        if self.exp_source >= len(self.data_output_options.expression_sources):
            self.exp_source = 0
        self.exp_source_combo.setCurrentIndex(self.exp_source)

        if not data_objects:
            self.Warning.no_expressions()
            return

        # Note: This here is to handle an edge case where we get
        #       different 'feature_type' data object in a collection.
        #       For now we raise a warning, but in the future we should
        #       discuss about how to properly handle different types of features.
        feature_types = {data.output['feature_type'] for data in data_objects}

        if len(feature_types) == 1 and 'gene' not in feature_types:
            self.Warning.unexpected_feature_type(feature_types.pop())
            # self.data_objects = []
            return

        if len(feature_types) > 1:
            self.Warning.multiple_feature_type()
            # self.data_objects = []
            return

        self.on_output_option_changed()

    def get_selected_row_data(self, column: int) -> Optional[str]:
        selection_model = self.table_view.selectionModel()
        rows = selection_model.selectedRows(column=column)
        if not rows:
            return

        return rows[0].data()

    def on_done(self, table: Table):
        if table:
            self.Outputs.table.send(table)

    def on_exception(self, ex):
        # if isinstance(ex, ResolweDataObjectsNotFound):
        #     self.Warning.no_data_objects()
        #     self.Outputs.table.send(None)
        #     self.data_table = None
        #     self.info.set_output_summary(StateInfo.NoOutput)
        #     self.update_info_box()
        # else:
        raise ex

    def on_partial_result(self, result: Any) -> None:
        pass

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()

    def sizeHint(self):
        return QSize(1280, 620)

    def runner(self, state: TaskState) -> Table:
        exp_type = self.data_output_options.expression_type[self.exp_type].type
        exp_source = self.data_output_options.expression_sources[self.exp_source]
        proc_slug = self.data_output_options.process[self.proc_slug].slug
        collection_id = self.selected_collection_id

        table = self.data_table
        progress_steps_download = iter(np.linspace(0, 50, 2))

        def callback(i: float, status=""):
            state.set_progress_value(i * 100)
            if status:
                state.set_status(status)
            if state.is_interruption_requested():
                raise Exception

        if not table:
            collection = self.res.get_collection_by_id(collection_id)
            coll_table = resdk.tables.RNATables(
                collection,
                expression_source=exp_source,
                expression_process_slug=proc_slug,
                progress_callable=wrap_callback(callback, end=0.5),
            )
            species = coll_table._data[0].output['species']
            sample = coll_table._samples[0]

            state.set_status('Downloading ...')
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            df_exp = coll_table.exp if exp_type != 'rc' else coll_table.rc
            df_exp = df_exp.rename(index=coll_table.readable_index)
            df_metas = coll_table.meta
            df_metas = df_metas.rename(index=coll_table.readable_index)
            df_qc = None
            if self.append_qc_data:
                # TODO: check if there is a way to detect if collection
                #       table contains QC data
                try:
                    df_qc = coll_table.qc
                    df_qc = df_qc.rename(index=coll_table.readable_index)
                except ValueError:
                    pass
            loop.close()

            state.set_status('To data table ...')

            duplicates = {
                item
                for item, count in Counter(
                    [
                        label.split('.')[1]
                        for label in df_metas.columns.to_list()
                        if '.' in label
                    ]
                ).items()
                if count > 1
            }

            # what happens if there is more nested sections?
            section_name_to_label = {
                section['name']: section['label']
                for section in sample.descriptor_schema.schema
            }

            column_labels = {}
            for field_schema, _, path in iterate_schema(
                sample.descriptor, sample.descriptor_schema.schema, path=''
            ):
                path = path[1:]  # this is ugly, but cant go around it
                if path not in df_metas.columns:
                    continue
                label = field_schema['label']
                section_name, field_name = path.split('.')
                column_labels[path] = (
                    label
                    if field_name not in duplicates
                    else f'{section_name_to_label[section_name]} - {label}'
                )

            df_exp = df_exp.reset_index(drop=True)
            df_metas = df_metas.astype('object')
            df_metas = df_metas.fillna(np.nan)
            df_metas = df_metas.replace('nan', np.nan)
            df_metas = df_metas.rename(columns=column_labels)
            if df_qc is not None:
                df_metas = pd.merge(df_metas, df_qc, left_index=True, right_index=True)

            xym, domain_metas = vars_from_df(df_metas)
            x, _, m = xym
            x_metas = np.hstack((x, m))
            attrs = [ContinuousVariable(col) for col in df_exp.columns]
            metas = domain_metas.attributes + domain_metas.metas
            domain = Domain(attrs, metas=metas)
            table = Table.from_numpy(domain, df_exp.to_numpy(), metas=x_metas)

            state.set_progress_value(next(progress_steps_download))

            state.set_status('Matching genes ...')
            progress_steps_gm = iter(np.linspace(50, 99, len(coll_table.gene_ids)))

            def gm_callback():
                state.set_progress_value(next(progress_steps_gm))

            tax_id = species_name_to_taxid(species)
            gm = GeneMatcher(tax_id, progress_callback=gm_callback)
            table = gm.match_table_attributes(table, rename=True)
            table.attributes[TableAnnotation.tax_id] = tax_id
            table.attributes[TableAnnotation.gene_as_attr_name] = True
            table.attributes[TableAnnotation.gene_id_attribute] = 'Entrez ID'
            self.data_table = table

        state.set_status('Normalizing ...')
        table = self.normalize(table)
        state.set_progress_value(100)

        return table


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWGenialisExpressions).run()
