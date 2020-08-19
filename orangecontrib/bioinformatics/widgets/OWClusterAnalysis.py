""" OWClusterAnalysis """
import sys
import itertools

import numpy as np
from scipy.stats import rankdata

from AnyQt.QtCore import Qt, QSize
from AnyQt.QtWidgets import QSplitter, QTableView, QHBoxLayout, QHeaderView, QListWidget

from Orange.data import Table, Domain, StringVariable, DiscreteVariable, ContinuousVariable
from Orange.widgets import settings
from Orange.widgets.gui import spin, vBox, comboBox, listView, widgetBox, doubleSpin, auto_commit, widgetLabel
from Orange.widgets.utils import itemmodels
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import Setting, ContextSetting, PerfectDomainContextHandler, vartype
from Orange.widgets.utils.signals import Input, Output

from orangecontrib.bioinformatics.geneset.utils import GeneSetException
from orangecontrib.bioinformatics.cluster_analysis import DISPLAY_GENE_SETS_COUNT, Cluster, ClusterModel
from orangecontrib.bioinformatics.ncbi.gene.config import ENTREZ_ID
from orangecontrib.bioinformatics.utils.statistics import score_hypergeometric_test
from orangecontrib.bioinformatics.widgets.utils.gui import HTMLDelegate, GeneScoringWidget, GeneSetsSelection
from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID,
    GENE_ID_COLUMN,
    GENE_ID_ATTRIBUTE,
    GENE_AS_ATTRIBUTE_NAME,
)


class ClusterAnalysisContextHandler(PerfectDomainContextHandler):
    def encode_setting(self, context, setting, value):
        if setting.name == 'cluster_indicators':
            value = [(var.name, 100 + vartype(var)) for var in value]
        return super().encode_setting(context, setting, value)

    def decode_setting(self, setting, value, domain=None):
        return (
            [domain[var[0]] for var in value]
            if setting.name == 'cluster_indicators'
            else super().decode_setting(setting, value, domain)
        )


class OWClusterAnalysis(OWWidget):
    name = "Cluster Analysis"
    description = (
        "The widget displays differentially expressed genes that characterize the cluster, "
        "and corresponding gene terms that describe differentially expressed genes"
    )
    icon = "../widgets/icons/OWClusterAnalysis.svg"
    priority = 110

    class Inputs:
        data_table = Input('Data', Table)
        custom_sets = Input('Custom Gene Sets', Table)

    class Outputs:
        selected_data = Output('Selected Data', Table)
        gene_scores = Output('Gene Scores', Table)
        gene_set_scores = Output('Gene Set Scores', Table)

    class Information(OWWidget.Information):
        pass

    class Warning(OWWidget.Warning):
        gene_enrichment = Msg('{}, {}.')
        no_selected_gene_sets = Msg('No gene set selected, select them from Gene Sets box.')

    class Error(OWWidget.Error):
        no_cluster_indicator = Msg('No cluster indicator in the input data')
        gene_as_attributes = Msg('Genes, in the input data, are expected as column names')
        organism_mismatch = Msg('Organism in input data and custom gene sets does not match')
        cluster_batch_conflict = Msg('Cluster and batch must not be the same variable')

    settingsHandler = ClusterAnalysisContextHandler()
    cluster_indicators = ContextSetting([])
    batch_indicator = ContextSetting(None)
    stored_gene_sets_selection = ContextSetting(())

    scoring_method_selection = ContextSetting(0)
    scoring_method_design = ContextSetting(0)
    scoring_test_type = ContextSetting(0)

    # genes filter
    max_gene_count = Setting(20)
    use_gene_count_filter = Setting(True)

    max_gene_p_value = Setting(0.1)
    use_gene_pval_filter = Setting(False)

    max_gene_fdr = Setting(0.1)
    use_gene_fdr_filter = Setting(True)

    # gene sets filter
    min_gs_count = Setting(5)
    use_gs_count_filter = Setting(True)

    max_gs_p_value = Setting(0.1)
    use_gs_pval_filter = Setting(False)

    max_gs_fdr = Setting(0.1)
    use_gs_max_fdr = Setting(True)

    # auto commit results
    auto_commit = settings.Setting(False)

    custom_gene_set_indicator = settings.Setting(None)

    def __init__(self):
        super().__init__()

        # widget attributes
        self.input_data = None
        self.store_input_domain = None
        self.input_genes_names = []
        self.input_genes_ids = []

        self.tax_id = None
        self.use_attr_names = None
        self.gene_id_attribute = None

        # custom gene set input
        self.feature_model = itemmodels.DomainModel(valid_types=(DiscreteVariable, StringVariable))
        self.custom_data = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None
        self.num_of_custom_sets = None

        self.rows_by_cluster = None
        self.rows_by_batch = None
        self.clusters = []
        self.new_cluster_profile = []

        # data model
        self.cluster_info_model = None

        # Info
        info_box = vBox(self.controlArea, 'Info')
        self.input_info = widgetLabel(info_box)

        # Cluster selection
        self.cluster_indicator_model = itemmodels.DomainModel(valid_types=(DiscreteVariable,), separators=False)
        self.cluster_indicator_box = widgetBox(self.controlArea, 'Cluster Indicator')

        self.cluster_indicator_view = listView(
            self.cluster_indicator_box,
            self,
            'cluster_indicators',
            model=self.cluster_indicator_model,
            selectionMode=QListWidget.MultiSelection,
            callback=self.invalidate,
            sizeHint=QSize(256, 70),
        )

        # Batch selection
        self.batch_indicator_model = itemmodels.DomainModel(
            valid_types=(DiscreteVariable,), separators=False, placeholder=""
        )
        box = widgetBox(self.controlArea, 'Batch Indicator')
        self.batch_indicator_combobox = comboBox(
            box,
            self,
            'batch_indicator',
            model=self.batch_indicator_model,
            sendSelectedValue=True,
            callback=self.batch_indicator_changed,
        )

        # Gene scoring
        box = widgetBox(self.controlArea, 'Gene Scoring')
        self.gene_scoring = GeneScoringWidget(box, self)
        self.gene_scoring.set_method_selection_area('scoring_method_selection')
        self.gene_scoring.set_method_design_area('scoring_method_design')
        self.gene_scoring.set_test_type('scoring_test_type')

        # Gene Sets widget
        gene_sets_box = widgetBox(self.controlArea, "Gene Sets")
        self.gs_widget = GeneSetsSelection(gene_sets_box, self, 'stored_gene_sets_selection')
        self.gs_widget.hierarchy_tree_widget.itemClicked.connect(self.__gene_sets_enrichment)

        # custom gene sets area
        box = vBox(self.controlArea, "Custom Gene Sets")

        if self.custom_gene_set_indicator not in self.feature_model:
            self.custom_gene_set_indicator = None

        self.gs_label_combobox = comboBox(
            box,
            self,
            "custom_gene_set_indicator",
            sendSelectedValue=True,
            model=self.feature_model,
            callback=self.handle_custom_gene_sets,
        )
        self.gs_label_combobox.setDisabled(True)

        # main area
        splitter = QSplitter(Qt.Horizontal, self.mainArea)
        self.mainArea.layout().addWidget(splitter)

        genes_filter = widgetBox(splitter, 'Filter Genes', orientation=QHBoxLayout())
        spin(
            genes_filter,
            self,
            'max_gene_count',
            0,
            10000,
            label='Count',
            tooltip='Minimum genes count',
            checked='use_gene_count_filter',
            callback=self.filter_genes,
            callbackOnReturn=True,
            checkCallback=self.filter_genes,
        )

        doubleSpin(
            genes_filter,
            self,
            'max_gene_p_value',
            0.0,
            1.0,
            0.0001,
            label='p-value',
            tooltip='Maximum p-value of the enrichment score',
            checked='use_gene_pval_filter',
            callback=self.filter_genes,
            callbackOnReturn=True,
            checkCallback=self.filter_genes,
        )

        doubleSpin(
            genes_filter,
            self,
            'max_gene_fdr',
            0.0,
            1.0,
            0.0001,
            label='FDR',
            tooltip='Maximum false discovery rate',
            checked='use_gene_fdr_filter',
            callback=self.filter_genes,
            callbackOnReturn=True,
            checkCallback=self.filter_genes,
        )

        gene_sets_filter = widgetBox(splitter, 'Filter Gene Sets', orientation=QHBoxLayout())
        spin(
            gene_sets_filter,
            self,
            'min_gs_count',
            0,
            DISPLAY_GENE_SETS_COUNT,
            label='Count',
            tooltip='Minimum genes count',
            checked='use_gs_count_filter',
            callback=self.filter_gene_sets,
            callbackOnReturn=True,
            checkCallback=self.filter_gene_sets,
        )

        doubleSpin(
            gene_sets_filter,
            self,
            'max_gs_p_value',
            0.0,
            1.0,
            0.0001,
            label='p-value',
            tooltip='Maximum p-value of the enrichment score',
            checked='use_gs_pval_filter',
            callback=self.filter_gene_sets,
            callbackOnReturn=True,
            checkCallback=self.filter_gene_sets,
        )

        doubleSpin(
            gene_sets_filter,
            self,
            'max_gs_fdr',
            0.0,
            1.0,
            0.0001,
            label='FDR',
            tooltip='Maximum false discovery rate',
            checked='use_gs_max_fdr',
            callback=self.filter_gene_sets,
            callbackOnReturn=True,
            checkCallback=self.filter_gene_sets,
        )

        self.cluster_info_view = QTableView()
        self.cluster_info_view.verticalHeader().setVisible(False)
        self.cluster_info_view.setItemDelegate(HTMLDelegate())
        self.cluster_info_view.horizontalHeader().hide()
        self.cluster_info_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

        self.mainArea.layout().addWidget(self.cluster_info_view)

    def sizeHint(self):
        return QSize(800, 600)

    def __update_info_box(self):
        info_string = ''
        if self.input_genes_ids:
            info_string += '{} samples, {} clusters\n'.format(
                self.input_data.X.shape[0], len(self.clusters) if self.clusters else '?'
            )
            info_string += '{:,d} unique genes\n'.format(len(self.input_genes_ids))
        else:
            info_string += 'No genes on input.\n'

        if self.custom_data:
            info_string += '{} marker genes in {} sets\n'.format(self.custom_data.X.shape[0], self.num_of_custom_sets)

        self.input_info.setText(info_string)

    def __set_cluster_info_model(self):
        self.cluster_info_view.setModel(None)

        self.cluster_info_model = ClusterModel(self)
        self.cluster_info_model.add_rows(self.clusters)

        # add model to the view
        self.cluster_info_view.setModel(self.cluster_info_model)
        # call sizeHint function
        self.cluster_info_view.resizeRowsToContents()
        self.cluster_info_view.selectionModel().selectionChanged.connect(self.commit)

    def __create_temp_class_var(self):
        """ See no evil !"""
        cluster_indicator_name = 'Cluster indicators'
        row_profile = None
        new_cluster_values = []
        var_index_lookup = {val: idx for var in self.cluster_indicators for idx, val in enumerate(var.values)}

        cart_prod = itertools.product(*[cluster.values for cluster in self.cluster_indicators])
        for comb in cart_prod:
            new_cluster_values.append(', '.join([val for val in comb]))
            self.new_cluster_profile.append([var_index_lookup[val] for val in comb])

        row_profile_lookup = {
            tuple(profile): indx for indx, (profile, _) in enumerate(zip(self.new_cluster_profile, new_cluster_values))
        }
        for var in self.cluster_indicators:
            if row_profile is None:
                row_profile = np.asarray(self.input_data.get_column_view(var)[0], dtype=int)
            else:
                row_profile = np.vstack((row_profile, np.asarray(self.input_data.get_column_view(var)[0], dtype=int)))

        ca_ind = DiscreteVariable.make(
            cluster_indicator_name, values=[val for val in new_cluster_values], ordered=True
        )

        domain = Domain(
            self.input_data.domain.attributes,
            self.input_data.domain.class_vars,
            self.input_data.domain.metas + (ca_ind,),
        )

        table = self.input_data.transform(domain)
        table[:, ca_ind] = np.array(
            [[row_profile_lookup[tuple(row_profile[:, i])]] for i in range(row_profile.shape[1])]
        )
        self.input_data = table
        return ca_ind

    def __set_clusters(self):
        self.clusters = []
        self.new_cluster_profile = []
        self.cluster_var = None

        if self.cluster_indicators and self.input_data:

            if isinstance(self.cluster_indicators, list) and len(self.cluster_indicators) > 1:
                self.cluster_var = self.__create_temp_class_var()
            else:
                self.cluster_var = self.cluster_indicators[0]

            self.rows_by_cluster = np.asarray(self.input_data.get_column_view(self.cluster_var)[0], dtype=int)
            for index, name in enumerate(self.cluster_var.values):
                cluster = Cluster(name, index)
                self.clusters.append(cluster)
                cluster.set_genes(self.input_genes_names, self.input_genes_ids)

    def __set_batch(self):
        self.Error.cluster_batch_conflict.clear()
        self.rows_by_batch = None

        if self.batch_indicator == self.cluster_var:
            self.Error.cluster_batch_conflict()
            return
        if self.batch_indicator and self.input_data:
            self.rows_by_batch = np.asarray(self.input_data.get_column_view(self.batch_indicator)[0], dtype=int)

    def __set_genes(self):
        self.input_genes_names = []
        self.input_genes_ids = []

        if self.use_attr_names:
            for variable in self.input_data.domain.attributes:
                self.input_genes_names.append(str(variable.name))
                self.input_genes_ids.append(str(variable.attributes.get(self.gene_id_attribute, np.nan)))

    def filter_genes(self):
        if self.cluster_info_model:
            # filter genes
            # note: after gene filter is applied, we need to recalculate gene set enrichment
            self.cluster_info_model.apply_gene_filters(
                self.max_gene_p_value if self.use_gene_pval_filter else None,
                self.max_gene_fdr if self.use_gene_fdr_filter else None,
                self.max_gene_count if self.use_gene_count_filter else None,
            )

            # recalculate gene set enrichment
            self.__gene_sets_enrichment()
            # call sizeHint function
            self.cluster_info_view.resizeRowsToContents()

            # commit changes after filter
            self.commit()

    def filter_gene_sets(self):
        if self.cluster_info_model:
            # filter gene sets
            self.cluster_info_model.apply_gene_sets_filters(
                self.max_gs_p_value if self.use_gs_pval_filter else None,
                self.max_gs_fdr if self.use_gs_max_fdr else None,
                self.min_gs_count if self.use_gs_count_filter else None,
            )

            # call sizeHint function
            self.cluster_info_view.resizeRowsToContents()

    def __gene_enrichment(self):
        design = bool(self.gene_scoring.get_selected_desig())  # if true cluster vs. cluster else cluster vs rest
        test_type = self.gene_scoring.get_selected_test_type()
        method = self.gene_scoring.get_selected_method()
        try:
            if method.score_function == score_hypergeometric_test:
                values = set(np.unique(self.input_data.X))
                if (0 not in values) or (len(values) != 2):
                    raise ValueError('Binary data expected (use Preprocess)')

            self.cluster_info_model.score_genes(
                design=design,
                table_x=self.input_data.X,
                rows_by_cluster=self.rows_by_cluster,
                rows_by_batch=self.rows_by_batch,
                method=method,
                alternative=test_type,
            )
        except ValueError as e:
            self.Warning.gene_enrichment(str(e), 'p-values are set to 1')

    def __gene_sets_enrichment(self):
        if self.input_data:
            self.Warning.no_selected_gene_sets.clear()
            all_sets = self.gs_widget.get_hierarchies()
            selected_sets = self.gs_widget.get_hierarchies(only_selected=True)

            if len(selected_sets) == 0 and len(all_sets) > 0:
                self.Warning.no_selected_gene_sets()

            # save setting on selected hierarchies
            self.stored_gene_sets_selection = tuple(selected_sets)
            ref_genes = set(self.input_genes_ids)

            try:
                self.cluster_info_model.gene_sets_enrichment(self.gs_widget.gs_object, selected_sets, ref_genes)
            except Exception as e:
                # TODO: possible exceptions?

                raise e

            self.filter_gene_sets()

    def invalidate(self, cluster_init=True):
        if self.input_data is not None and self.tax_id is not None:
            self.Warning.gene_enrichment.clear()

            if self.cluster_info_model is not None:
                self.cluster_info_model.cancel()

            self.__set_genes()
            if cluster_init:
                self.__set_clusters()
            self.__set_batch()
            self.__set_cluster_info_model()

            # note: when calling self.__gene_enrichment we calculate gse automatically.
            #       No need to call self.__gene_sets_enrichment here
            self.__gene_enrichment()
            self.__update_info_box()

    def batch_indicator_changed(self):
        self.invalidate(cluster_init=False)

    @Inputs.data_table
    def handle_input(self, data):
        self.closeContext()
        self.Warning.clear()
        self.Error.clear()

        self.input_data = None
        self.store_input_domain = None
        self.stored_gene_sets_selection = ()
        self.input_genes_names = []
        self.input_genes_ids = []
        self.tax_id = None
        self.use_attr_names = None
        self.gene_id_attribute = None
        self.clusters = None

        self.gs_widget.clear()
        self.gs_widget.clear_gene_sets()
        self.cluster_info_view.setModel(None)

        self.cluster_indicators = []
        self.cluster_var = None
        self.batch_indicator = None
        self.cluster_indicator_model.set_domain(None)
        self.batch_indicator_model.set_domain(None)

        self.__update_info_box()

        if data:
            self.input_data = data

            self.cluster_indicator_model.set_domain(self.input_data.domain)
            self.batch_indicator_model.set_domain(self.input_data.domain)

            # For Cluster Indicator do not use categorical variables that contain only one value.
            self.cluster_indicator_model.wrap([item for item in self.cluster_indicator_model if len(item.values) > 1])
            # First value in batch indicator model is a NoneType,
            # we can skip it when we validate categorical variables
            self.batch_indicator_model.wrap(
                self.batch_indicator_model[:1]
                + [item for item in self.batch_indicator_model[1:] if len(item.values) > 1]
            )

            self.tax_id = self.input_data.attributes.get(TAX_ID, None)
            self.use_attr_names = self.input_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.gene_id_attribute = self.input_data.attributes.get(GENE_ID_ATTRIBUTE, None)

            if not self.cluster_indicator_model:
                self.Error.no_cluster_indicator()
                return
            elif not self.use_attr_names:
                self.Error.gene_as_attributes()
                return

            self.openContext(self.input_data.domain)

            self.gs_widget.load_gene_sets(self.tax_id)
            if self.cluster_indicator_model and len(self.cluster_indicators) < 1:
                self.cluster_indicators = [self.cluster_indicator_model[0]]
            if self.batch_indicator_model and self.batch_indicator is None:
                self.batch_indicator = self.batch_indicator_model[0]

            self.invalidate()

            if self.custom_data:
                self.refresh_custom_gene_sets()
                self._handle_future_model()
                self.handle_custom_gene_sets()

    @Inputs.custom_sets
    def handle_custom_input(self, data):
        self.Error.clear()
        self.Warning.clear()
        self.closeContext()
        self.custom_data = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None
        self.num_of_custom_sets = None
        self.feature_model.set_domain(None)

        if data:
            self.custom_data = data
            self.feature_model.set_domain(self.custom_data.domain)
            self.custom_tax_id = str(self.custom_data.attributes.get(TAX_ID, None))
            self.custom_use_attr_names = self.custom_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.custom_gene_id_attribute = self.custom_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.custom_gene_id_column = self.custom_data.attributes.get(GENE_ID_COLUMN, None)

            self._handle_future_model()

        if self.input_data:
            self.openContext(self.input_data.domain)

        self.gs_label_combobox.setDisabled(True)
        self.refresh_custom_gene_sets()
        self.handle_custom_gene_sets(select_customs_flag=True)

    def __check_organism_mismatch(self):
        """ Check if organisms from different inputs match.

        :return: True if there is a mismatch
        """
        if self.tax_id is not None and self.custom_tax_id is not None:
            return self.tax_id != self.custom_tax_id
        return False

    def _handle_future_model(self):
        if self.custom_gene_set_indicator in self.feature_model:
            index = self.feature_model.indexOf(self.custom_gene_set_indicator)
            self.custom_gene_set_indicator = self.feature_model[index]
        else:
            if self.feature_model:
                self.custom_gene_set_indicator = self.feature_model[0]
            else:
                self.custom_gene_set_indicator = None

    def handle_custom_gene_sets(self, select_customs_flag=False):
        if self.custom_gene_set_indicator:
            if self.custom_data is not None and self.custom_gene_id_column is not None:

                if self.__check_organism_mismatch():
                    self.gs_label_combobox.setDisabled(True)
                    self.Error.organism_mismatch()
                    self.gs_widget.update_gs_hierarchy()
                    self.__gene_sets_enrichment()
                    return

                if isinstance(self.custom_gene_set_indicator, DiscreteVariable):
                    labels = self.custom_gene_set_indicator.values
                    gene_sets_names = [
                        labels[int(idx)] for idx in self.custom_data.get_column_view(self.custom_gene_set_indicator)[0]
                    ]
                else:
                    gene_sets_names, _ = self.custom_data.get_column_view(self.custom_gene_set_indicator)

                self.num_of_custom_sets = len(set(gene_sets_names))
                gene_names, _ = self.custom_data.get_column_view(self.custom_gene_id_column)
                hierarchy_title = (self.custom_data.name if self.custom_data.name else 'Custom sets',)
                try:
                    self.gs_widget.add_custom_sets(
                        gene_sets_names,
                        gene_names,
                        hierarchy_title=hierarchy_title,
                        select_customs_flag=select_customs_flag,
                    )
                except GeneSetException:
                    pass
                self.gs_label_combobox.setDisabled(False)
            else:
                self.gs_widget.update_gs_hierarchy()

        self.__gene_sets_enrichment()
        self.__update_info_box()

    def refresh_custom_gene_sets(self):
        self.gs_widget.clear_custom_sets()
        # self.gs_widget.update_gs_hierarchy()

    def gene_scores_output(self, selected_clusters):

        metas = [
            StringVariable('Gene'),
            StringVariable(ENTREZ_ID),
            StringVariable('Rank'),
            ContinuousVariable('Statistic score'),
            ContinuousVariable('P-value'),
            ContinuousVariable('FDR'),
        ]

        if len(self.new_cluster_profile):
            # note: order is important
            metas = self.cluster_indicators + metas

        domain = Domain([], metas=metas, class_vars=self.cluster_var)

        data = []
        for cluster in selected_clusters:
            num_of_genes = len(cluster.filtered_genes)

            scores = [gene.score for gene in cluster.filtered_genes]
            p_vals = [gene.p_val for gene in cluster.filtered_genes]
            fdr_vals = [gene.fdr for gene in cluster.filtered_genes]
            gene_names = [gene.input_identifier for gene in cluster.filtered_genes]
            gene_ids = [gene.gene_id for gene in cluster.filtered_genes]
            rank = rankdata(p_vals, method='min')

            if len(self.new_cluster_profile):
                profiles = [[cluster.index] * num_of_genes]
                [profiles.append([p] * num_of_genes) for p in self.new_cluster_profile[cluster.index]]
            else:
                profiles = [[cluster.index] * num_of_genes]

            for row in zip(*profiles, gene_names, gene_ids, rank, scores, p_vals, fdr_vals):
                data.append(list(row))

        out_data = Table(domain, data)
        out_data.attributes[TAX_ID] = self.tax_id
        out_data.attributes[GENE_AS_ATTRIBUTE_NAME] = False
        out_data.attributes[GENE_ID_COLUMN] = ENTREZ_ID
        self.Outputs.gene_scores.send(out_data)

    def gene_set_scores_output(self, selected_clusters):

        metas = [
            StringVariable('Term'),
            StringVariable('Term ID'),
            StringVariable('Rank'),
            ContinuousVariable('P-value'),
            ContinuousVariable('FDR'),
        ]

        if len(self.new_cluster_profile):
            # note: order is important
            metas = self.cluster_indicators + metas

        domain = Domain([], metas=metas, class_vars=self.cluster_var)

        data = []
        for cluster in selected_clusters:
            num_of_sets = len(cluster.filtered_gene_sets)

            p_vals = [gs.p_val for gs in cluster.filtered_gene_sets]
            fdr_vals = [gs.fdr for gs in cluster.filtered_gene_sets]
            gs_names = [gs.name for gs in cluster.filtered_gene_sets]
            gs_ids = [gs.gs_id for gs in cluster.filtered_gene_sets]
            rank = rankdata(p_vals, method='min')

            if len(self.new_cluster_profile):
                profiles = [[cluster.index] * num_of_sets]
                [profiles.append([p] * num_of_sets) for p in self.new_cluster_profile[cluster.index]]
            else:
                profiles = [[cluster.index] * num_of_sets]

            for row in zip(*profiles, gs_names, gs_ids, rank, p_vals, fdr_vals):
                data.append(list(row))

        self.Outputs.gene_set_scores.send(Table(domain, data))

    def commit(self):
        selection_model = self.cluster_info_view.selectionModel()
        selected_rows = selection_model.selectedRows()
        selected_clusters = []
        selected_cluster_indexes = set()
        selected_cluster_genes = set()

        if not self.input_data or not selected_rows:
            self.Outputs.selected_data.send(None)
            return

        for sel_row in selected_rows:
            cluster = sel_row.data()
            selected_clusters.append(cluster)
            selected_cluster_indexes.add(cluster.index)
            [selected_cluster_genes.add(gene.gene_id) for gene in cluster.filtered_genes]

        # get columns of selected clusters
        selected_columns = [
            column
            for column in self.input_data.domain.attributes
            if self.gene_id_attribute in column.attributes
            and str(column.attributes[self.gene_id_attribute]) in selected_cluster_genes
        ]

        domain = Domain(selected_columns, self.input_data.domain.class_vars, self.input_data.domain.metas)
        output_data = self.input_data.from_table(domain, self.input_data)

        # get rows of selected clusters
        selected_rows = [
            row_index
            for row_index, col_index in enumerate(self.rows_by_cluster)
            if col_index in selected_cluster_indexes
        ]

        # send to output signal
        self.Outputs.selected_data.send(output_data[selected_rows])
        self.gene_scores_output(selected_clusters)
        self.gene_set_scores_output(selected_clusters)


if __name__ == "__main__":

    def main(argv=None):
        from AnyQt.QtWidgets import QApplication

        app = QApplication(list(argv) if argv else [])

        w = OWClusterAnalysis()
        data = Table("https://datasets.orange.biolab.si/sc/aml-1k.pickle")
        w.show()
        w.handle_input(data)
        # w.cluster_indicator.initialize(data)
        rval = app.exec_()

        return rval

    sys.exit(main(sys.argv))
