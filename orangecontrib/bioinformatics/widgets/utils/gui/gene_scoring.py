""" Qt component for Gene scoring method selection """
from collections import namedtuple
from typing import Union

from AnyQt.QtWidgets import (
    QWidget
)

from Orange.widgets.gui import (
    radioButtons, QGroupBox, comboBox
)


from orangecontrib.bioinformatics.utils.statistics import score_t_test, score_mann_whitney, hypergeometric_test_vector


#: Selection types
LowTail, HighTail, TwoTail = 1, 2, 3
#: Test type - i.e a two sample (t-test, ...) or multi-sample (ANOVA) test
TwoSampleTest, VarSampleTest = 1, 2

gene_scoring_method = namedtuple('scoring_method', 'name, score_function, selection_type, test_type')


# TODO: Use this component in OWDifferentialExpression

class GeneScoringWidget(QWidget):

    # default methods, override to add custom scoring methods
    scores = [gene_scoring_method('T-test', score_t_test, TwoTail, TwoSampleTest),
              gene_scoring_method('Mann-Whitney', score_mann_whitney, LowTail, TwoSampleTest),
              gene_scoring_method('Hypergeometric Test', hypergeometric_test_vector, TwoTail, TwoSampleTest)]

    def __init__(self, box, parent,  **kwargs):
        # type: (Union[QGroupBox, QWidget], QWidget) -> None
        super().__init__(**kwargs)

        self.widget = box
        self.parent = parent

        # parent widget settings
        self.scoring_method_selection = None
        self.scoring_method_design = None

    def set_method_selection_area(self, settings_var):
        # type: (str) -> None
        self.scoring_method_selection = settings_var

        comboBox(self.widget, self.parent, self.scoring_method_selection,
                 items=[sm.name for sm in self.scores],
                 callback=self.on_method_selection_changed,
                 label='Method')

    def set_method_design_area(self, settings_var):
        # type: (str) -> None
        self.scoring_method_design = settings_var

        radioButtons(self.widget, self.parent, self.scoring_method_design,
                     ['Cluster vs. rest', 'Cluster vs. Cluster (max)'],
                     callback=self.on_design_selection_changed,
                     label='Design')

    def get_selected_method(self):
        # type: () -> gene_scoring_method
        return self.scores[self.parent.__getattribute__(self.scoring_method_selection)]

    def get_selected_desig(self):
        # type: () -> int
        return self.parent.__getattribute__(self.scoring_method_design)

    def on_method_selection_changed(self):
        """ Override this method if needed """
        self.parent.invalidate()

    def on_design_selection_changed(self):
        """ Override this method if needed """
        self.parent.invalidate()
