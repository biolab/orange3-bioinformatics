from Orange.data import Table
from Orange.widgets import settings

from orangecontrib.bioinformatics.widgets.components.gene_scoring import GeneScoringComponent


class GeneScoringComponentSettings(settings.ContextHandler):
    """
    Context settings for GeneScoring component. Setting are dependent only on
    column group candidates. See GeneScoringComponent.group_candidates for more details.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def new_context(self, data: Table):
        context = super().new_context()
        context.column_groups, _ = GeneScoringComponent.group_candidates(data)
        return context

    def match(self, context, data):

        column_groups, _ = GeneScoringComponent.group_candidates(data)

        if context.column_groups == column_groups:
            return self.PERFECT_MATCH

        return self.NO_MATCH


class SetContextHandler(settings.ContextHandler):
    def __init__(self, match_imperfect=False):
        super().__init__()
        self.match_imperfect = match_imperfect

    def match(self, context, items):
        items = set(items)

        if self.match_imperfect:
            intersection, union = items & context.items, items | context.items
            if len(union) > 0:
                return len(intersection) / len(union)
            else:
                return 0
        else:
            return 2 if items == context.items else 0

    def new_context(self, items):
        ctx = super().new_context()
        ctx.items = frozenset(items)
        return ctx

    def settings_to_widget(self, widget, *args):
        super().settings_to_widget(widget, *args)

        context = widget.current_context
        if context is None:
            return

        for setting, data, instance in self.provider.traverse_settings(context.values, widget):
            if isinstance(setting, settings.ContextSetting) and setting.name in data:
                value = self.decode_setting(settings, data[setting.name])
                setattr(instance, setting.name, value)


class OrganismContextHandler(settings.ContextHandler):
    def __init__(self):
        super().__init__()

    def match(self, context, taxonomy_id, *args):
        if not context.organism == taxonomy_id:
            return self.NO_MATCH

        return self.PERFECT_MATCH

    def new_context(self, tax_id):
        context = super().new_context()
        context.organism = tax_id
        return context

    def settings_from_widget(self, widget, *args):
        super().settings_from_widget(widget, *args)

        context = widget.current_context
        if context is None:
            return

        # get taxonomy id from the widget
        context.organism = widget.tax_id
