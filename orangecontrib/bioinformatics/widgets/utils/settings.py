from Orange.widgets import settings


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

    def settings_to_widget(self, widget):
        super().settings_to_widget(widget)

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
