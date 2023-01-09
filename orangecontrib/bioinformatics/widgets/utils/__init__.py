import itertools

from orangewidget.utils.signals import PartialSummary, summarize

from orangecontrib.bioinformatics.geneset import GeneSets


@summarize.register
def summarize_(gene_sets: GeneSets):  # pylint: disable=function-redefined
    n = len(gene_sets)
    if n == 0:
        details = 'empty gene sets'
    elif n <= 3:
        details = ', '.join(gs.name for gs in gene_sets)
    else:
        details = (
            ', '.join(gs.name for gs in itertools.islice(gene_sets, 3))
            + f' and {n - 3} others'
        )

    return PartialSummary(n, details)
