""" Utils """
from Orange.data import ContinuousVariable, StringVariable, TimeVariable, Domain, Table


def transpose_table(table):
    """ Transpose the rows and columns of the table.

    Args:
        table: Data in :obj:`Orange.data.Table`

    Returns:
         Transposed :obj:`Orange.data.Table`. (Genes as columns)
    """

    # TODO: remove this and use Orange.data.Table.transpose

    attrs = table.domain.attributes
    attr = [ContinuousVariable.make(ex['Gene'].value) for ex in table]
    #  Set metas
    new_metas = [StringVariable.make(name) if name is not 'Time' else TimeVariable.make(name)
                 for name in sorted(table.domain.variables[0].attributes.keys())]
    domain = Domain(attr, metas=new_metas)
    meta_values = [[exp.attributes[var.name] for var in domain.metas] for exp in attrs]

    return Table(domain, table.X.transpose(), metas=meta_values)


def etc_to_table(etc_json, time_var=False):
    """ Converts data from Json to :obj:`Orange.data.table`

    Args:
        etc_json (dict): Data in json like format from genesis
        time_var (bool): Create column of time points. Default is set to False.

    Returns:
        :obj:`Orange.data.Table`
    """

    variables = []
    time_point = 1
    for time in etc_json['etc']['timePoints']:
        var = ContinuousVariable('TP ' + str(time_point))
        var.attributes['Time'] = str(time)
        variables.append(var)
        time_point += 1

    meta_attr = StringVariable.make('Gene')
    domain = Domain(variables, metas=[meta_attr])

    table = []
    for row in etc_json['etc']['genes']:
        gene_expression = [exp for exp in etc_json['etc']['genes'][row]]
        gene_expression.append(row)
        table.append(gene_expression)

    orange_table = Table(domain, table)

    if time_var:
        orange_table = transpose_table(orange_table)

    return orange_table
