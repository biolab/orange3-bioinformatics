from typing import Dict, List, Optional
from functools import reduce
from collections import defaultdict

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi.gene import Gene


def _from_data_to_gene(line: str):
    gene = Gene()
    for attr, val in zip(['homology_group_id', 'tax_id', 'gene_id'], line.strip().split('\t')):
        setattr(gene, attr, val)
    return gene


class HomoloGene:
    """ Wrapper around NCBI HomoloGene database """

    def __init__(self):
        self.file_path: str = serverfiles.localpath_download('homologene', 'homologene.tab')

        with open(self.file_path, 'r') as fp:
            self._homologs: Dict[str, Gene] = {
                h.gene_id: h for h in [_from_data_to_gene(line) for line in fp.readlines()]
            }

        def _helper(groups, gene):
            groups[gene.homology_group_id].append(gene)
            return groups

        self._homologs_by_group: Dict[str, List[Gene]] = reduce(_helper, self._homologs.values(), defaultdict(list))

    def find_homolog(self, gene_id: str, organism: str) -> Optional[str]:
        """ Find homolog gene in organism. If the homolog does not exist, return None. """
        homology_group = self._homologs.get(str(gene_id), Gene()).homology_group_id
        homologs = [
            gene.gene_id for gene in self._homologs_by_group.get(homology_group, []) if gene.tax_id == organism
        ]
        if len(homologs) == 1:
            return homologs[0]
        else:
            # Is possible that find more then one gene?
            return None


if __name__ == "__main__":
    from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, load_gene_summary
    import Orange

    homology = HomoloGene()

    gm = GeneMatcher('4932')
    genes = Orange.data.Table("brown-selected")

    gm.genes = genes
    _homologs = [homology.find_homolog(str(gene.gene_id), '9606') for gene in gm.genes]
    _homologs = load_gene_summary('9606', _homologs)

    for gene, homolog in zip(gm.genes, _homologs):
        print(f'{gene} ----> {homolog}')
