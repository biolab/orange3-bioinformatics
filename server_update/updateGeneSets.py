""" update gene_sets """
import server_update

from server_update import SOURCE_SERVER
from server_update.tests.test_GeneSets import GeneSetsTest


import os
import io
from zipfile import ZipFile

from urllib.request import urlopen
from server_update import create_folder, sf_local, create_info_file
from orangecontrib.bioinformatics import go, kegg, omim
from orangecontrib.bioinformatics.dicty import phenotypes
from orangecontrib.bioinformatics.kegg import caching
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher
from orangecontrib.bioinformatics.geneset import (
    DOMAIN, GeneSet, GeneSets, GeneSetException,
    filename, GO_TERM_LINK, CYTOBAND_DOWNLOAD_LINK,
    REACTOME_DOWNLOAD_LINK, REACTOME_FILE_NAME, OMIM_LINK
)

DOMAIN_PATH = sf_local.localpath(DOMAIN)
create_folder(DOMAIN_PATH)


def go_gene_sets(org):
    """ Returns gene sets from GO.
    """

    ontology = go.Ontology()
    annotations = go.Annotations(org, ontology=ontology)

    gene_sets = []
    for termn, term in ontology.terms.items():
        genes = annotations.get_genes_by_go_term(termn)
        hier = ('GO', term.namespace)
        if len(genes) > 0:

            gs = GeneSet(gs_id=termn, name=term.name, genes=genes, hierarchy=hier,
                         organism=org, link=GO_TERM_LINK.format(termn))

            gene_sets.append(gs)

    return GeneSets(gene_sets)


def kegg_gene_sets(org):
    """ Returns gene sets from KEGG pathways.
    """
    caching.clear_cache()
    kegg_org = kegg.KEGGOrganism(taxonomy.name(org))
    ncbi_id_mapper = kegg_org.kegg_to_ncbi_mapper()
    genesets = []

    for id in kegg_org.pathways():
        pway = kegg.KEGGPathway(id)
        hier = ('KEGG', 'Pathways')

        if pway.pathway_attributes():
            kegg_names = kegg_org.get_genes_by_pathway(id)
            mapped_genes = []
            for gene in kegg_names:
                try:
                    mapped_genes.append(ncbi_id_mapper[gene.upper()])
                except KeyError:
                    # some kegg names can not be matched to ncbi ids
                    # they are included in geneset anyway
                    # remove prefix, that specifies kegg organism
                    # mapped_genes.append(gene.split(':')[-1])
                    pass

            gs = GeneSet(gs_id=id,
                         name=pway.title,
                         genes=mapped_genes,
                         hierarchy=hier,
                         organism=org,
                         link=pway.link)
            genesets.append(gs)

    return GeneSets(genesets)


def dicty_mutant_gene_sets(org):
    """ Return dicty mutant phenotype gene sets from Dictybase
    """
    if org == '352472':
        gene_sets = []
        gene_matcher = GeneMatcher('352472')

        for phenotype, mutants in phenotypes.phenotype_mutants().items():
            phenotype = phenotype.replace(",", " ")
            gene_symbols = [phenotypes.mutant_genes(mutant)[0] for mutant in mutants]
            gene_matcher.genes = gene_symbols
            gene_matcher.run_matcher()
            genes = []

            for gene in gene_matcher.genes:
                if gene.ncbi_id is not None:
                    genes.append(int(gene.ncbi_id))

            if len(gene_symbols) != len(genes):
                print(len(gene_symbols), len(genes))

            gs = GeneSet(gs_id=phenotype,
                         name=phenotype,
                         genes=genes,
                         hierarchy=('Dictybase', 'Phenotypes'),
                         organism='352472',
                         link='')

            gene_sets.append(gs)

        return GeneSets(gene_sets)


def cytoband_gene_sets(org):
    """ Create cytoband gene sets from Stanford Microarray Database
    """
    if org == '9606':
        gene_matcher = GeneMatcher('9606')

        with urlopen(CYTOBAND_DOWNLOAD_LINK) as stream:
            data = stream.read().splitlines()
            genesets = []

            for band in data:
                b = band.decode().split('\t')
                gene_symbols = b[2:]
                gene_matcher.genes = gene_symbols
                gene_matcher.run_matcher()

                genes = []
                for gene in gene_matcher.genes:
                    if gene.ncbi_id is not None:
                        genes.append(int(gene.ncbi_id))

                genesets.append(GeneSet(gs_id=b[0], name=b[1], genes=genes if b[2:] else [],
                                        hierarchy=('Cytobands',), organism='9606', link=''))

            return GeneSets(genesets)


def reactome_gene_sets(org):
    """ Prepare human pathways gene sets from reactome pathways
    """
    if org == '9606':
        gene_matcher = GeneMatcher('9606')

        with urlopen(REACTOME_DOWNLOAD_LINK) as url:
            memfile = io.BytesIO(url.read())

            with ZipFile(memfile, 'r') as myzip:
                f = myzip.open(REACTOME_FILE_NAME)
                content = f.read().decode().splitlines()
                genesets = []

                for path in content:
                    gene_symbols = path.split('\t')[2:] if path.split('\t')[2:] else []
                    gene_matcher.genes = gene_symbols
                    gene_matcher.run_matcher()
                    genes = []

                    for gene in gene_matcher.genes:
                        if gene.ncbi_id is not None:
                            genes.append(int(gene.ncbi_id))
                    pathway = path.split('\t')[0].replace(',', ' ')
                    gs = GeneSet(gs_id=pathway,
                                 name=pathway,
                                 genes=genes,
                                 hierarchy=('Reactome', 'pathways'),
                                 organism='9606', link='')

                    genesets.append(gs)

                return GeneSets(genesets)


def omim_gene_sets(org):
    """ Return gene sets from OMIM (Online Mendelian Inheritance in Man) diseses
    """
    if org == '9606':
        gene_matcher = GeneMatcher('9606')
        genesets = []

        for disease in omim.diseases():
            gene_symbols = omim.disease_genes(disease)
            gene_matcher.genes = gene_symbols
            gene_matcher.run_matcher()
            genes = []

            for gene in gene_matcher.genes:
                if gene.ncbi_id is not None:
                    genes.append(int(gene.ncbi_id))

            gs = GeneSet(gs_id=disease.id,
                         name=disease.name,
                         genes=genes,
                         hierarchy=('OMIM',),
                         organism='9606',
                         link=(OMIM_LINK.format(disease.id) if disease.id else None))
            genesets.append(gs)

        return GeneSets(genesets)


def register_serverfiles(genesets):
    """ Registers using the common hierarchy and organism. """
    org = genesets.common_org()
    hierarchy = genesets.common_hierarchy()
    fn = filename(hierarchy, org)

    file_path = os.path.join(DOMAIN_PATH, fn)

    if org is not None:
        taxname = taxonomy.name(org)
        title = "Gene sets: " + ", ".join(hierarchy) + ((" (" + taxname + ")") if org is not None else "")

        tags = list(hierarchy) + ["gene sets"] + ([taxname] if org is not None else []) + taxonomy.shortname(org)
        genesets.to_gmt_file_format(file_path)
        create_info_file(file_path,
                         domain=DOMAIN,
                         filename=fn,
                         source=SOURCE_SERVER,
                         title=title,
                         tags=tags)


def upload_genesets():
    """ Builds the default gene sets and
    """

    genesetsfn = [go_gene_sets,
                  kegg_gene_sets,
                  # omim_gene_sets,  # stop supporting OMIM. Did not update files since 2011
                  cytoband_gene_sets,
                  reactome_gene_sets,
                  dicty_mutant_gene_sets
                  ]

    organisms = taxonomy.common_taxids()
    for fn in genesetsfn:
        for org in organisms:
            try:
                # print("Uploading ORG {} {}".format(org, fn))
                try:
                    genesets = fn(org).split_by_hierarchy()
                except AttributeError as e:
                    # print(e)
                    # genesets = fn().split_by_hierarchy() print(e)
                    continue

                for gs in genesets:
                    # print("registering {}".format(str(gs.common_hierarchy())))
                    register_serverfiles(gs)  # server files register(gs)
            except taxonomy.UnknownSpeciesIdentifier:
                print("Organism ontology not available %s" % org)
            except GeneSetException:
                print("Empty gene sets. %s" % org)


upload_genesets()

helper = server_update.SyncHelper(DOMAIN, GeneSetsTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
