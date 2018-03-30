""" KEGG Pathway """
import sys
import gc
import numpy
import Orange
import webbrowser

from operator import add, itemgetter
from functools import reduce, partial
from contextlib import contextmanager
from collections import defaultdict

from AnyQt.QtWidgets import (
    QTreeWidget, QTreeWidgetItem, QSplitter, QAction, QMenu,
    QGraphicsView, QGraphicsScene, QGraphicsItem, QGraphicsPathItem,
    QGraphicsPixmapItem
)
from AnyQt.QtGui import (
    QBrush, QColor, QPen, QTransform, QPainter, QPainterPath, QPixmap,
    QKeySequence
)
from AnyQt.QtCore import Qt, QRectF, QSize, Slot, QItemSelectionModel

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils import itemmodels, concurrent

from orangecontrib.bioinformatics import kegg
from orangecontrib.bioinformatics import geneset
from orangecontrib.bioinformatics.utils import statistics
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, GENE_AS_ATTRIBUTE_NAME
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher


def relation_list_to_multimap(rellist, ncbi_gene_ids):
    """
    Convert a 'relations' list into a multimap and convert them to ncbi_ids on the fly

    Parameters
    ----------
    rellist : Iterable[Tuple[Hashable[K], Any]]
    ncbi_gene_ids :

    Returns
    -------
    multimap : Dict[K, List[Any]]

    Example
    -------
    relation_list_to_multimap([(1, "a"), (2, "c"), (1, "b")])
    {1: ['a', 'b'], 2: OWKEGGPathwayBrowser['c']}
    """
    mmap = defaultdict(list)
    for key, val in rellist:
        try:
            mmap[key].append(ncbi_gene_ids[val.upper()])
        except KeyError as e:
            # this means that gene id from paths do not exist in kegg gene databse, skip!
            pass
    return dict(mmap)


def split_and_strip(string, sep=None):
    return [s.strip() for s in string.split(sep)]


def path_from_graphics(graphics):
    """
    Return a constructed `QPainterPath` for a KEGG pathway graphics element.
    """
    path = QPainterPath()
    x, y, w, h = [int(graphics.get(c, 0)) for c in
                  ["x", "y", "width", "height"]]
    type = graphics.get("type", "rectangle")
    if type == "rectangle":
        path.addRect(QRectF(x - w / 2, y - h / 2, w, h))
    elif type == "roundrectangle":
        path.addRoundedRect(QRectF(x - w / 2, y - h / 2, w, h), 10, 10)
    elif type == "circle":
        path.addEllipse(QRectF(x - w / 2, y - h / 2, w, h))
    else:
        ValueError("Unknown graphics type %r." % type)
    return path


class EntryGraphicsItem(QGraphicsPathItem):
    """
    An Graphics Item with actions for an overlay of a KEGG pathway image.
    """
    def __init__(self, graphics, *args):
        QGraphicsPathItem.__init__(self, *args)
        path = path_from_graphics(graphics)
        self.setPath(path)
        self.setAcceptHoverEvents(True)
        self._actions = []
        self.link = None

    def hoverEnterEvent(self, event):
        self.setBrush(QBrush(QColor(0, 100, 0, 100)))

    def hoverLeaveEvent(self, event):
        self.setBrush(QBrush(Qt.NoBrush))

    def contextMenuEvent(self, event):
        if self._actions:
            self._menu = menu = QMenu()
            for action in self._actions:
                menu.addAction(action)
            menu.popup(event.screenPos())

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemSelectedHasChanged:
            self.setPen(QPen(Qt.red if self.isSelected() else Qt.blue, 2))

        return QGraphicsPathItem.itemChange(self, change, value)


class GraphicsPathwayItem(QGraphicsPixmapItem):
    """
    A Graphics Item displaying a KEGG Pathway image with optional
    marked objects.

    """
    def __init__(self, pathway, objects, *args, **kwargs):
        QGraphicsPixmapItem.__init__(self, *args)
        self.setTransformationMode(Qt.SmoothTransformation)
        self.setPathway(pathway)
        self.setMarkedObjects(objects,
                              name_mapper=kwargs.get("name_mapper", {}))

    def setPathway(self, pathway):
        """
        Set pathway
        """
        self.pathway = pathway
        if pathway:
            image_filename = pathway.get_image()
            self._pixmap = QPixmap(image_filename)
        else:
            self._pixmap = QPixmap()
        self.setPixmap(self._pixmap)

    def setMarkedObjects(self, objects, name_mapper={}):
        for entry in self.pathway.entries() if self.pathway else []:
            if entry.type == "group":
                continue
            graphics = entry.graphics
            contained_objects = [obj for obj in objects if obj in entry.name]
            item = EntryGraphicsItem(graphics, self)
            item.setToolTip(self.tooltip(entry, contained_objects,
                                         name_mapper))
            item._actions = self.actions(entry, contained_objects)
            item.marked_objects = contained_objects
            if contained_objects:
                item.setPen(QPen(Qt.blue, 2))
                item.setFlag(QGraphicsItem.ItemIsSelectable, True)

    def actions(self, entry, marked_objects=[]):
        actions = []
        type = entry.type
        if marked_objects:
            action = QAction("View genes on kegg website", None)
            org = set([s.split(":")[0] for s in marked_objects]).pop()
            genes = [s.split(":")[-1] for s in marked_objects]
            address = ("http://www.genome.jp/dbget-bin/www_bget?" +
                       "+".join([org] + genes))

            action.triggered.connect(partial(webbrowser.open, address))
            actions.append(action)
        elif hasattr(entry, "link"):
            action = QAction("View %s on KEGG website" % str(type), None)
            action.triggered.connect(partial(webbrowser.open, entry.link))
            actions.append(action)
        return actions

    def tooltip(self, entry, objects, name_mapper={}):
        names = [obj for obj in objects if obj in entry.name]
        names = [name_mapper.get(name, name) for name in names]
        text = entry.name[:16] + " ..." if len(entry.name) > 20 else entry.name
        text = "<p>%s</p>" % text
        if names:
            text += "<br>".join(names)
        return text

    def contextMenuEvent(self, event):
        self._menu = menu = QMenu()
        action = menu.addAction("View this pathway on KEGG website")
        address = ("http://www.kegg.jp/kegg-bin/show_pathway?%s%s" %
                   (self.pathway.org, self.pathway.number))

        action.triggered.connect(partial(webbrowser.open, address))
        menu.popup(event.screenPos())


class PathwayView(QGraphicsView):
    def __init__(self, master, *args):
        QGraphicsView.__init__(self, *args)
        self.master = master

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        self.setRenderHints(QPainter.Antialiasing)
        scene = QGraphicsScene(self)
        self.pixmapGraphicsItem = QGraphicsPixmapItem(None)
        scene.addItem(self.pixmapGraphicsItem)
        self.setScene(scene)

        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)

        self.setFocusPolicy(Qt.WheelFocus)

    def SetPathway(self, pathway=None, objects=[]):
        self.scene().clear()
        self.pathway = pathway
        self.objects = objects
        self.pathwayItem = GraphicsPathwayItem(
            pathway, objects, None,
            name_mapper=getattr(self.master, "uniqueGenesDict", {})
        )

        self.scene().addItem(self.pathwayItem)
        self.scene().setSceneRect(self.pathwayItem.boundingRect())
        self.updateTransform()

    def resizeEvent(self, event):
        self.updateTransform()
        return QGraphicsView.resizeEvent(self, event)

    def updateTransform(self):
        if self.master.autoResize:
            self.fitInView(self.scene().sceneRect().adjusted(-1, -1, 1, 1),
                           Qt.KeepAspectRatio)
        else:
            self.setTransform(QTransform())


class OWKEGGPathwayBrowser(widget.OWWidget):
    name = "KEGG Pathways"
    description = "Browse KEGG pathways that include an input set of genes."
    icon = "../widgets/icons/OWKEGGPathwayBrowser.svg"
    priority = 8

    inputs = [("Data", Orange.data.Table, "SetData", widget.Default),
              ("Reference", Orange.data.Table, "SetRefData")]
    outputs = [("Selected Data", Orange.data.Table, widget.Default),
               ("Unselected Data", Orange.data.Table)]

    settingsHandler = settings.DomainContextHandler()

    organismIndex = settings.ContextSetting(0)
    geneAttrIndex = settings.ContextSetting(0)
    useAttrNames = settings.ContextSetting(False)

    autoCommit = settings.Setting(False)
    autoResize = settings.Setting(True)
    useReference = settings.Setting(False)
    showOrthology = settings.Setting(True)

    Ready, Initializing, Running = 0, 1, 2

    def __init__(self, parent=None):
        super().__init__(parent)

        self.organismCodes = []
        self._changedFlag = False
        self.__invalidated = False
        self.__runstate = OWKEGGPathwayBrowser.Initializing
        self.__in_setProgress = False

        self.controlArea.setMaximumWidth(250)
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.widgetLabel(box, "No data on input\n")

        # Organism selection.
        box = gui.widgetBox(self.controlArea, "Organism")
        self.organismComboBox = gui.comboBox(
            box, self, "organismIndex",
            items=[],
            callback=self.Update,
            addSpace=True,
            tooltip="Select the organism of the input genes")

        # Selection of genes attribute
        box = gui.widgetBox(self.controlArea, "Gene attribute")
        self.geneAttrCandidates = itemmodels.VariableListModel(parent=self)
        self.geneAttrCombo = gui.comboBox(
            box, self, "geneAttrIndex", callback=self.Update)
        self.geneAttrCombo.setModel(self.geneAttrCandidates)

        gui.checkBox(box, self, "useAttrNames",
                    "Use variable names", disables=[(-1, self.geneAttrCombo)],
                    callback=self.Update)

        self.geneAttrCombo.setDisabled(bool(self.useAttrNames))

        gui.separator(self.controlArea)

        gui.checkBox(self.controlArea, self, "useReference",
                    "From signal", box="Reference", callback=self.Update)

        gui.separator(self.controlArea)

        gui.checkBox(self.controlArea, self, "showOrthology",
                     "Show pathways in full orthology", box="Orthology",
                     callback=self.UpdateListView)

        gui.checkBox(self.controlArea, self, "autoResize",
                     "Resize to fit", box="Image",
                     callback=self.UpdatePathwayViewTransform)

        box = gui.widgetBox(self.controlArea, "Cache Control")

        gui.button(box, self, "Clear cache",
                   callback=self.ClearCache,
                   tooltip="Clear all locally cached KEGG data.",
                   default=False, autoDefault=False)

        gui.separator(self.controlArea)

        gui.auto_commit(self.controlArea, self, "autoCommit", "Commit")

        gui.rubber(self.controlArea)

        spliter = QSplitter(Qt.Vertical, self.mainArea)
        self.pathwayView = PathwayView(self, spliter)
        self.pathwayView.scene().selectionChanged.connect(
            self._onSelectionChanged
        )
        self.mainArea.layout().addWidget(spliter)

        self.listView = QTreeWidget(
            allColumnsShowFocus=True,
            selectionMode=QTreeWidget.SingleSelection,
            sortingEnabled=True,
            maximumHeight=200)

        spliter.addWidget(self.listView)

        self.listView.setColumnCount(4)
        self.listView.setHeaderLabels(
            ["Pathway", "P value", "Genes", "Reference"])

        self.listView.itemSelectionChanged.connect(self.UpdatePathwayView)

        select = QAction(
            "Select All", self,
            shortcut=QKeySequence.SelectAll
        )
        select.triggered.connect(self.selectAll)
        self.addAction(select)

        self.data = None
        self.refData = None

        self._executor = concurrent.ThreadExecutor()
        self.setEnabled(False)
        self.setBlocking(True)
        progress = concurrent.methodinvoke(self, "setProgress", (float,))

        def get_genome():
            """Return a KEGGGenome with the common org entries precached."""
            genome = kegg.KEGGGenome()

            essential = genome.essential_organisms()
            common = genome.common_organisms()
            # Remove duplicates of essential from common.
            # (essential + common list as defined here will be used in the
            # GUI.)
            common = [c for c in common if c not in essential]

            # TODO: Add option to specify additional organisms not
            # in the common list.

            keys = list(map(genome.org_code_to_entry_key, essential + common))

            genome.pre_cache(keys, progress_callback=progress)
            return (keys, genome)

        self._genomeTask = task = concurrent.Task(function=get_genome)
        task.finished.connect(self.__initialize_finish)

        self.progressBarInit()
        self.infoLabel.setText("Fetching organism definitions\n")
        self._executor.submit(task)

    def __initialize_finish(self):
        if self.__runstate != OWKEGGPathwayBrowser.Initializing:
            return

        try:
            keys, genome = self._genomeTask.result()
        except Exception as err:
            self.error(0, str(err))
            raise

        self.progressBarFinished()
        self.setEnabled(True)
        self.setBlocking(False)

        entries = [genome[key] for key in keys]
        items = [entry.definition for entry in entries]
        codes = [entry.organism_code for entry in entries]

        self.organismCodes = codes
        self.organismComboBox.clear()
        self.organismComboBox.addItems(items)
        self.organismComboBox.setCurrentIndex(self.organismIndex)

        self.infoLabel.setText("No data on input\n")

    def Clear(self):
        """
        Clear the widget state.
        """
        self.queryGenes = []
        self.referenceGenes = []
        self.genes = {}
        self.uniqueGenesDict = {}
        self.revUniqueGenesDict = {}
        self.pathways = {}
        self.org = None
        self.geneAttrCandidates[:] = []

        self.infoLabel.setText("No data on input\n")
        self.listView.clear()
        self.pathwayView.SetPathway(None)

        self.send("Selected Data", None)
        self.send("Unselected Data", None)

    def SetData(self, data=None):
        if self.__runstate == OWKEGGPathwayBrowser.Initializing:
            self.__initialize_finish()

        self.data = data
        self.warning(0)
        self.error(0)
        self.information(0)

        if data is not None:
            vars = data.domain.variables + data.domain.metas
            vars = [var for var in vars
                    if isinstance(var, Orange.data.StringVariable)]
            self.geneAttrCandidates[:] = vars

            # Try to guess the gene name variable
            if vars:
                names_lower = [v.name.lower() for v in vars]
                scores = [(name == "gene", "gene" in name)
                          for name in names_lower]
                imax, _ = max(enumerate(scores), key=itemgetter(1))
            else:
                imax = -1

            self.geneAttrIndex = imax

            taxid = str(data.attributes.get(TAX_ID, ''))

            if taxid:
                try:
                    code = kegg.from_taxid(taxid)
                    self.organismIndex = self.organismCodes.index(code)
                except Exception as ex:
                    print(ex, taxid)

            self.useAttrNames = data.attributes.get(GENE_AS_ATTRIBUTE_NAME, self.useAttrNames)

            if len(self.geneAttrCandidates) == 0:
                self.useAttrNames = True
                self.geneAttrIndex = -1
            else:
                self.geneAttrIndex = min(self.geneAttrIndex,
                                         len(self.geneAttrCandidates) - 1)
        else:
            self.Clear()

        self.__invalidated = True

    def SetRefData(self, data=None):
        self.refData = data
        self.information(1)

        if data is not None and self.useReference:
            self.__invalidated = True

    def handleNewSignals(self):
        if self.__invalidated:
            self.Update()
            self.__invalidated = False

    def UpdateListView(self):
        self.bestPValueItem = None
        self.listView.clear()
        if not self.data:
            return

        allPathways = self.org.pathways()
        allRefPathways = kegg.pathways("map")

        items = []
        kegg_pathways = kegg.KEGGPathways()

        org_code = self.organismCodes[min(self.organismIndex,
                                          len(self.organismCodes) - 1)]

        if self.showOrthology:
            self.koOrthology = kegg.KEGGBrite("ko00001")
            self.listView.setRootIsDecorated(True)
            path_ids = set([s[-5:] for s in self.pathways.keys()])

            def _walkCollect(koEntry):
                num = koEntry.title[:5] if koEntry.title else None
                if num in path_ids:
                    return ([koEntry] +
                            reduce(lambda li, c: li + _walkCollect(c),
                                   [child for child in koEntry.entries],
                                   []))
                else:
                    c = reduce(lambda li, c: li + _walkCollect(c),
                               [child for child in koEntry.entries],
                               [])
                    return c + (c and [koEntry] or [])

            allClasses = reduce(lambda li1, li2: li1 + li2,
                                [_walkCollect(c) for c in self.koOrthology],
                                [])

            def _walkCreate(koEntry, lvItem):
                item = QTreeWidgetItem(lvItem)
                id = "path:" + org_code + koEntry.title[:5]

                if koEntry.title[:5] in path_ids:
                    p = kegg_pathways.get_entry(id)
                    if p is None:
                        # In case the genesets still have obsolete entries
                        name = koEntry.title
                    else:
                        name = p.name
                    genes, p_value, ref = self.pathways[id]
                    item.setText(0, name)
                    item.setText(1, "%.5f" % p_value)
                    item.setText(2, "%i of %i" % (len(genes), len(self.genes)))
                    item.setText(3, "%i of %i" % (ref, len(self.referenceGenes)))
                    item.pathway_id = id if p is not None else None
                else:
                    if id in allPathways:
                        text = kegg_pathways.get_entry(id).name
                    else:
                        text = koEntry.title
                    item.setText(0, text)

                    if id in allPathways:
                        item.pathway_id = id
                    elif "path:map" + koEntry.title[:5] in allRefPathways:
                        item.pathway_id = "path:map" + koEntry.title[:5]
                    else:
                        item.pathway_id = None

                for child in koEntry.entries:
                    if child in allClasses:
                        _walkCreate(child, item)

            for koEntry in self.koOrthology:
                if koEntry in allClasses:
                    _walkCreate(koEntry, self.listView)

            self.listView.update()
        else:
            self.listView.setRootIsDecorated(False)
            pathways = self.pathways.items()
            pathways = sorted(pathways, key=lambda item: item[1][1])

            for id, (genes, p_value, ref) in pathways:
                item = QTreeWidgetItem(self.listView)
                item.setText(0, kegg_pathways.get_entry(id).name)
                item.setText(1, "%.5f" % p_value)
                item.setText(2, "%i of %i" % (len(genes), len(self.genes)))
                item.setText(3, "%i of %i" % (ref, len(self.referenceGenes)))
                item.pathway_id = id
                items.append(item)

        self.bestPValueItem = items and items[0] or None
        self.listView.expandAll()
        for i in range(4):
            self.listView.resizeColumnToContents(i)

        if self.bestPValueItem:
            index = self.listView.indexFromItem(self.bestPValueItem)
            self.listView.selectionModel().select(
                index, QItemSelectionModel.ClearAndSelect
            )

    def UpdatePathwayView(self):
        items = self.listView.selectedItems()

        if len(items) > 0:
            item = items[0]
        else:
            item = None

        self.commit()
        item = item or self.bestPValueItem
        if not item or not item.pathway_id:
            self.pathwayView.SetPathway(None)
            return

        def get_kgml_and_image(pathway_id):
            """Return an initialized KEGGPathway with pre-cached data"""
            p = kegg.KEGGPathway(pathway_id)
            p._get_kgml()  # makes sure the kgml file is downloaded
            p._get_image_filename()  # makes sure the image is downloaded
            return (pathway_id, p)

        self.setEnabled(False)
        self._pathwayTask = concurrent.Task(
            function=lambda: get_kgml_and_image(item.pathway_id)
        )
        self._pathwayTask.finished.connect(self._onPathwayTaskFinshed)
        self._executor.submit(self._pathwayTask)

    def _onPathwayTaskFinshed(self):
        self.setEnabled(True)
        pathway_id, self.pathway = self._pathwayTask.result()
        self.pathwayView.SetPathway(
            self.pathway,
            self.pathways.get(pathway_id, [[]])[0]
        )

    def UpdatePathwayViewTransform(self):
        self.pathwayView.updateTransform()

    def Update(self):
        """
        Update (recompute enriched pathways) the widget state.
        """
        if not self.data:
            return

        self.error(0)
        self.information(0)

        # XXX: Check data in setData, do not even allow this to be executed if
        # data has no genes
        try:
            genes = self.GeneNamesFromData(self.data)
        except ValueError:
            self.error(0, "Cannot extract gene names from input.")
            genes = []

        if not self.useAttrNames and any("," in gene for gene in genes):
            genes = reduce(add, (split_and_strip(gene, ",")
                                 for gene in genes),
                           [])
            self.information(0,
                             "Separators detected in input gene names. "
                             "Assuming multiple genes per instance.")

        self.queryGenes = genes

        self.information(1)
        reference = None
        if self.useReference and self.refData:
            reference = self.GeneNamesFromData(self.refData)
            if not self.useAttrNames \
                    and any("," in gene for gene in reference):
                reference = reduce(add, (split_and_strip(gene, ",")
                                         for gene in reference),
                                   [])
                self.information(1,
                                 "Separators detected in reference gene "
                                 "names. Assuming multiple genes per "
                                 "instance.")

        org_code = self.SelectedOrganismCode()

        gm = GeneMatcher(kegg.to_taxid(org_code))
        gm.genes = genes
        gm.run_matcher()
        mapped_genes = {gene: str(ncbi_id) for gene, ncbi_id in gm.map_input_to_ncbi().items()}

        def run_enrichment(org_code, genes, reference=None, progress=None):
            org = kegg.KEGGOrganism(org_code)
            if reference is None:
                reference = org.get_ncbi_ids()

            # This is here just to keep widget working without any major changes.
            # map not needed, geneMatcher will not work on widget level.
            unique_genes = genes
            unique_ref_genes = dict([(gene, gene) for gene in set(reference)])

            taxid = kegg.to_taxid(org.org_code)
            # Map the taxid back to standard 'common' taxids
            # (as used by 'geneset') if applicable
            r_tax_map = dict((v, k) for k, v in
                             kegg.KEGGGenome.TAXID_MAP.items())
            if taxid in r_tax_map:
                taxid = r_tax_map[taxid]

            # We use the kegg pathway gene sets provided by 'geneset' for
            # the enrichment calculation.

            kegg_api = kegg.api.CachedKeggApi()
            linkmap = kegg_api.link(org.org_code, "pathway")
            converted_ids = kegg_api.conv(org.org_code, 'ncbi-geneid')
            kegg_sets = relation_list_to_multimap(linkmap, dict((gene.upper(), ncbi.split(':')[-1])
                                                                for ncbi, gene in converted_ids))

            kegg_sets = geneset.GeneSets(sets=[geneset.GeneSet(gs_id=ddi, genes=set(genes))
                                               for ddi, genes in kegg_sets.items()])

            pathways = pathway_enrichment(
                kegg_sets, unique_genes.values(),
                unique_ref_genes.keys(),
                callback=progress
            )
            # Ensure that pathway entries are pre-cached for later use in the
            # list/tree view
            kegg_pathways = kegg.KEGGPathways()
            kegg_pathways.pre_cache(
                pathways.keys(), progress_callback=progress
            )

            return pathways, org, unique_genes, unique_ref_genes

        self.progressBarInit()
        self.setEnabled(False)
        self.infoLabel.setText("Retrieving...\n")

        progress = concurrent.methodinvoke(self, "setProgress", (float,))

        self._enrichTask = concurrent.Task(
            function=lambda:
                run_enrichment(org_code, mapped_genes, reference, progress)
        )
        self._enrichTask.finished.connect(self._onEnrichTaskFinished)
        self._executor.submit(self._enrichTask)

    def _onEnrichTaskFinished(self):
        self.setEnabled(True)
        self.setBlocking(False)
        try:
            pathways, org, unique_genes, unique_ref_genes = \
                self._enrichTask.result()
        except Exception:
            raise

        self.progressBarFinished()

        self.org = org
        self.genes = unique_genes.keys()
        self.uniqueGenesDict = {ncbi_id: input_name for input_name, ncbi_id in unique_genes.items()}
        self.revUniqueGenesDict = dict([(val, key) for key, val in
                                        self.uniqueGenesDict.items()])
        self.referenceGenes = unique_ref_genes.keys()
        self.pathways = pathways

        if not self.pathways:
            self.warning(0, "No enriched pathways found.")
        else:
            self.warning(0)

        count = len(set(self.queryGenes))
        self.infoLabel.setText(
            "%i unique gene names on input\n"
            "%i (%.1f%%) genes names matched" %
            (count, len(unique_genes),
             100.0 * len(unique_genes) / count if count else 0.0)
        )

        self.UpdateListView()

    @Slot(float)
    def setProgress(self, value):
        if self.__in_setProgress:
            return

        self.__in_setProgress = True
        self.progressBarSet(value)
        self.__in_setProgress = False

    def GeneNamesFromData(self, data):
        """
        Extract and return gene names from `data`.
        """
        if self.useAttrNames:
            genes = [str(v.name).strip() for v in data.domain.attributes]
        elif self.geneAttrCandidates:
            assert 0 <= self.geneAttrIndex < len(self.geneAttrCandidates)
            geneAttr = self.geneAttrCandidates[self.geneAttrIndex]
            genes = [str(e[geneAttr]) for e in data
                     if not numpy.isnan(e[geneAttr])]
        else:
            raise ValueError("No gene names in data.")
        return genes

    def SelectedOrganismCode(self):
        """
        Return the selected organism code.
        """
        return self.organismCodes[min(self.organismIndex,
                                      len(self.organismCodes) - 1)]

    def selectAll(self):
        """
        Select all items in the pathway view.
        """
        changed = False
        scene = self.pathwayView.scene()
        with disconnected(scene.selectionChanged, self._onSelectionChanged):
            for item in scene.items():
                if item.flags() & QGraphicsItem.ItemIsSelectable and \
                        not item.isSelected():
                    item.setSelected(True)
                    changed = True
        if changed:
            self._onSelectionChanged()

    def _onSelectionChanged(self):
        # Item selection in the pathwayView/scene has changed
        self.commit()

    def commit(self):
        if self.data:
            selectedItems = self.pathwayView.scene().selectedItems()
            selectedGenes = reduce(set.union, [item.marked_objects
                                               for item in selectedItems],
                                   set())

            if self.useAttrNames:
                selected = [self.data.domain[self.uniqueGenesDict[gene]]
                            for gene in selectedGenes]
#                 newDomain = Orange.data.Domain(selectedVars, 0)
                data = self.data[:, selected]
#                 data = Orange.data.Table(newDomain, self.data)
                self.send("Selected Data", data)
            elif self.geneAttrCandidates:
                assert 0 <= self.geneAttrIndex < len(self.geneAttrCandidates)
                geneAttr = self.geneAttrCandidates[self.geneAttrIndex]
                selectedIndices = []
                otherIndices = []
                for i, ex in enumerate(self.data):
                    names = [self.revUniqueGenesDict.get(name, None)
                             for name in split_and_strip(str(ex[geneAttr]), ",")]
                    if any(name and name in selectedGenes for name in names):
                        selectedIndices.append(i)
                    else:
                        otherIndices.append(i)

                if selectedIndices:
                    selected = self.data[selectedIndices]
                else:
                    selected = None

                if otherIndices:
                    other = self.data[otherIndices]
                else:
                    other = None

                self.send("Selected Data", selected)
                self.send("Unselected Data", other)
        else:
            self.send("Selected Data", None)
            self.send("Unselected Data", None)

    def ClearCache(self):
        kegg.caching.clear_cache()

    def onDeleteWidget(self):
        """
        Called before the widget is removed from the canvas.
        """
        super().onDeleteWidget()

        self.org = None
        self._executor.shutdown(wait=False)
        gc.collect()  # Force collection (WHY?)

    def sizeHint(self):
        return QSize(1024, 720)


def pathway_enrichment(genesets, genes, reference, prob=None, callback=None):
    result_sets = []
    p_values = []
    if prob is None:
        prob = statistics.Hypergeometric()

    for i, gs in enumerate(genesets):
        cluster = gs.genes.intersection(genes)
        ref = gs.genes.intersection(reference)
        k = len(cluster)
        N = len(reference)
        m = len(ref)
        n = len(genes)
        if k:
            p_val = prob.p_value(k, N, m, n)
            result_sets.append((gs.gs_id, cluster, ref))
            p_values.append(p_val)
        if callback is not None:
            callback(100.0 * i / len(genesets))

    # FDR correction
    p_values = statistics.FDR(p_values)

    return dict([(id, (genes, p_val, len(ref)))
                 for (id, genes, ref), p_val in zip(result_sets, p_values)])


@contextmanager
def disconnected(signal, slot):
    """
    A context manager disconnecting a slot from a signal.
    ::

        with disconnected(scene.selectionChanged, self.onSelectionChanged):
            # Can change item selection in a scene without
            # onSelectionChanged being invoked.
            do_something()

    """
    signal.disconnect(slot)
    try:
        yield
    finally:
        signal.connect(slot)


def test_main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(sys.argv)

    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = "brown-selected"
    data = Orange.data.Table(filename)

    w = OWKEGGPathwayBrowser()
    w.show()
    w.raise_()

    w.SetData(data)
    w.handleNewSignals()

    r = app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
    return r


if __name__ == "__main__":
    sys.exit(test_main())
