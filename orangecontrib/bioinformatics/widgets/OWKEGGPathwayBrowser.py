""" KEGG Pathway """
import gc
import sys
import webbrowser
from typing import Dict, List, Tuple, TypeVar, Iterable, Sequence
from functools import reduce, partial
from itertools import chain
from contextlib import contextmanager
from collections import defaultdict

from AnyQt.QtGui import QPen, QBrush, QColor, QPixmap, QPainter, QTransform, QKeySequence, QPainterPath
from AnyQt.QtCore import Qt, Slot, QSize, QRectF, QItemSelectionModel
from AnyQt.QtWidgets import (
    QMenu,
    QAction,
    QSplitter,
    QTreeWidget,
    QGraphicsItem,
    QGraphicsView,
    QGraphicsScene,
    QTreeWidgetItem,
    QGraphicsPathItem,
    QGraphicsPixmapItem,
)

import Orange
from Orange.widgets import gui, widget, settings
from Orange.widgets.utils import concurrent

from orangecontrib.bioinformatics import kegg, geneset
from orangecontrib.bioinformatics.utils import statistics
from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID,
    GENE_ID_COLUMN,
    GENE_ID_ATTRIBUTE,
    GENE_AS_ATTRIBUTE_NAME,
    ERROR_ON_MISSING_TAX_ID,
    ERROR_ON_MISSING_GENE_ID,
    ERROR_ON_MISSING_ANNOTATION,
)

A = TypeVar("A")
B = TypeVar("B")
C = TypeVar("C")
D = TypeVar("D")


def relation_list_to_multimap(rellist: Sequence[Tuple[A, B]]) -> Dict[A, Sequence[B]]:
    """
    Convert a 'relations' list into a multimap

    Parameters
    ----------
    rellist : Iterable[Tuple[Hashable[K], Any]]

    Returns
    -------
    multimap : Dict[K, List[Any]]

    Example
    -------
    >>> relation_list_to_multimap([(1, "a"), (2, "c"), (1, "b")])
    {1: ['a', 'b'], 2: ['c']}
    """
    mmap = defaultdict(list)
    for key, val in rellist:
        mmap[key].append(val)
    return dict(mmap)


def relation_join(rela: Sequence[Tuple[A, B]], relb: Sequence[Tuple[C, D]]) -> Sequence[Tuple[A, D]]:
    res = []
    relb_mmap = relation_list_to_multimap(relb)
    for a, b in rela:
        res.extend((a, d) for d in relb_mmap.get(b, ()))
    return res


def relation_map(relation: Sequence[Tuple[A, B]], values: Sequence[A]) -> Sequence[Tuple[A, Sequence[B]]]:
    rel_mmap = relation_list_to_multimap(relation)
    return [rel_mmap.get(v, []) for v in values]


def flatten(seq: Iterable[Iterable[A]]) -> Sequence[A]:
    return list(chain.from_iterable(seq))


def split_and_strip(string, sep=None):
    return [s.strip() for s in string.split(sep)]


def path_from_graphics(graphics):
    """
    Return a constructed `QPainterPath` for a KEGG pathway graphics element.
    """
    path = QPainterPath()
    x, y, w, h = [int(graphics.get(c, 0)) for c in ["x", "y", "width", "height"]]
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
        self.setMarkedObjects(objects)

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

    def setMarkedObjects(self, objects):
        for entry in self.pathway.entries() if self.pathway else []:
            if entry.type == "group":
                continue
            graphics = entry.graphics
            contained_objects = [obj for obj in objects if obj in entry.name]
            item = EntryGraphicsItem(graphics, self)
            item.setToolTip(self.tooltip(entry, contained_objects))
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
            org = {s.split(":")[0] for s in marked_objects}.pop()
            genes = [s.split(":")[-1] for s in marked_objects]
            address = "http://www.genome.jp/dbget-bin/www_bget?" + "+".join([org] + genes)

            action.triggered.connect(partial(webbrowser.open, address))
            actions.append(action)
        elif hasattr(entry, "link"):
            action = QAction("View %s on KEGG website" % str(type), None)
            action.triggered.connect(partial(webbrowser.open, entry.link))
            actions.append(action)
        return actions

    def tooltip(self, entry, objects):
        names = [obj for obj in objects if obj in entry.name]
        text = entry.name[:16] + " ..." if len(entry.name) > 20 else entry.name
        text = "<p>%s</p>" % text
        if names:
            text += "<br>".join(names)
        return text

    def contextMenuEvent(self, event):
        self._menu = menu = QMenu()
        action = menu.addAction("View this pathway on KEGG website")
        address = "http://www.kegg.jp/kegg-bin/show_pathway?%s%s" % (self.pathway.org, self.pathway.number)

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
        self.pathwayItem = GraphicsPathwayItem(pathway, objects, None)

        self.scene().addItem(self.pathwayItem)
        self.scene().setSceneRect(self.pathwayItem.boundingRect())
        self.updateTransform()

    def resizeEvent(self, event):
        self.updateTransform()
        return QGraphicsView.resizeEvent(self, event)

    def updateTransform(self):
        if self.master.autoResize:
            self.fitInView(self.scene().sceneRect().adjusted(-1, -1, 1, 1), Qt.KeepAspectRatio)
        else:
            self.setTransform(QTransform())


class OWKEGGPathwayBrowser(widget.OWWidget):
    name = "KEGG Pathways"
    description = "Browse KEGG pathways that include an input set of genes."
    icon = "../widgets/icons/OWKEGGPathwayBrowser.svg"
    priority = 70

    inputs = [("Data", Orange.data.Table, "SetData", widget.Default), ("Reference", Orange.data.Table, "SetRefData")]
    outputs = [("Selected Data", Orange.data.Table, widget.Default), ("Unselected Data", Orange.data.Table)]

    autoCommit = settings.Setting(False)
    autoResize = settings.Setting(True)
    useReference = settings.Setting(False)
    showOrthology = settings.Setting(True)

    Ready, Initializing, Running = 0, 1, 2

    class Error(widget.OWWidget.Error):
        missing_annotation = widget.Msg(ERROR_ON_MISSING_ANNOTATION)
        missing_gene_id = widget.Msg(ERROR_ON_MISSING_GENE_ID)
        missing_tax_id = widget.Msg(ERROR_ON_MISSING_TAX_ID)

    def __init__(self, parent=None):
        super().__init__(parent)

        self._changedFlag = False
        self.__invalidated = False
        self.__runstate = OWKEGGPathwayBrowser.Initializing
        self.__in_setProgress = False

        self.controlArea.setMaximumWidth(250)
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.widgetLabel(box, "No data on input\n")

        gui.separator(self.controlArea)

        gui.checkBox(self.controlArea, self, "useReference", "From signal", box="Reference", callback=self.Update)

        gui.separator(self.controlArea)

        gui.checkBox(
            self.controlArea,
            self,
            "showOrthology",
            "Show pathways in full orthology",
            box="Orthology",
            callback=self.UpdateListView,
        )

        gui.checkBox(
            self.controlArea,
            self,
            "autoResize",
            "Resize to fit",
            box="Image",
            callback=self.UpdatePathwayViewTransform,
        )

        box = gui.widgetBox(self.controlArea, "Cache Control")

        gui.button(
            box,
            self,
            "Clear cache",
            callback=self.ClearCache,
            tooltip="Clear all locally cached KEGG data.",
            default=False,
            autoDefault=False,
        )

        gui.separator(self.controlArea)

        gui.auto_commit(self.controlArea, self, "autoCommit", "Commit")

        gui.rubber(self.controlArea)

        spliter = QSplitter(Qt.Vertical, self.mainArea)
        self.pathwayView = PathwayView(self, spliter)
        self.pathwayView.scene().selectionChanged.connect(self._onSelectionChanged)
        self.mainArea.layout().addWidget(spliter)

        self.listView = QTreeWidget(
            allColumnsShowFocus=True, selectionMode=QTreeWidget.SingleSelection, sortingEnabled=True, maximumHeight=200
        )

        spliter.addWidget(self.listView)

        self.listView.setColumnCount(4)
        self.listView.setHeaderLabels(["Pathway", "P value", "Genes", "Reference"])

        self.listView.itemSelectionChanged.connect(self.UpdatePathwayView)

        select = QAction("Select All", self, shortcut=QKeySequence.SelectAll)
        select.triggered.connect(self.selectAll)
        self.addAction(select)

        self.data = None
        self.input_genes = []
        self.tax_id = None
        self.use_attr_names = None
        self.gene_id_attribute = None
        self.gene_id_column = None

        self.ref_data = None
        self.ref_genes = []
        self.ref_tax_id = None
        self.ref_use_attr_names = None
        self.ref_gene_id_attribute = None
        self.ref_gene_id_column = None

        self.pathways = {}
        self.ncbi_gene_map = []
        self.org = None

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

        self.infoLabel.setText("No data on input\n")

    def clear(self):
        """
        Clear the widget state.
        """
        self.pathways = {}
        self.org = None

        self.infoLabel.setText("No data on input\n")
        self.listView.clear()
        self.pathwayView.SetPathway(None)

        self.send("Selected Data", None)
        self.send("Unselected Data", None)

    def SetData(self, data=None):
        if self.__runstate == OWKEGGPathwayBrowser.Initializing:
            self.__initialize_finish()

        self.Error.clear()
        if data:
            self.data = data
            self.tax_id = str(self.data.attributes.get(TAX_ID, None))
            self.use_attr_names = self.data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.gene_id_attribute = self.data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.gene_id_column = self.data.attributes.get(GENE_ID_COLUMN, None)

            if not (
                self.use_attr_names is not None and ((self.gene_id_attribute is None) ^ (self.gene_id_column is None))
            ):

                if self.tax_id is None:
                    self.Error.missing_annotation()
                    return

                self.Error.missing_gene_id()
                return

            elif self.tax_id is None:
                self.Error.missing_tax_id()
                return

            self.warning(0)
            self.error(0)
            self.information(0)

            self.__invalidated = True
        else:
            self.clear()

    def SetRefData(self, data=None):
        self.information(1)

        if data is not None and self.useReference:
            self.ref_data = data
            self.ref_tax_id = str(self.ref_data.attributes.get(TAX_ID, None))
            self.ref_use_attr_names = self.ref_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.ref_gene_id_attribute = self.ref_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.ref_gene_id_column = self.ref_data.attributes.get(GENE_ID_COLUMN, None)

            if not (
                self.ref_use_attr_names is not None
                and ((self.ref_gene_id_attribute is None) ^ (self.ref_gene_id_column is None))
            ):

                if self.ref_tax_id is None:
                    self.Error.missing_annotation()
                    return

                self.Error.missing_gene_id()
                return

            elif self.ref_tax_id is None:
                self.Error.missing_tax_id()
                return

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

        org_code = self.org.org_code

        if self.showOrthology:
            self.koOrthology = kegg.KEGGBrite("ko00001")
            self.listView.setRootIsDecorated(True)
            path_ids = {s[-5:] for s in self.pathways.keys()}

            def _walkCollect(koEntry):
                num = koEntry.title[:5] if koEntry.title else None
                if num in path_ids:
                    return [koEntry] + reduce(
                        lambda li, c: li + _walkCollect(c), [child for child in koEntry.entries], []
                    )
                else:
                    c = reduce(lambda li, c: li + _walkCollect(c), [child for child in koEntry.entries], [])
                    return c + (c and [koEntry] or [])

            allClasses = reduce(lambda li1, li2: li1 + li2, [_walkCollect(c) for c in self.koOrthology], [])

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
                    item.setText(2, "%i of %i" % (len(genes), len(self.input_genes)))
                    item.setText(3, "%i of %i" % (ref, len(self.ref_genes)))
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
                item.setText(2, "%i of %i" % (len(genes), len(self.input_genes)))
                item.setText(3, "%i of %i" % (ref, len(self.ref_genes)))
                item.pathway_id = id
                items.append(item)

        self.bestPValueItem = items and items[0] or None
        self.listView.expandAll()
        for i in range(4):
            self.listView.resizeColumnToContents(i)

        if self.bestPValueItem:
            index = self.listView.indexFromItem(self.bestPValueItem)
            self.listView.selectionModel().select(index, QItemSelectionModel.ClearAndSelect)

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
        self._pathwayTask = concurrent.Task(function=lambda: get_kgml_and_image(item.pathway_id))
        self._pathwayTask.finished.connect(self._onPathwayTaskFinshed)
        self._executor.submit(self._pathwayTask)

    def _onPathwayTaskFinshed(self):
        self.setEnabled(True)
        pathway_id, self.pathway = self._pathwayTask.result()
        objects = self.pathways.get(pathway_id, [[]])[0]  # [ncbi_gene_id]
        # map ncbi_gene_id to keg_id for display
        objects = flatten(relation_map(self.ncbi_gene_map, objects))
        self.pathwayView.SetPathway(self.pathway, objects)

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
            self.__get_input_genes()
            self.input_genes = set(self.input_genes)
        except ValueError:
            self.error(0, "Cannot extract gene names from input.")

        self.information(1)

        self.org = kegg.KEGGOrganism(kegg.from_taxid(self.tax_id))

        if self.useReference and self.ref_data:
            self.__get_ref_genes()
            self.ref_genes = set(self.ref_genes)
        else:
            self.ref_genes = self.org.get_ncbi_ids()

        def run_enrichment(org_code, genes, reference, progress=None):
            # We use the kegg pathway gene sets provided by 'geneset' for
            # the enrichment calculation.

            kegg_api = kegg.api.CachedKeggApi()
            link_map = kegg_api.link(org_code, "pathway")  # [(pathway_id, kegg_gene_id)]
            ncbi_gene_map = kegg_api.conv(org_code, 'ncbi-geneid')  # [(ncbi_gene_id, kegg_gene_id)]
            ncbi_gene_map = [(_1.split(":", 1)[1], _2) for _1, _2 in ncbi_gene_map]
            link_map = relation_join(link_map, [(_2, _1) for _1, _2 in ncbi_gene_map])  # [(pathway_id, ncbi_gene_id)]
            kegg_sets = relation_list_to_multimap(link_map)  #  {pathway_id -> [ncbi_gene_ids]}
            # map kegg gene ids to ncbi_gene_ids.
            kegg_sets = geneset.GeneSets(
                sets=[geneset.GeneSet(gs_id=ddi, genes=set(genes)) for ddi, genes in kegg_sets.items()]
            )
            pathways = pathway_enrichment(kegg_sets, genes, reference, callback=progress)
            # Ensure that pathway entries are pre-cached for later use in the
            # list/tree view
            kegg_pathways = kegg.KEGGPathways()
            kegg_pathways.pre_cache(pathways.keys(), progress_callback=progress)

            return pathways, ncbi_gene_map

        self.progressBarInit()
        self.setEnabled(False)
        self.infoLabel.setText("Retrieving...\n")

        progress = concurrent.methodinvoke(self, "setProgress", (float,))

        run_func = partial(run_enrichment, self.org.org_code, self.input_genes, self.ref_genes, progress)
        self._enrichTask = concurrent.Task(function=run_func)
        self._enrichTask.finished.connect(self._onEnrichTaskFinished)
        self._executor.submit(self._enrichTask)

    def _onEnrichTaskFinished(self):
        self.setEnabled(True)
        self.setBlocking(False)
        try:
            pathways, ncbi_gene_map = self._enrichTask.result()
        except Exception:
            raise

        self.progressBarFinished()

        self.pathways = pathways
        self.ncbi_gene_map = ncbi_gene_map

        if not self.pathways:
            self.warning(0, "No enriched pathways found.")
        else:
            self.warning(0)

        self.infoLabel.setText("{} unique gene names on input\n".format(len(set(self.input_genes))))

        self.UpdateListView()

    @Slot(float)
    def setProgress(self, value):
        if self.__in_setProgress:
            return

        self.__in_setProgress = True
        self.progressBarSet(value)
        self.__in_setProgress = False

    def __get_input_genes(self):
        """
        Extract and return gene names from `data`.
        """
        self.input_genes = []

        if self.use_attr_names:
            for variable in self.data.domain.attributes:
                self.input_genes.append(str(variable.attributes.get(self.gene_id_attribute, '?')))
        else:
            genes, _ = self.data.get_column_view(self.gene_id_column)
            self.input_genes = [str(g) for g in genes]

        if len(self.input_genes) <= 0:
            raise ValueError("No gene names in data.")

    def __get_ref_genes(self):
        """
        Extract and return gene names from `data`.
        """
        self.ref_genes = []

        if self.ref_use_attr_names:
            for variable in self.ref_data.domain.attributes:
                self.ref_genes.append(str(variable.attributes.get(self.ref_gene_id_attribute, '?')))
        else:
            genes, _ = self.ref_data.get_column_view(self.ref_gene_id_column)
            self.ref_genes = [str(g) for g in genes]

    def selectAll(self):
        """
        Select all items in the pathway view.
        """
        changed = False
        scene = self.pathwayView.scene()
        with disconnected(scene.selectionChanged, self._onSelectionChanged):
            for item in scene.items():
                if item.flags() & QGraphicsItem.ItemIsSelectable and not item.isSelected():
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
            selectedGenes = reduce(set.union, [item.marked_objects for item in selectedItems], set())
            # map kegg_ids back to ncbi_gene_id
            backmap = [(_2, _1) for _1, _2 in self.ncbi_gene_map]
            selectedGenes = set(flatten(relation_map(backmap, selectedGenes)))
            if self.use_attr_names:
                selected = [
                    column
                    for column in self.data.domain.attributes
                    if self.gene_id_attribute in column.attributes
                    and str(column.attributes[self.gene_id_attribute]) in selectedGenes
                ]
                data = self.data[:, selected]
                self.send("Selected Data", data)
            else:
                selected_indices = []
                other_indices = []

                for row_index, row in enumerate(self.data):
                    gene_in_row = str(row[self.gene_id_column])
                    if gene_in_row in self.input_genes and gene_in_row in selectedGenes:
                        selected_indices.append(row_index)
                    else:
                        other_indices.append(row_index)

                if selected_indices:
                    selected = self.data[selected_indices]
                else:
                    selected = None

                if other_indices:
                    other = self.data[other_indices]
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

    return {_id: (genes, p_val, len(ref)) for (_id, genes, ref), p_val in zip(result_sets, p_values)}


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


if __name__ == "__main__":

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

    sys.exit(test_main())
