# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.8.2] - 2024-03-20
### Changed
- [#319](https://github.com/biolab/orange3-bioinformatics/pull/319): GEO Datasets: Move GDSInfo initialization in a worker thread

## [4.8.1] - 2024-01-10
### Changed
- [#345](https://github.com/biolab/orange3-bioinformatics/pull/345): Raise the required resdk package version
### Fixed
- [#346](https://github.com/biolab/orange3-bioinformatics/pull/346): Fix OWDifferentialExpression data inputs and outputs


## [4.8.0] - 2023-09-01
### Changed
- [#338](https://github.com/biolab/orange3-bioinformatics/pull/338): Bioinformatics add-on now has improved documentation.
- [#344](https://github.com/biolab/orange3-bioinformatics/pull/344): GEO Data Sets now outputs metas as discrete and not as string variables.
- [#343](https://github.com/biolab/orange3-bioinformatics/pull/344): OWdictyExpress now outputs 'Time' as continuous and not as string variable.

### Fixed
- [#340](https://github.com/biolab/orange3-bioinformatics/pull/340): Delaunay attribute vertices was deprecated since Scipy version 0.12.0 in favour of simplices and it was finally removed in version 1.11.0.
- [#344](https://github.com/biolab/orange3-bioinformatics/pull/344): fixes sorting that did not work when you clicked on the column in Gene Set Enrichment widget.

 
## [4.7.2] - 2023-07-06
### Changed
- [#335](https://github.com/biolab/orange3-bioinformatics/pull/335):  Time is now column name for Dicty express data

## [4.7.1] - 2023-02-06

## [4.7.0] - 2023-02-06

## [4.6.0] - 2022-11-24
### Changed
- [#311](https://github.com/biolab/orange3-bioinformatics/pull/311):  Add QC data to Genialis Expression widget
- [#322](https://github.com/biolab/orange3-bioinformatics/pull/322):  (OWGeneSetEnrichment) Export gene sets enrichment report
- [#311](https://github.com/biolab/orange3-bioinformatics/pull/311):  Add QC data to Genialis Expression widget
- 
### Fixed
- [#308](https://github.com/biolab/orange3-bioinformatics/pull/308):  Annotate Projection: Fix commit invoke
- [#313](https://github.com/biolab/orange3-bioinformatics/pull/313):  Fix compatibility with Python 3.10
- [#314](https://github.com/biolab/orange3-bioinformatics/pull/314):  Replace colorpalette and other necessary fixes
- [#316](https://github.com/biolab/orange3-bioinformatics/pull/316):  kegg: Fix kegg.from_taxid mapping due to change in api search results
- [#320](https://github.com/biolab/orange3-bioinformatics/pull/320):  QtCore.Qt.ItemIsTristate is deprecated. Use ItemIsUserTristate instead
- [#323](https://github.com/biolab/orange3-bioinformatics/pull/323):  (OWGenialisExpressions) fix sign in dialog, remove info box, use defered commit


## [4.5.0] - 2021-10-07

### Changed
- [#290](https://github.com/biolab/orange3-bioinformatics/pull/290): (OWAnnotateProjection) follow widget UI style guidelines
- [#289](https://github.com/biolab/orange3-bioinformatics/pull/289): Refactor GUI components of resolwe widgets
- [#286](https://github.com/biolab/orange3-bioinformatics/pull/286): Change the documentation theme

### Fixed
- [#289](https://github.com/biolab/orange3-bioinformatics/pull/289): Fix for breaking changes in requests-cache package
- [#288](https://github.com/biolab/orange3-bioinformatics/pull/288): Fix for issues that occurred after changes in widget auto-summary were introduced


## [4.4.0] - 2021-06-08

### Changed
- [#282](https://github.com/biolab/orange3-bioinformatics/pull/285): (OWGenialisExpressions) use CollectionTables to fetch data

## [4.3.1] - 2021-03-01

### Fixed
- [#281](https://github.com/biolab/orange3-bioinformatics/pull/281): (OWKEGGPathwayBrowser) Fix selection.

### Changed
- [#282](https://github.com/biolab/orange3-bioinformatics/pull/282): (OWdictyExpress) use CredentialManager and ConcurrentWidgetMixin.


## [4.3.0] - 2020-12-23

### Added
- [#278](https://github.com/biolab/orange3-bioinformatics/pull/278): (OWGenialisExpressions) Metadata now includes Relations and optional Clinical data.


## [4.2.0] - 2020-10-27

### Fixed
- [#267](https://github.com/biolab/orange3-bioinformatics/pull/267): (OWGenialisExpressions) Fixed animation toggle.

### Changed
- [#269](https://github.com/biolab/orange3-bioinformatics/pull/269): (OWGenialisExpressions) Changed default window size.
- [#274](https://github.com/biolab/orange3-bioinformatics/pull/274): (OWGenialisExpressions) Use number of samples/genes as quantiles in Quantile transformation.

### Added
- [#268](https://github.com/biolab/orange3-bioinformatics/pull/268): (OWGenialisExpressions) Use Orange.widgets.credentials.CredentialManager to securely store user password.
- [#270](https://github.com/biolab/orange3-bioinformatics/pull/270): (OWGenialisExpressions) Additionally, filter data objects by process type and input annotations.
- [#272](https://github.com/biolab/orange3-bioinformatics/pull/272), [#276](https://github.com/biolab/orange3-bioinformatics/pull/276): (OWGenialisExpressions) Rename table attributes when matching genes and display report.

### Other
The release includes some other minor bug fixes and project maintenance.

## [4.1.0] - 2020-08-27

### Fixed
- [#214](https://github.com/biolab/orange3-bioinformatics/pull/214): Fixed a problem in OWGEODatasets with filter component, [related isue](https://github.com/biolab/orange3-bioinformatics/issues/210). 
- [#233](https://github.com/biolab/orange3-bioinformatics/pull/233): Fixed output summary and stop showing links if they are not available, [related isue](https://github.com/biolab/orange3-bioinformatics/issues/228).  
- [#231](https://github.com/biolab/orange3-bioinformatics/pull/231): Properly use gene matcher in dictyExpress widget.
- [#237](https://github.com/biolab/orange3-bioinformatics/pull/237): Fixed a problem with enrichment. Results from Gene Set Enrichemnt and GOBrowser widgets are now consistent, [related isue](https://github.com/biolab/orange3-bioinformatics/issues/234).
- [#241](https://github.com/biolab/orange3-bioinformatics/pull/241): OWGEODatasets: Remove `wrap` reimplementation.
- [#242](https://github.com/biolab/orange3-bioinformatics/pull/242): Marker genes settings migration.

### Changed
- [#215](https://github.com/biolab/orange3-bioinformatics/pull/215): Bump minimum version of point-annotator module.
- [#216](https://github.com/biolab/orange3-bioinformatics/pull/216): Update icon for Gene Set Enrichment widget.
- [#221](https://github.com/biolab/orange3-bioinformatics/pull/221): Update icons for Volcano plot and Annotate Projection widgets.
- [#260](https://github.com/biolab/orange3-bioinformatics/pull/260): Update icon for dictyExpress widget.
- [#219](https://github.com/biolab/orange3-bioinformatics/pull/219): Complete redesign of Marker Genes widget.
- [#244](https://github.com/biolab/orange3-bioinformatics/pull/244): Unify add-on descriptions, [related isue](https://github.com/biolab/orange3/issues/4850).  

### Added
- [#215](https://github.com/biolab/orange3-bioinformatics/pull/215): Add support for sparse data in Annotation module.
- [#245](https://github.com/biolab/orange3-bioinformatics/pull/245): setup Github Actions and tox.
- [#252](https://github.com/biolab/orange3-bioinformatics/pull/251): Add GenialisExpressions widget, [project status](https://github.com/biolab/orange3-bioinformatics/projects/2).

## [4.0.0] - 2019-09-05

### Fixed
- [#203](https://github.com/biolab/orange3-bioinformatics/pull/203): Correct column data in GeneSets widget.
- [#205](https://github.com/biolab/orange3-bioinformatics/pull/205): Fix issue where GeneSets widget crashed on network error. Display friendly warning that cached files are used.
- [#206](https://github.com/biolab/orange3-bioinformatics/pull/206): Fix issue where GEODatasets widget crashed on network error. Display friendly warning that cached files are used.

### Changed
- [#168](https://github.com/biolab/orange3-bioinformatics/pull/168): 
Overhaul of gene module. We no longer use pickled dictionaries to map gene names to their corresponding Entrez Ids. 
We use SQLite tables that are indexed using FTS5 SQLite extension.
- [#204](https://github.com/biolab/orange3-bioinformatics/pull/204):
Package cellannotation replaced with point-annotator 

### Added
- [13b54e7](https://github.com/biolab/orange3-bioinformatics/pull/168/commits/13b54e7d93e09283ca5edfae4f11468fc2c0b12b): HomoloGene module (Wrapper around NCBI HomoloGene database)
- [#178](https://github.com/biolab/orange3-bioinformatics/pull/178): add Homology widget.

### Removed
- [a713abe](https://github.com/biolab/orange3-bioinformatics/commit/a713abe3b799efcfbd40f50aa724a4924fcf6df8):
removed legacy serverfiles update scripts in favor of [bioinformatics-serverfiles.](https://github.com/JakaKokosar/bioinformatics-serverfiles)

## [3.5.0] - 2019-07-19
Did not keep a changelog...
