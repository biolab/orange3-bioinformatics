# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
