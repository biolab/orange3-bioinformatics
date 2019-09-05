# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
