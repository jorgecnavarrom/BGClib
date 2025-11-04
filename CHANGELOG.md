# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)

## [Unreleased]

### Changed

- Changed contact email
- Improved BGC.load style
- BGCtoolkit: correctly leave a gap if a BGC isn't found. Fixed code for 'ribbon' shape
- Domain colors:
  * Added colors for all Pfam 37.1
  * Added hex colors in the color file `domain_color_file_IDh.tsv`
  * Domains colors (except for CBP domains) were re-colored so that domains belonging to the same clan (`Pfam-A.clans.tsv`) were colored similarly
- HMM_DB
  * hmmpress and extracting domain metadata now using PyHmmer
- ArrowerOps
  * More thorough checking of values when loading
- Improved readability:
  * Changed to fstrings
- BGC.load: Fixed bug where some proteins in a biosynthetic region would not be annotated as such
- SVG
  * Correctly center a single SVG vertically
  * `horizontal_margin` is now called `topbottom_margin`
  * New option for shape: "Ribbon" for narrower figures
  * Correctly applied `topbottom_margin` between each figure (stacked)
  * Fixed correct width in stacked figures
- add option 'Ribbon' for the gene shape

### Planned/TODO: 

- separate BGClib into subfiles (e.g. one for collections, one for proteins, one for a BGC, etc)
- separate arrower logic (different folder?). Make a class to produce svgs with some sort of interface. A demo script would take a gbk or a pair of fasta/gff and that class to quickly produce a figure. 
- SVG arrower: main arrow class should be 'CDS', not 'exon'
- SVG arrower: optimize stacked multi-svgs as possible
- What to do with all the print statements? These are usually warnings. Is there a better way (logging?). Also, better to be super descriptive? (e.g. print also the name of the class/method triggering the warning)
- `BGCCollection.predict_domains`: use cpus to paralellize calculate_domain_sets
- Threads or cpus? read https://newvick.com/python-concurrency/
- Create functions to (un)pickle collections and bgcs inside each class
- `BGCLocus.gene_coordinates`: check note about splicinc causing mismatches to the coordinates
- Turn pickled files (.bgc, .bgccase) into SQLite files