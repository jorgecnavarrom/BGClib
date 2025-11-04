# BGC lib

**Biosynthetic Gene Cluster Library** is a library with classes and methods that allow easy manipulation and analysis of BGCs. In a nutshell, it stores protein data from genomic loci that can be annotated with domain information.

Some of its classes:

## HMM_DB

Stores information about hmm databases to be used with [HMMER](http://hmmer.org/)
- Add hmm databases
- Store aliases for known domains (e.g. `ketoacyl-synt` -> `KS`)
- Store domain colors

## BGCCollection

Intended for making whole-bgc set operations like:
- Domain prediction (can also be done individually at the protein level)
- Distance calculation (between e.g. two organism's BGCs, a new set of BGCs against a precomputed database, etc.)

## BGC

Main class. Will have BGCProtein objects associated to them (ideally through another level of organization: an object of the BGCLocus class) and supporting information to *classify* these clusters.

## BGCProtein

Currently the most worked on class. Some of its properties:
* An object of this class can have three identifiers: the official `identifier` attribute (intended for a string of the type BGCName:CDS#:Ref_Accession), the original `accession` and a de-replicated `ref_accession` provided by [NCBI's IPG database](https://www.ncbi.nlm.nih.gov/ipg)
