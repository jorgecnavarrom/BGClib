# BGC lib

This is the official repository of the **Biosynthetic Gene Cluster Library**, a Python library with classes and methods that facilitates manipulation and analysis of BGCs.

For more information on BGClib, see [here](./BGClib/Readme.md).

:warning: This code is under heavy development!

# BGC Toolkit

## What is it

It is a tool that facilitates domain annotation of protein data from GenBank and fasta files. For this, it harnesses the capabilities of BGClib

## Features overview

What can it do:

* Annotate domains
* ~~Extract protein sequences annotated with requested domains~~
* ~~Extract the subsequence of requested domains~~
* Generate high quality SVG files

Additionally, if working with biosynthetic gene clusters you can:

* Use antiSMASH results as input to label proteins as **core biosynthetic proteins** (CBPs)
* Extract protein (sub)sequences of CBPs
* Print metadata of the files per BGC and CBP (e.g. protein id, CBP content, internal CBP identifiers)
        
And finally, if you're working with fungal biosynthetic gene clusters in particular, you also can:

* Classify CBPs into more specific types using included hmm models (e.g. non-reducing PKS)

Other features:

* GenBank, fasta and binary [^1] input 
* Input filters (default values tailored for antiSMASH output; based on file/BGC name)
* Domain annotation using multiple hmm libraries.
* Save as binary files [^1]

[^1]: Currently uses Python's Pickle module for serialization of `BGC` or `BGCCollection` objects from BGClib

# Requirements

Installation of the required libraries to make BGClib work through (mini)conda is recommended. Here's a list of what you'll need. Version in parenthesis is the one that is known to work, but newer versions should work as well.

* biopython (1.78)
* hmmer (3.3.2)
* lxml (4.5.0)

# Documentation

Please see the **[BGCtoolkit wiki](https://github.com/jorgecnavarrom/BGClib/wiki)** for the documentation

# Results overview

A quick overview of results using the [example output](https://fungismash.secondarymetabolites.org/upload/fungal-example/index.html) from fungiSMASH (on _Aspergillus fumigatus_ Af293).

BGGtoolkit can produce metadata from all the input files, such as:

* [BGC Metadata](examples/AfumigatusAf293.metadata.BGCs.tsv):

| BGC | antiSMASH products | Core Biosynthetic Protein content | Core Biosynthetic Protein IDs | Core Biosynthetic Protein Identifiers | Metabolites |
| --- | ------------------ | --------------------------------- | ----------------------------- | ------------------------------------- | ----------- |
| CM000169.1.region001 | T1PKS | rPKS | EAL87813.1 | CM000169.1.region001~L0+CDS7 |
| CM000169.1.region002 | NRPS | NRPS | EAL90366.1 | CM000169.1.region002~L0+CDS8 |
| CM000169.1.region003 | terpene |  |  |  |
| CM000169.1.region004 | betalactone | NRPS-like | EAL90832.1 | CM000169.1.region004~L0+CDS5 |


* [CBP Metadata](examples/AfumigatusAf293.metadata.CBPs.tsv):

| BGC | Core Biosynthetic Protein type | Protein identifier | Protein Id | Gene Id | Domain architecture |
| --- | ------------------------------ | ------------------ | ---------- | ------- | ------------------- |
| CM000169.1.region001 | rPKS | CM000169.1.region001~L0+CDS7 | EAL87813.1 | | `KS | KS_C | KS_Ce | AT | DH | CMeT | KR | T/ACP >` |
|CM000169.1.region002 |NRPS | CM000169.1.region002~L0+CDS8 | EAL90366.1 | | `A | T/ACP | C | C | A | C | A | T/ACP | C | A | T/ACP | C | C | T/ACP | C | T/ACP >` |
| CM000169.1.region004 | NRPS-like | CM000169.1.region004~L0+CDS5 | EAL90832.1 | | `DMAP_binding | A | A >` |
| CM000169.1.region005 | NRPS-like | CM000169.1.region005~L0+CDS9 | EAL91049.1 | | `A | A_C >` |

* And a [summary](examples/AfumigatusAf293.metadata.summary.txt) of all the regions found, and their CBP compositions:

```
AfumigatusAf293 summary file

This collection contains
* 37 BGCs
* 0 Proteins

Core Biosynthetic Composition count (BGCs):
7       NRPS
6       NRPS-like
5       rPKS
...
```

fungiSMASH results report a 80% similarity between region `CM000175.1.region001` and neosartoricin B (cluster `BGC00001144` in [MIBiG](https://mibig.secondarymetabolites.org/repository/BGC0001144/index.html#r1c1)). Here's a comparison of both regions (top: MIBiG entry; bottom, `CM000175.1.region001`):

![neosartoricin](examples/neosartoricin_gcf.svg)

Here, colored boxes represent genomic regions that will code for predicted domains. Introns are drawn by default.

