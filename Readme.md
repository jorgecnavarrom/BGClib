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

* Metadata at BGC level

![Metadata BGC](examples/AfumigatusAf293.metadata.BGCs.tsv)

* Metadata at core biosynthetic protein (CBP) level

![Metadata CBP](examples/AfumigatusAf293.metadata.CBPs.tsv)

* And a summary of all the regions found, and their CBP compositions

![Summary](examples/AfumigatusAf293.metadata.summary.txt)

fungiSMASH results report an 80% similarity between region `CM000175.1.region001` and neosartoricin (cluster `BGC00001144` in [MIBiG](https://mibig.secondarymetabolites.org/repository/BGC0001144/index.html#r1c1)). Here's a comparison of both regions (top: MIBiG entry, bottom, `CM000175.1.region001`):

![neosartoricin](examples/neosartoricin_gcf.svg)



