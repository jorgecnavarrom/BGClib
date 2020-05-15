# BGC lib

This is the official repository of the **Biosynthetic Gene Cluster Library**, a Python library with classes and methods that facilitates manipulation and analysis of BGCs.

For more information on BGClib, see [here](./BGClib/Readme.md).

:warning: This code is under heavy development!

As an application of what it can do, two companion scripts are included:

# BGC SVG generator

This script will generate high-quality figures from GenBank files. It is not limited to biosynthetic gene clusters.

## Features

* Input: GenBank or .bgc/.bgccase files (see [here](#bgc_toolkit)). GenBank files (.gb/.gbk) must have the protein translation.
* Data filters: There are three types of filters
  - Include strings. Only filenames (or BGC identifiers, in the case of .bgc/.bgccase files) with these strings will be included. Default: 'region', 'cluster'.
  - Exclude strings. Filenames or identifiers with this string(s) will be ignored. Default: 'final'. 
  
  The default values for these two filters ensure that results from antiSMASH (v4, v5) will not include the whole genome.
  
  - BGC list. A tab-separated text file. First column: cluster/identifier name. Second (optional) column: Protein Id. This file can be used to 1) further filter the input data, 2) define an order (in case all BGCs are included in a single figure) and 3) to align.
  
* There are some options to change the styling of the figure through an external configuration file.
* Figures can be produced individually or in a single SVG.
* Additionally, the genes can be decorated by drawing regions inside them that depict domains in their coresponding translations. 

## The whole options

```
Input:
  -i INPUTFOLDERS [INPUTFOLDERS ...], --inputfolders INPUTFOLDERS [INPUTFOLDERS ...]
                        Folder(s) to look for .gb and .gbk files (note: output folder will not preserve the structure of input folders).
  -f FILES [FILES ...], --files FILES [FILES ...]
                        File(s) used to draw the figures (accepted: .gb .gbk, .bgc, .bgccase). Note: for .bgc and .bgccase files, inclusion rules by --include, --exclude and --bgclist will be applied to the internal BGC identifier, not the name of the file
  --hmm HMM [HMM ...]   Location of .hmm file(s). This will also enable internal hmm models. Note that if the SVG style options have 'draw_domains=False', no domain prediction will be made, even if .hmm files are specified
  -l BGCLIST, --bgclist BGCLIST
                        A file containing a list of BGC identifiers (i.e. filename without .gb or .gbk extension). If specified, use it to filter all the BGCs found with --inputfolders or --files. If --stacked is used, this list will determine the order (and
                        filename). An optional second column (separated by a tab) with the Protein ID can be specified. If this column is present, the BGC will be mirrored if needed such that the gene encoding the Protein ID in the second column is in the forward
                        strand. Additionally, if --stacked is used and all Protein IDs are present, the corresponding gene will also be used for horizontal alignment. Any extra columns or rows starting with the '#' character will be ignored. The BGC identifiers in
                        this file are treated case-sensitive.
  --include [INCLUDE [INCLUDE ...]]
                        Specify string(s) to filter which BGCs will be accepted. In the case of .gb or .gbk files, the filter is applied to the filename. For data stored as .bgc or .bgccase files, the filter is applied to the BGC(s) identifier. If the argument is
                        present but no parameters are given, the filter will be ignored. If the argument is not present, the default is to use the strings 'region' and 'cluster')
  --exclude [EXCLUDE [EXCLUDE ...]]
                        Specify string(s) to filter which BGCs will be rejected. Similar rules are applied as with --include. If the argument is not present, the default is to use 'final'.

Processing options:
  --cfg CFG             Configuration file with SVG style. Default: 'SVG_arrow_options.cfg'
  -m, --mirror          Toggle to mirror each BGC figure. Ignored with --stacked or --bgclist
  --override            Use domain prediction in .bgc and .bgccase files, even if they already contain domain data (does not overwrite input files).
  -c CPUS, --cpus CPUS  Number of CPUs used for domain prdiction. Default: all available

Output:
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Folder where results will be put (default='output')
  -s, --stacked         If used, all BGCs will be put in the same figure. Default: each BGC has its own SVG.
  -g, --gaps            If --stacked is used, toggle this option to leave gaps when a particular BGC is not found in the input data
```

## Examples

Here is a simple example with the cladofulvin BGC:
```
./BGC_SVG_generator.py --files ./examples/data/cladofulvin/Cladofulvin_Final.gbk --outputfolder examples/output/cladofulvin --include
```

It has been renamed from a fungiSMASH result, so it doesn't include the usual 'region' string and we have to disable that filter by using `--include`. This command will produce

![Cladofulvin_Final](examples/output/cladofulvin/Cladofulvin_Final.svg "Cladofulvin. color_mode=random-pastel")

mibig 1.4

sophorolipid single
./BGC_SVG_generator.py --files examples/data/BGC0001274.1.gbk --hmm ~/Databases/pfam/32_top_JGI/top.hmm --outputfolder examples/output/sophorolipid --include
./BGC_SVG_generator.py --files examples/data/BGC0001274.1.gbk --hmm ~/Databases/pfam/32_top_JGI/top.hmm --outputfolder examples/output/sophorolipid --include --mirror

cercosporin

mibig 2.0
mycophenolic acid (single, stacked, aligned)

different styles

# BGC Toolkit

# Requirements

Installation of the required libraries to make BGClib work through (mini)conda is recommended. Here's a list of what you'll need. Version in parenthesis is the one that is known to work, but newer versions should work as well.

* biopython (1.76)
* hmmer (3.3)
* lxml (4.5.0)
