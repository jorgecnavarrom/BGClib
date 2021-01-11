#!/usr/bin/env python

"""
Handle Biosynthetic Gene Cluster objects
- Input:
    * .gb/.gbk
    * .bgc/.bgccase
    * .fasta
- Filter:
    * Include string
    * Exclude string
    * The user's own list (bgclist)
- Processing:
    * Domain prediction
    * Merge split domains
- Output:
    * Figures (individual or stacked)
    * BGC objects (.bgc or .bgccase)
    * Metadata (list of protein ids)
    * Sequences (core biosynthetic proteins + selected domains)

TODO: 
    * Convert fasta files into .bgcprotein or .bgcproteincase
"""

import os
import sys
from shutil import copyfile
from typing import Tuple
import json
import argparse
import pickle
from collections import defaultdict
from itertools import chain
# import shelve
from lxml import etree
from pathlib import Path
from multiprocessing import cpu_count
from BGClib import HMM_DB, BGC, BGCLocus, BGCCollection, ProteinCollection, \
    ArrowerOpts, BGCDomain, valid_CBP_types, valid_CBP_types_antiSMASH, \
    valid_CBP_types_fungal
from datetime import datetime

__author__ = "Jorge Navarro"
__version__ = "0.3"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"


# Prepare arguments for the script
def CMD_parser():
    parser = argparse.ArgumentParser(description="BGC toolkit v{}. Tools to \
        facilitate biosynthetic gene cluster handling and \
        visualization.".format(__version__))
    
    group_input = parser.add_argument_group("Input")
    
    # TODO: change to --inputfoldersgbk? (and add --inputfoldersbgc?)
    group_input.add_argument("-i", "--inputfolders", nargs='+', type=Path, \
        help="Folder(s) to search (recursively) for .gb, .gbk and .gbff files.\
        If the file starts with generic terms ('scaffold', 'contig', \
        'Supercontig'), the parent folder's name will be used in the internal \
        BGC name for that file.")
    group_input.add_argument("-f", "--files", nargs='+', type=Path, \
        help="Input individual files (accepted: .gb .gbk, .bgc, .bgccase, \
        .fasta, .proteincase).")

    group_filtering = parser.add_argument_group(title="Filtering",
        description="Include or exclude BGCs based on their names. Note:\
        for .bgc and .bgccase files, inclusion rules by --include, --exclude\
        and --bgclist will be applied to the internal BGC identifier, not to \
        the name of the file.")
    group_filtering.add_argument("--include", nargs='*', default=['region', 'cluster'], 
        help="Specify string(s) to filter which BGCs will be included. In the \
        case of .gb or .gbk files, the filter is applied to the filename. For \
        data stored as .bgc or .bgccase files, the filter is applied to the \
        BGC(s) identifier. If the option is present but no arguments are \
        given, the filter will be ignored (all BGCs will be included). If the \
        argument is not present, the default is to use the strings 'region' \
        and 'cluster').")
    group_filtering.add_argument("--exclude", nargs='*', default=['final'], 
        help="Specify string(s) to filter which BGCs will be rejected. \
        Similar rules are applied as with --include. If the argument is not \
        present, the default is to use 'final'.")
    group_filtering.add_argument("-l", "--bgclist", help="A tab-separated file \
        containing a list of BGC names (first column) and protein names (second\
        column). The BGC names filter input from --files and --inputfolder. If \
        only the protein column is present, it will be used to filter input from \
        .proteincase or .fasta files. If used with --svg and --stacked, the \
        contents of this file also define the order (and protein to align) in \
        the final SVG figure. Any extra columns or rows starting with the '#' \
        character will be ignored.", type=Path)   
        
    group_processing = parser.add_argument_group("Domain options")
    group_processing.add_argument("--hmms", nargs='*', help="List of paths to \
        .hmm file(s). This will also enable internal hmm models (use without \
        arguments to only use internal models).")
    group_processing.add_argument("--update-domains", help="Use domain prediction \
        on input, even if they were marked as already having domain prediction. \
        If new and old hits overlap, only the best-scoring hit will be kept.", \
        default=False, action="store_true")
    group_processing.add_argument("--clear_domains", help="Clean all domain data \
        from input before domain prediction.", default=False, action="store_true")
    group_processing.add_argument("-c", "--cpus", type=int, default=cpu_count(), 
        help="Number of CPUs used for domain prediction. Default: all available.")
    internal_alias_file = Path(__file__).parent / "BGClib" / "data" / "CBP_domains.tsv"
    group_processing.add_argument("--alias", type=Path, help="Argument is a \
        tab-separated file with domain alias. First column: domain ID, Second \
        column: alias. If domain prediction is used, this information will also \
        be stored (e.g. to use afterwards when making SVG figures). Default \
        aliases are stored in the special file {}".format(internal_alias_file), \
        dest='alias_file')

    group_post_proc = parser.add_argument_group("Post-processing options")

    group_post_proc.add_argument("--merge", default=False, action="store_true",
        help="Try to fix successive domains of the same type that have been \
        split. This could be useful for phylogenetic analysis but will not \
        fix incorrect gene predictions. Use with care.")
    # group_post_proc.add_argument("--jsonfolders", nargs='+', type=Path, \
    #     help="Use MIBiG-style JSON files to annotate metadata (organism, \
    #     metabolites, literature)")
    
    group_output = parser.add_argument_group(title="Output", 
        description="Basic output options")
    
    default_output = Path("./") / "output"
    group_output.add_argument("-o", "--outputfolder", 
        default=default_output, help="Base folder where \
        results will be put (default='{}').".format(default_output))
    group_output.add_argument("--metadata", type=str, help="Writes \
        information files at three levels: a summary of the whole collection, \
        a summary per BGC (CBP content) and a summary per core protein (domain \
         organization). Argument is the basename of these files (no extension).\
        Activated by default if --bgccase is used.")

    group_svg = parser.add_argument_group("SVG options")

    group_svg.add_argument("--svg", default=False, 
        action="store_true", help="Toggle to enable SVG output for each BGC.")
    group_svg.add_argument("--svgcfg", type=Path,
        default=(Path(__file__).parent/"SVG_arrow_options.cfg"),
        help="Configuration file with SVG style. Default: \
        'SVG_arrow_options.cfg'")
    group_svg.add_argument("-m", "--mirror", default=False, 
        action="store_true", help="Toggle to mirror each BGC figure. Ignored \
        with --stacked or --bgclist")
    group_svg.add_argument("-s", "--stacked", type=str, help="If used with \
        --svg, all BGC SVGs will be put in the same figure. The argument of \
        this parameter is the filename (no extension).")
    # NOTE: if we have more than one comparison method, use 
    # "choices=['rock', 'paper', 'scissors']"
    group_svg.add_argument("--comparative", default=False, action="store_true",
        help="If --stacked and --bgclist are used, calculate protein similarity \
        between each pair of BGCs and represent it as connecting bands.")
    group_svg.add_argument("--gaps", default=False, action="store_true",
        help="If --stacked is used, toggle this option to leave gaps\
        when a particular BGC or protein is not found in the input data \
        (--bgclist).")
    
    group_organize = parser.add_argument_group(title="Organize biosynthetic output", \
        description="(Optional) Use these parameters with the remaining output \
        options to separate data into different sub-folders according to \
        selected Core Biosynthetic Type(s) (CBTs). If any of these options are \
        used, the final biosynthetic types used will be those common between \
        a) the types in the dataset, b) the ones in the cbt file and c) the\
        ones in the cbt-include list, minus those from the cbt-exclude option.\
        Currently supported types from antiSMASH: {}. Currently supported \
        fungal types (these will need domain prediction): {}.".format(\
        ", ".join(["'{}'".format(x) for x in sorted(valid_CBP_types_antiSMASH.keys())]), \
        ", ".join(["'{}'".format(x) for x in sorted(valid_CBP_types_fungal.keys())])))
    group_organize.add_argument("--cbt-file", type=Path, nargs='?', 
        const=(Path(__file__).parent/"CBP_output_types.cfg"), \
        help="A file with valid core biosynthetic types. If no argument \
        is given, the file 'CBP_output_types.cfg' will be read.")
    # TODO explain alias
    group_organize.add_argument("--cbt-include", type=str, nargs='*', \
        help="A space-separated list of core biosynthetic types to use when \
        organizing output. If 'all' or not argument is given, all posible CBTs \
        will be used. Additionally for this parameter, it is possible to \
        search for specific domains (e.g. 'all:C' or 'nrPKS:SAT'). Domains \
        may be the model's ID, the stored alias, or the alias given in the \
        file specified with --alias")
    group_organize.add_argument("--cbt-exclude", type=str, nargs='+', \
        help="A space-separated list of all the CBTs to exclude.")

    group_file_output = parser.add_argument_group(title="File output", \
        description="These options will trigger metadata output")
    group_file_output.add_argument("--bgc", default=False,
        action="store_true", help="Toggle to save binary (.bgc) files of each \
        BGC (from .bgc, .bgccase and GenBank files). Data from fasta input \
        files will not be included")
    group_file_output.add_argument("--bgccase", type=str, help="Output a single \
        binary file with all BGCs (.bgccase). The argument of this parameter \
        is the filename (no extension). Data from fasta input files will not \
        be included.")
    group_file_output.add_argument("-p", "--proteincase", type=str, \
        help="Output .proteincase files for each type of Core Biosynthetic\
        Protein type specified in the file 'CBP_output_types.cfg'. The \
        argument of this parameter is the basename for all .proteincase \
        files")
    group_file_output.add_argument("--genbank", default=False, \
        action="store_true", help="Toggle to create GenBank files. Currently, \
        only works for gbk input files. Use with 'cbt' options to separate by\
        biosynthetic type.")
    group_organize.add_argument("--cbt-fasta", default=False, action="store_true",
        help="Toggle to output sequences of (specialized metabolite's) Core \
        Biosynthetic Proteins for each CBT defined in the \
        'CBP_output_types.cfg' file. Only available when cbt* options are \
        enabled.")

    return parser.parse_args()


def sanitize(filename):
    """
    Remove illegal characters of strings that are intended to be used as filenames
    """
    illegal = r"/|\:<>*"
    sanitized_filename = filename
    for c in illegal:
        sanitized_filename = sanitized_filename.replace(c, "_")

    return sanitized_filename


# Basic validation of input files
def check_parameters(args):
    """
    Checkes whether parameters are valid and paths exist
    """

    inputfiles = args.files
    inputfolders = args.inputfolders

    # Stop if no input was specified
    if not inputfiles and not inputfolders:
        sys.exit("\nStop: no input data. See options using -h")

    # Stop if no output was specified
    if not (args.svg or args.bgc or args.bgccase or args.proteincase \
            or args.metadata or args.cbt_fasta or args.genbank):
        sys.exit("\nStop: no output options were selected (use --svg, --bgc, etc.)")

    if args.cpus < 1:
        sys.exit("Error: parameter --cpus must be a positive integer")

    if args.cbt_file:
        if not args.cbt_file.is_file():
            sys.exit("Error (--cbt-file): cannot open specified file")

    if args.alias_file:
        if not args.alias_file.is_file():
            sys.exit("Error (--alias): cannot open specified file")

    if inputfiles:
        for file_ in inputfiles:
            f = Path(file_)
            if not f.is_file():
                sys.exit("Error (--files): {} is not a file".format(f))
            # if f.suffix == ".bgccase":
            #     try:
            #         s = shelve.open(str(f), 'r')
            #     except:
            #         sys.exit("Error (--files): {} not a valid bgccase database".format(f))
            #     else:
            #         s.close()
            # else:
            #     if not f.is_file():
            #         sys.exit("Error (--files): {} is not a file".format(f))
                
    if inputfolders:
        for folder in inputfolders:
            if not Path(folder).is_dir():
                sys.exit("Error (--inputfolders): {} is not a folder".format(folder))
            
    if args.hmms:
        for hmm in args.hmms:
            hmm_file = Path(hmm)
            if not hmm_file.is_file():
                sys.exit("Error (--hmm): {} is not a file".format(hmm_file))
            if hmm_file.suffix not in {".hmm", ".HMM"}:
                sys.exit("Error (--hmm): {} does not have a .hmm extension".format(hmm_file))

    if args.bgclist:
        if not Path(args.bgclist).is_file():
            sys.exit("Error: (--bgclist): {} is not a file".format(args.bgclist))

    # Style options
    if args.svg:
        if not args.svgcfg.is_file():
            sys.exit("Error (--svgcfg): {} is not a file".format(args.svgcfg))
        
    if args.metadata and args.bgccase:
        print("Warning: --metadata and --bgccase. The argument of --metadata will be ignored")


# Input stage
def valid_name(name, include, exclude, filter_bgc):
    """
    Checks whether a filename is valid and should be included in the analysis
    based on the allowed strings (args.include) or strings to be avoided
    (args.exclude) as well as the criterium that the BGC is included in the
    filter list.
    It is expected that the parameter 'name' is the name of the file \
    without extension
    """
    
    if include:
        if not any([word in name for word in include]):
            return False
        
    if exclude:
        if exclude != [] and any(word in name for word in exclude):
            return False
    
    if filter_bgc:
        if name not in filter_bgc:
            return False
    
    return True


def read_bgc_list(bgclist):
    """
    Reads the bgclist file
    * Tab-separated file
    * First column: bgc id. Can be empty
    * Second column: protein identifier or protein_id
    * Third column: comment (or any line starting with #)
    """
    filterlist = list()
    with open(args.bgclist) as f:
        for line in f:
            if line[0] == "#" or line.strip() == "":
                continue
            
            xline = line.split("\t")

            if len(xline) > 1:
                # assume second column is protein_id
                bgc_id = xline[0]
                protein_id = xline[1]
                filterlist.append((bgc_id.strip(), protein_id.strip()))
            else:
                filterlist.append((xline[0], ""))
    
    if len(filterlist) == 0:
        sys.exit("Error: filter BGC list given but the file is empty...")

    return(filterlist)


def read_alias_file(alias_file: Path):
    external_alias = dict()
    with open(alias_file) as a:
        for line in a:
            if line[0] == "#" or line.strip() == "":
                continue
            xline = line.strip().split("\t")
            try:
                ID, alias = xline[0], xline[1]
            except IndexError:
                # there should be at least two columns, separated by tabs
                continue
            else:
                external_alias[ID] = alias

    return external_alias


def get_files(args, filter_bgc_prot:list) -> Tuple[BGCCollection, ProteinCollection]:
    """
    Collect all data from the input, using filters defined by user to produce
     the collection of data to work with.
    """

    inputfolders = args.inputfolders
    files = args.files
    include = args.include
    exclude = args.exclude

    # If these words are found at the beginning of the filename (which can
    # typically happen in generic antiSMASH output), append the parent folder's
    # name to the BGC identifier
    forbidden_words = {'scaffold', 'contig', 'Supercontig', 'c00'}

    bgc_col = BGCCollection()
    prot_col = ProteinCollection()
    gbk_files = dict()
    all_bgcs = defaultdict(list) # includes repeated bgc ids

    # Build filter sets for both BGCs and proteins
    filter_bgc_set = set()
    filter_prot = set()
    for bgc_id, protein_id in filter_bgc_prot:
        if bgc_id != "":
            filter_bgc_set.add(bgc_id)
        else:
            filter_prot.add(protein_id)

    if inputfolders:
        for folder in inputfolders:
            for gbk_file in chain(folder.glob("**/*.gb"), folder.glob("**/*.gbk")):
                if valid_name(gbk_file.stem, include, exclude, filter_bgc_set):
                    identifier = gbk_file.stem
                    if any([identifier.startswith(fword) for fword in forbidden_words]):
                        identifier = "{}.{}".format(gbk_file.parts[-2], identifier)
                    gbk_files[identifier] = gbk_file
                    all_bgcs[identifier].append(gbk_file)

    if files:
        for f in files:
            if f.suffix.lower() in {".gb", ".gbk", ".gbff"}:
                if valid_name(f.stem, include, exclude, filter_bgc_set):
                    identifier = f.stem
                    if any([identifier.startswith(fword) for fword in forbidden_words]):
                        identifier = "{}.{}".format(f.resolve().parts[-2], identifier)
                    gbk_files[identifier] = f
                    all_bgcs[identifier].append(f)

            elif f.suffix.lower() == ".bgc":
                with open(f, "rb") as dc:
                    bgc = pickle.load(dc)
                    bgc_id = bgc.identifier
                if valid_name(bgc_id, include, exclude, filter_bgc_set):
                    bgc_col.bgcs[bgc_id] = bgc
                    all_bgcs[bgc_id].append(f)

                    if bgc_id in gbk_files:
                        print("\tWarning: BGC {} has been found in a .gbk/.gb file and a .bgc file (using the latter)".format(bgc_id))
                        del gbk_files[bgc_id]

            elif f.suffix.lower() == ".bgccase":
                # TODO: fix shelves to always use the same backend or find
                # another solution
                # with shelve.open(str(f), flag='r') as col:
                with open(f, "rb") as dc:
                    col = pickle.load(dc)

                # if we've got a filter list, it's faster to use it
                if filter_bgc_set:
                    for bgc_id in filter_bgc_set:
                        try:
                            bgc = col.bgcs[bgc_id]
                        except KeyError:
                            pass
                        else:
                            # still have to check include and exclude filters
                            if valid_name(bgc_id, include, exclude, set()):
                                bgc_col.bgcs[bgc_id] = bgc
                                all_bgcs[bgc_id].append(f)
                                # if we already had it in a gbk file, use this version instead
                                if bgc_id in gbk_files:
                                    print("\tWarning: BGC {} has been found in a .gbk/.gb file and a .bgccase collection (using the latter)".format(bgc_id))
                                    del gbk_files[bgc_id]
                else:
                    for bgc_id, bgc in col.bgcs.items():
                        if valid_name(bgc_id, include, exclude, filter_bgc_set):
                            bgc_col.bgcs[bgc_id] = bgc
                            all_bgcs[bgc_id].append(f)

                            # if we already had it in a gbk file, use this version instead
                            if bgc_id in gbk_files:
                                print("\tWarning: BGC {} has been found in a .gbk/.gb file and a .bgccase collection (using the latter)".format(bgc_id))
                                del gbk_files[bgc_id]

            elif f.suffix.lower() == ".fasta":
                prot_col.fasta_load(f)

            elif f.suffix.lower() == ".proteincase":
                with open(f, 'rb') as pc:
                    temp_prot_col = pickle.load(pc)

                if temp_prot_col.proteins.keys() & prot_col.proteins.keys():
                    print("\tWarning! Over-writing proteins with the same header. Make sure you have unique headers in all your protein inputs. This will probably result in unintended results")
                    prot_col.proteins.update(temp_prot_col)
                prot_col.proteins.update(temp_prot_col.proteins)

            else:
                print("\tWarning: unknown format ({})".format(f))
    
    # Applying filtering to the protein collection. This is tricky because
    # the filtering may include protein_ids, so we have to traverse the 
    # whole collection. Also, more practical to do it with all the proteins
    # already in a collection (as opposed to filtering-while-adding-them as
    # with BGCs)
    if filter_prot:
        final_prot_col = ProteinCollection()
        for protein in prot_col.proteins.values():
            if protein.identifier in filter_prot or protein.protein_id in filter_prot:
                final_prot_col.proteins[protein.identifier] = protein
    else:
        final_prot_col = prot_col

    # Report repeated BGCs
    print_warning_repeated = False
    for bgc_id, file_list in all_bgcs.items():
        if len(file_list) > 1:
            print_warning_repeated = True
            break
    if print_warning_repeated:
        rbl_file = o / "repeated_bgc_locations.txt"
        print("\tWarning: multiple files for the same BGC were found. See file {}".format(rbl_file))
        with open(rbl_file, "w") as rbl:
            for bgc_id, file_list in all_bgcs.items():
                if len(file_list) > 1:
                    rbl.write("{}\n".format(bgc_id))
                    for f in file_list:
                        rbl.write("\t{}\n".format(f))
                    rbl.write("\n")

    # Finally, add the gbk files to the BGC collection
    for identifier, gbk_path in gbk_files.items():
        bgc_col.add_gbk(gbk_path, identifier)

    return bgc_col, final_prot_col, gbk_files


def predict_domains(args, hmmdbs, bgc_collection, protein_collection):
    """
    Runs domain prediction on BGCs and Proteins depending or parameters
     Note that we're using their 'attempted_domain_prediction' attribute
     as opposed to using the ones that currently don't have any domain 
     predicted
    """

    bgc_collection_need_dom_pred = BGCCollection()
    protein_collection_need_dom_pred = ProteinCollection()

    # if clear_domains: clear all data first, re-predict domains
    # Notice that the protein role and type attributes are NOT cleared
    if args.clear_domains or args.update_domains:
        # Need to convert generator to list b/c we're removing items in the loop
        for bgc_id in list(bgc_collection.bgcs.keys()):
            if args.clear_domains:
                bgc_collection.bgcs[bgc_id].clear_domain_predictions()
            bgc_collection_need_dom_pred.bgcs[bgc_id] = bgc_collection.bgcs.pop(bgc_id)
        
        for protein_identifier in list(protein_collection.proteins.keys()):
            if args.clear_domains:
                protein_collection.proteins[protein_identifier].clear_domain_predictions()
            protein_collection_need_dom_pred.proteins[protein_identifier] = protein_collection.proteins.pop(protein_identifier)

    # nothing selected: only predict domains for "clean" objects
    else:
        for bgc_id in list(bgc_collection.bgcs.keys()):
            if not bgc_collection.bgcs[bgc_id].attempted_domain_prediction:
                bgc_collection_need_dom_pred.bgcs[bgc_id] = bgc_collection.bgcs.pop(bgc_id)

        for protein_identifier in list(protein_collection.proteins.keys()):
            if not protein_collection.proteins[protein_identifier].attempted_domain_prediction:
                protein_collection_need_dom_pred.proteins[protein_identifier] = protein_collection.proteins.pop(protein_identifier)

    print("\tPredicting domains...")
    bgc_collection_need_dom_pred.predict_domains(hmmdbs, cpus=hmmdbs.cores)
    protein_collection_need_dom_pred.predict_domains(hmmdbs, cpus=hmmdbs.cores)
    print("\t...done!\n")

    print("\nApplying classification rules")
    bgc_collection_need_dom_pred.classify_proteins(cpus=hmmdbs.cores)
    protein_collection_need_dom_pred.classify_proteins(cpus=hmmdbs.cores)
    print("\t...done!\n")

    # Finally, the data we'll be processing
    bgc_collection.bgcs.update(bgc_collection_need_dom_pred.bgcs)
    protein_collection.proteins.update(protein_collection_need_dom_pred.proteins)

    return


# Post-processing step
def fix_core_split_domains(bgc_col, prot_col):
    """
    Finds repeated domains and tries to merge them in a single one
    Currently for domains found in core biosynthetic proteins
    """

    # key: BGCProtein. 
    # Value: list whose elements are themselves lists of consecutive domain IDs
    protein_to_broken_doms = defaultdict(list) 

    # collect data
    for bgc in bgc_col.bgcs.values():
        for cbp_type in bgc.CBPcontent:
            for protein in bgc.CBPcontent[cbp_type]:
                if not protein.domain_list:
                    continue
                
                consecutive_domain_list = [] # consecutive dom objs relative to hmm
                current_id = ""
                for current_dom_num in range(len(protein.domain_list)-1):
                    current_dom = protein.domain_list[current_dom_num]
                    if current_id != current_dom.ID:
                        if consecutive_domain_list:
                            # a new type of dom was found; append list of fragments
                            protein_to_broken_doms[protein].append(consecutive_domain_list)
                        
                        consecutive_domain_list = []
                        
                        # look ahead to see if we need to start a new list
                        next_dom = protein.domain_list[current_dom_num+1]
                        if current_dom.ID == next_dom.ID and \
                                current_dom.hmm_to < next_dom.hmm_from:
                            consecutive_domain_list.append(current_dom)

                        current_id = current_dom.ID

                    # Current and previous domains are the same type. Check coords.
                    elif consecutive_domain_list:
                        if consecutive_domain_list[-1].hmm_to < current_dom.hmm_from:
                            consecutive_domain_list.append(current_dom)
                # check last domain
                current_dom = protein.domain_list[-1]
                if current_id != current_dom.ID:
                    if consecutive_domain_list:
                        protein_to_broken_doms[protein].append(consecutive_domain_list)
                # could have same ID but hmm coordinates don't indicate splitting
                elif consecutive_domain_list:
                    if consecutive_domain_list[-1].hmm_to < current_dom.hmm_from:
                        consecutive_domain_list.append(current_dom)
                        protein_to_broken_doms[protein].append(consecutive_domain_list)

    # check now proteins from Protein Collection
    for protein in prot_col.proteins.values():
        if not protein.role == "biosynthetic":
            continue

        consecutive_domain_list = [] # consecutive dom objs relative to hmm
        current_id = ""
        for current_dom_num in range(len(protein.domain_list)-1):
            current_dom = protein.domain_list[current_dom_num]
            if current_id != current_dom.ID:
                if consecutive_domain_list:
                    # a new type of dom was found; append list of fragments
                    protein_to_broken_doms[protein].append(consecutive_domain_list)
                
                consecutive_domain_list = []
                
                # look ahead to see if we need to start a new list
                next_dom = protein.domain_list[current_dom_num+1]
                if current_dom.ID == next_dom.ID and \
                        current_dom.hmm_to < next_dom.hmm_from:
                    consecutive_domain_list.append(current_dom)

                current_id = current_dom.ID

            # Current and previous domains are the same type. Check coords.
            elif consecutive_domain_list:
                if consecutive_domain_list[-1].hmm_to < current_dom.hmm_from:
                    consecutive_domain_list.append(current_dom)
        # check last domain
        current_dom = protein.domain_list[-1]
        if current_id != current_dom.ID:
            if consecutive_domain_list:
                protein_to_broken_doms[protein].append(consecutive_domain_list)
        # could have same ID but hmm coordinates don't indicate splitting
        elif consecutive_domain_list:
            if consecutive_domain_list[-1].hmm_to < current_dom.hmm_from:
                consecutive_domain_list.append(current_dom)
                protein_to_broken_doms[protein].append(consecutive_domain_list)

    # try to stitch fragments
    for protein in protein_to_broken_doms:
        # print(protein.identifier)
        for consecutive_domain_list in protein_to_broken_doms[protein]:
            # for d in consecutive_domain_list:
                # print("{}\t{}\t{}\t{}\t{}\t{}".format(protein.protein_id, \
                #     d.ID, d.hmm_from, d.hmm_to, d.ali_from, d.ali_to))
                # print(d.get_sequence())
            # print("")

            first_dom = consecutive_domain_list[0]
            last_dom = consecutive_domain_list[-1]
            # Warning! Score of the new domain is naively asssigned to be the
            # simple sum of the scores of its parts
            merged_domain = BGCDomain(protein, first_dom.ID, \
                first_dom.AC, first_dom.DE, first_dom.alias, \
                first_dom.ali_from, last_dom.ali_to, first_dom.hmm_from, \
                last_dom.hmm_to, 0, \
                sum(d.score for d in consecutive_domain_list), 0, "")
            for d in consecutive_domain_list:
                protein.domain_list.remove(d)
            protein.domain_list.append(merged_domain)
            protein.domain_list.sort(key=lambda x:x.ali_from)


# Output: svg options
def draw_svg_individual(args, o: Path, bgc_col: BGCCollection, \
    prot_col: ProteinCollection, hmmdbs: HMM_DB, filter_bgc_prot) -> None:
    """
    Draw individual SVGs per BGC or Protein
    """

    mirror = args.mirror
    svgopts = ArrowerOpts(args.svgcfg)

    filter_bgc = dict()
    filter_prot = set()
    for bgc_id, protein_id in filter_bgc_prot:
        if bgc_id != "":
            filter_bgc[bgc_id] = protein_id
        else:
            filter_prot.add(protein_id)

    # BGCs
    filtered_bgc_col = BGCCollection()
    # Use filter_bgc (if needed) to build the collection with bgcs we need to visualize
    if filter_bgc:
        for bgc_id in (filter_bgc.keys() & bgc_col.bgcs.keys()):
            filtered_bgc_col.bgcs[bgc_id] = bgc_col.bgcs[bgc_id]
    else:
        filtered_bgc_col = bgc_col

    for bgc_id, bgc in filtered_bgc_col.bgcs.items():
        m = mirror
        # see if we need to mirror the cluster
        if bgc_id in filter_bgc:
            pid = filter_bgc[bgc_id]
            
            if pid != "":
                for p in bgc.protein_list:
                    if pid == p.protein_id or pid == p.identifier:
                        if not p.forward:
                            m = True
                            break
        
        coregenearch = ""
        if len(bgc.CBPtypes) > 0:
            coregenearch = "_[{}]".format(",".join(bgc.CBPtypes))
                
        m_info = ""
        if m:
            m_info = "_m"
            
        bgc_name = o / "{}{}{}.svg".format(bgc_id, coregenearch, m_info)
        bgc.BGC_SVG(bgc_name, hmmdbs, svgopts, mirror=m)

    # Proteins. Id can be protein_id or identifier
    filtered_prot_col = ProteinCollection()
    if filter_prot:
        for protein in prot_col.proteins.values:
            if protein.identifier in filter_prot or protein.protein_id in filter_prot:
                filtered_prot_col.proteins[protein.identifier] = protein
    else:
        filtered_prot_col = prot_col

    for pid, protein in filtered_prot_col.proteins.items():
        cbp_arch = ""
        if protein.role == "biosynthetic":
            cbp_arch = "_[{}]".format(protein.protein_type)
            
        protein_name = o / "{}{}.svg".format(sanitize(pid), cbp_arch)
        protein.arrow_SVG(protein_name, hmmdbs, svgopts)

    return


def draw_svg_stacked(args, o: Path, bgc_col: BGCCollection, \
    prot_col: ProteinCollection, hmmdbs: HMM_DB, filter_bgc_prot):
    """
    Draws a stacked BGC figure
    """

    # Discard "trivial" case of not having an ordered list of BGCs/proteins
    if len(filter_bgc_prot) == 0:
        return draw_svg_stacked_simple(args, o, bgc_col, prot_col, hmmdbs)

    filename = o / "{}.svg".format(args.stacked)
    svgopts = ArrowerOpts(args.svgcfg)
    gaps = args.gaps
    thickness = svgopts.gene_contour_thickness
    
    # Read extra info for comparative fig.
    if args.comparative:
        comparative_figure_spacing = 1.0
        with open(args.svgcfg) as f:
            for line in f:
                if line.startswith("comparative_figure_spacing"):
                    try:
                        comparative_figure_spacing = float(line.split("=")[1].strip())
                    except:
                        pass
                    else:
                        if comparative_figure_spacing < 0:
                            comparative_figure_spacing = 1.0
                    break

    # Read original filtering options and make two lists: one with BGCs and 
    # the other with the internal protein identifier of reference for 
    # horizontal alignment. Stand-alone proteins will be converted to 
    # pseudo-BGCs for the figure
    draw_order_bgcs = list()
    scaling = svgopts.scaling
    H = svgopts.arrow_height # used for loci spacing
    bgc_lengths = dict()
    needs_mirroring = dict()
    bgc_distance_to_target = dict()
    start_to_target_max_offset = 0
    for bgc_id, pid in filter_bgc_prot:
        if bgc_id != "":
            if bgc_id not in bgc_col.bgcs:
                print("\tSVG (stacked, bgclist): --bgclist used but {} not found in BGC data".format(bgc_id))
                draw_order_bgcs.append(None)
                continue

            bgc = bgc_col.bgcs[bgc_id]
            draw_order_bgcs.append(bgc)
            # second term in the addition: inter-loci spacing element
            L = sum([locus.length/scaling for locus in bgc.loci]) \
                + H * (len(bgc.loci)-1) \
                + thickness
            bgc_lengths[bgc_id] = L

            if pid == "":
                print("\tSVG (stacked, bgclist): Warning, {} has not reference Protein Id".format(bgc_id))
                needs_mirroring[bgc_id] = False
                bgc_distance_to_target[bgc_id] = -1
                continue

            protein = None
            for locus_num, locus in enumerate(bgc.loci):
                for candidate_protein in locus.protein_list:
                    if pid == candidate_protein.protein_id \
                            or pid == candidate_protein.identifier:
                        needs_mirroring[bgc_id] = not candidate_protein.forward
                        protein = candidate_protein

                        # get start_to_target distance
                        if protein.forward:
                            # lenghts of each loci until target + 
                            # inter-loci spacing + distance of current \
                            # locus to target_protein
                            start_to_target = \
                                sum([locus.length/scaling for locus in bgc.loci[:locus_num]]) \
                                    + H * locus_num \
                                    + protein.cds_regions[0][0]/scaling 
                            bgc_distance_to_target[bgc_id] = start_to_target
                        else:
                            start_to_target = \
                                sum([locus.length/scaling for locus in bgc.loci[locus_num+1:]]) \
                                    + H * (len(bgc.loci) - locus_num - 1) \
                                    + (locus.length - protein.cds_regions[-1][1])/scaling
                            bgc_distance_to_target[bgc_id] = start_to_target

                        if start_to_target > start_to_target_max_offset:
                            start_to_target_max_offset = start_to_target

                        break

            # Couldn't find protein specified by user; typo?
            if protein is None:
                print("\tSVG (stacked, bgclist): Warning, cannot find reference Protein Id [{}] for {}".format(pid, bgc_id))
                needs_mirroring[bgc_id] = False
                bgc_distance_to_target[bgc_id] = -1
                continue

        # got stand-alone protein, convert to BGC
        else:
            protein = None
            try:
                # Easy case: given pid is internal identifier
                protein = prot_col.proteins[pid]
            except KeyError:
                # then pid come from annotations? look for protein_id
                for protein_id, candidate_protein in prot_col.proteins.items():
                    if protein_id == pid:
                        protein = candidate_protein
                        break
                if protein is None:
                    print("\tSVG (stacked, bgclist): Warning, stand-alone protein {} not found".format(pid))
                    draw_order_bgcs.append(None)
                    continue

            L = (protein.cds_regions[-1][1] - protein.cds_regions[0][0])/scaling + thickness

            # Create fake BGC
            locus = BGCLocus()
            locus.length = L
            locus.protein_list.append(protein)
            locus.gene_coordinates.append((protein.cds_regions[0][0], protein.cds_regions[-1][1]))

            bgc = BGC()
            bgc.identifier = protein.identifier
            bgc.loci.append(locus)
            bgc.protein_list.append(protein)
            bgc.proteins[protein.identifier] = protein

            draw_order_bgcs.append(bgc)
            bgc_lengths[bgc.identifier] = L
            needs_mirroring[bgc.identifier] = False
            bgc_distance_to_target[bgc.identifier] = 0


    # obtain max_L considering all the starting offsets
    max_L = 0
    for bgc_id in bgc_distance_to_target:
        if bgc_distance_to_target[bgc_id] == -1:
            max_L = max(max_L, bgc_lengths[bgc_id])
        else:
            max_L = max(max_L, bgc_lengths[bgc_id] \
                + start_to_target_max_offset \
                - bgc_distance_to_target[bgc_id])

    # Start SVG internal structure
    row_height = 2*svgopts.arrow_height # one for the arrow, 0.5 + 0.5 for the head height
    inner_row_height = row_height + thickness
    base_attribs = {"version":"1.1", 
                    "baseProfile":"full", 
                    "width":str(int(max_L))
                    }
    root = etree.Element("svg", attrib=base_attribs, nsmap={None:'http://www.w3.org/2000/svg'})
    
    # Add each figure
    Yoffset = 0
    rows = 0
    for item in draw_order_bgcs:
        Yoffset = rows * inner_row_height

        if item is None:
            if gaps:
                rows += 1
            continue
        else:
            bgc = item
            bgc_id = bgc.identifier
        
        # Marked BGCs with no reference Protein Id won't have offset
        Xoffset = 0
        if bgc_distance_to_target[bgc_id] != -1:
            Xoffset = start_to_target_max_offset - bgc_distance_to_target[bgc_id]
        root.append(bgc.xml_BGC(Xoffset, Yoffset, hmmdbs, svgopts, needs_mirroring[bgc_id]))
        rows += 1
            
        #Yoffset = rows * inner_bgc_height
    
    # Now that we now how many BGCs were actually drawn, add height property
    root.attrib["height"] = str(int(thickness + rows*(row_height+thickness)))
    
    # Write SVG
    with open(filename, "bw") as f:
        f.write(etree.tostring(root, pretty_print=True))
        

def draw_svg_stacked_simple(args, o, bgc_col, prot_col, hmmdbs) -> None:
    """
    Draws a stacked BGC figure. 
    In this simple version, no order or alignment are given; all BGCs, then 
    proteins, will be drawn in arbitrary order
    """
    filename = o / "{}.svg".format(args.stacked)
    svgopts = ArrowerOpts(args.svgcfg)
    thickness = svgopts.gene_contour_thickness
    
    # Now obtain a list of protein objects that will be used for alignment
    # Also get information about mirroring and offset distancing
    scaling = svgopts.scaling
    H = svgopts.arrow_height # used for loci spacing

    # find max_L
    max_L = 0
    for bgc in bgc_col.bgcs.values():
        L = sum([locus.length/scaling for locus in bgc.loci]) \
            + H * (len(bgc.loci)-1) \
            + thickness
        if L > max_L:
            max_L = L
    for protein in prot_col.proteins.values():
        L = (protein.cds_regions[-1][1] - protein.cds_regions[0][0])/scaling + thickness
        if L > max_L:
            max_L = L

    # Start SVG internal structure
    row_height = 2*svgopts.arrow_height # one for the arrow, 0.5 + 0.5 for the head height
    inner_row_height = row_height + thickness
    base_attribs = {"version":"1.1", 
                    "baseProfile":"full", 
                    "width":str(int(max_L))
                    }
    root = etree.Element("svg", attrib=base_attribs, nsmap={None:'http://www.w3.org/2000/svg'})
    
    # Add each figure
    Yoffset = 0
    rows = 0
    for bgc in bgc_col.bgcs.values():
        Yoffset = rows * inner_row_height
        root.append(bgc.xml_BGC(0, Yoffset, hmmdbs, svgopts, False))
        rows += 1
    for protein in prot_col.proteins.values():
        Yoffset = rows * inner_row_height
        root.append(protein.xml_arrow(hmmdbs, svgopts, 0, Yoffset))
        rows += 1

    # Now that we now how many BGCs were actually drawn, add height property
    root.attrib["height"] = str(int(thickness + rows*(row_height+thickness)))
    
    # Write SVG
    with open(filename, "bw") as f:
        f.write(etree.tostring(root, pretty_print=True))
        
    return


# Output: metadata
def write_metadata(o: Path, metadata_base: str, bgc_col: BGCCollection, \
    prot_col: ProteinCollection, alias: dict):
    """
    Writes information files at three levels:
    - A summary of the whole collection
    - A summary per BGC (core biosynthetic protein content)
    - A summary per core protein (domain organization)
    """

    # Whole collection summary
    with open(o / "{}.metadata.summary.txt".format(metadata_base), "w") as s:
        s.write("{} summary file\n\n".format(metadata_base))
        s.write("This collection contains\n")
        s.write("* {} BGCs\n".format(len(bgc_col.bgcs)))
        s.write("* {} Proteins\n".format(len(prot_col.proteins)))

        if len(bgc_col.bgcs) > 0:
            cbp_type_histogram = defaultdict(int)
            # cbp_type_bgcs = defaultdict(list)
            for bgc in bgc_col.bgcs.values():
                cbp_comp = ", ".join(sorted(bgc.CBPtypes_set))
                if cbp_comp == "":
                    cbp_comp = "Other"
                cbp_type_histogram[cbp_comp] += 1
                # cbp_type_bgcs[cbp_comp].append(bgc_id) 
            s.write("\nCore Biosynthetic Composition count (BGCs):\n")
            for cbp_comp in sorted(cbp_type_histogram, key=cbp_type_histogram.get, reverse=True):
                # s.write("{}\t{}\t{}\n".format(cbp_type_histogram[cbp_comp], \
                #     cbp_comp, ", ".join(cbp_type_bgcs[cbp_comp])))
                s.write("{}\t{}\n".format(cbp_type_histogram[cbp_comp], cbp_comp))

        if len(prot_col.proteins) > 0:
            cbp_type_histogram = defaultdict(int)
            # cbp_type_prots = defaultdict(list)
            for protein in prot_col.proteins.values():
                if protein.role == "biosynthetic":
                    cbp_comp = protein.protein_type
                    cbp_type_histogram[cbp_comp] += 1
                    # cbp_type_prots[cbp_comp].append(protein.identifier)
            s.write("\nCore Biosynthetic Composition count (Proteins):\n")
            for cbp_comp in sorted(cbp_type_histogram, key=cbp_type_histogram.get, reverse=True):
                # s.write("{}\t{}\t{}\n".format(cbp_type_histogram[cbp_comp], \
                #     cbp_comp, ", ".join(cbp_type_prots[cbp_comp])))
                s.write("{}\t{}\n".format(cbp_type_histogram[cbp_comp], cbp_comp))

    # BGC summary
    if len(bgc_col.bgcs) > 0:
        with open(o / "{}.metadata.BGCs.tsv".format(metadata_base), "w") as m:
            header = "BGC\tantiSMASH products\tCore Biosynthetic Protein content\tCore Biosynthetic Protein IDs\tCore Biosynthetic Protein Identifiers\tMetabolites\n"
            m.write(header)
            for bgc_id in sorted(bgc_col.bgcs):
                bgc = bgc_col.bgcs[bgc_id]
                list_core_types = []
                list_protein_ids = []
                list_protein_identifiers = []
                all_domains = set() # TODO remove this. Just for testing
                for protein in bgc.protein_list:
                    all_domains.update(protein.domain_set)
                    if protein.role == "biosynthetic":
                        list_core_types.append(protein.protein_type)
                        # TODO currently doesn't print anything if protein_id is
                        # missing, but we could leave the comma-separated blank 
                        # spaces 
                        if protein.protein_id != "":
                            list_protein_ids.append(protein.protein_id)
                        list_protein_identifiers.append(protein.identifier)
                metabolites = ", ".join(m.name for m in bgc.metabolites)
                m.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bgc.identifier, 
                    ", ".join(bgc.products), 
                    ", ".join(list_core_types), 
                    ", ".join(list_protein_ids), 
                    ", ".join(list_protein_identifiers),
                    metabolites))

    # CBP summary
    with open(o / "{}.metadata.CBPs.tsv".format(metadata_base), "w") as c:
        header = "BGC\tCore Biosynthetic Protein type\tProtein identifier\tProtein Id\tGene Id\tDomain architecture\n"
        c.write(header)
        for bgc_id in sorted(bgc_col.bgcs):
            bgc = bgc_col.bgcs[bgc_id]
            for cbp_type in bgc.CBPcontent:
                for protein in bgc.CBPcontent[cbp_type]:
                    c.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bgc_id, \
                        cbp_type, protein.identifier, protein.protein_id, \
                        protein.gene, protein.domain_string(alias)))

        for identifier in sorted(prot_col.proteins):
            protein = prot_col.proteins[identifier]
            if protein.role != "biosynthetic":
                continue
            c.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(protein.parent_cluster_id, \
                protein.protein_type, identifier, protein.protein_id, \
                protein.gene, protein.domain_string(alias)))


# Output: organize data
def read_cbp_cfg(cbt_file):
    cbp_types = set()

    with open(cbt_file) as f:
        for line in f:
            if line[0] == "#":
                continue
            bgctype, decision = line.strip().split("=")[:2]
            if decision.strip().lower() == "true":
                cbp_types.add(bgctype.strip())
    return cbp_types


def get_cbt_types(args, bgc_col, prot_col) -> Tuple[set, dict]:
    """
    Reads input parameters from the "Organize biosynthetic output" section
    and produces a set that contains the valid core biosynthetic types
    """

    cbt_types = set()
    cbt_domains = dict() # key: cbt, value: domain ID list

    if not (args.cbt_file or args.cbt_include is not None or args.cbt_exclude):
        return cbt_types, cbt_domains

    cbt_file = args.cbt_file
    include_set = set()
    if args.cbt_include:
        include_set = set(args.cbt_include) # Can also contain domain options        
    cbt_include = set() # define but populate later as domain options need separation
    cbt_exclude = set()
    if args.cbt_exclude:
        cbt_exclude = set(args.cbt_exclude)

    # Check which types we actually have in the input data
    types_in_dataset = set()
    for bgc in bgc_col.bgcs.values():
        types_in_dataset.update(bgc.CBPtypes_set)
    for protein in prot_col.proteins.values():
        types_in_dataset.update(protein.domain_set)

    # If cbt-include was empty or someone explicitly put "all", activate
    # subfolder organization for all known types
    if "all" in set(include_set) or args.cbt_include == []:
        cbt_include = valid_CBP_types
    else:
        # separate domain options
        for item in list(include_set):
            if ":" in item:
                cbt, ID = item.split(":")
                include_set.remove(item)
                try:
                    cbt_domains[cbt].add(ID)
                except KeyError:
                    cbt_domains[cbt] = set([ID])
            else:
                cbt_include.add(item)

    if cbt_file:
        types_from_file = read_cbp_cfg(cbt_file)
    else:
        types_from_file = valid_CBP_types

    # combination of * what is valid, + what is actually in data, 
    # + what is in the file, + what is included, - what is excluded
    cbt_types = valid_CBP_types \
        & types_in_dataset \
        & types_from_file \
        & cbt_include \
        - cbt_exclude

    # Check domains
    # if cbt_domains and not (types_in_dataset & cbt_domains.keys() - {"all", "ALL"})):
    #     print("Warning: Domain output was specified but the input data does not contain the parent type")
    #     cbt_domains = dict()

    return cbt_types, cbt_domains


def create_folder(folder_path: Path) -> None:
    """
    Created a folder if necessary
    """
    if not folder_path.is_dir():
        os.makedirs(folder_path, exist_ok=True)

    return


# Output: files
def save_bgc_output(args, o: Path, cbt_types: set, cbt_domains: dict, \
    bgc_col: BGCCollection, alias: dict) -> None:
    """
    Save .bgc or .bgccase files to the output (sub)folder(s)
    """

    # Could happen if we only have protein input
    if not(bgc_col.bgcs):
        print("Warning, BGC output selected but no BGCs in input data")
        return

    # general output; all BGCs in the same folder
    if not (cbt_types or cbt_domains):
        if args.bgc:
            print("Saving individual BGCs...")
            for bgc_id, bgc in bgc_col.bgcs.items():
                target = o / "{}.bgc".format(sanitize(bgc_id))
                with open(target, "wb") as b:
                    pickle.dump(bgc, b)

        if args.bgccase:
            print("Saving BGC collection...")
            target = o / "{}.bgccase".format(sanitize(args.bgccase))
            with open(target, "wb") as c:
                pickle.dump(bgc_col, c)
            write_metadata(o, sanitize(args.bgccase), bgc_col, ProteinCollection(), alias)
        return

    #  For output data organization:
    folder_to_collection = dict()
    folder_to_suffix = dict()

    # Strategy is put every BGC in a collection linked to the target folder
    for bgc_id, bgc in bgc_col.bgcs.items():
        # target by type
        for cbt in (bgc.CBPtypes_set & cbt_types):
            target_folder = o / cbt
            folder_to_suffix[target_folder] = "_{}".format(cbt)
            try:
                folder_to_collection[target_folder].bgcs[bgc_id] = bgc
            except KeyError:
                folder_to_collection[target_folder] = BGCCollection()
                folder_to_collection[target_folder].bgcs[bgc_id] = bgc

        # Extend the set of the current BGC to hold all valid aliases
        internal_aliases = set() # work out internal aliases first
        for p in bgc.protein_list:
            for d in p.domain_list:
                if d.alias != "":
                    internal_aliases.add(d.alias)
        extended_bgc_domain_set = bgc.domain_set | \
            internal_aliases | \
            set(alias[d] for d in bgc.domain_set if d in alias)

        # target by domain content
        for cbt in (bgc.CBPtypes_set & cbt_domains.keys()):
            # dom_ID will contain the alias so target names will use it
            for dom_ID in (extended_bgc_domain_set & cbt_domains[cbt]):
                target_folder = o / cbt / "{}_{}".format(cbt, dom_ID)
                folder_to_suffix[target_folder] = "_{}_{}".format(cbt, dom_ID)
                try:
                    folder_to_collection[target_folder].bgcs[bgc_id] = bgc
                except KeyError:
                    folder_to_collection[target_folder] = BGCCollection()
                    folder_to_collection[target_folder].bgcs[bgc_id] = bgc

        if "all" in cbt_domains and (extended_bgc_domain_set & cbt_domains['all']):
            for dom_ID in extended_bgc_domain_set & cbt_domains['all']:
                target_folder = o / "domains_{}".format(dom_ID)
                folder_to_suffix[target_folder] = "_domains_{}".format(dom_ID)
                try:
                    folder_to_collection[target_folder].bgcs[bgc_id] = bgc
                except KeyError:
                    folder_to_collection[target_folder] = BGCCollection()
                    folder_to_collection[target_folder].bgcs[bgc_id] = bgc
    
    # Everything should be separated. Dump all data
    for target_folder, bgc_sub_col in folder_to_collection.items():
        create_folder(target_folder)

        if args.bgc:
            for bgc_id, bgc in bgc_sub_col.bgcs.items():
                target = target_folder / "{}.bgc".format(sanitize(bgc_id))
                with open(target, "wb") as b:
                    pickle.dump(bgc, b)

        if args.bgccase:
            suffix = folder_to_suffix[target_folder]
            target_base_name = "{}{}".format(sanitize(args.bgccase), suffix)
            target = target_folder / "{}.bgccase".format(target_base_name)
            with open(target, "wb") as c:
                pickle.dump(bgc_sub_col, c)
            # Output metadata by default
            write_metadata(target_folder, target_base_name, bgc_sub_col, ProteinCollection(), alias)
    
    return


def save_protein_output(args, o: Path, cbt_types: set, cbt_domains: dict, \
    bgc_col: BGCCollection, prot_col: ProteinCollection, alias: dict) -> None:
    """
    Save ProteinCollection objects to the output (sub)folder(s)
    """

    # Make a single collection with all CBPs
    all_proteins = ProteinCollection()
    all_proteins.proteins.update(prot_col.proteins)
    throw_warning = False
    for bgc in bgc_col.bgcs.values():
        for protein_list in bgc.CBPcontent.values():
            for protein in protein_list:
                if protein.identifier in all_proteins.proteins:
                    throw_warning = True
                all_proteins.proteins[protein.identifier] = protein
    if throw_warning:
        print("Protein collection output: warning, data from BGC input has overwritten protein data (same identifier)")

    # General output: all proteins in the same folder
    if not (cbt_types or cbt_domains):
        print("Saving Protein collection...")
        basename = sanitize(args.proteincase)
        target = o / "{}.proteincase".format(basename)
        with open(target, "wb") as c:
            pickle.dump(all_proteins, c)
        write_metadata(o, basename, BGCCollection(), all_proteins, alias)
        return


    # For output data organization:
    folder_to_collection = dict()
    folder_to_suffix = dict()

    # Same as BGC collection: make a sub collection linked to each target folder
    for p_id, protein in all_proteins.proteins.items():
        ptype = protein.protein_type
        # target by type
        if ptype in cbt_types:
            target_folder = o / ptype
            folder_to_suffix[target_folder] = "_{}".format(ptype)
            try:
                folder_to_collection[target_folder].proteins[p_id] = protein
            except KeyError:
                folder_to_collection[target_folder] = ProteinCollection()
                folder_to_collection[target_folder].proteins[p_id] = protein

        # extend the protein domains to include their aliases
        extended_prot_domain_set = protein.domain_set | \
            set(d.alias for d in protein.domain_list if d.alias != "") | \
            set(alias[d] for d in protein.domain_set if d in alias)

        # target by type + domain content
        if ptype in cbt_domains:
            for dom_ID in (extended_prot_domain_set & cbt_domains[ptype]):
                target_folder = o / ptype / "{}_{}".format(ptype, dom_ID)
                folder_to_suffix[target_folder] = "_{}_{}".format(ptype, dom_ID)
                try:
                    folder_to_collection[target_folder].proteins[p_id] = protein
                except KeyError:
                    folder_to_collection[target_folder] = ProteinCollection()
                    folder_to_collection[target_folder].proteins[p_id] = protein

        # target by domain
        if "all" in cbt_domains and (extended_prot_domain_set & cbt_domains["all"]):
            for dom_ID in (extended_prot_domain_set & cbt_domains["all"]):
                target_folder = o / "domains_{}".format(dom_ID)
                folder_to_suffix[target_folder] = "_domains_{}".format(dom_ID)
            try:
                folder_to_collection[target_folder].proteins[p_id] = protein
            except KeyError:
                folder_to_collection[target_folder] = ProteinCollection()
                folder_to_collection[target_folder].proteins[p_id] = protein
    
    # Dump data
    for target_folder, prot_sub_col in folder_to_collection.items():
        create_folder(target_folder)
        suffix = folder_to_suffix[target_folder]
        target_base_name = "{}{}".format(sanitize(args.proteincase), suffix)
        target = target_folder / "{}.proteincase".format(target_base_name)
        with open(target, "wb") as c:
            pickle.dump(prot_sub_col, c)
        write_metadata(target_folder, target_base_name, BGCCollection(), prot_sub_col, alias)

    return


def save_fasta(o: Path, cbt_types: set, cbt_domains: dict, bgc_col: BGCCollection, \
            prot_col: ProteinCollection, alias: dict) -> None:
    """
    Save fasta to the output (sub)folder(s)
    """

    # Need specific output
    if not (cbt_types or cbt_domains):
        print("Fasta output selected, but no specific CBT or domain selected...")
        return

    # Make a single collection with all CBPs
    all_proteins = ProteinCollection()
    all_proteins.proteins.update(prot_col.proteins)
    throw_warning = False
    for bgc in bgc_col.bgcs.values():
        for protein_list in bgc.CBPcontent.values():
            for protein in protein_list:
                if protein.identifier in all_proteins.proteins:
                    throw_warning = True
                all_proteins.proteins[protein.identifier] = protein
    if throw_warning:
        print("Fasta output: warning, data from BGC input has overwritten protein data (same identifier)")

    # For output data organization:
    folder_to_col_cbt = dict()
    folder_to_col_dom = dict() # these collections need to have a
    folder_to_suffix = dict()
    folder_to_dom = dict() # need to know which domain to extract from sequence

    # Same as BGC collection: make a sub collection linked to each target folder
    for p_id, protein in all_proteins.proteins.items():
        ptype = protein.protein_type
        # target by type
        if ptype in cbt_types:
            target_folder = o / ptype
            folder_to_suffix[target_folder] = ptype
            try:
                folder_to_col_cbt[target_folder].proteins[p_id] = protein
            except KeyError:
                folder_to_col_cbt[target_folder] = ProteinCollection()
                folder_to_col_cbt[target_folder].proteins[p_id] = protein

        # extend the protein domains to include their aliases
        extended_prot_domain_set = protein.domain_set | \
            set(d.alias for d in protein.domain_list if d.alias != "") | \
            set(alias[d] for d in protein.domain_set if d in alias)
        
        # target by type + domain content
        if ptype in cbt_domains:
            for dom_ID in (extended_prot_domain_set & cbt_domains[ptype]):
                target_folder = o / ptype / "{}_{}".format(ptype, dom_ID)
                folder_to_suffix[target_folder] = "{}_{}".format(ptype, dom_ID)
                folder_to_dom[target_folder] = dom_ID
                try:
                    folder_to_col_dom[target_folder].proteins[p_id] = protein
                except KeyError:
                    folder_to_col_dom[target_folder] = ProteinCollection()
                    folder_to_col_dom[target_folder].proteins[p_id] = protein

        # target by domain
        if "all" in cbt_domains and (extended_prot_domain_set & cbt_domains["all"]):
            for dom_ID in (extended_prot_domain_set & cbt_domains["all"]):
                target_folder = o / "domains_{}".format(dom_ID)
                folder_to_suffix[target_folder] = "all_{}_domains".format(dom_ID)
                folder_to_dom[target_folder] = dom_ID
            try:
                folder_to_col_dom[target_folder].proteins[p_id] = protein
            except KeyError:
                folder_to_col_dom[target_folder] = ProteinCollection()
                folder_to_col_dom[target_folder].proteins[p_id] = protein
    
    # Dump data by CBT
    for target_folder, prot_sub_col in folder_to_col_cbt.items():
        create_folder(target_folder)
        suffix = folder_to_suffix[target_folder]
        target = target_folder / "{}.fasta".format(suffix)

        with open(target, "w") as f:
            f.write(prot_sub_col.get_fasta())
        write_metadata(target_folder, "{}.fasta".format(suffix), BGCCollection(), prot_sub_col, alias)

    for target_folder, prot_sub_col in folder_to_col_dom.items():
        create_folder(target_folder)
        suffix = folder_to_suffix[target_folder]
        target = target_folder / "{}.fasta".format(suffix)
        target_meta = target_folder / "{}.fasta.metadata.tsv".format(suffix)
        dom_ID = folder_to_dom[target_folder]
        with open(target, "w") as f, open(target_meta, "w") as m:
            m.write("Protein Identifier\tProtein Id\tDomain\tCBT\n")
            for p_id, protein in prot_sub_col.proteins.items():
                # This compression gets a bit tricky because we want all domains 
                # from the protein's domain list whose ID is the one we're looking
                # for. But because the target ID can be an alias, we need to check
                # that too
                for d_num, d in enumerate([dom \
                        for dom in protein.domain_list \
                        if dom.ID == dom_ID or (dom.ID in alias and alias[dom.ID] == dom_ID)]):
                    header = "{}_{}{} ProteinId:{} GeneId:{}".format(protein.identifier, \
                        dom_ID, d_num+1, protein.protein_id, protein.gene)
                    f.write(">{}\n{}".format(header, protein.sequence80(d.ali_from, d.ali_to)))
                    m.write("{}\t{}\t{}\t{}\n".format(protein.identifier, \
                        protein.protein_id, d_num+1, protein.protein_type))

    return


def save_genbank(o: Path, cbt_types: set, cbt_domains: dict, \
        bgc_col: BGCCollection, gbk_files: dict, alias: dict) -> None:
    """
    Saves the original gbk files that were used as input. Currently, this works
    mostly to filter and organize files
    """

    if not cbt_types or cbt_domains:
        for bgc_id in bgc_col.bgcs.keys():
            try:
                original_gbk = gbk_files[bgc_id]
            except KeyError:
                # that bgc came from .bgccase or .bgc files
                continue

            if original_gbk.stem != bgc_id:
                # if a different bgc_id was given, it's probably b/c it had a generic name
                copyfile(original_gbk, o / "{}.gbk".format(bgc_id))
            else:
                copyfile(original_gbk, o / original_gbk.name)
        return

    # For output organization:
    folder_to_list = defaultdict(list)

    for bgc_id, bgc in bgc_col.bgcs.items():
        # discard bgc if not linked to physical .gbk file
        if bgc_id not in gbk_files:
            continue

        # target by type
        for cbt in (bgc.CBPtypes_set & cbt_types):
            target_folder = o / cbt
            folder_to_list[target_folder].append(bgc_id)

        # Extend set of domains of this BGC with its aliases
        internal_aliases = set() # work out internal aliases first
        for p in bgc.protein_list:
            for d in p.domain_list:
                if d.alias != "":
                    internal_aliases.add(d.alias)
        extended_bgc_domain_set = bgc.domain_set | \
            internal_aliases | \
            set(alias[d] for d in bgc.domain_set if d in alias)

        # target by domain content
        for cbt in (bgc.CBPtypes_set & cbt_domains.keys()):
            # dom_target can be a proper dom.ID or a dom.alias
            for dom_target in (extended_bgc_domain_set & cbt_domains[cbt]):
                target_folder = o / cbt / "{}_{}".format(cbt, dom_target)
                folder_to_list[target_folder].append(bgc_id)

        # special case, "all"
        if "all" in cbt_domains and (extended_bgc_domain_set & cbt_domains['all']):
            for dom_target in extended_bgc_domain_set & cbt_domains['all']:
                target_folder = o / "domains_{}".format(dom_target)
                folder_to_list[target_folder].append(bgc_id)

    # Everything classified, copy:
    for target_folder, bgc_list in folder_to_list.items():
        create_folder(target_folder)
        
        for bgc_id in bgc_list:
            original_gbk = gbk_files[bgc_id]
            if original_gbk.stem != bgc_id:
                # Probably b/c this BGC had a generic name (e.g. scaffoldXX.region...)
                copyfile(original_gbk, target_folder / "{}.gbk".format(bgc_id))
            else:
                copyfile(original_gbk, target_folder / original_gbk.name)
    return


if __name__ == "__main__":
    args = CMD_parser()
    start_time = datetime.now()

    # Verify parameters
    check_parameters(args)
    
    # Read filter list
    filter_bgc_prot = list()   # each item is a tuple (bgc, protein)
    if args.bgclist:
        filter_bgc_prot = read_bgc_list(args.bgclist)
    
    # Read alias file
    external_alias = dict()
    if args.alias_file:
        external_alias = read_alias_file(args.alias_file)

    # Add hmm databases
    hmmdbs = HMM_DB()
    hmmdbs.cores = args.cpus
    if args.hmms is not None:
        hmmdbs.add_included_database()
        for hmm in args.hmms:
            hmmdbs.add_database(Path(hmm))
    # use domain aliases given by user, if present.
    if external_alias:
        hmmdbs.alias.update(external_alias)

    # Read annotation files
    # if args.jsonfolders:
    #     print(args.jsonfolders)
    #     sys.exit("wip")

    # Output folder. Even if no files to work with are found, it needs
    #  to be created now in case there are files with repeated names/ids 
    if args.outputfolder:
        o = Path(args.outputfolder)
        if not o.is_dir():
            print("Trying to create output folder")
            os.makedirs(o, exist_ok=True) # recursive folder creation
    else:
        # There is a default value for this parameter so we should actually not
        # have this case...
        o = (Path(__file__).parent / "output")
        if not o.is_dir():
            os.makedirs(o, exist_ok=True)

    # Read input data:
    print("\nCollecting data and filtering")
    if len(args.include) > 0:
        print(" - Including only BGCs with the following:")
        print("\t{}".format("\n\t".join(args.include)))
    if len(args.exclude) > 0:
        print(" - Excluding all BGCs with the following:")
        print("\t{}".format("\n\t".join(args.exclude)))
    
    bgc_collection, protein_collection, gbk_files = get_files(args, filter_bgc_prot)
    print(" ...done")
    if len(bgc_collection.bgcs) + len(protein_collection.proteins) == 0:
        sys.exit("Stop: no valid BGCs or proteins found (check filters)")
    else:
        print("\nWorking with {} BGC(s) and {} Protein(s)".format(len(bgc_collection.bgcs), len(protein_collection.proteins)))    

    # Domain prediction stage
    if args.hmms is not None:
        print("Domain prediction stage")
        print("\tPreparing data")

        predict_domains(args, hmmdbs, bgc_collection, protein_collection)

    # Post-processing stage
    if args.merge:
        print("\nTrying to fix split domains (--merge)")
        fix_core_split_domains(bgc_collection, protein_collection)
        print("\tdone!")


    # OUTPUT stage

    # SVG output
    if args.svg:
        if args.stacked:
            print("SVG: Generating stacked figure")
            # Note: the HMM_DB object provides color data
            draw_svg_stacked(args, o, bgc_collection, protein_collection, \
                hmmdbs, filter_bgc_prot)
        else:
            print("SVG: Generating individual figures")
            draw_svg_individual(args, o, bgc_collection, protein_collection, \
                hmmdbs, filter_bgc_prot)

    # Generic output
    if args.metadata and not (args.bgccase or args.proteincase):
        write_metadata(o, args.metadata, bgc_collection, protein_collection, \
            hmmdbs.alias)

    # For the last options, check for organized output
    cbt_types, cbt_domains = get_cbt_types(args, bgc_collection, protein_collection)
    
    # Output stage
    if args.bgc or args.bgccase:
        save_bgc_output(args, o, cbt_types, cbt_domains, bgc_collection, \
            hmmdbs.alias)

    if args.proteincase:
        save_protein_output(args, o, cbt_types, cbt_domains, bgc_collection, \
            protein_collection, hmmdbs.alias)

    if args.cbt_fasta:
        save_fasta(o, cbt_types, cbt_domains, bgc_collection, \
            protein_collection, hmmdbs.alias)
        
    if args.genbank:
        save_genbank(o, cbt_types, cbt_domains, bgc_collection, gbk_files, hmmdbs.alias)
        

    end_time = datetime.now()
    print('Time: ', end_time - start_time) 
