#!/usr/bin/env python

"""
Handle Biosynthetic Gene Cluster objects
- Input:
    * .gb/.gbk files
    * .bgc/.bgccase files
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
import json
import argparse
import pickle
from collections import defaultdict
# import shelve
from lxml import etree
from pathlib import Path
from multiprocessing import cpu_count
from BGClib import HMM_DB, BGC, BGCCollection, ArrowerOpts, valid_CBP_types, \
    ProteinCollection, BGCDomain

__author__ = "Jorge Navarro"
__version__ = "1"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"


# Prepare arguments for the script
def CMD_parser():
    parser = argparse.ArgumentParser(description="BGC toolkit.\
        Store protein and domain information from GenBank files and generate \
        SVG figures. [v{}]".format(__version__))
    
    group_input = parser.add_argument_group("Input")
    
    group_input.add_argument("-i", "--inputfolders", nargs='+', help="Folder(s)\
        to look for .gb and .gbk files (note: output folder will not preserve \
        the structure of input folders).")
    group_input.add_argument("-f", "--files", nargs='+', help="File(s) used to \
        draw the figures (accepted: .gb .gbk, .bgc, .bgccase). Note: for\
        .bgc and .bgccase files, inclusion rules by --include, --exclude\
        and --bgclist will be applied to the internal BGC identifier, not the\
        name of the file")

    group_filtering = parser.add_argument_group("Filtering")
    group_filtering.add_argument("-l", "--bgclist", help="A file containing a list \
        of BGC identifiers (i.e. filename without .gb or .gbk extension).\
        If specified, use it to filter all the BGCs found with --inputfolders \
        or --files. \
        If --stacked is used, this list will determine the order (and filename).\
        An optional second column (separated by a tab) with the Protein ID can\
        be specified. If this column is present, the BGC will be mirrored if \
        needed such that the gene encoding the Protein ID in the second column\
        is in the forward strand. Additionally, if --stacked is used and all \
        Protein IDs are present, the corresponding gene will also be used for \
        horizontal alignment.\
        Any extra columns or rows starting with the '#' character will be \
        ignored.\
        The BGC identifiers in this file are treated case-sensitive.",
        type=Path)
    group_filtering.add_argument("--include", nargs='*', default=['region', 'cluster'], 
        help="Specify string(s) to filter which BGCs will be accepted. In the \
        case of .gb or .gbk files, the filter is applied to the filename. For \
        data stored as .bgc or .bgccase files, the filter is applied to the \
        BGC(s) identifier. If the argument is present but no parameters are \
        given, the filter will be ignored. If the argument is not present, \
        the default is to use the strings 'region' and 'cluster')")
    group_filtering.add_argument("--exclude", nargs='*', default=['final'], 
        help="Specify string(s) to filter which BGCs will be rejected. \
        Similar rules are applied as with --include. If the argument is not \
        present, the default is to use 'final'.")
    
    group_processing = parser.add_argument_group("Processing options")
    
    group_processing.add_argument("--hmm", nargs='*', help="List of paths to \
        .hmm file(s). This will also enable internal hmm models (use without \
        arguments to only use internal models).")
    group_processing.add_argument("--override", help="Use domain prediction in \
        .bgc and .bgccase files, even if they already contain domain \
        data.", default=False, action="store_true")
    group_processing.add_argument("--jsonfolders", nargs='+', type=Path, \
        help="Use MIBiG-style JSON files to annotate metadata (organism, \
        metabolites, literature)")
    group_processing.add_argument("--merge", default=False, action="store_true",
        help="Try to fix successive domains that have been split ")
    group_processing.add_argument("-c", "--cpus", type=int, default=cpu_count(), 
        help="Number of CPUs used for domain prdiction. Default: all available")
    
    group_output = parser.add_argument_group("Output")
    
    group_output.add_argument("-o", "--outputfolder", 
        default=(Path(__file__).parent/"output"), help="Folder where results \
        will be put (default='output')")
    group_output.add_argument("--svg", default=False, 
        action="store_true", help="Toggle to enable SVG output")
    group_output.add_argument("-s", "--stacked", type=str, help="If used with \
        --svg, all BGC SVGs will be put in the same figure. The argument of \
        this parameter is the filename (no extension)")
    group_output.add_argument("--gaps", default=False, action="store_true",
        help="If --stacked is used, toggle this option to leave gaps\
        when a particular BGC is not found in the input data (--bgclist).")
    group_output.add_argument("--svgcfg", 
        default=(Path(__file__).parent/"SVG_arrow_options.cfg"),
        help="Configuration file with SVG style. Default: \
        'SVG_arrow_options.cfg'")
    group_output.add_argument("-m", "--mirror", default=False, 
        action="store_true", help="Toggle to mirror each BGC figure. Ignored \
        with --stacked or --bgclist")
    group_output.add_argument("--bgc", default=False,
        action="store_true", help="Toggle to save binary files with the \
            content of each BGC (.bgc)")
    group_output.add_argument("-g", "--group", type=str, help="If used with \
        --bgc, all BGCs will be put in the same binary file (.bgccase). The \
        argument of this parameter is the filename (no extension)")
    group_output.add_argument("--metadata", type=str, help="Output a \
        tab-separated-file with metadata of the core biosynthetic proteins \
        contained in the input. Argument is the basename of the file (no \
        extension)")
    group_output.add_argument("--fasta", default=False, action="store_true",
        help="Toggle to output sequences of Core Biosynthetic Proteins \
        defined in the 'Core_Biosynthetic_Protein_fasta_options.cfg' file")
    
    return parser.parse_args()


def check_input_data(inputfiles, inputfolders, hmms, bgclist):
    """
    Checkes whether all paths to files and folders are valid
    """
    
    if not inputfiles and not inputfolders:
        sys.exit("Error: no input data. See options using -h")
        
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
            
    if hmms:
        for hmm in hmms:
            hmm_file = Path(hmm)
            if not hmm_file.is_file():
                sys.exit("Error (--hmm): {} is not a file".format(hmm_file))
            if hmm_file.suffix not in {".hmm", ".HMM"}:
                sys.exit("Error (--hmm): {} does not have a .hmm extension".format(hmm_file))

    if bgclist:
        if not Path(bgclist).is_file():
            sys.exit("Error: (--bgclist): {} is not a file".format(bgclist))


def valid_name(name, include, exclude, filter_bgc):
    """
    Checks whether a filename is valid and should be included in the analysis
    based on the allowed strings (args.include) or strings to be avoided
    (args.exclude) as well as the criterium that the BGC is included in the
    filter list.
    It is expected that the parameter 'name' is the name of the file \
    without extension
    """

    if len(include) > 0:
        if not any([word in name for word in include]):
            return False
        
    if len(exclude) > 0:
        if exclude != [] and any(word in name for word in exclude):
            return False
    
    if len(filter_bgc) > 0:
        if name not in filter_bgc:
            return False
    
    return True


def get_bgc_files(inputfolders, files, include, exclude, filter_bgc):
    """
    Reads various data sources, applies filters (strings, list)
    """
    
    input_bgc_files = dict() # Keeps a record of location of files
    collection_working = BGCCollection() # BGCs that need domain prediction
    collection_external = BGCCollection() # From .bgc or .bgccase files
                                        # They may not need domain prediction
    
    if inputfolders:
        for x in inputfolders:
            inputfolder = Path(x)
            for gb_file in inputfolder.glob("**/*.gb"):
                if valid_name(gb_file.stem, include, exclude, filter_bgc):
                    input_bgc_files[gb_file.stem] = gb_file
                    collection_working.bgcs[gb_file.stem] = BGC(gb_file)
            
            for gbk_file in inputfolder.glob("**/*.gbk"):
                if valid_name(gbk_file.stem, include, exclude, filter_bgc):
                    # If we have duplicate cluster names, use the .gbk
                    if gbk_file.stem in input_bgc_files:
                        print("Warning: substituting {} \n\twith {}".format(
                            input_bgc_files[gbk_file.stem], gbk_file))
                        input_bgc_files[gbk_file.stem] = gbk_file
                    collection_working.bgcs[gbk_file.stem] = BGC(gbk_file)
    
    if files:
        for str_file in files:
            f = Path(str_file)
            if f.suffix.lower() in {".gb", ".gbk"}:
                if valid_name(f.stem, include, exclude, filter_bgc):
                    if f.stem in input_bgc_files:
                        print("Warning: substituting {} \n\twith {}".format(input_bgc_files[f.stem], f))
                    input_bgc_files[f.stem] = f
                    collection_working.bgcs[f.stem] = BGC(f)
                    
            elif f.suffix.lower() == ".bgc":
                with open(f, "rb") as dc:
                    bgc = pickle.load(dc)
                    bgc_id = bgc.identifier
                if valid_name(bgc_id, include, exclude, filter_bgc):
                    if bgc_id in input_bgc_files:
                        print("Warning: substituting {} \n\twith {}\
                                ".format(input_bgc_files[bgc_id], f))
                        if not override:
                            del collection_working.bgcs[bgc_id]
                    
                    if override:
                        # flag this bgc to re-predict domains
                        collection_working.bgcs[bgc_id] = bgc
                    else:
                        collection_external.bgcs[bgc_id] = bgc
                    
                    input_bgc_files[bgc_id] = f
                        
            elif f.suffix.lower() == ".bgccase":
                # TODO: fix shelves to always use the same backend or find
                # another solution
                # with shelve.open(str(f), flag='r') as col:
                with open(f, "rb") as dc:
                    col = pickle.load(dc)

                    # if we've got a filter list, use it
                    if len(filter_bgc) > 0:
                        for bgc_id in filter_bgc:
                            try:
                                bgc = col[bgc_id]
                            except KeyError:
                                continue
                            
                            if valid_name(bgc_id, include, exclude, \
                                    filter_bgc):
                                if bgc_id in input_bgc_files:
                                    print("Warning: substituting {} \n\t with \
                                          bgc in collection {}".format(\
                                              input_bgc_files[bgc_id], f))
                                    if not override:
                                        del collection_working.bgcs[bgc_id]
                                
                                if override:
                                    collection_working.bgcs[bgc_id] = bgc
                                else:
                                    collection_external.bgcs[bgc_id] = bgc
                    # no filter list. Check all content in the collection
                    else:
                        for bgc_id in col.bgcs:
                            bgc = col.bgcs[bgc_id]
                            if valid_name(bgc_id, include, exclude, \
                                filter_bgc):
                                if bgc_id in input_bgc_files:
                                    print("Warning: substituting {} \n\t with \
                                            bgc in collection {}".format(\
                                                input_bgc_files[bgc_id], f))
                                    if not override:
                                        del collection_working.bgcs[bgc_id]
                                
                                if override:
                                    collection_working.bgcs[bgc_id] = bgc
                                else:
                                    collection_external.bgcs[bgc_id] = bgc
            else:
                print("Warning: unknown format ({})".format(f))

    return collection_working, collection_external


def fix_core_split_domains(output_collection):
    """
    Finds repeated domains and tries to merge them in a single one
    Currently for domains found in core biosynthetic proteins
    """

    # key: BGCProtein. 
    # Value: list whose elements are themselves lists of consecutive domain IDs
    protein_to_broken_doms = defaultdict(list) 

    # collect data
    for bgc_id in output_collection.bgcs:
        bgc = output_collection.bgcs[bgc_id]
        for cbp_type in bgc.CBPcontent:
            for protein in bgc.CBPcontent[cbp_type]:
                consecutive_domain_list = [] # consecutive doms relative to hmm
                current_id = ""
                for current_dom_num in range(len(protein.domain_list)-1):
                    current_dom = protein.domain_list[current_dom_num]
                    if current_id != current_dom.ID:
                        if consecutive_domain_list:
                            protein_to_broken_doms[protein].append(consecutive_domain_list)
                        
                        consecutive_domain_list = []
                        
                        # look ahead to see if we need to start a new list
                        next_dom = protein.domain_list[current_dom_num+1]
                        if current_dom.ID == next_dom.ID and \
                                current_dom.hmm_to < next_dom.hmm_from:
                            consecutive_domain_list.append(current_dom)

                        current_id = current_dom.ID
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
                first_dom.ali_from, last_dom.ali_to, first_dom.hmm_from, \
                last_dom.hmm_to, 0, \
                sum(d.score for d in consecutive_domain_list), 0, "")
            for d in consecutive_domain_list:
                protein.domain_list.remove(d)
            protein.domain_list.append(merged_domain)
            protein.domain_list.sort(key=lambda x:x.ali_from)


def draw_svg_individual(o, svg_collection, svgopts, hmmdbs, mirror, filter_bgc):
    for bgc_id in svg_collection.bgcs:
        bgc = svg_collection.bgcs[bgc_id]
        
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


def draw_svg_stacked(filename, svg_collection, svgopts, hmmdbs, gaps, filter_bgc, filter_bgc_order):
    """
    Draws a stacked BGC figure
    """
    thickness = svgopts.gene_contour_thickness
    
    # Obtain a final list of bgc_ids in the order that they must be printed
    if len(filter_bgc_order) > 0:
        # provided list of BGCs. Not all of them were necessarily found
        draw_order = filter_bgc_order
    else:
        # No list provided. Include everything in the order it was read
        draw_order = [*svg_collection.bgcs] # convert iterable dict's keys into list
    
    # Now obtain a list of protein objects that will be used for alignment
    # Also get information about mirroring and offset distancing
    scaling = svgopts.scaling
    H = svgopts.arrow_height # used for loci spacing
    needs_mirroring = dict()
    bgc_lengths = dict()
    bgc_distance_to_target = dict()
    target_to_start_max_offset = 0
    for bgc_id in draw_order:
        try:
            bgc = svg_collection.bgcs[bgc_id]
        except KeyError:
            print(" Warning: Cannot find BGC {} in input data".format(bgc_id))
            continue
        
        # second term in the addition: inter-loci spacing element
        L = sum([locus.length/scaling for locus in bgc.loci]) \
            + H * (len(bgc.loci)-1) \
            + thickness
        bgc_lengths[bgc_id] = L
        
        if bgc_id in filter_bgc:
            if filter_bgc[bgc_id] == "":
                print(" Warning (--bgclist): {} has not reference Protein Id".format(bgc_id))
                needs_mirroring[bgc_id] = False
                bgc_distance_to_target[bgc_id] = -1
            else:
                # try to find protein with given protein_id
                pid = filter_bgc[bgc_id]
                target_protein = None
                for locus_num in range(len(bgc.loci)):
                    locus = bgc.loci[locus_num]
                    for protein in locus.protein_list:
                        if pid == protein.protein_id or pid == protein.identifier:
                            needs_mirroring[bgc_id] = not protein.forward
                            target_protein = protein
                            if protein.forward:
                                # lenghts of each loci until target + 
                                # inter-loci spacing + distance of current \
                                # locus to target_protein
                                target_to_start = \
                                    sum([locus.length/scaling for locus in bgc.loci[:locus_num]]) \
                                    + H * locus_num \
                                    + protein.cds_regions[0][0]/scaling 
                                bgc_distance_to_target[bgc_id] = target_to_start
                            else:
                                target_to_start = \
                                    sum([locus.length/scaling for locus in bgc.loci[locus_num+1:]]) \
                                    + H * (len(bgc.loci) - locus_num - 1) \
                                    + (locus.length - protein.cds_regions[-1][1])/scaling
                                bgc_distance_to_target[bgc_id] = target_to_start
                                
                            if target_to_start > target_to_start_max_offset:
                                target_to_start_max_offset = target_to_start
                            break
                        
                # typo in Protein Id?
                if target_protein == None:
                    print(" Warning (--bgclist): cannot find reference Protein Id {} for {}".format(pid, bgc_id))
                    needs_mirroring[bgc_id] = False
                    bgc_distance_to_target[bgc_id] = -1
        else:
            needs_mirroring[bgc_id] = False
            bgc_distance_to_target[bgc_id] = -1
    
    # obtain max_L considering all the starting offsets
    max_L = 0
    for bgc_id in bgc_distance_to_target:
        if bgc_distance_to_target[bgc_id] == -1:
            max_L = max(max_L, bgc_lengths[bgc_id])
        else:
            max_L = max(max_L, bgc_lengths[bgc_id] \
                + target_to_start_max_offset \
                - bgc_distance_to_target[bgc_id])

    # Start SVG internal structure
    bgc_height = 2*svgopts.arrow_height # one for the arrow, 0.5 + 0.5 for the head height
    inner_bgc_height = bgc_height + thickness
    base_attribs = {"version":"1.1", 
                    "baseProfile":"full", 
                    "width":str(int(max_L))
                    }
    root = etree.Element("svg", attrib=base_attribs, nsmap={None:'http://www.w3.org/2000/svg'})
    
    # Add each figure
    Yoffset = 0
    rows = 0
    for bgc_id in draw_order:
        Yoffset = rows * inner_bgc_height
        try:
            bgc = svg_collection.bgcs[bgc_id]
        except KeyError:
            if gaps:
                rows += 1
            continue
        
        # Marked BGCs with no reference Protein Id won't have offset
        if bgc_distance_to_target[bgc_id] == -1:
            Xoffset = 0
        else:
            Xoffset = target_to_start_max_offset - bgc_distance_to_target[bgc_id]
        root.append(bgc.xml_BGC(Xoffset, Yoffset, hmmdbs, svgopts, needs_mirroring[bgc_id]))
        rows += 1
            
        #Yoffset = rows * inner_bgc_height
    
    # Now that we now how many BGCs were actually drawn, add height property
    root.attrib["height"] = str(int(thickness + rows*(bgc_height+thickness)))
    
    # Write SVG
    with open(filename, "bw") as f:
        f.write(etree.tostring(root, pretty_print=True))
        

def read_cbp_cfg(cfg_file):
    cbp_types = set()

    with open(cfg_file) as f:
        for line in f:
            if line[0] == "#":
                continue
            bgctype, decision = line.strip().split("=")[:2]
            if decision.strip().lower() == "true":
                cbp_types.add(bgctype.strip())
    return cbp_types


def make_fasta_files(o, requested_cbp_types, output_collection, dom_alias):
    """
    Based on options found in Core_Biosynthetic_Protein_fasta_options.cfg,
    this function extracts protein sequences of core biosynthetic proteins.
    It optionally also extracts selected domains (A, C and KS)
    """
    
    extract_A_domain = False
    if "NRPS-A-Domains" in requested_cbp_types:
        extract_A_domain = True
        proteins_w_A_domains = ProteinCollection()
        requested_cbp_types.remove("NRPS-A-Domains")
    extract_C_domain = False
    if "NRPS-C-Domains" in requested_cbp_types:
        extract_C_domain = True
        proteins_w_C_domains = ProteinCollection()
        requested_cbp_types.remove("NRPS-C-Domains")
    extract_KS_domain = False
    if "PKS-KS-Domain" in requested_cbp_types:
        extract_KS_domain = True
        proteins_w_KS_domains = ProteinCollection()
        requested_cbp_types.remove("PKS-KS-Domain")

    types_with_A_domain = {"NRPS", "NRPS-like", "PKS-NRPS_hybrid", \
        "PKS-mmNRPS_hybrid", "NRPS-PKS_hybrid"}
    types_with_C_domain = {"NRPS", "PKS-NRPS_hybrid", "PKS-mmNRPS_hybrid"}
    types_with_KS_domain = {"nrPKS", "rPKS", "other_PKS", \
        "PKS-NRPS_hybrid", "PKS-mmNRPS_hybrid", "NRPS-PKS_hybrid"}

    # Make sure user didn't change types in cfg file
    working_cbp_types = requested_cbp_types & set(valid_CBP_types)
    
    # populate protein collections for each target type
    # Basically
    # valid_CBP_types >= required_cbp_types >= types actually in data
    proteins_by_type = dict()
    for bgc_id in output_collection.bgcs:
        bgc = output_collection.bgcs[bgc_id]
        for cbp_type in (bgc.CBPtypes_set & working_cbp_types):
            # p = protein object with target type
            for p in bgc.CBPcontent[cbp_type]:
                try:
                    proteins_by_type[cbp_type].proteins[p.identifier] = p
                except KeyError:
                    proteins_by_type[cbp_type] = ProteinCollection()
                    proteins_by_type[cbp_type].proteins[p.identifier] = p

                # TODO: optimize this
                if extract_A_domain and cbp_type in types_with_A_domain:
                    proteins_w_A_domains.proteins[p.identifier] = p
                if extract_C_domain and cbp_type in types_with_C_domain:
                    proteins_w_C_domains.proteins[p.identifier] = p
                if extract_KS_domain and cbp_type in types_with_KS_domain:
                    proteins_w_KS_domains.proteins[p.identifier] = p

    # write sequences
    for cbp_type in working_cbp_types:
        try:
            prot_col = proteins_by_type[cbp_type]
        except KeyError:
            continue
        
        folder = o / cbp_type
        if not folder.is_dir():
            os.makedirs(folder, exist_ok=True)

        fasta_file = folder / (cbp_type + ".fasta")
        with open(fasta_file, "w") as f:
            f.write(prot_col.get_fasta())

        # write domain composition file
        domain_composition_file = folder / "{}_domains.tsv".format(cbp_type)
        with open(domain_composition_file, "w") as f:
            for p_id in sorted(prot_col.proteins):
                protein = prot_col.proteins[p_id]
                organism = ""
                if protein.parent_cluster.organism is not None:
                    organism = protein.parent_cluster.organism.fullname
                f.write("{}\t{}\t{}\t{}\t{}\n".format(protein.identifier, \
                    protein.protein_id, protein.gene, \
                    protein.domain_string(dom_alias), organism))

    # write domain subsequences
    if len(proteins_w_A_domains.proteins) > 0:
        folder = o / "All_A_domains"
        if not folder.is_dir():
            os.makedirs(folder, exist_ok=True)
        
        subsequences = []
        metadata = []
        for pid in proteins_w_A_domains.proteins:
            p = proteins_w_A_domains.proteins[pid]
            d_num = 1
            for d in p.domain_list:
                if d.ID == "AMP-binding":
                    header = "{}_A{} ProteinId:{} GeneId:{}".format(p.identifier, \
                        d_num, p.protein_id, p.gene)
                    subsequences.append(">{}\n{}".format(header, \
                        p.sequence80(d.ali_from, d.ali_to)))
                    metadata.append((p.identifier, p.protein_id, \
                        str(d_num), p.protein_type))
                    d_num += 1
        with open(folder / "All_A_domains.fasta", "w") as f:
            f.write("".join(subsequences))
        with open(folder / "All_A_domains.metadata.tsv", "w") as f:
            f.write("Internal protein identifier\tAnnotated protein id\tDomain number\tProtein type\n")
            for m in metadata:
                f.write("{}\n".format("\t".join(m)))

    if len(proteins_w_C_domains.proteins) > 0:
        folder = o / "All_C_domains"
        if not folder.is_dir():
            os.makedirs(folder, exist_ok=True)

        subsequences = []
        metadata = []
        for pid in proteins_w_C_domains.proteins:
            p = proteins_w_C_domains.proteins[pid]
            d_num = 1
            for d in p.domain_list:
                if d.ID == "Condensation":
                    header = "{}_C{} ProteinId:{} GeneId:{}".format(p.identifier, \
                        d_num, p.protein_id, p.gene)
                    subsequences.append(">{}\n{}".format(header, \
                        p.sequence80(d.ali_from, d.ali_to)))
                    metadata.append((p.identifier, p.protein_id, \
                        str(d_num), p.protein_type))
                    d_num += 1
        with open(folder / "All_C_domains.fasta", "w") as f:
            f.write("".join(subsequences))
        with open(folder / "All_C_domains.metadata.tsv", "w") as f:
            f.write("Internal protein identifier\tAnnotated protein id\tDomain number\tProtein type\n")
            for m in metadata:
                f.write("{}\n".format("\t".join(m)))

    if len(proteins_w_KS_domains.proteins) > 0:
        folder = o / "All_KS_domains"
        if not folder.is_dir():
            os.makedirs(folder, exist_ok=True)
        
            subsequences = []
            metadata = []
            for pid in proteins_w_KS_domains.proteins:
                p = proteins_w_KS_domains.proteins[pid]
                # d_num = 1 # should not have more than 1 KS domain
                for d in p.domain_list:
                    if d.ID == "ketoacyl-synt":
                        header = "{}_KS ProteinId:{} GeneId:{}".format(p.identifier, \
                            p.protein_id, p.gene)
                        subsequences.append(">{}\n{}".format(header, \
                            p.sequence80(d.ali_from, d.ali_to)))
                        # d_num += 1
                        metadata.append((p.identifier, p.protein_id, \
                            p.protein_type))
        with open(folder / "All_KS_domains.fasta", "w") as f:
            f.write("".join(subsequences))
        with open(folder / "All_KS_domains.metadata.tsv", "w") as f:
            f.write("Internal protein identifier\tAnnotated protein id\tProtein type\n")
            for m in metadata:
                f.write("{}\n".format("\t".join(m)))



if __name__ == "__main__":
    args = CMD_parser()

    # Verify user typed paths correctly
    check_input_data(args.files, args.inputfolders, args.hmm, args.bgclist)

    # Get parameters
    outputfolder = args.outputfolder
    override = args.override
    mirror = args.mirror
    gaps = args.gaps
    
    # Read filter list
    filter_bgc = dict()
    filter_bgc_order = list() # Dictionaries should keep order in Python
                            # 3.something but let's make sure
    if args.bgclist:
        with open(args.bgclist) as f:
            for line in f:
                if line[0] == "#":
                    continue
                
                xline = line.strip().split("\t")
                filter_bgc_order.append(xline[0])
                if len(xline) > 1:
                    # assume second column is protein_id
                    filter_bgc[xline[0]] = xline[1]
                else:
                    filter_bgc[xline[0]] = ""
        
        if len(filter_bgc) == 0:
            sys.exit("Error: filter BGC list given but the file is empty...")

    # Style options
    if args.svg:
        svgopts = ArrowerOpts(args.svgcfg)
    
    # Add hmm databases
    hmmdbs = HMM_DB()
    if (args.svg and svgopts.draw_domains and args.hmm) or \
        (args.bgc and args.hmm is not None):
        hmmdbs = HMM_DB()
        hmmdbs.cores = args.cpus
        hmmdbs.add_included_database()
        
        for hmm in args.hmm:
            hmmdbs.add_database(Path(hmm))
            
    # Read annotation files
    if args.jsonfolders:
        print(args.jsonfolders)
        sys.exit("wip")

    # Read input data:
    print("Collecting data")
    if len(args.include) > 0:
        print(" - Including only BGCs with the following:")
        print("\t{}".format("\n\t".join(args.include)))
    if len(args.exclude) > 0:
        print(" - Excluding all BGCs with the following:")
        print("\t{}".format("\n\t".join(args.exclude)))
    print("")
    collection_working, collection_external = get_bgc_files(args.inputfolders, 
                                                            args.files, 
                                                            args.include, 
                                                            args.exclude, 
                                                            filter_bgc)
    total_bgcs = len(collection_working.bgcs) + len(collection_external.bgcs)
    if total_bgcs == 0:
        sys.exit("No valid files were found")
    else:
        print("Working with {} BGC(s)".format(total_bgcs))
    
    #  Output folder
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
            
    # Ready to start
    if args.hmm is not None:
        print("Predicting domains...")
        collection_working.predict_domains(hmmdbs, cpus=hmmdbs.cores)
        print("\tdone!")
    
    output_collection = BGCCollection() # BGCs that will be rendered/processed
    output_collection.bgcs.update(collection_working.bgcs)
    output_collection.bgcs.update(collection_external.bgcs)
    print("Applying classification rules")
    output_collection.classify_proteins(args.cpus)
    print("\tdone!")

    if args.merge:
        print("\nTrying to fix split domains (--merge)")
        fix_core_split_domains(output_collection)
        print("\tdone!")

    if args.bgc:
        if args.group:
            filename = o / "{}.bgccase".format(args.group)
            with open(filename, "wb") as c:
                pickle.dump(output_collection, c)
            # with shelve.open(str(filename)) as bgc_shelf:
            #     for bgc_id in output_collection.bgcs:
            #         bgc_shelf[bgc_id] = output_collection.bgcs[bgc_id]
        else:
            for bgc_id in output_collection.bgcs:
                with open(o / "{}.bgc".format(bgc_id), "wb") as b:
                    pickle.dump(output_collection.bgcs[bgc_id], b)

    if args.metadata:
        with open(o / "{}.metadata.tsv".format(args.metadata), "w") as m:
            m.write("BGC\tDefinition\tantiSMASH products\tCore Biosynthetic Protein content\tCore Biosynthetic Protein IDs\tCore Biosynthetic Protein Identifiers\n")
            for bgc_id in sorted(output_collection.bgcs):
                bgc = output_collection.bgcs[bgc_id]
                list_core_types = []
                list_protein_ids = []
                list_protein_identifiers = []
                all_domains = set() # TODO remove this. Just for testing
                for protein in bgc.protein_list:
                    all_domains.update(protein.domain_set)
                    if protein.role == "biosynthetic":
                        list_core_types.append(protein.protein_type)
                        list_protein_ids.append(protein.protein_id)
                        list_protein_identifiers.append(protein.identifier)
                m.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bgc.identifier, 
                    bgc.definition, ", ".join(bgc.products), 
                    ", ".join(list_core_types), 
                    ", ".join(list_protein_ids), 
                    ", ".join(list_protein_identifiers)))

    if args.fasta:
        cfg_file = Path(__file__).parent/"CBP_fasta_options.cfg"
        if not cfg_file.is_file():
            print("Error (--fasta): mising configuration file for fasta extraction")
            print("(Core_Biosynthetic_Protein_fasta_options.cfg)")
        else:
            requested_cbp_types = read_cbp_cfg(cfg_file)
            make_fasta_files(o, requested_cbp_types, output_collection, hmmdbs.alias)

    if args.svg:
        if args.stacked:
            print("Generating stacked figure")
            filename = o / "{}.svg".format(args.stacked)
            draw_svg_stacked(filename, output_collection, svgopts, hmmdbs, gaps, filter_bgc, filter_bgc_order)
        else:
            print("Generating individual figures")
            draw_svg_individual(o, output_collection, svgopts, hmmdbs, mirror, filter_bgc)
        

