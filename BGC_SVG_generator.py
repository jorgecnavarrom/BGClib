#!/usr/bin/env python

"""
This script will use the BGC library to predict protein domains of a GenBank
file and create svg figures of the whole cluster as well as for every protein
"""

from BGClib import *
import os
from multiprocessing import cpu_count
import pickle
import argparse


# Prepare arguments for the script
def CMD_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pfamfolder", help="Location of the Pfam-A.hmm file. \
                        If not present, location of this script will be assumed",
                        default=Path(__file__).parent)
    parser.add_argument("-i", "--inputfolder", nargs='+', help="Folder(s) to \
                        look for .gbk files (note: outputfolder does\
                        not preserve inputfolder's structure). Not recursive")
    parser.add_argument("-f", "--files", nargs='+', help="File(s) to analyze and\
                        produce figures (can be a list).")
    parser.add_argument("-o", "--outputfolder", help="Folder where results will be \
                        put (default='output')", default=(Path(__file__).parent/"output"))
    parser.add_argument("--override", help="Override previous domain data stored in\
                        the Cache folder. Default: tries to read filtered domain \
                        info from './Cache/Domain lists'", default=False, 
                        action="store_true")
    parser.add_argument("--no_domains", action="store_true", default=False, 
                        help="Don't use domain information; only make gene-arrow\
                        figures (protein-box figure will not be produced even if \
                        dbox is activated)")
    parser.add_argument("-b", "--dbox", action="store_true", default=False, 
                        help="Toggle to create individual protein domain-box\
                        type of figures")
    parser.add_argument("-a", "--arrows", action="store_true", default=False,
                        help="Toggle to create individual protein arrow figures")
    parser.add_argument("-m", "--mirror", action="store_true", default=False,
                        help="Toggle to mirror the BGC figure (does not work with\
                        the protein-box or protein-arrows figures)")
    return parser.parse_args()


if __name__ == "__main__":
    args = CMD_parser()

    override = args.override
    use_domains = not args.no_domains
    dbox = args.dbox
    arrows = args.arrows
    mirror_bgc = args.mirror
        
    # Check input parameters:
    input_files = []
    # Does not check if they are even gbk files. But biopython will surely complain
    if args.files:
        for f in args.files:
            pf = Path(f)
            if not pf.is_file():
                print("Skipping file {} (Not a file?)".format(f))
                continue
            if pf.suffix != ".gbk":
                print("Skipping file {} (not GenBank file?)".format(str(f)))
                continue
            input_files.append(Path(f))
    
    input_folders = []
    if args.inputfolder:
        for i in args.inputfolder:
            pi = Path(i)
            if not pi.is_dir():
                print("Skipping folder {} (not a folder?)".format(i))
                continue
            
            for gbk in pi.glob("*.gbk"):
                input_files.append(gbk)
    
    if len(input_files) == 0:
        sys.exit("No input files found...")
    
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
    
    # Prepare database of hmm models
    hmmdbs = HMM_DB()
    hmmdbs.cores = cpu_count()
    if use_domains:
        pfam_path = Path(args.pfamfolder)
        if not pfam_path.is_dir():
            sys.exit("Error: --pfamfolder does not point to a folder")
        pfam_path = pfam_path / "Pfam-A.hmm"        
        hmmdbs.add_database(pfam_path)
        
        hmmdbs.add_included_database()
    
    # Options for whole cluster
    svgopts = ArrowerOpts("./SVG_arrow_options.cfg")
        
    # Options for box-plot domain content figure
    if dbox:
        box_svgopts = ArrowerOpts()
        box_svgopts.scaling = 4
        box_svgopts.arrow_height = 15
        box_svgopts.stripe_thickness = 2


    # Three BGC collection objects:    
    cached_collection = BGCCollection() # From previous results. Stored in pickled file
    working_collection = BGCCollection() # BGCs that need domain prediction
    svg_collection = BGCCollection() # BGCs that will be rendered
    
    
    # always load cached data. It's easier
    data_cache = (Path(__file__).parent / "Cached_BGC_collection.pickle")
    if data_cache.exists():
        with open(data_cache, "rb") as dc:
            cached_collection = pickle.load(dc)
    for gbk in input_files:
        identifier = gbk.stem
        if identifier in cached_collection.bgcs and not override:
            svg_collection.bgcs[identifier] = cached_collection.bgcs[identifier]
        else:
            # New BGCs that we hadn't in the cached file
            bgc = BGC(gbk)
            working_collection.bgcs[bgc.identifier] = bgc
    
    if use_domains and len(working_collection.bgcs) > 0:
        print("Predicting domains for {} BGC(s)".format(len(working_collection.bgcs)))
        working_collection.predict_domains(hmmdbs, cpus=hmmdbs.cores)
        
    for b in working_collection.bgcs:
        bgc = working_collection.bgcs[b]
        for p in bgc.protein_list:
            p.classify_sequence(hmmdbs)
        
        svg_collection.bgcs[bgc.identifier] = bgc
        
        if use_domains:
            cached_collection.bgcs[bgc.identifier] = bgc
        
            
    # Update cache only if domains were predicted
    with open(data_cache, "wb") as dc:
        pickle.dump(cached_collection, dc)
    

    # Render figures
    for b in svg_collection.bgcs:
        bgc = svg_collection.bgcs[b]
        for protein in bgc.protein_list:
            protein.classify_sequence(hmmdbs)
            
        bgc_name = o / (bgc.identifier + ".svg")
        if mirror_bgc:
            bgc_name = o / (bgc.identifier + ".m.svg")
        print("Saving {}".format(bgc_name.name))
        bgc.BGC_SVG(bgc_name, hmmdb=hmmdbs, svg_options=svgopts, mirror=mirror_bgc)
        
        if dbox or arrows:
            for protein in bgc.protein_list:
                if dbox and use_domains:
                    protein.domain_SVG(o / "{}_{}_box.svg".format(bgc.identifier, protein.identifier), hmmdbs, box_svgopts)
                    
                if arrows:
                    protein.arrow_SVG(o / "{}_{}_arrow.svg".format(bgc.identifier, protein.identifier), hmmdbs, svgopts)
