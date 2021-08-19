#!/usr/bin/env python

"""
BGC library

Classes for data storage and analysis of Biosynthetic Gene Clusters

"""

import sys
from pathlib import Path
from subprocess import PIPE, Popen
from multiprocessing import Pool, cpu_count
from random import uniform
from colorsys import hsv_to_rgb
from colorsys import rgb_to_hsv
from lxml import etree
from operator import itemgetter
from copy import deepcopy
from io import StringIO
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

__author__ = "Jorge Navarro"
__version__ = "0.7.4"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"


# Static data

# TODO: make a single variable BGCcpus that will be used across all methods
global BGC_cpus


# Don't include Acyl_transf_1 domain as some NRPS have it
PKS_domains = {"SAT", "ketoacyl-synt", "Ketoacyl-synt_C", "KAsynt_C_assoc",
                "TIGR04532"}
reducing_domains = {"PKS_ER", "KR", "PS-DH"}
PKS3_domains = {"Chal_sti_synt_N", "Chal_sti_synt_C"}
NRPS_domains = {"Condensation", "AMP-binding", "AMP-binding_C"}
NRPS_Independent_Siderophore_domains = {"IucA_IucC"}
Terpene_meroterpenoid_domains = {"mero_tc"}
Terpene_diterpene_domains = {"diterpene_tc"}
Terpene_triterpene_domains = {"SQHop_cyclase_N", "SQHop_cyclase_C"}
Terpene_sesquiterpene_domains = {"TRI5", "Terpene_syn_C_2"}
Terpene_sesquiterpene_bifunc_domains = {"Terpene_syn_C_2", "polyprenyl_synt"}
Terpene_squalene_domains = {"SQS_PSY"}
Terpene_carotenoid_domains = {"TIGR03462"}
Terpene_UbiA_domains = {"UbiA"}
Other_terpene_domains = {"Terpene_synth", "Terpene_synth_C", "Lycopene_cycl",
    "Prenyltrans"}
Terpene_domains =  Terpene_meroterpenoid_domains | Terpene_diterpene_domains \
    | Terpene_triterpene_domains | Terpene_sesquiterpene_domains \
    | Terpene_sesquiterpene_bifunc_domains | Terpene_squalene_domains \
    | Terpene_UbiA_domains | Other_terpene_domains
                   
DMATS_domain = {"Trp_DMAT"}

# Precursors
FAS_domains_A = {"Fas_alpha_ACP" ,"FAS_I_H", "ACPS"}
FAS_domains_B = {"DUF1729", "FAS_meander", "MaoC_dehydrat_N", "MaoC_dehydratas"}
precursor_domains = FAS_domains_A | FAS_domains_B

# TODO: add PT domain in this set until we have a better model?
hmmdbs_without_tc = {"FNP_terpene_models"}

# TODO: choose colors for the last ones
valid_CBP_types_fungal = {\
        "nrPKS": "#76b7f4", # blue
        "rPKS": "#2c9cdc", # darker blue
        "t3PKS": "#3368c8", #007dfb", # a bit more dark blue
        "NRPS": "#ffc755", # orange
        "NRPS-like": "#ffeb87", # light orange / yellow
        "other_PKS": "#80d8eb", # another blue. Intense-ish 
        "unknown_PKS": "#aaffff", # lighter version of previous
        "PKS-NRPS_hybrid": "#aa007f", # purple
        "PKS-mmNRPS_hybrid": "#9a116f", # darkish purple
        "NRPS-PKS_hybrid": "#a25fe6", # violet
        "NIS": "#c20c28", # red blood
        "Meroterpenoid_synthase": "#b42d2d", # shades of brown/dark orange
        "Diterpene_synthase": "#ba5454", # shades of brown/dark orange
        "Triterpene_synthase": "#e04d37", # shades of brown/dark orange
        "Sesquiterpene_synthase": "#995358", # shades of brown/dark orange
        "Sesquiterpene_bifunctional_synthase": "#f4f4fc", # shades of brown/dark orange
        "Carotenoid_synthase": "#f4815c", # shades of brown/dark orange
        "Squalene_synthase": "#be480a", # shades of brown/dark orange
        "UbiA-type_terpene": "#da746e", # shades of brown/dark orange
        "Terpene_other": "#ca8a73", # shades of brown/dark orange
        "DMATS": "#c5c55d" # "gold"
    }

# antiSMASH list retrieved 2020-08-19
valid_CBP_types_antiSMASH = {\
        'T1PKS': '#34b4eb',
        'T2PKS': '#999999',
        'T3PKS': '#3368c8',
        'transAT-PKS': '#999999',
        'transAT-PKS-like': '#999999', # "Trans-AT PKS fragment, with trans-AT domain not found"
        'PpyS-KS': '#999999',
        'hglE-KS': '#999999',
        'CDPS': '#cccdac', # "tRNA-dependent cyclodipeptide synthases"
        'PKS-like': '#999999',
        'arylpolyene': '#cdc5c0',
        'resorcinol': '#999999',
        'ladderane': '#cdc69f',
        'PUFA': '#999999',
        'nrps': '#999999',
        'nrps-like': '#999999',
        'thioamide-NRP': '#999999',
        'terpene': '#b44c3a',
        'lanthipeptide': '#5b9950', # Obsolete: split into subclasses in antiSMASH6
        'lipolanthine': '#999999',
        'bacteriocin': '#70b598',
        'betalactone': '#b8cdc5',
        'thiopeptide': '#90b442',
        'linaridin': '#999999',
        'cyanobactin': '#999999',
        'glycocin': '#999999',
        'LAP': '#999999',
        'lassopeptide': '#999999',
        'sactipeptide': '#999999',
        'bottromycin': '#999999',
        'head_to_tail': '#999999',
        'microviridin': '#31b568',
        'proteusin': '#999999',
        'blactam': '#999999',
        'amglyccycl': '#999999',
        'aminocoumarin': '#999999',
        'siderophore': '#c20c28',
        'ectoine': '#999999',
        'butyrolactone': '#b1b0bc',
        'indole': '#999999',
        'nucleoside': '#999999',
        'phosphoglycolipid': '#999999',
        'melanin': '#999999',
        'oligosaccharide': '#999999',
        'furan': '#999999',
        'hserlactone': '#bca8ad',
        'phenazine': '#999999',
        'phosphonate': '#f6f9d4',
        'fused': '#999999',
        'PBDE': '#999999',
        'acyl_amino_acids': '#999999',
        'tropodithietic-acid': '#999999',
        'NAGGN': '#999999',
        'RaS-RiPP': '#999999',
        'fungal-RiPP': '#3cb5a1',
        'TfuA-related': '#999999',
        'other': '#999999',
        'saccharide': '#999999',
        'fatty_acid': '#999999',
        'halogenated': '#999999'
    }
valid_CBP_types_antiSMASH_set = set(valid_CBP_types_antiSMASH.keys())


valid_CBP_types_colors = valid_CBP_types_fungal
valid_CBP_types_colors.update(valid_CBP_types_antiSMASH)

valid_CBP_types = valid_CBP_types_fungal.keys() | valid_CBP_types_antiSMASH
# other colors
# "unknown": "#f4f4fc", # super light lilac
# "other": "#fcdcdc", # very light pink
# no_domains": "#ffffff", # white. For stuff like RiPPs (TODO)

role_colors = {"biosynthetic":"#f06c6e", # red, rgb(240, 108, 110)
               "tailoring":"#8fc889", # green, rgb(143, 200, 137)
               "transport":"#f0d963", # yellow, rgb(240, 217, 99)
               "regulatory":"#33c1f0", # blue, rgb(51, 193, 240)
               "other":"#eff0f1", # light gray, rgb(239, 240, 241)
               "precursor":"#9797dc", # lilac, rgb(151,151,220)
               "unknown":"#eff0f1", # gray, rgb(220, 220, 220) dcdcdc
               "resistance":"#f0a1ac", # light red, rgb(240, 161, 172) 
               "biosynthetic-additional":"#f0986b"} # orange rgb(240, 152, 107)

# Auxiliary functions
def random_color_tuple(h_, s_, v_):
    """
    returns a random color in hex, with the hue, saturation and value being
    (possibly) bound by input parameters
    
    Additional info:
    https://en.wikipedia.org/wiki/HSL_and_HSV
    and http://stackoverflow.com/a/1586291
    
    input:
    h_, s_ and v_ are tuples of lower/upper numbers for hue/saturation/value
    """
    
    h = uniform(h_[0], h_[1])
    s = uniform(s_[0], s_[1])
    v = uniform(v_[0], v_[1])
    
    return "#{}".format("".join( hex(int(c * 255))[2:] for c in hsv_to_rgb(h, s, v) ))



# Classes definition
class HMM_DB:
    """
    This class keeps information about HMM databases to be used by the other
    classes
    """
    
    def __init__(self):
        self.db_list = []           # list of paths to hmm databases
        self.alias = dict()         # ID to alias
        self.colors = dict()        # ID to tuple(r,g,b)
        self.colors_hex = dict()
        self.color_outline = dict() # Darker version of self.colors. Precalculated
        self.colors_outline_hex = dict()
        self.cores = 0              # for hmmer. Remember that it always uses an
                                    # extra core for reading the database
                                    
        self.ID_to_AC = {}          # ID: short name; AC: accession number (w/
        self.ID_to_DE = {}          # version); DE: description
        self.ID_to_LENG = {}        
        
        self.ID_to_role = {}        # manually assigned roles (see role_colors)
        self.domain_signature_to_protein_type = {}  # curated using annotated 
                                                    # proteins from MIBiG
                                                    # Domain Signature = tilde-
                                                    # separated domain_IDs
        
        self.read_domain_colors(Path(__file__).parent / "data/domain_color_file_ID.tsv")
        self.read_domain_roles(Path(__file__).parent / "data/SM_domain_roles.tsv")
        self.read_protein_types(Path(__file__).parent / "data/protein_types.tsv")
        self.read_domain_alias_file(Path(__file__).parent / "data/domains_alias.tsv")
        
        return
    

    def add_included_database(self):
        """
        Reads hmm profiles included in the library
        """
        
        for hmm in (Path(__file__).parent).glob("data/Domain_models/*.hmm"):
            self.add_database(hmm)
        for hmm in (Path(__file__).parent).glob("data/Domain_models/*.HMM"):
            self.add_database(hmm)
            
        return
    
    
    def add_database(self, db_path):
        #TODO rename db_path to db_hmm
        """
        Adds a database of hmm models. It also stores information about each
        individual domain (linking its ID to Accession and Description)
        
        For the latter, a local file will be kept in the same location as the 
        hmm user database. This has the effect that, in order to work on previous
        data (e.g. a pickled file with domains already predicted), and ensure
        that the data hasn't changed (e.g. the pfam accession numbers), she will
        have to point to the exact same hmm database.
        
        
        input:
            db_path: Path object that points to .hmm file
        """
        
        if not db_path.is_file():
            print("Not able to add hmm database (not a file. Wrong path?): {}".format(str(db_path)))
            return False
        elif db_path.suffix.lower() != ".hmm":
            print("Not able to add hmm datase (not a .hmm file)")
            return False
        # make sure database is already "pressed" for hmmscan
        try:
            assert Path(db_path.parent / (db_path.name + ".h3i")).is_file()
        except AssertionError:
            command = ["hmmpress", str(db_path)]
            print(" Pressing hmm file with command {}".format(" ".join(command)))
            proc_hmmpress = Popen(command, shell=False)
            proc_hmmpress.wait()
            
        try:
            assert Path(db_path.parent / (db_path.name + ".h3i")).is_file()
        except AssertionError:
            sys.exit("Not able to hmmpress the database file")
            
        """
        check if there is a file in the same location as the database that we
        have previously stored (*.domain_info.tsv) 
        if there isn't:
        check if in the hmmdb folder there is a .dat file (e.g. Pfam-A.hmm.dat)
        if there is, 
            read it for AC, DE
        else
            read the whole Pfam-A.hmm file
        
        """
        
        # save the domain info file in the same place as the database so it's 
        # accessible to the user
        domain_info_file = db_path.parent / (db_path.name + ".domain_info.tsv")
        db_ID_to_AC = {}
        db_ID_to_DE = {}
        
        if domain_info_file.is_file():
            with open(domain_info_file, "r") as dif:
                for line in dif:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    ID, AC, DE = line.strip().split("\t")
                    
                    db_ID_to_AC[ID] = AC
                    db_ID_to_DE[ID] = DE
        else:
            dat_file = db_path.parent / (db_path.name + ".dat")
            # If the .dat file is found, it's much faster read
            if dat_file.is_file():
                print("\tFound {} file".format(db_path.name + ".dat"))
                with open(dat_file, "r") as dat:
                    putindict = False
                    
                    for line in dat:
                        if line[5:7] == "ID":
                            ID = line.strip()[10:]
                        if line[5:7] == "AC":
                            AC = line.strip()[10:]
                        if line[5:7] == "DE":
                            DE = line.strip()[10:]
                            putindict = True
                            
                        if putindict:
                            putindict = False
                            db_ID_to_AC[ID] = AC
                            db_ID_to_DE[ID] = DE
                            
            # Have to read the complete file. This will take a few seconds...
            else:
                with open(db_path, "r") as pfam:
                    print("\tReading domain info from {} file".format(db_path.name))
                    putindict = False
                    # assuming that the order of the information never changes
                    for line in pfam:
                        if line[:4] == "NAME":
                            ID = line.strip()[6:]
                        if line[:3] == "ACC":
                            #AC = line.strip()[6:].split(".")[0]
                            AC = line.strip()[6:]
                        if line[:4] == "DESC":
                            DE = line.strip()[6:]
                            putindict = True
                            
                        if putindict:
                            putindict = False
                            db_ID_to_AC[ID] = AC
                            db_ID_to_DE[ID] = DE

            # TODO domain_info_file is empty...?
            with open(domain_info_file, "w") as dif:
                for ID in db_ID_to_AC:
                    AC = db_ID_to_AC[ID]
                    DE = db_ID_to_DE[ID]
                    dif.write("{}\t{}\t{}\n".format(ID, AC, DE))
        
        self.ID_to_AC.update(db_ID_to_AC)
        self.ID_to_DE.update(db_ID_to_DE)
        self.db_list.append(db_path)
        
        return True


    def read_domain_alias_file(self, alias_file):
        try:
            with open(alias_file, "r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    line = line.split("\t")
                    hmm_ID = line[0]
                    hmm_alias = line[1]
                    
                    self.alias[hmm_ID] = hmm_alias
        except FileNotFoundError:
            print("Could not open domain alias file ({})".format(alias_file))

        
    def read_domain_colors(self, colors_file, append=True):
        """
        Expects the input file to be in a tsv with two columns:
        hmm_ID \t r,g,b
        
        input:
            colors_file: a file with domain ID \t rgb colors
            append: toggle to false to reset all colors
        """
        if not append:
            self.colors.clear()
            self.color_outline.clear()
        
        try:
            with open(colors_file, "r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    hmm_ID, colors = line.strip().split("\t")
                    rgb = [int(c) for c in colors.split(",")]
                    r, g, b = rgb
                    self.colors[hmm_ID] = tuple(rgb)
                    self.colors_hex[hmm_ID] = f"#{r:02x}{g:02x}{b:02x}"
                    
                    h, s, v = rgb_to_hsv(r/255.0, g/255.0, b/255.0)
                    self.color_outline[hmm_ID] = tuple(str(int(round(c * 255))) for c in hsv_to_rgb(h, s, 0.8*v))
                    rgb_out = tuple(int(round(c * 255)) for c in hsv_to_rgb(h, s, 0.8*v))
                    self.colors_outline_hex[hmm_ID] = f"#{rgb_out[0]:02x}{rgb_out[1]:02x}{rgb_out[2]:02x}"

        except FileNotFoundError:
            print("Could not open domain colors file ({})".format(str(colors_file)))
        
        return


    def read_domain_roles(self, domain_roles_file, append=True):
        """
        input:
            domain_roles_file. A tsv file with two columns:
            hmm_ID \t role
            Where role is one of biosynthetic, tailoring, transport, 
            regulation, other, unknown
            
        append:
            If False, the ID_to_role will be erased. Otherwise values will only
            be overwritten for previously assigned IDs
        """
        if not append:
            self.ID_to_role.clear()
        
        try:
            with open(domain_roles_file, "r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    ID, role = line.strip().split("\t")
                    
                    if role not in role_colors:
                        print("Warning. Found unknown role when loading domain role file:\n{}\t{}".format(ID, role))
                    
                    self.ID_to_role[ID] = role
        except FileNotFoundError:
            print("Could not open domain role file({})".format(str(domain_roles_file)))
            
        return
    
    
    # TODO is it possible to use and integrate interpro rules?
    def read_protein_types(self, protein_types_file, append=True):
        """
        input:
            protein_types_file: A tsv file with two columns:
                Domain Signature (tilde separated hmm_IDs) \t protein_type
            
            protein_type comes from Fungal MIBiG BGCs with annotations
        append:
            If False, the domain_signature_to_protein_type dictionary will be 
            erased. Otherwise values will only be overwritten for previously
            assigned values
        """
        if not append:
            self.ID_to_role.clear()
        
        try:
            with open(protein_types_file, "r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    stripped_line = line.strip().split("\t")
                    domain_signature = stripped_line[0]
                    protein_type = stripped_line[1]
                    self.domain_signature_to_protein_type[domain_signature] = protein_type
        except FileNotFoundError:
            print("Could not open protein types file({})".format(str(protein_types_file)))
            
        return



class ArrowerOpts:
    """
    Options for Arrower-like figures. Colors are linked with domains so see the
    HMM_DB class for them.
    """
    
    def __init__(self, cfg=""):
        self.scaling = 20                    # px per bp
        
        self.arrow_height = 40               # H: Height of arrows' body
        #self.arrow_head_height = 15          # h: Additional width of arrows' head
                                        #  (i.e. total arrow head = H + 2*aH)
                                        # internally, h = H/2 to preserve 45° angle
                                        
        #self.arrow_head_length = 30        # Distance from start of head till end.
        
        self.gene_contour_thickness = 2      # note: thickness grows outwards
        
        self.internal_domain_margin = 3
        self.domain_contour_thickness = 1
        
        self.stripe_thickness = 3
                
        self._color_mode = "white"
        self.valid_color_modes = {"white", "gray", "random-pastel", "random-dark", 
                                  "random", "roles", "domains", "none"}
        
        self.outline = True
        self.show_domains = True
        self.show_introns = True # Show the actual gap where the intron is

        self.original_orientation = True # If false, all arrows will point forward
        
        if cfg != "":
            self.load_options(cfg)

        return

    @property
    def color_mode(self):
        return self._color_mode
    @color_mode.setter
    def color_mode(self, cm):
        if cm in self.valid_color_modes:
            self._color_mode = cm
        else:
            print("Color mode not supported; defaulting to 'white'")
            print("Valid color modes are ['{}']".format("', '".join(self.valid_color_modes)))
            self._color_mode = "white"


    def load_options(self, cfg):
        cfgf = Path(cfg)
        
        if not cfgf.is_file():
            print("Error: Not possible to load options file for ArrowerOpts ({})".format(str(cfg)))
            return False
        
        for line in open(cfgf, "r"):
            if line[0] == "#" or line.strip() == "":
                continue
            
            try:
                option, value = line.strip().split("=")
            except ValueError:
                print("ArrowerOpts configuration file: ignoring bad line ({})".format(line.strip()))
                continue
            else:
                option = option.strip().lower()
                value = value.strip()
                
                truefalse_errors = []
                
                if option == "scaling":
                    self.scaling = int(value)
                elif option == "arrow_height":
                    self.arrow_height = int(value)
                elif option == "gene_contour_thickness":
                    self.gene_contour_thickness = int(value)
                elif option == "internal_domain_margin":
                    self.internal_domain_margin = int(value)
                elif option == "domain_contour_thickness":
                    self.domain_contour_thickness = int(value)
                elif option == "stripe_thickness":
                    self.stripe_thickness = int(value)
                elif option == "color_mode":
                    value = value.replace("'", "").replace('"', '').lower()
                    if value in self.valid_color_modes:
                        self.color_mode = value
                elif option in {"outline", "show_domains", "show_introns",
                        "original_orientation"}:
                    value = value.capitalize()
                    if value in {"True","False"}:
                        if value == "True":
                            value = True
                        else:
                            value = False
                            
                        if option == "outline":
                            self.outline = value
                        elif option == "show_domains":
                            self.show_domains = value
                        elif option == "original_orientation":
                            self.original_orientation = value
                        elif option == "show_introns":
                            self.show_introns = value
                    else:
                        truefalse_errors.append(option)
                else:
                    pass
                    # print("ArrowerOpts configuration file: unknown option {}".format(option))
                    
                if len(truefalse_errors) > 0:
                    print("ArrowerOpts configuration file: the following options must have a value of True or False {}".format(option))
                    
        return True
    


class BGCCollection:
    """
    This class will allow implementation of collection-wide functions such as 
    single-step prediction of domains and comparisons between collections
    """
    
    def __init__(self):
        self.bgcs = {}
        self.name = ""


    def __len__(self):
        return len(self.bgcs)
    

    def add_gbk(self, gbk, identifier=""):
        """
        Creates a BGC from a gbk file and adds it to the collection
        """
        
        bgc = BGC(gbk, identifier)
        if bgc.identifier in self.bgcs:
            print("Warning: {} already in BGC Collection. Skipping".format(bgc.identifier))
        else:
            self.bgcs[bgc.identifier] = bgc

    
    def predict_domains(self, hmmdb, domtblout_path="", cpus=1, tc=True, filterdoms=True):
        """
        Compile the protein sequences of the whole collection and apply hmmscan
        on them to predict domains. Assign predicted domains to each protein
        and filter
        """
        assert(isinstance(hmmdb, HMM_DB))
        if domtblout_path != "":
            assert(isinstance(domtblout_path, Path))
        
        pc = ProteinCollection()
        
        missing_identifier = False
        for bgc in self.bgcs.values():
            # bgc = self.bgcs[b]
            for protein in bgc.protein_list:
                if protein.identifier != "":
                    pc.proteins[protein.identifier] = protein
                else:
                    missing_identifier = True
        
        if missing_identifier:
            print("Warning: one or more proteins don't have a valid identifier")
            
        pc.predict_domains(hmmdb, domtblout_path, cpus, tc, filterdoms)
        
        for bgc in self.bgcs.values():
            bgc.attempted_domain_prediction = True
            bgc.calculate_domain_sets()

        return


    def classify_proteins(self, cpus=1):
        """
        Once domains have been predicted, use them to assign protein type
        and protein role for all proteins in the collection
        """
        
        with Pool(cpus) as pool:
            for bgc in self.bgcs.values():
                pool.apply_async(bgc.classify_proteins())
            pool.close()
            pool.join()
        return


    def clear_domain_predictions(self, cpus=1):
        """
        Remove all domain annotations from the proteins in this collection
        """

        with Pool(cpus) as pool:
            for bgc in self.bgcs.values():
                pool.apply_async(bgc.clear_domain_predictions())
            pool.close()
            pool.join()
    

    def clear_protein_roles(self):
        """
        Removes all protein role annotations from the collection.
        Note that this may affect annotations imported from the original
        GenBank files
        """


        for bgc in self.bgcs.values():
            for protein in bgc.protein_list:
                protein.role = "unknown"



class BGC:
    def __init__(self, gbkfile=None, identifier=""):
        self.identifier = identifier        # usually the file name
        
        self.CBPtypes = []          # Core Biosynthetic Protein List: simple
                                    #  5'-3' ordered list of biosynthetic types
        self.CBPtypes_set = set()
        self.CBPcontent = {}        # Core Biosynthetic Protein Content. Every 
                                    #  biosynthetic-type returns a list of 
                                    #  pointers to each Protein object
        
        self.products = set()       # All different 'product' qualifiers  
                                    #  annotated by antiSMASH
        self.contig_edge = False    # antiSMASH v4+ was not able to fully 
                                    #  complete the extension phase. BGC *might*
                                    #  be fragmented
    
        self.protein_list = []      # should also be present in loci
        self.proteins = {}          # direct access with protein identifier

        self.domain_set = set()
        self.domain_set_core = set()
        self.domain_set_complement = set()

        self.loci = []

        # These will be linked to Metadata Objects
        self.metabolites = []
        self.organism = None

        self.attempted_domain_prediction = False
        
        # try to initialize object with GenBank content
        if gbkfile is not None:
            self.load(gbkfile, identifier)
        
        
    def load(self, _gbk, _id):
        """
        Initializes the object by reading all CDS fields from a GenBank file
        
        Input is intended to be a Path() object (but it should work with a string)
        """
        
        if isinstance(_gbk, Path):
            gbk = _gbk
        else:
            gbk = Path(_gbk)
        
        if self.identifier == "":
            self.identifier = gbk.stem
        
        try:
            records = list(SeqIO.parse(str(gbk), "genbank"))
        except ValueError as e:
            print("Error, not able to parse file {}: {}".format(str(gbk), str(e)))
        else:
            self.accession = records[0].id
            self.definition = records[0].description

            if records[0].annotations["organism"] != "":
                self.organism = Organism()
                self.organism.fullname = records[0].annotations["organism"]
                self.organism.lineage = records[0].annotations["taxonomy"]
            
            cds_list = []
            
            # traverse all possible records in the file. There's usually only 1
            locus_num = 0
            for record in records:
                locus = BGCLocus()
                locus.length = len(record.seq)
                current_biosynthetic_region = FeatureLocation(locus.length+1, locus.length+2)
                
                cds_num = 0

                for feature in record.features:
                    # antiSMASH <= 4
                    if feature.type == "cluster":
                        if "product" in feature.qualifiers:
                            for product in feature.qualifiers["product"]:
                                for p in product.replace(" ","").split("-"):
                                    self.products.add(p)
                                
                        if "contig_edge" in feature.qualifiers:
                            if feature.qualifiers["contig_edge"][0] == "True":
                                self.contig_edge = True
                                
                        current_biosynthetic_region = feature.location
                        continue
                    
                    # antiSMASH = 5
                    if feature.type == "region":
                        if "product" in feature.qualifiers:
                            for product in feature.qualifiers["product"]:
                                self.products.add(product)
                                
                        if "contig_edge" in feature.qualifiers:
                            # there might be mixed contig_edge annotations
                            # in multi-record files. Turn on contig_edge when
                            # there's at least one annotation
                            if feature.qualifiers["contig_edge"][0] == "True":
                                self.contig_edge = True

                        current_biosynthetic_region = feature.location
                        continue
                                
                    if feature.type == "CDS":
                        cds_num += 1
                        
                        CDS = feature
                        
                        cds_start = max(0, int(CDS.location.start))
                        cds_end = max(0, int(CDS.location.end))

                        identifier = "{}~L{}+CDS{}".format(self.identifier, locus_num, cds_num)

                        product = ""
                        if "product" in CDS.qualifiers:
                            product = ", ".join(CDS.qualifiers["product"])
                            
                        # NOTE: JGI annotations have a non-standard qualifier:
                        # "proteinId" or "proteinID", which is just a number (sometimes). 
                        protein_id = ""
                        # found in NCBI:
                        if "protein_id" in CDS.qualifiers:
                            protein_id = CDS.qualifiers["protein_id"][0]
                            
                        # found in JGI
                        elif "proteinID" in CDS.qualifiers:
                            protein_id = CDS.qualifiers["proteinID"][0]
                        elif "proteinId" in CDS.qualifiers:
                            protein_id = CDS.qualifiers["proteinId"][0]

                        gene = ""
                        if "gene" in CDS.qualifiers:
                            gene = CDS.qualifiers["gene"][0]

                        role = ""
                        if "gene_kind" in CDS.qualifiers:
                            role = CDS.qualifiers["gene_kind"][0]

                        # TODO: test this. It's a bit tricky with antiSMASH 5 now with regions
                        protein_type = ""
                        if role == "biosynthetic" and "gene_functions" in CDS.qualifiers:
                            protein_types = []
                            for x in CDS.qualifiers["gene_functions"]:
                                if x.startswith("biosynthetic (rule-based-clusters)"):
                                    protein_types.append(x.split("biosynthetic (rule-based-clusters) ")[1].split(":")[0])
                            if len(set(protein_types)) == 1:
                                protein_type = protein_types[0]
                            elif "NRPS_PKS" in CDS.qualifiers:
                                protein_type = "PKS-NRPS_hybrid"
                                # TODO delete these
                                #for x in CDS.qualifiers["NRPS_PKS"]:
                                    #if x.startswith("type:"):
                                        #protein_type = x.split("type: ")[1]
                                # as they're not informative. e.g.: 
                                # /NRPS_PKS="type: Type I Iterative PKS"

                        protein = BGCProtein()
                        
                        protein.identifier = identifier
                        protein.product = product
                        protein.protein_id = protein_id
                        protein.gene = gene
                        protein.role = role
                        protein.protein_type = protein_type

                        if cds_start in current_biosynthetic_region and cds_end in current_biosynthetic_region:
                            protein.in_biosynthetic_region = True
                        
                        # TODO: check if the 'translation' annotation is really 
                        # there. If not, try to manually translate from dna
                        if "translation" not in CDS.qualifiers:
                            print(" Warning. Skipping CDS without 'translation' qualifier: {}".format(identifier))
                            continue
                        protein.sequence = CDS.qualifiers["translation"][0]
                        # TODO: verify sequence doesn't contain stop codons at the end (anywhere?)

                        # TODO: verify that len(CDS.location.parts) == len(protein.sequence)
                        
                        # NOTE what happens if there is no strand info (strand=None)?
                        if CDS.location.strand != 1:
                            protein.forward = False
                        
                        protein.parent_cluster = self
                        protein.parent_cluster_id = self.identifier
                        protein.parent_locus = locus
                        
                        # TODO: detect overlapping CDS (e.g. gene annotation 
                        #  errors, multiple splicing)
                        del cds_list[:]
                        for x in CDS.location.parts:
                            cds_list.append((int(x.start), int(x.end)))
                        protein.cds_regions = tuple(sorted(cds_list, key=itemgetter(0)))

                        self.protein_list.append(protein)
                        self.proteins[protein.identifier] = protein
                        locus.identifier = "{}~{}".format(self.identifier, locus_num)
                        locus.protein_list.append(protein)
                        locus.gene_coordinates.append((cds_start,cds_end))

                self.loci.append(locus)
                locus_num += 1
    
        self.set_CBP_content()

        return

    
    # TODO: finish
    def load_fasta(self, fasta):
        """
        Creates faux BGCs from a fasta file.
        Each entry in the fasta file will be treated as a separate locus.
        Assumes that the fasta file contains a protein sequence.
        """
        pass
    
    
    def inter_loci_element(self, xoffset, yoffset, svg_options=ArrowerOpts()):
        """
        Draws the following SVG figure as xml code to link loci genomic stripes:
        __/ /__
         / /
        
        and returns the corresponding etree element
        """
        
        H = svg_options.arrow_height
        
        inter_loci_main = etree.Element("g")
        
        inter_loci_attribs = {
            "stroke": "#464646",
            "stroke-width": str(svg_options.stripe_thickness)
            }
            
        # left
        inter_loci_attribs["d"] = "M{} {} L {} {}".format(int(xoffset), int(yoffset + 0.5*H), int(xoffset+0.375*H), int(yoffset + 0.5*H))
        inter_loci_main.append(etree.Element("path", attrib=inter_loci_attribs))
        
        # right
        inter_loci_attribs["d"] = "M{} {} L {} {}".format(int(xoffset+0.625*H), int(yoffset + 0.5*H), int(xoffset+H), int(yoffset + 0.5*H))
        inter_loci_main.append(etree.Element("path", attrib=inter_loci_attribs))
        
        # slash1
        inter_loci_attribs["d"] = "M{} {} L {} {}".format(int(xoffset + 0.25*H), int(yoffset + 0.75*H), int(xoffset + 0.5*H), int(yoffset + 0.25*H))
        inter_loci_main.append(etree.Element("path", attrib=inter_loci_attribs))
        
        # slash2
        inter_loci_attribs["d"] = "M{} {} L {} {}".format(int(xoffset + 0.5*H), int(yoffset + 0.75*H), int(xoffset + 0.75*H), int(yoffset + 0.25*H))
        inter_loci_main.append(etree.Element("path", attrib=inter_loci_attribs))
        
        return inter_loci_main
    

    # TODO: implement the extra_label parameter. The main issue is how to 
    # calculate the width of the text element to take into account for the SVG width
    def BGC_SVG(self, filename, hmmdb=HMM_DB(), svg_options=ArrowerOpts(), extra_label="", xoffset=0, yoffset=0, mirror = False):
        """
        Writes an SVG figure for the cluster
        """
        
        # Arrows need to be represented in their original orientation. 
        # Be careful as not to modify input svg options though
        if svg_options.original_orientation:
            bgc_svg_options = svg_options
        else:
            bgc_svg_options = deepcopy(svg_options)
            bgc_svg_options.original_orientation = True
        
        H = bgc_svg_options.arrow_height
        
        # length of the SVG figure
        # second term accounts for the inter-loci figure
        L = sum([locus.length/bgc_svg_options.scaling for locus in self.loci]) + (H*(len(self.loci)-1))

        # main node
        # Note: we are in dna space (so amino acid sequence must be scaled)
        # width, height, xoffset and yoffset must be corrected or the corners
        # with sharp angles will be chopped off
        thickness = bgc_svg_options.gene_contour_thickness
        Xoffset = xoffset + thickness
        Yoffset = yoffset + thickness
        base_attribs = {"version":"1.1", 
                        "baseProfile":"full",
                        "width":str(int(Xoffset + L + thickness)),
                        "height":str(int(Yoffset + 2*H + thickness))}
        root = etree.Element("svg", attrib=base_attribs, nsmap={None:'http://www.w3.org/2000/svg'})
        
        inner_tree = self.xml_BGC(Xoffset, Yoffset, hmmdb, bgc_svg_options, mirror=mirror)
        root.append(inner_tree)
        
        with open(filename, "bw") as f:
            f.write(etree.tostring(root, pretty_print=True))


    def xml_BGC(self, xoffset, yoffset, hmmdb=HMM_DB(), svg_options=ArrowerOpts(), mirror = False):
        """
        Creates a BGC SVG figure (xml structure)
        """
        main_group = etree.Element("g")
        
        bgc_title = etree.Element("title")
        bgc_title.text = self.identifier
        main_group.append(bgc_title)
        
        H = svg_options.arrow_height
        h = 0.5*H
        
        l = 1
        Xoffset = xoffset # start of each locus
        loci_list = self.loci
        if mirror:
            loci_list = reversed(self.loci)
        for locus in loci_list:
            L = locus.length/svg_options.scaling
            
            line_attribs = {
                "x1": str(int(Xoffset)),
                "y1": str(int(yoffset + h + 0.5*H)),
                "x2": str(int(Xoffset+L)),
                "y2": str(int(yoffset + h + 0.5*H)),
                "stroke": "#464646",
                "stroke-width": str(int(svg_options.stripe_thickness))
                }
            line = etree.Element("line", attrib=line_attribs)
            main_group.append(line)
            
            p = 0
            intron_offset = 0.0
            for protein in locus.protein_list:
                gene_start, gene_end = locus.gene_coordinates[p]
                
                gene_start = gene_start / svg_options.scaling
                gene_position = gene_start-intron_offset
                if mirror:
                    if svg_options.intron_regions:
                        gene_length = (protein.cds_regions[-1][1] - protein.cds_regions[0][0] - 3)/svg_options.scaling
                    else:
                        gene_length = (protein.length*3)/svg_options.scaling
                    gene_position = L - gene_position - gene_length
                inner_tree = protein.xml_arrow(hmmdb, svg_options, \
                    Xoffset + gene_position, yoffset, mirror=mirror)
                main_group.append(inner_tree)
                
                p += 1
                
            Xoffset += L
            
            if l < len(self.loci):
                # draw inter-locus fig of length H
                # NOTE: the mirror parameter will flip each individual locus, not
                # the whole collection of them. This means that the inter-loci
                # element is still drawn at the end of each locus
                main_group.append(self.inter_loci_element(Xoffset, yoffset + 0.5*H, svg_options))
            
                Xoffset += H
            l += 1
        
        return main_group


    def classify_proteins(self):
        """
        Uses predicted domains to classify (Secondary Metabolism) role and 
        protein_type of all its proteins.
        """
        
        for p in self.protein_list:
            p.classify_sequence()
        
        self.set_CBP_content()

        return

    
    def set_CBP_content(self):
        del self.CBPtypes[:]        # ordered list
        self.CBPtypes_set.clear()   # a set of all CBTs
        self.CBPcontent.clear()     # CBT to list of BGCProtein
        for protein in self.protein_list:
            if protein.role == "biosynthetic":
                self.CBPtypes.append(protein.protein_type)
                try:
                    self.CBPcontent[protein.protein_type].append(protein)
                except KeyError:
                    self.CBPcontent[protein.protein_type] = [protein]
        self.CBPtypes_set = set(self.CBPtypes)


    def predict_domains(self, hmmdb, domtblout_path="", cpus=1, tc=True, filterdoms=True):
        """
        Compile the protein sequences of the whole collection and apply hmmscan
        on them to predict domains. Assign predicted domains to each protein
        and filter
        """
        assert(isinstance(hmmdb, HMM_DB))
        if domtblout_path != "":
            assert(isinstance(domtblout_path, Path))
        
        pc = ProteinCollection()
        
        missing_identifier = False

        for protein in self.protein_list:
            if protein.identifier != "":
                pc.proteins[protein.identifier] = protein
            else:
                missing_identifier = True
        
        if missing_identifier:
            print("Warning: one or more proteins don't have a valid identifier")
            
        pc.predict_domains(hmmdb, domtblout_path, cpus, tc, filterdoms)
        
        self.attempted_domain_prediction = True
        self.calculate_domain_sets()

        return
    

    def calculate_domain_sets(self):
        for protein in self.protein_list:
            self.domain_set.update(protein.domain_set)
            if protein.role == "biosynthetic":
                self.domain_set_core.update(protein.domain_set)
            else:
                self.domain_set_complement.update(protein.domain_set)
    

    def clear_domain_predictions(self):
        """
        Remove all domain information for all the proteins in this BGC
        """

        for protein in self.protein_list:
            protein.domain_list = []
            protein.domain_set = set()
        
        self.attempted_domain_prediction = False



class BGCLocus:
    """
    BGC might be divided in more than one locus (because of sequencing, or, as
    is more often the case in fungi, because of genomics)
    
    This class can be used to organize proteins
    """
    def __init__(self):
        self.identifier = ""        # e.g. "bgc~L01"
        self.protein_list = []
        self.gene_coordinates = []  # a list of tuples. end-start might not match
                                    # the corresponding protein length due to 
                                    # splicing, but will be useful for inter-
                                    # protein region length when making arrower
                                    # figure.
        self.length = 0



class ProteinCollection:
    """
    Intended for mass-prediction of domains
    """
    
    def __init__(self, fastafile=None):
        self.name = ""              # Name of the collection
        self.proteins = {}          # key = identifier
        self.proteins_by_pid = {}   # key = protein id. Because we can't ensure
                                    #  uniqueness here, only non-empty, 
                                    #  non-repeated keys will be used, so this
                                    #  will be a subset of self.proteins
        
        if fastafile is not None:
            self.fasta_load(fastafile)
    
    
    def __len__(self):
        return len(self.proteins)

            
    def fasta_load(self, fastafile):
        """
        Tries to load a fasta file
        """
        
        # initial string
        random_string = "l5aRx]3DahVr9bnIGVVaMHwpFhOu6iML#DWDbY7DJxTEhc9F4Vw_ZoFRdcA4"
        
        with open(fastafile) as f:
            header = random_string
            sequence = ""
            for line in f:
                if line.strip() == "":
                    continue
                
                if line[0] == ">":
                    if header == "":
                        print("Warning: skipping headerless sequence")
                        print(" ({}...)".format(sequence[:50]))
                        sequence = ""
                        header = line.strip()[1:]
                        continue
                    
                    if header != random_string:
                        if header in self.proteins:
                            print("Warning: {} already in Protein Collection. \
                                  Skipping...".format(header))
                            print(" ({}...)".format(sequence[:50]))
                        else:
                            protein = BGCProtein()
                            protein.identifier = header
                            protein.sequence = sequence
                            self.proteins[header] = protein
                    sequence = ""
                    header = line.strip()[1:]
                else:
                    sequence += line.strip()
                    
        if header == "":
            print("Warning: skipping headerless sequence")
            print(" ({}...)".format(sequence[:50]))
        else:
            if header in self.proteins:
                print("Warning: {} already in Protein Collection. Skipping...".format(header))
                print(" ({})".format(sequence[:50]))
            else:
                protein = BGCProtein()
                protein.identifier = header
                protein.sequence = sequence
                self.proteins[header] = protein
        
        
    # TODO: break work on sets of 3 cpus
    # TODO: evaluate whether hmmsearch is better than hmmscan
    # TODO: TEST!
    def predict_domains(self, hmmdb, domtblout_path="", cpus=1, tc=True, filterdoms=True):
        """
        Uses hmmscan to search the protein sequences for hmm models specified
        
        input:
            hmmdb: contains a list to all hmm databases
            domtblout_path: where to store HMMER's output file (optional)
            cpus: number of cpus to pass to HMMER. Remember it uses cpu + 1
            tc: if true, use the model's Trusted Cutoff score (lower score of all
                true positives)
        """
        
        if domtblout_path != "":
            try:
                assert not domtblout_path.is_file()
            except AssertionError:
                sys.exit("BGClib.ProteinCollection.predict_domains: domtblout_path should not be a file")
            
        protein_list = []
        for protein_id in self.proteins:
            protein = self.proteins[protein_id]
            protein_list.append(">{}\n{}".format(protein.identifier, protein.sequence))
                    
        if len(protein_list) == 0:
            return
        
        # TODO: 
        #  * merge all hmmdbs and remove this loop
        #  * Divide the protein collection into groups where each group is cpus // 3
        #  * do "with Pool(groups) as pool"
        for db in hmmdb.db_list:
            #command = ['hmmscan', '--cpu', str(cpus), '--noali', '--notextw']
            command = ['hmmscan', '--cpu', str(cpus), '--notextw', '--qformat', 'fasta']
            if tc and db.stem not in hmmdbs_without_tc:
                command.append('--cut_tc')
            if domtblout_path != "":
                path = str(domtblout_path / ("output_" + db.stem + ".domtable"))
                command.extend(['--domtblout', path ])
            dbpath = str(db)
            command.extend([ dbpath, '-'])
            print(" ".join(command))

            proc_hmmscan = Popen(command, shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = proc_hmmscan.communicate(input="\n".join(protein_list).encode("utf-8"))
            
            # "SearchIO parses a search output file's contents into a hierarchy of four 
            # nested objects: QueryResult, Hit, HSP, and HSPFragment"
            # https://biopython.org/DIST/docs/api/Bio.SearchIO-module.html
            # https://biopython.org/DIST/docs/api/Bio.SearchIO.HmmerIO-module.html
            results = SearchIO.parse(StringIO(out.decode()), 'hmmer3-text')
            #TODO: parallelize this loop
            for qresult in results:
                # if the original header of the fasta sequence contains a space,
                # this thing takes everything after the space as the "description"
                # and the first part as the "id" but probably this isn't what 
                # we want
                if qresult.description != "<unknown description>":
                    seq_identifier = "{} {}".format(qresult.id, qresult.description)
                else:
                    seq_identifier = qresult.id
                                
                for hit in qresult:
                    #print(hit.description)
                    for hsp in hit:
                        hspf = hsp[0] # access to HSPFragment

                        #seq_identifier = qresult.id
                        hmm_id = hit.id
                        try:
                            AC = hmmdb.ID_to_AC[hmm_id]
                        except KeyError:
                            AC = ""

                        try:
                            DE = hmmdb.ID_to_DE[hmm_id]
                        except KeyError:
                            DE = ""
                        
                        try:
                            alias = hmmdb.alias[hmm_id]
                        except KeyError:
                            alias = ""
                        
                        # NOTE: these are 1-based indices but biopython takes
                        # care of this and shifts automatically (e.g. hmmscan
                        # position 1 transforms to array index 0)
                        ali_from = hspf.query_start
                        ali_to = hspf.query_end
                        hmm_from = hspf.hit_start
                        hmm_to = hspf.hit_end
                        #env_from = hsp.env_start
                        #env_to = hsp.env_end
                        
                        # NOTE hspf.query is a SeqRecord object
                        algn_seq = hspf.query.seq
                        
                        # NOTE: The alignment of the query sequence vs the hmm
                        # model can be obtained directly from the hmmscan output.
                        # In order to get same-sized alignments (allowed amino 
                        # acids + deletions in the query), we need to remove all
                        # insertions (of the query, with the model as reference. 
                        # marked as lowercase symbols in the query).
                        # For this we could use the '#=GC RF' line from a proper
                        # Stockholm alignment (e.g. using hmmalign) as a mask. 
                        # hmmscan does not report that line as part of its output
                        # but we can use the "hit" sequence to find the positions
                        # of the insertions (from the perspective of the hmm with
                        # the query as reference, these will be viewed as deletions
                        # which are marked with periods in the hmm sequence).
                        # A second option would be to filter the lowercase
                        # symbols in the sequence
                        #print(hspf.query.seq)
                        #print(hspf.hit.seq)
                        #print(hspf.aln_span, len(hspf.hit.seq))
                        
                        hit_seq = hspf.hit.seq
                        ## this is probably super expensive...
                        #hit_seq_split = hit_seq.split(".")
                        
                        #x = "xxx.x.xxxx.x...xx".split(".")
                        #print([len(a) for a in x])
                        #sys.exit()
                        
                        
                        hmm_size = len(hit_seq)
                        
                        Evalue = hsp.evalue
                        score = hsp.bitscore
                        
                        domain = BGCDomain(self.proteins[seq_identifier], \
                                           hmm_id, AC, DE, alias, \
                                            ali_from, ali_to, 
                                           hmm_from, hmm_to, hmm_size, \
                                           score, Evalue, algn_seq)
                        
                        self.proteins[seq_identifier].domain_list.append(domain)
        
        with Pool(cpus) as pool:
            for p in self.proteins:
                protein = self.proteins[p]
                if filterdoms:
                    pool.apply_async(protein.filter_domains())
                protein.attempted_domain_prediction = True
            pool.close()
            pool.join()
    
        return


    # # TODO: finish this. Challenge is target now is the query
    # def predict_domains_search(self, hmmdb, domtblout_path="", cpus=1, tc=True):
    #     """
    #     Uses hmmsearch to search the protein sequences for hmm models specified
        
    #     input:
    #         hmmdb: contains a list to all hmm databases
    #         domtblout_path: where to store HMMER's output file (optional)
    #         cpus: number of cpus to pass to HMMER. Remember it uses cpu + 1
    #         tc: if true, use the model's Trusted Cutoff score (lower score of all
    #             true positives)
    #     """
        
    #     if domtblout_path != "":
    #         try:
    #             assert not domtblout_path.is_file()
    #         except AssertionError:
    #             sys.exit("BGClib.ProteinCollection.predict_domains: domtblout_path should not be a file")
            
    #     protein_list = []
    #     for protein_id in self.proteins:
    #         protein = self.proteins[protein_id]
    #         protein_list.append(">{}\n{}".format(protein.identifier, protein.sequence))
                    
    #     if len(protein_list) == 0:
    #         return
        
    #     for db in hmmdb.db_list:
    #         command = ['hmmsearch', '--cpu', str(cpus), '--notextw']
    #         if tc:
    #             command.append('--cut_tc')
    #         if domtblout_path != "":
    #             path = str(domtblout_path / ("output_" + db.stem + ".domtable"))
    #             command.extend(['--domtblout', path ])
    #         dbpath = str(db)
    #         command.extend([ dbpath, '-'])
            
    #         proc_hmmscan = Popen(command, shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    #         out, err = proc_hmmscan.communicate(input="\n".join(protein_list).encode("utf-8"))
            
    #         # "SearchIO parses a search output file's contents into a hierarchy of four 
    #         # nested objects: QueryResult, Hit, HSP, and HSPFragment"
    #         # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html
    #         results = SearchIO.parse(StringIO(out.decode()), 'hmmer3-text')
    #         for qresult in results:
    #             #print(qresult)
    #             for hit in qresult:
    #                 for hsp in hit:
    #                     hspf = hsp[0] # access to HSPFragment

    #                     seq_identifier = qresult.id
    #                     hmm_id = hit.id
    #                     ali_from = hspf.query_start
    #                     ali_to = hspf.query_end
    #                     hmm_from = hspf.hit_start
    #                     hmm_to = hspf.hit_end
    #                     env_from = hsp.env_start
    #                     env_to = hsp.env_end
                        
    #                     Evalue = hsp.evalue
    #                     score = hsp.bitscore
                        
    #                     domain = BGCDomain(self.proteins[seq_identifier], 
    #                                        hmm_id, env_from, env_to, ali_from, 
    #                                        ali_to, hmm_from, hmm_to, score, 
    #                                        Evalue, algn_seq)
                        
    #                     self.proteins[seq_identifier].domain_list.append(domain)
        
    #     with Pool(cpus) as pool:
    #         for p in self.proteins:
    #             protein = self.proteins[p]
    #             pool.apply_async(protein.filter_domains())
    #             protein.attempted_domain_prediction = True
    #         pool.close()
    #         pool.join()
    
    #     return
    

    def classify_proteins(self, cpus=1):
        """
        Once domains have been predicted, classify proteins from the whole 
        collection
        
        Note that directly classifying a protein will not update BGC.CBPcontent,
         BGC.CBPtypes and BGC.CBPtypes_set of its parent BGC
        """
        
        with Pool(cpus) as pool:
            for protein in self.proteins.values():
                pool.apply_async(protein.classify_sequence())
            pool.close()
            pool.join()
        return
    
    
    def get_fasta(self, sort=True):
        """
        Returns a string with the fasta sequences of all proteins in the
        collection. Header will be the protein accession, if present. Otherwise
        the protein's identifier will be used. 
        """
        sequences = []
        if sort:
            seq_list = sorted(self.proteins.keys())
        else:
            seq_list = self.proteins.keys()
        for p_id in seq_list:
            p = self.proteins[p_id]
            header = "{} ProteinId:{} GeneId:{}".format(p.identifier, \
                p.protein_id, p.gene)

            sequences.append(">{}\n{}".format(header, p.sequence80()))

        return "".join(sequences)



class BGCProtein:
    """
    Information about a Protein encoded by a gene which is part of the gene cluster
    """
        
    def __init__(self):
        self.parent_cluster = None  # Should point to an object of the BGC class
        self.parent_cluster_id = ""
        self.parent_locus = None    # Should point to an object of the BGCLocus class
        
        self.identifier = ""        # Internal identifier. Should be unique:
                                    # BGCname~L#+CDS#

        self.protein_id = ""        # As annotated in GenBank
        self.gene = ""              # As annotated in GenBank
        
        #self._accession = ""        # original NCBI accession if possible
        self._ref_accession = ""    # "RefSeq Selected Product" from NCBI's IPG DB
        self.ncbi_id = ""           # NCBI internal id in the Protein database
        
        self.ncbi_ipg_id = ""       # NCBI internal id in the Identical Protein 
                                    # Group database. Should be unique
        
        # NOTE: other NCBI annotations:
        # /db_xref="GeneID:1095902"
        # /gene_synonym="SCF76.19c"
        # /locus_tag="SCO0479"
        
        self._sequence = ""
        self.length = 0             # automatically calculated when sequence is 
                                    # added. Calculated in amino acid space
        
        self.protein_type = ""      # Protein type, e.g.: nrPKS, FAS_A
        self.product = ""           # From NCBI's GenBank annotations
        self.role = "unknown"       # e.g. [biosynthetic, transporter, tailoring, 
                                    #   resistance, unknown]. See role_colors
        self.in_biosynthetic_region = False # From e.g. antiSMASH predictions
        
        
        # self.compound_family = ""   # Compound family e.g. "emodin-like"
        # self.compound = ""
        # self.source = "unknown"     # e.g. MIBiG, curated document etc.
        # self.organism = ""
        # self.TaxId = ""
        
        self.forward = True         # orientation of the gene encoding this protein
        
        self.cds_regions = tuple()  # collection of start/stop CDS regions (to 
                                    # draw intron position).
                                    # These come from the original GenBank file
                                    # so are relative to the start of the locus.
                                    # They are in nucleotide space.
        
        self.domain_list = []       # list of BGCDomain objects, ordered by ali_from
        self.domain_set = set()     # set of unique domain IDs
        self.attempted_domain_prediction = False
        
        self.metabolites = []       # scaffold for PKSs, NRPSs

        return
    
           
    @property
    def ref_accession(self):
        return self._ref_accession
    @ref_accession.setter
    def ref_accession(self, ra):
        self._ref_accession = ra
        if self.identifier == "":
            self.identifier = ra

    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, seq):
        # TODO: check if contains valid amino acid characters
        self._sequence = seq.replace("\n", "").strip().replace(" ", "").upper()
        self.length = len(seq)
        
        # cds_regions not yet imported from genbank. Do this now to be able to
        # draw standalone SVG arrow figures
        if len(self.cds_regions) == 0:
            self.cds_regions = ([0, 3*self.length], )
    

    def sequence80(self, start=0, end=None):
        """
        Returns a prettified sequence (broken at 80 chars.)
        0-based coordinates (i.e. 'start' coordinate is included but 'end' is
        not)
        """

        if end == None: # NOTE parameter can't be 'end=self.length' (check)
            end = self.length
            
        if start >= end or start < 0 or end > self.length:
            return ""
        
        length = end - start
        seq = self._sequence[start:end]
        
        part_one = "\n".join([seq[row*80:(row+1)*80] for row in \
            range(length // 80)])
        
        part_two = ""
        remainder = length % 80
        if remainder > 0:
            part_two = "{}\n".format(seq[-remainder:])
        return "{}\n{}".format(part_one, part_two)
        
        
    # TODO rename to "get_fasta"?
    def fasta(self, start=0, end=None):
        """
        Returns a fasta-formatted data
        Tries to use reference accession for header. If not present, falls back
        to original accession. If either are absent, use identifier (should
        always be present)
        """
        try:
            assert(type(start) == int)
        except AssertionError:
            return ""
        
        if end == None:
            end = self.length
        
        molecules = ""
        if self.metabolites:
            metabolite_list = []
            for metabolite in self.metabolites:
                if metabolite.alias:
                    metabolite_list.append(metabolite.alias[0])
                else:
                    metabolite_list.append(metabolite.name)

            molecules = ",".join(m.replace("/", "") for m in metabolite_list)
            
        return ">{}|{}|{}|{}\n{}".format(self.identifier, self.protein_id, \
            self.gene, molecules, self.sequence80(start, end))


    def domain_string(self, domain_alias, original_orientation=False, simple=True):
        """Returns a basic domain-organization string
        
        Hierarchy for each domain's label:
        - If present, use alias from external domain_alias. If not,
        - If present, use the original alias recorded in the domain object. If not,
        - Uses the hmm model ID (not accession number)
        
        Does NOT check if the domains are already ordered by start coordinate
        
        in:
            domain_alias (optional): external dictionary for domain ID aliases
        out: 
            a string
        """

        # get names; try to convert them to alias
        # First try to get internal alias. Then try external alias dictionary
        domain_alias_list = []
        for d in self.domain_list:
            label = d.ID

            try:
                label = domain_alias[d.ID]
            except KeyError:
                if d.alias != "":
                    label = d.alias
                else:
                    label = d.ID

            domain_alias_list.append(label)
        
        # tag = self.identifier
        # if self.organism != "":
        #     tag += "_[{}]".format(self.organism)
        # if self.compound != "":
        #     tag += "_[{}]".format(self.compound)

        if original_orientation and not self.forward:
            if simple:
                return "< {}".format(" | ".join(reversed(domain_alias_list)))
            else:
                return "<--[{}]--|".format("]-[".join(reversed(domain_alias_list)))
        else:
            if simple:
                return "{} >".format(" | ".join(domain_alias_list))
            else:
                return "|--[{}]-->".format("]-[".join(domain_alias_list))
        
        
    def recursive_interval(self, interval_list):
        """every interval in interval_list is a tuple of 
        (original index, score, start, end)
        
        Returns: a list of the original indices that need to be deleted using the 
        score as defining criteria
        
        Inspired by the Interval Tree algorithm
        
        Sort by score, descending.
        Take strongest member (first). It will not be deleted.
        Populate 3 lists:
        - elements at the right of strongest:   right_list
        - elements at the left of strongest:    left_list
        - elements overlapping with strongest:  center_list
        - analyze center_list against strongest:
            If weaker contains strongest, or strongest contains weaker (complete
                overlap), mark for deletion
            If overlap is small, (less than or equal to 5% of the shortest 
                between weaker and strongest), allow it to survive. Put in either
                left_list or right_list
        - Repeat until you only have a single element
        """
        if len(interval_list) < 2:
            return []
        
        interval_list.sort(key=lambda x:x[1], reverse=True)
        
        left_list = []
        right_list = []
        center_list = []
        
        # as the list is already sorted, the first interval should have the best
        # score
        strongest = interval_list[0]
        for interval in interval_list[1:]:
            if interval[3] <= strongest[2]:
                left_list.append(interval)
            elif interval[2] >= strongest[3]:
                right_list.append(interval)
            else:
                # some overlap
                center_list.append(interval)
        
        delete_indices = []
        # analyze overlapping list. They are all deleted unless the overlap is small
        for interval in center_list:
            # start is at the left
            if interval[2] <= strongest[2]:
                if interval[3] >= strongest[3]:
                    # interval is longer or equal than strongest on both sides
                    delete_indices.append(interval[0])
                else:
                    # check if overlap < 5% of the shortest interval
                    if interval[3]-strongest[2] > 0.05*min(interval[3]-interval[2], strongest[3]-strongest[2]):
                        delete_indices.append(interval[0])
                    else:
                        # has to compete with other intervals at the left
                        left_list.append(interval)
            else:
                if interval[3] <= strongest[3]:
                    # interval is embedded in strongest, delete
                    delete_indices.append(interval[0])
                else:
                    # check if overlap < 5% of the shortest interval
                    if strongest[3]-interval[2] > 0.05*min(interval[3]-interval[2], strongest[3]-strongest[2]):
                        delete_indices.append(interval[0])
                    else:
                        # has to compete with other intervals at the right
                        right_list.append(interval)
                        
        delete_indices.extend(self.recursive_interval(left_list))
        delete_indices.extend(self.recursive_interval(right_list))
            
        return delete_indices
        
        
    # TODO merge splitted domains
    def filter_domains(self):
        """
        Filters domains which overlap with other, higher-scoring ones
        Also re-orders self.domain_list according to hmm hit _alignment_ 
        coordinates
        """
        
        # No need to filter anything if we only have one or no domains
        if len(self.domain_list) == 0:
            return
            
        if len(self.domain_list) == 1:
            self.domain_set = {self.domain_list[0].ID}
            return
        
        interval_list = []
        for idx in range(len(self.domain_list)):
            domain = self.domain_list[idx]
            
            interval_list.append((idx, domain.score, domain.ali_from, domain.ali_to))
            
        deletion_list = sorted(self.recursive_interval(interval_list), reverse=True)
    
        for d in deletion_list:
            del self.domain_list[d]
            
        # sort by env_from
        self.domain_list.sort(key=lambda x:x.ali_from)
        
        # finally, set up the domain set
        self.domain_set = set([x.ID for x in self.domain_list])
        
        return
    
    
    # TODO substitute for a call to ProteinCollection.predict_domains
    def predict_domains(self, hmmdb, domtblout_path="", cpus=1, tc=True, filterdoms=True):
        """
        domtblout_path is the path where the domtableout file will be deposited. 
        If present, a Path-like object is expected
        """

        assert(isinstance(hmmdb, HMM_DB))
        if domtblout_path != "":
            assert(isinstance(domtblout_path, Path))
        
        pc = ProteinCollection()
        pc.proteins[self.identifier] = self
        pc.predict_domains(hmmdb, domtblout_path, cpus, tc, filterdoms)
        
        self.classify_sequence(hmmdb)
        
        return
    
    
    def classify_sequence(self, hmmdb = HMM_DB()):
        """Classifies a sequence based on its predicted domains according to a 
        set of rules.
        
        This will be a work in progress as new rules get incorporated.
        
        """

        # Immediately return to preserve potential pre-existing annotations
        # (e.g. from antiSMASH) # TODO: use antiSMASH's predictions
        if len(self.domain_set) == 0:
            return

        self.protein_type = "unknown"
        self.role = "unknown"
        
        cbp_type = ""
        
        # pks/nrps hybrid
        if self.domain_set & precursor_domains:
            self.role = "precursor"
            
            if self.domain_set & FAS_domains_A & FAS_domains_B:
                self.protein_type = "Fatty Acid Synthase"
            elif self.domain_set & FAS_domains_A:
                self.protein_type = "Fatty Acid Synthase A"
            elif self.domain_set & FAS_domains_B:
                self.protein_type = "Fatty Acid Synthase B"
            else:
                # shouldn't really happen
                self.protein_type = "unknown_precursor"
            return
        elif (self.domain_set & PKS_domains) and (self.domain_set & NRPS_domains):
            for d in self.domain_list:
                if d.ID in PKS_domains:
                    # Detect PKS/multi-modular NRPS. Count both C and A domains 
                    # just in case there is a broken domain counted twice
                    A_domains = len([a.ID for a in self.domain_list if a.ID == "AMP-binding"])
                    C_domains = len([c.ID for c in self.domain_list if c.ID == "Condensation"])
                    if A_domains > 1 and C_domains > 1:
                        self.protein_type = "PKS-mmNRPS_hybrid"
                    else:
                        self.protein_type = "PKS-NRPS_hybrid"
                    self.role = "biosynthetic"
                    
                    return
                elif d.ID in NRPS_domains:
                    self.protein_type = "NRPS-PKS_hybrid"
                    self.role = "biosynthetic"
                    return
                else:
                    pass
        # nrPKS or (h/p)rPKS
        elif self.domain_set & PKS_domains:
            # try to set appart FAS-like sequences of e.g. azaphilone
            if "ketoacyl-synt" not in self.domain_set:
                cbp_type = "unknown_PKS"
                self.role = "unknown"
                
            elif self.domain_set & {"TIGR04532","SAT"}:
                # assume that having a SAT domain is enough for labeling as nrPKS
                # but note that in this category, there seem to be three cases for PT:
                # - PT detected by TIGR04532
                # - DH detected by pfam PF14765 (TIGR doesn't pass --cut_tc)
                # - no detection
                # There are either several types of PT domains or not a general model
                # to detect all three+ cases
                cbp_type = "nrPKS"
                self.role = "biosynthetic"
            elif self.domain_set & reducing_domains:
                cbp_type = "rPKS"
                self.role = "biosynthetic"
            else:
                cbp_type = "other_PKS"
                self.role = "biosynthetic"
        elif self.domain_set & PKS3_domains:
            cbp_type = "t3PKS"
            self.role = "biosynthetic"
        elif self.domain_set & NRPS_domains:
            if "Condensation" in self.domain_set:
                cbp_type = "NRPS"
            else:
                cbp_type = "NRPS-like"
            self.role = "biosynthetic"
        elif self.domain_set & NRPS_Independent_Siderophore_domains:
            cbp_type = "NIS"
            self.role = "biosynthetic"
        elif self.domain_set & Terpene_domains:
            self.role = "biosynthetic"
            if self.domain_set & Terpene_meroterpenoid_domains:
                cbp_type = "Meroterpenoid_synthase"
            elif self.domain_set & Terpene_diterpene_domains:
                cbp_type = "Diterpene_synthase"
            elif self.domain_set & Terpene_triterpene_domains:
                cbp_type = "Triterpene_synthase"
            elif self.domain_set & Terpene_sesquiterpene_domains:
                cbp_type = "Sesquiterpene_synthase"
            elif self.domain_set & Terpene_sesquiterpene_bifunc_domains:
                if Terpene_sesquiterpene_bifunc_domains <= self.domain_set:
                    cbp_type = "Sesquiterpene_bifunctional_synthase"
                elif "Terpene_syn_C_2" in self.domain_set:
                    cbp_type = "Terpene_other"
                else:
                    self.role = "precursor"
                    cbp_type = "Polyprenyl transferase"
                    return
            elif self.domain_set & Terpene_squalene_domains:
                if self.domain_set & Terpene_carotenoid_domains:
                    cbp_type = "Carotenoid_synthase"
                else:
                    cbp_type = "Squalene_synthase"
            elif self.domain_set & Terpene_UbiA_domains:
                cbp_type = "UbiA-type_terpene"
            else:
                cbp_type = "Terpene_other"
        elif self.domain_set & DMATS_domain:
            cbp_type = "DMATS"
            self.role = "biosynthetic"
        else:
            domain_signature = "~".join(x.ID for x in self.domain_list)
            try:
                cbp_type = hmmdb.domain_signature_to_protein_type[domain_signature]
            except KeyError:
                cbp_type = "other"
            
            role_set = set()
            for d in self.domain_set:
                if d in hmmdb.ID_to_role:
                    role_set.add(hmmdb.ID_to_role[d])
            
            if len(role_set) == 1:
                for r in role_set:
                    self.role = r
            else:
                self.role = "unknown"

        self.protein_type = cbp_type
        
        return


    def clear_domain_predictions(self):
        self.domain_list = []
        self.domain_set = set()
        
        self.attempted_domain_prediction = False


    # TODO FIX
    def domain_SVG(self, filename, hmmdb, svg_options=ArrowerOpts()):
        """Creates an SVG figure with the domains as boxes
    
        in:
            name: name of the file including path
            hmmdb: contains alias and color information
            introns: make a dashed line signaling intron position (if there's any)
            
        out:
            an svg file
        """
    
        H = svg_options.arrow_height
        h = H/2
        stripe_thickness = svg_options.stripe_thickness
        curve_radius = int(H/10)
        scaling = svg_options.scaling
        intron_regions = svg_options.show_introns
        intron_break = False

        # main node
        base_attribs = {"version":"1.1", 
                        "baseProfile":"full",
                        "width":str(int(self.length/scaling)),
                        "height":str(int(2*h + H))}
        root = etree.Element("svg", attrib=base_attribs, nsmap={None:'http://www.w3.org/2000/svg'})
        main_group = etree.SubElement(root, "g")
        
        # add title
        title = etree.Element("title")
        title.text = self.identifier
        main_group.append(title)
        
        # genomic locus line
        genomic_locus_attribs = {"x1":"0", 
                                 "y1":str(int(h+H/2)),
                                 "x2":str(int(self.length/scaling)),
                                 "y2":str(int(h+H/2)),
                                 "style":"stroke:#0a0a0a; stroke-width:{}".format(stripe_thickness)}
        genomic_locus = etree.Element("line", attrib=genomic_locus_attribs)
        main_group.append(genomic_locus)


        for domain in self.domain_list:
            # General properties of the domain: color and title
            try:
                color = ",".join([str(x) for x in hmmdb.colors[domain.ID]])
            except KeyError:
                color = "150,150,150"
                
            title = domain.ID
            try:
                title = hmmdb.alias[domain.ID]
            except KeyError:
                if domain.alias != "":
                    title = domain.alias

            domain_title = etree.Element("title")
            domain_title.text = title
            
            domain_box_attribs = {"x":str(int(domain.ali_from/scaling)),
                                  "y":"0",
                                  "rx":str(curve_radius),
                                  "ry":str(curve_radius),
                                  "width":str(int((domain.ali_to - domain.ali_from)/scaling)),
                                  "height":str(int(2*H)),
                                  "fill":"rgb({})".format(color)}
            domain_box = etree.Element("rect", attrib = domain_box_attribs)
            
            domain_node_main = etree.Element("g")
            domain_node_main.append(domain_title)
            domain_node_main.append(domain_box)
            main_group.append(domain_node_main)
            

        if intron_break and len(self.cds_regions) > 1:
            l = 0
            if not self.forward:
                l = int(self.length / scaling)
            
            start = int(self.cds_regions[0][0]/(3*scaling))
            end = int(self.cds_regions[0][1]/(3*scaling))
            
            offset = start
            current_x = end
            for i in self.cds_regions[1:]:
                start = int(i[0]/(3*scaling))
                end = int(i[1]/(3*scaling))
                
                intron_main = etree.Element("g")
                
                # intron background
                if self.forward:
                    intron_attribs = {"d":"M {:d} 0 L {:d} {:d}".format(current_x-offset, current_x-offset, 2*H),
                                    "stroke":"#ffffff",
                                    "stroke-width":"2"}
                else:
                    # if reverse, mirror position of introns
                    intron_attribs = {"d":"M {:d} 0 L {:d} {:d}".format(l-current_x+offset, l-current_x+offset, 2*H),
                                    "stroke":"#ffffff",
                                    "stroke-width":"2"}
                intron_main.append(etree.Element("path", attrib=intron_attribs))
                
                # intron dashed line
                intron_attribs["stroke"] = "#787878"
                intron_attribs["stroke-linecap"] = "round"
                intron_attribs["stroke-dasharray"] = "4,4"
                intron_main.append(etree.Element("path", attrib=intron_attribs))
                
                # add the element to the main group
                main_group.append(intron_main)
                
                current_x += end-start
           

        with open(filename, "bw") as f:
            f.write(etree.tostring(root, standalone=True, pretty_print=True))
    

    def arrow_colors(self, mode, hmmdb=HMM_DB()):
        if mode == "white":
            return "#ffffff"
        elif mode == "random-pastel":
            s_ = (0.3, 0.5)
            v_ = (0.9, 0.95)
            return random_color_tuple((0.0, 1.0), s_, v_)
        elif mode == "random-dark":
            s_ = (0.5, 0.75)
            v_ = (0.4, 0.6)
            return random_color_tuple((0.0, 1.0), s_, v_)
        elif mode == "random":
            s_ = (0.4, 0.6)
            v_ = (0.5, 0.8)
            return random_color_tuple((0.0, 1.0), s_, v_)
        elif mode == "gray":
            return "#bababa"
        elif mode == "none":
            return "none"
        elif mode == "roles":
            try:
                return role_colors[self.role]
            except KeyError:
                return "#eff0f1"
        elif mode == "domains":
            # colors based on protein domains. 
            # If role is 'biosynthetic' or 'precursor' use color or role
            # If protein has 1 domain, use its color
            # If protein has more than one domain, use color from first one
            
            # If it's a CBP, it has a defined color already
            if self.role == "biosynthetic":
                try:
                    color = valid_CBP_types_colors[self.protein_type]
                except KeyError:
                    # unknown core gene type??
                    color = "#d59d7c"
                return color
            
            if len(self.domain_list) > 0:
                domain = self.domain_list[0]
                try:
                    color = "#{}".format("".join( hex(int(c))[2:] for c in hmmdb.colors[domain.ID] ))
                except KeyError:
                    color = "#eff0f1"
            else:
                color = "#eff0f1"
                
            return color
            
        else:
            # "white"
            return "#ffffff"
        

    def arrow_SVG(self, filename, hmmdb=HMM_DB(), svg_options=ArrowerOpts(), xoffset=0, yoffset=0):
        """Creates an SVG figure with the domains as boxes
    
        in:
            name: name of the file including path
            hmmdb: contains alias and color information
            svg_options: contains size of various elements in the figure
            x/y offset: upper left corner of the arrow. In nucleotide space
            
        out:
            an svg file
        """
        
        H = svg_options.arrow_height
        
        thickness = svg_options.gene_contour_thickness
        
        # length of the SVG figure
        L = self.length*3/svg_options.scaling
        if svg_options.intron_regions:
            L = (self.cds_regions[-1][1] - self.cds_regions[0][0])/svg_options.scaling
        
        # main node
        # Note: we are in dna space (so amino acid sequence must be scaled)
        # width, height, xoffset and yoffset must be corrected or the corners
        # with sharp angles will be chopped off
        correction = 2*thickness
        base_attribs = {"version":"1.1", 
                        "baseProfile":"full",
                        "width":str(int(xoffset + L + correction)),
                        "height":str(int(yoffset + 2*H + correction))}
        root = etree.Element("svg", attrib=base_attribs, nsmap={None:'http://www.w3.org/2000/svg'})

        inner_tree = self.xml_arrow(hmmdb, svg_options, xoffset+thickness, yoffset+thickness, mirror=False)
        root.append(inner_tree)
        
        with open(filename, "bw") as f:
            f.write(etree.tostring(root, pretty_print=True))
        
    
    def xml_arrow(self, hmmdb=HMM_DB(), svg_options=ArrowerOpts(), xoffset=0, yoffset=0, mirror=False, needs_class_data=True):
        """Creates an arrow SVG figure (xml structure)
        
        Starting from upper left corner (X,Y), the vertices of the arrow are
        (clockwise)
        A, B, C, D (tip of the arrow), E, F, G
        (not considering vertices created by intron regions)
        But note that here the first coordinate to be annotated will be the one
        in the bottom left corner (G or E, depending on whether the arrow is 
        complete or will be a squeezed version)
        
        Domains are similar, except that b and c follow the slope of the arrow
        head. Vertices:
        a, b, c, d (if domain end matches the arrow's), e, f, g
    
        in:
            options through an ArrowerOpts object
                L: total arrow length
                l: length of arrow head (45° angle)
                H: height of arrow's body (A-G, B-F)
                h: height of top/bottom arrow edge (B-C, E-F)
                
            an object of the HMM_DB class which contains domain information
            
            xoffset /yoffset: upper left corner of the canvas
            
            (svg_options.original_orientation): draw arrow according to strand information
            mirror: draw arrow with the opposite direction of the strand orientation
        out:
            an etree tree
            
        X,Y: coordinate of the upper left corner of the arrow's tail (if it's
        pointing forward). A = (X,Y) = (xoffset, yoffset+h)
        
        Note that scaling only happens on the x axis

        needs_class_data: to be implemented. Highest level call (protein, BGC, 
            BGC list) should define style classes for genes, introns and domains 
            to save even more space
        """
        
        fill_color = self.arrow_colors(svg_options.color_mode, hmmdb)
        
        if mirror:
            flip = self.forward
        else:
            flip = not self.forward and svg_options.original_orientation
        
        L = self.cds_regions[-1][1] - self.cds_regions[0][0]
        
        scaling = svg_options.scaling

        H = svg_options.arrow_height
        h = H/2
        l = H * scaling # needs to be scaled along the x-axis
        Y = h # X = 0 (xoffset)
        
        idm = svg_options.internal_domain_margin
        show_introns = svg_options.show_introns
        show_domains = svg_options.show_domains
        domain_contour_thickness = svg_options.domain_contour_thickness

        main_group = etree.Element("g")
        
        # add title
        arrow_title = etree.Element("title")
        
        arrow_identifier = self.identifier
        if self.protein_id:
            arrow_identifier = f"{self.identifier} [{self.protein_id}]"
        
        # If only one domain and color_mode == 'domains', don't draw domains even
        # if requested (instead, whole arrow will be filled with domain's color)
        core_type = ""
        if self.role == "biosynthetic":
            core_type = f"\n{self.protein_type}"

        d_info = ""
        if svg_options.color_mode == "domains":
            show_domains = False
                
            try:
                d_info_str = " + ".join([d.DE for d in self.domain_list])
            except KeyError:
                d_info_str = " + ".join([d.ID for d in self.domain_list])
            d_info = f"\n{d_info_str}"
        
        arrow_title.text = f"{arrow_identifier}{core_type}{d_info}"        
        main_group.append(arrow_title)

        # precalculate a couple of numbers
        vertices = []
        head_start = L - l # In local 'arrow coordinates'
        if head_start <= 0:
            head_start = 0
        center = Y + h
        Hl = H/l # alpha = (H*(L-start)/l) distance from center to colission with head
        HL = H/L # for shorter arrows

        #
        # DRAW MAIN GENE AREA
        vertices = list()
        if head_start == 0:
            # short arrow
            C = [0, 0]
            D = [L, center]
            E = [0, Y + H + h]
            vertices = [C, D, E]
        else:
            # full arrow
            A = [0, Y]
            B = [head_start, Y]
            C = [head_start, 0]
            D = [L, center]
            E = [head_start, Y + H + h]
            F = [head_start, Y + H]
            G = [0, Y + H]
            vertices = [A, B, C, D, E, F, G]

        # flip if requested and gene is in the reverse strand
        if flip:
            for v in vertices:
                v[0] = L - v[0]
        
        # rescale everything
        for v in vertices:
            v[0] = v[0]/scaling
        
        # shape properties
        string_vertices = []
        for v in vertices:
            string_vertices.append(f"{round(v[0]+xoffset, 2)},{round(v[1]+yoffset, 2)}")
        
        region_attribs = dict()
        region_attribs["class"] = "exon"
        region_attribs["points"] = " ".join(string_vertices)
        region_attribs["fill"] = fill_color

        if needs_class_data:
            region_attribs["stroke-width"] = str(svg_options.gene_contour_thickness)
            if svg_options.outline:
                region_attribs["stroke"] = "#0a0a0a"
                
        gene = etree.Element("polygon", attrib = region_attribs)
        main_group.append(gene)

        #
        # DRAW INTRONS
        intron_regions = list()
        intron_elements = list()
        if show_introns and len(self.cds_regions) > 1:
            # calculate regions.
            # Coordinates are relative to the start of the gene-arrow
            if self.forward:
                cds_start = self.cds_regions[0][0]
                for cds_num in range(len(self.cds_regions) - 1):
                    intron_regions.append((self.cds_regions[cds_num][1] + 1 - cds_start, 
                        self.cds_regions[cds_num+1][0] - 1 - cds_start))
            else:
                cds_end = self.cds_regions[-1][1] + 1
                for cds_num in range(len(self.cds_regions)-1, 0, -1):
                    intron_regions.append((cds_end - self.cds_regions[cds_num][0], 
                        cds_end - self.cds_regions[cds_num-1][1] - 1))

            for intron in intron_regions:
                start = intron[0]
                end = intron[1]
                vertices = list()

                if start < head_start:
                    G = [start, Y + H]
                    vertices.append(G)
                    A = [start, Y]
                    vertices.append(A)
                else:
                    if head_start == 0:
                        alpha = (HL*(L-start))
                    else:
                        alpha = (Hl*(L-start))
                    Eprime = [start, center + alpha]
                    vertices.append(Eprime)
                    Cprime = [start, center - alpha]
                    vertices.append(Cprime)

                # draw right-most part of region
                if end <= head_start:
                    Aprime = [end, Y]
                    vertices.append(Aprime)
                    Gprime = [end, Y + H]
                    vertices.append(Gprime)
                else:
                    if start < head_start:
                        B = [head_start, Y]
                        vertices.append(B)
                        C = [head_start, Y-h]
                        vertices.append(C)
                    
                    if head_start == 0:
                        alpha = (HL*(L-end))
                    else:
                        alpha = (Hl*(L-end))
                    Cprimeprime = [end, center - alpha]
                    vertices.append(Cprimeprime)
                    Eprimeprime = [end, center + alpha]
                    vertices.append(Eprimeprime)
                        
                    if start < head_start:
                        E = [head_start, Y + H + h]
                        vertices.append(E)
                        F = [head_start, Y + H]
                        vertices.append(F)

                # flip if requested and gene is in the reverse strand
                if flip:
                    for v in vertices:
                        v[0] = L - v[0]
                
                # rescale everything
                for v in vertices:
                    v[0] = v[0]/scaling
                
                # shape properties
                region_attribs = dict()
                region_attribs["stroke-width"] = str(svg_options.gene_contour_thickness)
                region_attribs["stroke-opacity"] = "0.8"
                region_attribs["class"] = "intron"
                region_attribs["fill"] = "#d8d8d8"
                region_attribs["fill-opacity"] = "0.5"
                # region_attribs["stroke-dasharray"] = "3,2"
                if svg_options.outline:
                    region_attribs["stroke"] = "#919191"

                string_vertices = []
                for v in vertices:
                    string_vertices.append(f"{round(v[0]+xoffset, 2)},{round(v[1]+yoffset, 2)}")
                region_attribs["points"] = " ".join(string_vertices)
                intron_element = etree.Element("polygon", attrib=region_attribs)

                intron_elements.append(intron_element)

        #
        # DOMAINS 
        domain_elements = list()
        if show_domains and len(self.domain_list) > 0:
            # split each domain into dna regions
            offset = 0

            # In contrast with intron_regions, elements here will
            # be a dictionary with two items: regions and ID
            domain_regions = list()

            # calculate regions.
            # Coordinates are relative to the start of the gene-arrow
            # Another option would be to use the intron regions, with the
            # caveat that if the user doesn't choose to show_introns, those
            # regions won't be calculated
            forward_regions = list()
            offset = self.cds_regions[0][0]
            if self.forward:
                for s, e in self.cds_regions:
                    forward_regions.append(tuple([s-offset, e-offset]))
            else:
                end = self.cds_regions[-1][1]
                for e, s in self.cds_regions[::-1]:
                    forward_regions.append(tuple([end-s, end-e]))

            # prepare domain regions
            for domain in self.domain_list:
                # get cDNA coordinates
                dstart = domain.ali_from * 3    # inclusive dstart
                dend = (domain.ali_to + 1) * 3  # exclusive dend

                # get genomic coordinates
                offset = 0
                regions = list()
                for current_cds, cds in enumerate(forward_regions):
                    # CDS is before domain. Push domain forward
                    if cds[1] < dstart:
                        offset = forward_regions[current_cds+1][0] - cds[1] + 1
                        dstart += offset
                        dend += offset
                    # CDS is after domain. Finish
                    elif cds[0] >= dend:
                        regions.append(tuple([dstart, dend]))
                        break
                    # CDS end splits domain region
                    elif cds[1] < dend:
                        regions.append(tuple([dstart, cds[1]+1]))
                        if current_cds < len(forward_regions) - 1:
                            offset = forward_regions[current_cds+1][0] - cds[1] + 1
                            dstart = forward_regions[current_cds+1][0]
                            dend += offset
                    # CDS overlaps completely with domains
                    else:
                        regions.append(tuple([dstart, dend]))
                        break
                domain_regions.append({
                    "ID": domain.ID,
                    "DE": domain.DE,
                    "AC": domain.AC,
                    "alias": domain.alias,
                    "regions": regions})

            if head_start == 0:
                arrow_collision = L - ( (h-idm)/HL )
            else:
                arrow_collision = head_start + (h+idm)/Hl
            
            # "linker half height"
            lhh = 0.125 * H
            
            # Work every domain
            for domain in domain_regions:
                ID = domain["ID"]
                DE = domain["DE"]
                alias = domain["alias"]
                AC = domain["AC"]
                regions = domain["regions"]

                vertices = list()
                vertices_top = list()
                vertices_bottom = list()

                # draw each domain segment
                for r_num, (dbox_start, dbox_end) in enumerate(regions):
                    # Case i) Full rectangle    
                    if dbox_end <= arrow_collision:
                        # Points for linker. Left
                        if r_num > 0:
                            x = [dbox_start, center - lhh]
                            vertices_top.append(x)

                            y = [dbox_start, center + lhh]
                            vertices_bottom.append(y)

                        A = [dbox_start, Y + idm]
                        aprime = [dbox_end, Y + idm]
                        vertices_top.extend([A, aprime])

                        G = [dbox_start, Y + H - idm]
                        gprime = [dbox_end, Y + H - idm]
                        vertices_bottom.extend([G, gprime])
                        
                        # Points for linker. Right
                        if r_num < len(regions)-1:
                            xprime = [dbox_end, center - lhh]
                            vertices_top.append(xprime)

                            yprime = [dbox_end, center + lhh]
                            vertices_bottom.append(yprime)
                    
                    # Case ii) Rectangle + trepezoid/triangle
                    elif dbox_start < arrow_collision and dbox_end > arrow_collision:
                        # Points for linker. Left
                        if r_num > 0:
                            x = [dbox_start, center - lhh]
                            y = [dbox_start, center + lhh]
                            vertices_top.append(x)
                            vertices_bottom.append(y)

                        a = [dbox_start, Y + idm]
                        b = [arrow_collision, Y + idm]
                        vertices_top.append(a)
                        vertices_top.append(b)

                        g = [dbox_start, Y + H - idm]
                        f = [arrow_collision, center + h - idm]
                        vertices_bottom.append(g)
                        vertices_bottom.append(f)

                        if dbox_end == L:
                            d = [L, center]
                            vertices_top.append(d)
                        else:
                            if head_start == 0:
                                alpha = int(HL*(L-dbox_end))
                            else:
                                alpha = int(Hl*(L-dbox_end))
                        
                            c = [dbox_end, center - alpha]
                            vertices_top.append(c)

                            e = [dbox_end, center + alpha]
                            vertices_bottom.append(e)

                            # Points for linker. Right
                            if r_num < len(regions)-1:
                                # peek into the next domain region to see how
                                # thin the linker should be here
                                dbox_start_prime = regions[r_num+1][0]
                                if head_start == 0:
                                    alpha1 = HL*(L-dbox_start_prime)
                                else:
                                    alpha1 = Hl*(L-dbox_start_prime)
                                lhh_prime = min(lhh, alpha1)
                                xprime = [dbox_end, center - lhh_prime]
                                vertices_top.append(xprime)

                                yprime = [dbox_end, center + lhh_prime]
                                vertices_bottom.append(yprime)
                        
                    # Case iii) trapezoid
                    else:
                        if head_start == 0:
                            alpha1 = HL*(L-dbox_start)
                        else:
                            alpha1 = Hl*(L-dbox_start)

                        # Points for linker. Left
                        if r_num > 0:
                            lhh_prime = min(lhh, alpha1)
                            x = [dbox_start, center - lhh_prime]
                            vertices_top.append(x)

                            y = [dbox_start, center + lhh_prime]
                            vertices_bottom.append(y)

                        bprime = [dbox_start, center - alpha1]
                        vertices_top.append(bprime)

                        fprime = [dbox_start, center + alpha1]
                        vertices_bottom.append(fprime)

                        if dbox_end == L:
                            d = [L, center]
                            vertices_top.append(d)
                        else:
                            if head_start == 0:
                                alpha2 = HL*(L-dbox_end)
                            else:
                                alpha2 = Hl*(L-dbox_end)

                            c = [dbox_end, center - alpha2]
                            vertices_top.append(c)

                            e = [dbox_end, center + alpha2]
                            vertices_bottom.append(e)

                            # Points for linker. Right
                            if r_num < len(regions)-1:
                                # peek into the next domain region to see how
                                # thin the linker should be here
                                dbox_start_prime = regions[r_num+1][0]
                                if head_start == 0:
                                    alpha1 = HL*(L-dbox_start_prime)
                                else:
                                    alpha1 = Hl*(L-dbox_start_prime)
                                lhh_prime = min(lhh, alpha1)
                                xprime = [dbox_end, center - lhh_prime]
                                vertices_top.append(xprime)

                                yprime = [dbox_end, center + lhh_prime]
                                vertices_bottom.append(yprime)

                # ready to draw domains
                vertices += vertices_top
                vertices += list(reversed(vertices_bottom))

                # Flip if needed
                if flip:
                    for v in vertices:
                        v[0] = L - v[0]

                # Rescale
                for v in vertices:
                    v[0] = v[0]/scaling

                # Make element
                # General properties of the current domain: color and title
                try:
                    color = hmmdb.colors_hex[ID]
                    color_outline = hmmdb.colors_outline_hex[ID]
                except KeyError:
                    color = "#969696"
                    color_outline = "#d2d2d2"
                
                title = ""
                title = DE
                try:
                    title = f"{title} [{hmmdb.alias[ID]}]"
                except KeyError:
                    if alias != "":
                        title = f"{title} [{alias}]"
                title = f"{title}\n"
                if AC != "":
                    title = f"{title} {AC} - "
                title = f"{title}{ID}"
                domain_title = etree.Element("title")
                domain_title.text = title
                
                domain_attribs = {"class": f"domain,{ID}"}
                domain_node_main = etree.Element("g", attrib=domain_attribs)
                domain_node_main.append(domain_title)

                string_vertices = []
                for v in vertices:
                    string_vertices.append(f"{round(v[0]+xoffset, 2)},{round(v[1]+yoffset, 2)}")
                domain_inner_attribs = {
                    "points": " ".join(string_vertices),
                    "fill": color,
                    "stroke": color_outline,
                    "stroke-width": str(domain_contour_thickness),
                    "stroke-linejoin":"round"
                    }
                domain_inner = etree.Element("polygon", attrib=domain_inner_attribs)
                domain_node_main.append(domain_inner)
                domain_elements.append(domain_node_main)

        # Draw intron elements
        for intron_element in intron_elements:
            main_group.append(intron_element)

        # Draw domains on top
        for domain_element in domain_elements:
            main_group.append(domain_element)

        

        return main_group



class BGCDomain:
    def __init__(self, protein, ID, AC, DE, alias, ali_from, ali_to, \
            hmm_from, hmm_to, hmm_size, score, Evalue, algn_seq):
        self.protein = protein
        self.ID = ID                # domain short name e.g. 'ketoacyl-synt'
        self.AC = AC                # domain accession. e.g. 'PF00109.27'
        self.DE = DE                # domain description. e.g. 'Beta-ketoacyl synthase, N-terminal domain'
        self.alias = alias          # domain alias e.g. 'KS'
        self.ali_from = ali_from    # Position in target sequence at which the 
        self.ali_to = ali_to        #   hit starts/ends
        self.hmm_from = hmm_from    # Position in the hmm at which the hit 
        self.hmm_to = hmm_to        #   starts/ends
        self.hmm_size = hmm_size    # Size of the hmm model
        self.score = score          # Only depends on profile HMM and target seq.
        self.Evalue = Evalue        # Based on score and DB size
        self.algn_seq = algn_seq    # Sequence of hit, aligned to the hmm profile
        
    def get_sequence(self):
        return self.protein.sequence[self.ali_from:self.ali_to]

    def get_aligned_sequence(self):
        return "{}{}{}".format("-"*self.hmm_from, self.algn_seq, "-"*(self.hmm_size-self.hmm_to))



class Organism:
    def __init__(self):
        self.shortname = ""
        self.fullname = ""
        self.taxid = ""                 # NCBI tax Id
        self.lineage = []



class Metabolite:
    def __init__(self):
        self.name = ""                  # e.g. "1,3,6,8-tetrahydroxynaphthalene"
        self.alias = []                 # Alternative names e.g. "THN"
        self.database = {}              # Key=database (PubChem, 
                                        #  ChemSpider, etc.). Item=id
        self.SMILES = ""
        self.compound_family = ""       # semi-arbitrary e.g. "azaphilone
        self.role = ""                  # defense, communication etc.
        self.properties = []            # e.g. "volatile"
        
        self.characterized = False      # toggle to True only if supported by 
                                        # publication
        self.comment = ""
