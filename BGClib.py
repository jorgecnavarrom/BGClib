#!/usr/bin/env python

"""
BGC library

Classes for data storage and analysis of Biosynthetic Gene Clusters

"""

import sys
from pathlib import Path
from subprocess import PIPE, Popen
import warnings
from multiprocessing import Pool, cpu_count
from random import uniform
from colorsys import hsv_to_rgb
from colorsys import rgb_to_hsv
from collections import defaultdict
from math import sin, pi, atan2
from lxml import etree
from operator import itemgetter
from copy import deepcopy

try:
    from io import StringIO
    from Bio import BiopythonExperimentalWarning
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonExperimentalWarning)
        from Bio import SearchIO
    from Bio import SeqIO
except ModuleNotFoundError:
    sys.exit("BGC lib did not find all needed dependencies")

__author__ = "Jorge Navarro"
__version__ = "0.4.0"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@westerdijkinstitute.nl"


# Static data
valid_CBP_types = {"nrPKS", "rPKS", "NRPS", "t3PKS", "unknown", "other",
        "PKS-NRPS_hybrid", "PKS-mmNRPS_hybrid", "NRPS-PKS_hybrid", "other_PKS", 
        "unknown_PKS", "no_domains", "NIS"}

# this should all have an alias in CBP_domains.tsv
PKS_domains = {"SAT", "ketoacyl-synt", "Ketoacyl-synt_C", "KAsynt_C_assoc",
                "Acyl_transf_1", "TIGR04532"}
PKS3_domains = {"Chal_sti_synt_N", "Chal_sti_synt_C"}
NRPS_domains = {"Condensation", "AMP-binding", "AMP-binding_C"}
reducing_domains = {"PKS_ER_names_mod", "KR", "PS-DH"}
NRPS_Independent_Siderophore_domains = {"IucA_IucC"}
#Terpene_domains = set({"Terpene_synth", "Terpene_synth_C"})
#Squalene_domains = set({"SQHop_cyclase_N", "SQHop_cyclase_C"})

FAS_domains_A = {"Fas_alpha_ACP" ,"FAS_I_H", "ACPS"}
FAS_domains_B = {"DUF1729", "FAS_meander", "MaoC_dehydrat_N", "MaoC_dehydratas"}

role_colors = {"biosynthetic":"#f06c6e", # red, rgb(240, 108, 110)
               "tailoring":"#8fc889", # green, rgb(143, 200, 137)
               "transporter":"#f0d963", # yellow, rgb(240, 217, 99)
               "transcription factor":"#33c1f0", # blue, rgb(51, 193, 240)
               "other":"#f0d963", # light gray, rgb(239, 240, 241)
               "unknown":"#dcdcdc"} # gray, rgb(220, 220, 220)


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
    
    return "#{}".join( hex(int(c * 255)) for c in hsv_to_rgb(h, s, v) )


# Classes definition
class HMM_DB:
    """
    This class keeps information about HMM databases to be used by the other
    classes
    """
    
    def __init__(self):
        self.db_list = []           # list of paths to hmm databases
        self.alias = {}             # ID to alias
        self.colors = {}            # ID to tuple(r,g,b)
        self.color_outline = {}
        self.cores = 0              # for hmmer. Remember that it always uses an
                                    # extra core for reading the database
                                    
        self.ID_to_acc = {}   # TODO: load these data...
        self.ID_to_desc = {}
        
        return
    
    def add_database(self, db_path):
        if not db_path.is_file():
            print("Not able to add hmm database (not a file. Wrong path?)")
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
        
        self.db_list.append(db_path)
        return True

    def read_alias_file(self, alias_file):
        try:
            with open(alias_file, "r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    line = line.split("\t")
                    hmm_alias = line[2]
                    hmm_ID = line[3]
                    
                    self.alias[hmm_ID] = hmm_alias
        except FileNotFoundError:
            print("Could not open domain alias file ({})".format(alias_file))

    def read_domain_colors(self, colors_file):
        """
        Expects the input file to be in a tsv with two columns:
        hmm_ID \t r,g,b
        """
        try:
            with open(colors_file, "r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    hmm_ID, colors = line.strip().split("\t")
                    self.colors[hmm_ID] = tuple(colors.split(","))
                    
                    r, g, b = colors.split(",")
                    h, s, v = rgb_to_hsv(int(r)/255.0, int(g)/255.0, int(b)/255.0)
                    self.color_outline[hmm_ID] = tuple(str(int(round(c * 255))) for c in hsv_to_rgb(h, s, 0.8*v))
        except FileNotFoundError:
            print("Could not open domain colors file ({})".format(colors_file))
        pass


class ArrowerOpts():
    """
    Options for Arrower-like figures. Colors will be elsewhere.
    """
    def __init__(self):
        self.scaling = 30                    # px per bp
        
        self.arrow_height = 30               # H: Height of arrows' body
        #self.arrow_head_height = 15          # h: Additional width of arrows' head
                                        #  (i.e. total arrow head = H + 2*aH)
                                        # internally, h = H/2 to preserve 45° angle
        self.arrow_head_length = 30
        
        self.gene_contour_thickness = 2      # note: thickness grows outwards
        
        self.internal_domain_margin = 3
        self.domain_contour_thickness = 1
        
        self.stripe_thickness = 3
        self.stripe_color = (100, 100, 100)
        
        self.fontsize = 30
        self.label_arrow_margin = 10
        
        self.write_id = True
        self._color_mode = "white"
        self.valid_color_modes = {"white", "gray", "random-pastel", "random-dark", 
                                  "random", "roles"}
        
        self.outline = True
        self.draw_domains = True # If false, all arrows will point forward
        
        self.original_orientation = False
        
        self.intron_break = False # Show a dashed line where the intron breaks gene
                                # This means that this is not an actual 
        self.intron_regions = False # Show the actual gap where the intron is

    @property
    def color_mode(self):
        return self._color_mode
    @color_mode.setter
    def color_mode(self, cm):
        if cm in self.valid_color_modes:
            self._color_mode = cm
        else:
            print("Color mode not supported; defaulting to 'white'")
            self._color_mode = "white"


class BGCCollection:
    """
    This class will allow implementation of collection-wide functions such as 
    single-step prediction of domains and comparisons between collections
    """
    
    def __init__(self):
        self.bgcs = {}
        self.name = ""
    
    
    def predict_domains(self, hmmdb, domtblout_path="", cpus=1):
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
        for b in self.bgcs:
            bgc = self.bgcs[b]
            for protein in bgc.protein_list:
                if protein.identifier != "":
                    pc.proteins[protein.identifier] = protein
                else:
                    missing_identifier = True
        
        if missing_identifier:
            print("Warning: one or more proteins don't have a valid identifier")
            
        pc.predict_domains(hmmdb, domtblout_path, cpus)
        
        return


class BGC:
    def __init__(self, gbkfile=None):
        self.identifier = ""    # usually the file name
        self.gca = []           # Gene Cluster Architecture. List of CBP types 
                                #  in the cluster
        self.gcf = []           # should be a list of CBP signatures (numbers?)
        self.products = set()   # All different 'product' qualifiers annotated 
                                #  by antiSMASH
        self.contig_edge = False    # antiSMASH v4+ was not able to fully 
                                    #  complete the extension phase. BGC *might*
                                    #  be fragmented
    
        self.protein_list = []      # should also be present in loci
        self.loci = []

        self.accession = ""
        self.definition = ""
        self.organism = ""
        self.TaxId = ""
        
        # try to initialize object with GenBank content
        if gbkfile is not None:
            self.load(gbkfile)
        
        
    def load(self, _gbk):
        """
        Initializes the object by reading all CDS fields from a GenBank file
        
        Input is intended to be a Path() object (but it should work with a string)
        """
        
        if isinstance(_gbk, Path):
            gbk = _gbk
        else:
            gbk = Path(_gbk)
        clusterName = gbk.stem
        
        self.identifier = clusterName
        
        try:
            records = list(SeqIO.parse(str(gbk), "genbank"))
        except ValueError as e:
            print("Error, not able to parse file {}: {}".format(str(gbk), str(e)))
        else:
            self.accession = records[0].id
            self.definition = records[0].description
            self.organism = records[0].annotations["organism"]
            
            cds_list = []
            
            # traverse all possible records in the file. There's usually only 1
            for record in records:
                locus = BGCLocus()
                locus.length = len(record.seq)
                
                cds_num = 0
                product = ""

                for feature in record.features:
                    if feature.type == "cluster":
                        if "product" in feature.qualifiers:
                            # I don't think this'll every happen...
                            if len(feature.qualifiers["product"]) > 1:
                                print("  WARNING: more than product annotated in {}".format(gbk.name))
                            for product in feature.qualifiers["product"]:
                                for p in product.split("-"):
                                    self.products.add(p.replace(" ",""))
                                
                        if "contig_edge" in feature.qualifiers:
                            if feature.qualifiers["contig_edge"][0] == "True":
                                self.contig_edge = True
                                
                        continue
                                
                    if feature.type == "CDS":
                        cds_num += 1
                        
                        CDS = feature
                        
                        cds_start = max(0, int(CDS.location.start))
                        cds_end = max(0, int(CDS.location.end))
                        
                        
                        accession = ""
                        if "protein_id" in CDS.qualifiers:
                            accession = CDS.qualifiers["protein_id"][0]
                        
                            
                        identifier = clusterName + "~" + "CDS" + str(cds_num)
                        if accession != "":
                            identifier += "~" + accession
                            
                        protein = BGCProtein()
                        protein.identifier = identifier
                        protein.accession = accession
                        protein.sequence = CDS.qualifiers["translation"][0]
                        
                        if CDS.location.strand != 1:
                            protein.forward = False
                        
                        # should only point to BGC objects but we're still
                        # not outputting self.identifier into the annotation file
                        protein.parent_cluster = clusterName
                        
                        del cds_list[:]
                        for x in CDS.location.parts:
                            cds_list.append((int(x.start), int(x.end)))
                        protein.cds_regions = tuple(sorted(cds_list, key=itemgetter(0)))
                        
                        self.protein_list.append(protein)
                        locus.protein_list.append(protein)
                        locus.gene_coordinates.append((cds_start,cds_end))
              
                self.loci.append(locus)
  
    
    def inter_loci_element(self, xoffset, yoffset, svg_options=ArrowerOpts()):
        """
        Draws the following SVG figure to link loci genomic stripes:
        __/ /__
         / /
        
        And returns the corresponding etree element
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
        # substract all non-coding regions (but intergenic space remains
        # unchanged). As we are mixing generic dna regions with only-coding regions
        # this is just a _representation_ of the BGC; not meant to be accurate
        intron_correction = 0.0
        if not bgc_svg_options.intron_regions:
            for locus in self.loci:
                for p in locus.protein_list:
                    intron_correction += ((p.cds_regions[-1][1] - p.cds_regions[0][0] - 3) - p.length*3) / bgc_svg_options.scaling
                    
        
        # main node
        # Note: we are in dna space (so amino acid sequence must be scaled)
        # width, height, xoffset and yoffset must be corrected or the corners
        # with sharp angles will be chopped off
        thickness = bgc_svg_options.gene_contour_thickness
        correction = 2*thickness
        Xoffset = xoffset + thickness
        Yoffset = yoffset + thickness
        base_attribs = {"version":"1.1", 
                        "baseProfile":"full",
                        "width":str(int(Xoffset + L + thickness - intron_correction)),
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
        
        H = svg_options.arrow_height
        h = 0.5*H
        
        l = 1
        Xoffset = xoffset # start of each locus
        for locus in self.loci:
            L = locus.length/svg_options.scaling
            
            # substract all non-coding regions (but intergenic space remains
            # unchanged). As we are mixing generic dna regions with only-coding regions
            # this is just a _representation_ of the BGC; not meant to be accurate
            intron_correction = 0.0
            if not svg_options.intron_regions:
                for p in locus.protein_list:
                    intron_correction += ((p.cds_regions[-1][1] - p.cds_regions[0][0] - 3) - p.length*3) / svg_options.scaling
                L -= intron_correction
            line_attribs = {
                "x1": str(int(Xoffset)),
                "y1": str(int(yoffset + h + 0.5*H)),
                "x2": str(int(Xoffset+L)),
                "y2": str(int(yoffset + h + 0.5*H)),
                "stroke": "#464646",
                "stroke-width": "2"
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
                inner_tree = protein.xml_arrow(hmmdb, svg_options, Xoffset + gene_position, yoffset, mirror=mirror)
                main_group.append(inner_tree)
                
                p += 1
                if not svg_options.intron_regions:
                    intron_offset += ((protein.cds_regions[-1][1] - protein.cds_regions[0][0] - 3) - protein.length*3) / svg_options.scaling
        
            Xoffset += L
            
            if l < len(self.loci):
                # draw inter-locus fig of length H
                # NOTE: the mirror parameter will flip each individual locus, not
                # the whole collection of them. This means that the inter-loci
                # element is still drawn at the end of each locus
                main_group.append(inter_loci_element(Xoffset, yoffset + 0.5*H, svg_options))
            
                Xoffset += H
            l += 1
        
        return main_group

    

class BGCLocus:
    """
    BGC might be divided in more than one locus (because of sequencing, or, as
    is more often the case in fungi, because of genomics)
    
    This class can be used to organize proteins
    """
    def __init__(self):
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
    
    def __init__(self):
        self.proteins = {}    # key = identifier
        
        
    # TODO: break work on sets of 4 cpus
    # TODO: evaluate whether hmmsearch is better than hmmscan
    # TODO: TEST!
    def predict_domains(self, hmmdb, domtblout_path="", cpus=1):
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
        
        for db in hmmdb.db_list:
            command = ['hmmscan', '--cpu', str(cpus), '--cut_tc', '--noali', '--notextw']
            if domtblout_path != "":
                path = str(domtblout_path / ("output_" + db.stem + ".domtable"))
                command.extend(['--domtblout', path ])
            dbpath = str(db)
            command.extend([ dbpath, '-'])
            
            proc_hmmscan = Popen(command, shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = proc_hmmscan.communicate(input="\n".join(protein_list).encode("utf-8"))
            
            # "SearchIO parses a search output file's contents into a hierarchy of four 
            # nested objects: QueryResult, Hit, HSP, and HSPFragment"
            # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html
            results = SearchIO.parse(StringIO(out.decode()), 'hmmer3-text')
            for qresult in results:
                for hit in qresult:
                    for hsp in hit:
                        hspf = hsp[0] # access to HSPFragment

                        seq_identifier = qresult.id
                        hmm_id = hit.id
                        ali_from = hspf.query_start
                        ali_to = hspf.query_end
                        hmm_from = hspf.hit_start
                        hmm_to = hspf.hit_end
                        env_from = hsp.env_start
                        env_to = hsp.env_end
                        
                        Evalue = hsp.evalue
                        score = hsp.bitscore
                        
                        domain = BGCDomain(self.proteins[seq_identifier], 
                                           hmm_id, env_from, env_to, ali_from, 
                                           ali_to, hmm_from, hmm_to, score, Evalue)
                        
                        self.proteins[seq_identifier].domain_list.append(domain)
        
        with Pool(cpus) as pool:
            for p in self.proteins:
                protein = self.proteins[p]
                pool.apply_async(protein.filter_domains())
                protein.attempted_domain_prediction = True
            pool.close()
            pool.join()
    
        return
    

class BGCProtein:
    """
    Information about a Protein encoded by a gene which is part of the gene cluster
    """
        
    def __init__(self):
        self.parent_cluster = None  # Should point to an object of the BGC class
        
        self.identifier = ""        # Should be unique. In principle:
                                    # BGCname:CDS#:ref_accession
                                    # (use _accession if ref_accession not available)

        
        self._accession = ""        # original NCBI accession if possible
        self._ref_accession = ""    # "RefSeq Selected Product" from NCBI's IPG DB
        self.ncbi_id = ""           # NCBI internal id in the Protein database
        
        # NCBI internal id in the Identical Protein Group database. Should be unique
        self.ncbi_ipg_id = ""
        
        self._sequence = ""
        self.length = 0             # automatically calculated when sequence is added
                                    # Amino acid space
        
        self.role = "unknown"       # e.g. [biosynthetic, transporter, tailoring, 
                                    #   resistance, unknown]
        self._CBP_type = "unknown"  # should be one from "valid_CBP_types"
        self.s_cbp = ""             # Sub classes of CBP e.g.: Group V
        self.ss_cbp = ""            # e.g.: Group V2
        
        # e.g.: nrPKS, nrPKS+hrPKS etc.
        # TODO: redefine as "Gene Cluster Architecture"?
        # then GCF would be defined by the list of CBP-signatures
        self._gca = []
        
        
        self.compound_family = ""   # Compound family e.g. "emodin-like"
        self.compound = ""
        self.source = "unknown"     # e.g. MIBiG, curated document etc.
        self.organism = ""
        self.TaxId = ""
        
        self.forward = True
        
        self.cds_regions = tuple()  # collection of start/stop CDS regions (to 
                                    # draw intron position).
                                    # These come from the original GenBank file
                                    # so are relative to the start of the locus
                                    # They are in nucleotide space.
        
        self.domain_list = []       # list of BGCDomain objects
        self.domain_set = set()     # set of unique domain IDs
        self.attempted_domain_prediction = False
        
        return
    
    # always try to have an identifier. Either the reference accession or the 
    # original accession
    @property
    def accession(self):
        return self._accession
    @accession.setter
    def accession(self, acc):
        self._accession = acc
        if self.identifier == "":
            self.identifier = acc
            
    @property
    def ref_accession(self):
        return self._ref_accession
    @ref_accession.setter
    def ref_accession(self, ra):
        self._ref_accession = ra
        if self.identifier == "":
            self.identifier = ra
    
    @property
    def CBP_type(self):
        return self._CBP_type    
    @CBP_type.setter
    def CBP_type(self, cbptype):
        try:
            assert cbptype in valid_CBP_types
        except AssertionError:
            print("{} is not a valid CBP type".format(cbptype))
            self._CBP_type = "unknown"
        else:
            self._CBP_type = cbptype
    
    @property
    def gca(self):
        return "+".join(self._gca)
    @gca.setter
    def gca(self, gca_in):
        if isinstance(gca_in, list):
            self._gca = gca_in
        else:
            self._gca = gca_in.split("+")
    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, seq):
        # note: does not check if contains valid amino acid characters yet
        self._sequence = seq.replace("\n", "").strip().replace(" ", "").upper()
        self.length = len(seq)
    
    
    def get_annotations(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.parent_cluster, self.gca, self.identifier, self.accession, 
            self.ref_accession, self.ncbi_id, self.ncbi_ipg_id, self.CBP_type, 
            self.s_cbp, self.ss_cbp, self.compound_family, self.compound, 
            self.source, self.organism, self.TaxId)
    
    
    def sequence80(self):
        return "\n".join([self._sequence[i:i+80] for i in range(0, self.length, 80)])
        
        
    def fasta(self):
        """
        Returns a fasta-formatted data
        Tries to use reference accession for header. If not present, falls back
        to original accession. If either are absent, use identifier (should
        always be present)
        """
        compound = ""
        if self.compound != "":
            compound = " {}".format(self.compound.replace("/",""))
            
        #if self.ref_accession != "":
            #return ">{}{}\n{}".format(self.ref_accession, compound, self.sequence80())
        #elif self.accession != "":
            #return ">{}{}\n{}".format(self.accession, compound, self.sequence80())
        #else:
        return ">{}{}\n{}".format(self.identifier, compound, self.sequence80())


    def domain_string(self, domain_alias):
        """Returns a basic domain-organization string
        
        - It uses the hmm model name (not accession number)
        - Does NOT check if the domains are already ordered by start coordinate
        
        in:
            domain_alias: extracted from CBP_domains.tsv
            domains: list of tuples with domain hit information
            seq_id: identifier of the sequence
        out: 
            a string
        """
        
        # get names
        domain_name_list = [x.ID for x in self.domain_list]
        
        # try to convert names to alias
        domain_alias_list = []
        for d in domain_name_list:
            if d in domain_alias:
                domain_alias_list.append(domain_alias[d])
            else:
                domain_alias_list.append(d)
        
        tag = self.identifier
        if self.organism != "":
            tag += "_[{}]".format(self.organism)
        if self.compound != "":
            tag += "_[{}]".format(self.compound)
        
        if self.forward:
            return "|--[{}]-->\t{}".format("]-[".join(domain_alias_list), tag)
        else:
            return "<--[{}]--|\t{}".format("]-[".join(list(reversed(domain_alias_list))), tag)
        
        
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
    
    
    def predict_domains(self, hmmdb, domtblout_path="", cpus=0):
        """
        domtblout_path is the path where the domtableout file will be deposited. If
        present, a Path-like object is expected
        """
        assert(isinstance(hmmdb, HMM_DB))
        if domtblout_path != "":
            assert(isinstance(domtblout_path, Path))
        
        for db in hmmdb.db_list:
            command = ['hmmscan', '--cpu', str(cpus), '--cut_tc', '--noali', '--notextw']
            if domtblout_path != "":
                path = str(domtblout_path / (self.identifier + "_" + db.stem + ".domtable"))
                command.extend(['--domtblout', path ])
            dbpath = str(db)
            command.extend([ dbpath, '-'])
            
            proc_hmmscan = Popen(command, shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = proc_hmmscan.communicate(input=self.fasta().encode("utf-8"))
            
            # "SearchIO parses a search output file's contents into a hierarchy of four 
            # nested objects: QueryResult, Hit, HSP, and HSPFragment"
            # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html
            results = SearchIO.parse(StringIO(out.decode()), 'hmmer3-text')
            for qresult in results:
                for hit in qresult:
                    for hsp in hit:
                        hspf = hsp[0] # access to HSPFragment

                        seq_id = qresult.id
                        hmm_id = hit.id
                        ali_from = hspf.query_start
                        ali_to = hspf.query_end
                        hmm_from = hspf.hit_start
                        hmm_to = hspf.hit_end
                        env_from = hsp.env_start
                        env_to = hsp.env_end
                        
                        Evalue = hsp.evalue
                        score = hsp.bitscore
                        
                        domain = BGCDomain(self, hmm_id, env_from, env_to, ali_from, 
                                           ali_to, hmm_from, hmm_to, score, Evalue)
                        
                        self.domain_list.append(domain)
        
        self.filter_domains()
        self.attempted_domain_prediction = True
        
        return
    
    
    def classify_sequence(self):
        """Classifies a sequence based on its predicted domains according to a 
        set of rules.
        
        This will be a work in progress as new rules are discovered
        """
        
        self._CBP_type = "unknown"
        self.role = "unknown"
        
        # Basic filtering
        #if len(self.domain_set) == 1:
            #self._CBP_type = "unknown"
            #self.role = "unknown"
            #return
            #if "Trp_DMAT" in domain_set:
                #return "DMAT"
            #else:
                #return "other"
        
        sequence_type = ""

        # pks/nrps hybrid
        if len(self.domain_set & PKS_domains) > 0 and len(self.domain_set & NRPS_domains) > 0:
            for d in self.domain_list:
                if d.ID in PKS_domains:
                    # Detect PKS/multi-modular NRPS. Count both C and A domains 
                    # just in case there is a broken domain counted twice
                    A_domains = len([a.ID for a in self.domain_list if a.ID == "AMP-binding"])
                    C_domains = len([c.ID for c in self.domain_list if c.ID == "Condensation"])
                    if A_domains > 1 and C_domains > 1:
                        self._CBP_type = "PKS-mmNRPS_hybrid"
                    else:
                        self._CBP_type = "PKS-NRPS_hybrid"
                    self.role = "biosynthetic"
                    return
                elif d.ID in NRPS_domains:
                    self._CBP_type = "NRPS-PKS_hybrid"
                    self.role = "biosynthetic"
                    return
                else:
                    pass
            
        # nrPKS or (h/p)rPKS
        elif len(self.domain_set & PKS_domains) > 0:
            # try to set appart FAS-like sequences of e.g. azaphilone
            if "ketoacyl-synt" not in self.domain_set:
                sequence_type = "unknown_PKS"
                self.role = "unknown"
                
            elif len(self.domain_set & {"TIGR04532","SAT"}) > 0:
                # assume that having a SAT domain is enough for labeling as nrPKS
                # but note that in this category, there seem to be three cases for PT:
                # - PT detected by TIGR04532
                # - DH detected by pfam PF14765 (TIGR doesn't pass --cut_tc)
                # - no detection
                # There are either several types of PT domains or not a general model
                # to detect all three+ cases
                sequence_type = "nrPKS"
                self.role = "biosynthetic"
            elif len(self.domain_set & reducing_domains) > 0:
                sequence_type = "rPKS"
                self.role = "biosynthetic"
            else:
                sequence_type = "other_PKS"
                self.role = "biosynthetic"
        elif len(self.domain_set & PKS3_domains) > 0:
            sequence_type = "t3PKS"
            self.role = "biosynthetic"
        elif len(self.domain_set & NRPS_domains) > 0:
            sequence_type = "NRPS"
            self.role = "biosynthetic"
            
        elif len(self.domain_set & NRPS_Independent_Siderophore_domains) > 0:
            sequence_type = "NIS"
            self.role = "biosynthetic"
            
        else:
            sequence_type = "other"
            self.role = "unknown"

        self._CBP_type = sequence_type
        
        return


    # TODO polish code like in arrow_SVG
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
        intron_break = svg_options.intron_break
        intron_regions = svg_options.intron_regions

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
                
            try:
                title = hmmdb.alias[domain.ID]
            except KeyError:
                title = domain.ID
            
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
    

    def arrow_colors(self, mode):
        if mode == "white":
            return "#ffffff"
        elif mode == "random-pastel":
            s_ = (0.15, 0.4)
            v_ = (0.9, 0.95)
            return self.random_color_tuple((0.0, 1.0), s_, v_)
        elif mode == "random-dark":
            s_ = (0.5, 0.75)
            v_ = (0.4, 0.5)
            return self.random_color_tuple((0.0, 1.0), s_, v_)
        elif mode == "random":
            s_ = (0.6, 0.75)
            v_ = (0.65, 0.85)
            return self.random_color_tuple((0.0, 1.0), s_, v_)
        elif mode == "gray":
            return (180, 180, 180)
        elif mode == "role":
            return role_colors[self.role]
        else:
            # "white"
            return "#ffffff"
        

    def arrow_SVG(self, filename, hmmdb, svg_options=ArrowerOpts(), xoffset=0, yoffset=0):
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
        
    
    def xml_arrow(self, hmmdb=HMM_DB(), svg_options=ArrowerOpts(), xoffset=0, yoffset=0, mirror=False):
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
                L: total arrow length. Will increase if considering introns
                l: length of arrow head
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
        """
        original_orientation = svg_options.original_orientation
        
        if mirror:
            flip = self.forward
        else:
            flip = not self.forward and original_orientation
        
        L = self.length*3/svg_options.scaling
        if svg_options.intron_regions:
            # if using nucleotides to calculate length, remember to remove three
            # that correspond to the stop codon
            L = (self.cds_regions[-1][1] - self.cds_regions[0][0] - 3)/svg_options.scaling
        l = svg_options.arrow_head_length
        
        H = svg_options.arrow_height
        h = H/2
        
        #X = 0
        Y = h
        
        scaling = svg_options.scaling
        stripe_thickness = svg_options.stripe_thickness
        
        idm = svg_options.internal_domain_margin
        
        intron_break = svg_options.intron_break
        intron_regions = svg_options.intron_regions

        main_group = etree.Element("g")
        
        # add title
        title = etree.Element("title")
        title.text = self.identifier
        main_group.append(title)
        
        exon_elements = []
        intron_elements = []
        domain_elements = []
        
        # Use original list of CDS regions to make full list of tuples of
        # (start, end, type) where type 0 = exons and type 1 = introns
        regions = []
        if intron_regions and len(self.cds_regions) > 1:
            start = int(self.cds_regions[0][0]/scaling)
            end = int(self.cds_regions[0][1]/scaling)
            
            offset = start
            
            regions.append((0, end-offset, 0))
            for cds in self.cds_regions[1:]:
                new_start = int(cds[0]/scaling)
                new_end = int(cds[1]/scaling)
                
                regions.append((end - offset + 1, new_start - offset - 1, 1))
                regions.append((new_start - offset, new_end - offset, 0))
                
                start = new_start
                end = new_end
        else:
            regions.append((0, L, 0))
            
        
        # precalculate a couple of numbers
        vertices = []
        head_start = L - l # In local 'arrow coordinates'
        if head_start <= 0:
            head_start = 0
        center = Y + h
        Hl = H/l # alpha = (H*(L-start)/l) distance from center to colission with head
        HL = H/L # for shorter arrows
        
        # Keep track of current domain index and other info.
        # Domains should be sorted by position
        if svg_options.draw_domains and len(self.domain_list) > 0:
            current_domain = 0
            
            domain = self.domain_list[current_domain]
            # General properties of the current domain: color and title
            try:
                color = ",".join([str(c) for c in hmmdb.colors[domain.ID]])
                color_outline = ",".join([str(c) for c in hmmdb.color_outline[domain.ID]])
            except KeyError:
                color = "150,150,150"
                color_outline = "210,210,210"
            
            try:
                title = hmmdb.alias[domain.ID]
            except KeyError:
                title = domain.ID
            domain_title = etree.Element("title")
            domain_title.text = title
            
            domain_attribs = {"class": "domain,{}".format(domain.ID)}
            domain_node_main = etree.Element("g", attrib=domain_attribs)
            domain_node_main.append(domain_title)
            
            # get coordinates in dna space
            dstart = (domain.ali_from*3)/svg_options.scaling
            dend = (domain.ali_to*3)/svg_options.scaling
        else:
            # this will prevent trying to draw domains later on
            current_domain = len(self.domain_list) + 1
        
        intron_offset = 0
        # keep track of region number so we can draw the tip on the last region
        r = 0
        # Traverse all exon/intron regions
        for region in regions:
            r += 1

            del vertices[:]
            start = region[0]
            end = region[1]
            region_type = region[2]
            
            # draw left-most part of region
            if r == 1:
                # short arrow
                if head_start == 0:
                    E = [start, Y + H + h]
                    vertices.append(E)
                    C = [start, Y - h]
                    vertices.append(C)
                else:
                    G = [start, Y + H]
                    vertices.append(G)
                    A = [start, Y]
                    vertices.append(A)
            else:
                if start < head_start:
                    G = [start, Y + H]
                    vertices.append(G)
                    A = [start, Y]
                    vertices.append(A)
                else:
                    if head_start == 0:
                        alpha = int(HL*(L-start))
                    else:
                        alpha = int(Hl*(L-start))
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
                
                # if last region, only one vertex
                if r == len(regions):
                    D = [end, center]
                    vertices.append(D)
                else:
                    if head_start == 0:
                        alpha = int(HL*(L-end))
                    else:
                        alpha = int(Hl*(L-end))
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
            
            # shape properties
            string_vertices = []
            for v in vertices:
                string_vertices.append("{},{}".format(int(v[0]+xoffset), int(v[1])+yoffset))
            
            region_attribs = {
                "points": " ".join(string_vertices),
                "stroke-width": str(svg_options.gene_contour_thickness)
                }
            
            if region_type == 0:
                region_attribs["class"] = "exon"
                region_attribs["fill"] = self.arrow_colors(svg_options.color_mode)
                if svg_options.outline:
                    region_attribs["stroke"] = "#0a0a0a"
            else:
                region_attribs["class"] = "intron"
                region_attribs["fill"] = "#bbbbbb"
                region_attribs["stroke-dasharray"] = "2,2"
                if svg_options.outline:
                    region_attribs["stroke"] = "#acacac"
                
            region_element = etree.Element("polygon", attrib = region_attribs)
            
            
            # DOMAINS
            if region_type == 0 and current_domain <= len(self.domain_list) - 1:
                # get domain coordinates considering intron spaces
                dstarti = int(dstart + intron_offset)
                dendi = int(dend + intron_offset)
                
                # indicate we cannot draw another domain, go to next region
                move_on = False

                # Repeat for all domains that we can draw in the current region
                while ((dstarti >= start and dstarti < end) or (dstarti < start and dendi > end) or (dendi > start and dendi <=end)) and (current_domain < len(self.domain_list)) and not move_on:
                    # remember that head_start and arrow_collision are in local
                    # 'arrow coordinates'
                    if head_start == 0:
                        arrow_collision = head_start + HL*(h+idm)
                    else:
                        arrow_collision = head_start + Hl*(h+idm)
                    
                    # get domain coordinates considering intron spaces (for current
                    # domain)
                    dstarti = int(dstart + intron_offset)
                    dendi = int(dend + intron_offset)
                    
                    # Domain viewboxes:
                    dbox_start = dstarti
                    dbox_end = dendi
                    
                    ## Find if there is some portion of the domain in this region
                    #draw = False
                    
                    # Indicate if the domain is split by an intron on any side. 
                    # This will change the configuration of points used on the 
                    # outline path
                    outline_left = True
                    outline_right = True
                    
                    next_domain = False # indicate whether we need to go to the next domain
                    
                    # start of domain is in this region
                    if dstarti >= start and dstarti < end:
                        # ...but the end of it is outside region
                        if dendi > end:
                            dbox_end = end
                            outline_right = False
                            move_on = True
                        else:
                            next_domain = True
                    # domain continues in this region
                    elif dstarti < start and dendi > end:
                        dbox_start = start
                        dbox_end = end
                        outline_left = False
                        outline_right = False
                        move_on = True
                    # domain didn't start in this region, but it ends here
                    else:
                        # This is 'if dendi > start and dendi <= end'
                        dbox_start = start
                        outline_left = False
                        next_domain = True
                    
                    del vertices[:]
                    
                    # - Inner domain vertices -
                    
                    # case i) rectangle
                    if dbox_end <= arrow_collision:
                        g = [dbox_start, Y + H - idm]
                        vertices.append(g)
                        a = [dbox_start, Y + idm]
                        vertices.append(a)
                        aprime = [dbox_end, Y + idm]
                        vertices.append(aprime)
                        gprime = [dbox_end, Y + H - idm]
                        vertices.append(gprime)
                        
                    # case ii) rectangle + trapezoid/triangle
                    elif dbox_start < arrow_collision and dbox_end > arrow_collision:
                        g = [dbox_start, Y + H - idm]
                        vertices.append(g)
                        a = [dbox_start, Y + idm]
                        vertices.append(a)
                        
                        b = [arrow_collision, Y + idm]
                        vertices.append(b)
                        
                        if int(dbox_end) == int(L):
                            d = [L, center]
                            vertices.append(d)
                        else:
                            if head_start == 0:
                                alpha = int(HL*(L-dbox_end))
                            else:
                                alpha = int(Hl*(L-dbox_end))
                            c = [dbox_end, center - alpha]
                            vertices.append(c)
                            e = [dbox_end, center + alpha]
                            vertices.append(e)
                            
                        f = [arrow_collision, center + h - idm]
                        vertices.append(f)
                        
                    # case iii) trapezoid
                    else:
                        if head_start == 0:
                            alpha1 = int(HL*(L-dbox_start))
                        else:
                            alpha1 = int(Hl*(L-dbox_start))
                        bprime = [dbox_start, center - alpha1]
                        vertices.append(bprime)
                        
                        if int(dbox_end) == int(L):
                            d = [L, center]
                        else:
                            if head_start == 0:
                                alpha2 = int(HL*(L-dbox_end))
                            else:
                                alpha2 = int(Hl*(L-dbox_end))
                            c = [dbox_end, center - alpha2]
                            vertices.append(c)
                            e = [dbox_end, center + alpha2]
                            vertices.append(e)
                    
                        fprime = [dbox_start, center + alpha1]
                        vertices.append(fprime)
                    
                    
                    # - Domain outline path -
                    vertices_outline = []
                    vertices_outline_string = []
                    
                    # Not necessary
                    #if flip and (not outline_left or not outline_right):
                        #outline_left = not outline_left
                        #outline_right = not outline_right
                    
                    # general attributes
                    domain_outline_attribs = {
                        "stroke": "rgb({})".format(color_outline),
                        "stroke-linejoin": "round",
                        "stroke-width": str(svg_options.gene_contour_thickness)
                        }
                    
                    # two paths
                    if not outline_left and not outline_right:
                        vertices_outline_top = []
                        vertices_outline_bottom = []
                        
                        # case i) rectangle
                        if dbox_end <= arrow_collision:
                            vertices_outline_top = [a, aprime]
                            vertices_outline_bottom = [g, gprime]
                            
                        # case ii) rectangle + trapezoid
                        elif dbox_start < arrow_collision and dbox_end > arrow_collision:
                            vertices_outline_top = [a, b, c]
                            vertices_outline_bottom = [g, f, e]
                        # case iii) trapezoid
                        else:
                            vertices_outline_top = [bprime, c]
                            vertices_outline_bottom = [fprime, e]
                        
                        if flip:
                            vertices_outline_top = [(L - v[0], v[1]) for v in vertices_outline_top]
                            vertices_outline_bottom = [(L - v[0], v[1]) for v in vertices_outline_bottom]
                            
                        vertices_outline_string = []
                        vertices_outline_string.append("M{} {}".format(int(xoffset+vertices_outline_top[0][0]), int(yoffset+vertices_outline_top[0][1])))
                        vertices_outline_string.append("L")
                        vertices_outline_string.append(", ".join("{} {}".format(int(xoffset+v[0]),int(yoffset+v[1])) for v in vertices_outline_top[1:]))
                        domain_outline_attribs["d"] = " ".join(vertices_outline_string)
                        domain_outline_top = etree.Element("path", attrib=domain_outline_attribs)
                        domain_node_main.append(domain_outline_top)
                        
                        vertices_outline_string = []
                        vertices_outline_string.append("M{} {}".format(int(xoffset+vertices_outline_bottom[0][0]), int(yoffset+vertices_outline_bottom[0][1])))
                        vertices_outline_string.append("L")
                        vertices_outline_string.append(", ".join("{} {}".format(int(xoffset+v[0]),int(yoffset+v[1])) for v in vertices_outline_bottom[1:]))
                        domain_outline_attribs["d"] = " ".join(vertices_outline_string)
                        domain_outline_bottom = etree.Element("path", attrib=domain_outline_attribs)
                        domain_node_main.append(domain_outline_bottom)
                        
                    else:
                        # case i) rectangle
                        if dbox_end <= arrow_collision:
                            if outline_left:
                                vertices_outline = [aprime, a, g, gprime]
                            else:
                                vertices_outline = [a, aprime, gprime, g]
                        # case ii) rectangle + trapezoid/triangle
                        elif dbox_start < arrow_collision and dbox_end > arrow_collision:
                            if outline_right:
                                vertices_outline = [a, b]
                                if int(dbox_end) == int(L):
                                    vertices_outline.append(d)
                                else:
                                    vertices_outline.append(c)
                                    vertices_outline.append(e)
                                vertices_outline.append(f)
                                vertices_outline.append(g)
                            else:
                                vertices_outline = [c, b, a, g, f, e]
                        # case iii) trapezoid/triangle
                        else:
                            if outline_left:
                                vertices_outline = [c, bprime, fprime, e]
                            else:
                                if int(dbox_end) == int(L):
                                    vertices_outline = [bprime, d, fprime]
                                else:
                                    vertices_outline = [bprime, c, e, fprime]
                        
                        if flip:
                            vertices_outline = [(L - v[0], v[1]) for v in vertices_outline]
                            
                        vertices_outline_string = []
                        vertices_outline_string.append("M{} {}".format(int(xoffset+vertices_outline[0][0]), int(yoffset+vertices_outline[0][1])))
                        vertices_outline_string.append("L")
                        vertices_outline_string.append(", ".join("{} {}".format(int(xoffset+v[0]),int(yoffset+v[1])) for v in vertices_outline[1:]))
                        if outline_left and outline_right:
                            vertices_outline_string.append("Z")
                        
                        domain_outline_attribs["d"] = " ".join(vertices_outline_string)
                        
                        domain_outline_attribs["fill"] = "none"
                        domain_outline = etree.Element("path", attrib=domain_outline_attribs)
                        domain_node_main.append(domain_outline)
                    

                    # - Render domain _inner_ block -
                        
                    # flip if requested and gene is in the reverse strand
                    if flip:
                        for v in vertices:
                            v[0] = L - v[0]
                    # shape properties
                    string_vertices = []
                    for v in vertices:
                        string_vertices.append("{},{}".format(int(xoffset+v[0]), int(yoffset+v[1])))
                    
                    domain_inner_attribs = {
                        "points": " ".join(string_vertices),
                        "fill": "rgb({})".format(color),
                        "stroke-linejoin":"round"
                        }
                    domain_inner = etree.Element("polygon", attrib=domain_inner_attribs)
                    domain_node_main.append(domain_inner)
                    
                    #domain_node_main.append(domain_outer)
                    domain_elements.append(domain_node_main)


                    if next_domain:
                        current_domain += 1
                        
                        # is there yet another domain?
                        if current_domain < len(self.domain_list):
                            domain = self.domain_list[current_domain]
                            # General properties of the current domain: color and title
                            try:
                                color = ",".join([str(c) for c in hmmdb.colors[domain.ID]])
                                color_outline = ",".join([str(c) for c in hmmdb.color_outline[domain.ID]])
                            except KeyError:
                                color = "150,150,150"
                                color_outline = "210,210,210"
                            
                            try:
                                title = hmmdb.alias[domain.ID]
                            except KeyError:
                                title = domain.ID
                            domain_title = etree.Element("title")
                            domain_title.text = title
                            
                            domain_attribs = {"class": "domain,{}".format(domain.ID)}
                            domain_node_main = etree.Element("g", attrib=domain_attribs)
                            domain_node_main.append(domain_title)
                            
                            # get coordinates in dna space
                            dstart = (domain.ali_from*3)/svg_options.scaling
                            dend = (domain.ali_to*3)/svg_options.scaling
                            
                            dstarti = int(dstart + intron_offset)
                            dendi = int(dend + intron_offset)
            
            
            # DOMAIN LINKER
            # if the domain runs across two exons, draw a linker line in the intron
            if region_type == 1 and current_domain <= len(self.domain_list) - 1:
                if dstarti < start and dendi + end > end:
                    # +/- 1 because linker goes beyond intron region
                    if flip:
                        x1 = L - start + 1
                        x2 = L - end - 1
                    else:
                        x1 = start - 1
                        x2 = end + 1
                    domain_linker_attribs = {
                        "x1": str(int(xoffset+x1)),
                        "y1": str(int(yoffset+center)),
                        "x2": str(int(xoffset+x2)),
                        "y2": str(int(yoffset+center)),
                        "stroke": "rgb({})".format(color),
                        "stroke-width": str(int(H/4))
                        }
                    domain_linker = etree.Element("line", attrib=domain_linker_attribs)
                    domain_node_main.append(domain_linker)

        
            # Finish this region by appending data
            if region_type == 0:
                exon_elements.append(region_element)
            else:
                intron_offset += end - start
                intron_elements.append(region_element)
            
        
        # INTRON BREAKS
        # Draw intron breaks if intron regions were not requested
        if intron_break and len(self.cds_regions) > 1 and not intron_regions:
            # prepare regions
            regions = []
            start = int(self.cds_regions[0][0]/scaling)
            end = int(self.cds_regions[0][1]/scaling)
            
            offset = start
            
            regions.append((0, end-offset, 0))
            for cds in self.cds_regions[1:]:
                new_start = int(cds[0]/scaling)
                new_end = int(cds[1]/scaling)
                
                regions.append((end - offset + 1, new_start - offset - 1, 1))
                regions.append((new_start - offset, new_end - offset, 0))
                
                start = new_start
                end = new_end
            
            intron_offset = 0
            vertices = []
            vertices_string = []
            r = 0
            for region in regions:
                r += 1
                
                start, end, region_type = region
                
                if region_type == 1:
                    intron_offset += end - start
                else:
                    intron_break = end - intron_offset
            
                    # end of last region doesn't have an intron break
                    if r < len(regions):
                        del vertices[:]
                        del vertices_string[:]
                        
                        if intron_break < head_start:
                            a = [intron_break, Y]
                            vertices.append(a)
                            b = [intron_break, Y+H]
                            vertices.append(b)
                        else:
                            alpha = int(Hl*(L-intron_break))
                            
                            a = [intron_break, center - alpha]
                            vertices.append(a)
                            b = [intron_break, center + alpha]
                            vertices.append(b)
                            
                        if flip:
                            for v in vertices:
                                v[0] = L - v[0]
                                
                        vertices = [[str(int(xoffset+v[0])), str(int(yoffset+v[1]))] for v in vertices]
                            
                        intron_attribs = {"class":"intron"}
                        intron_main = etree.Element("g", attrib=intron_attribs)
                        
                        # intron background
                        intron_attribs = {
                            "d":"M {} {} L {} {}".format(vertices[0][0], vertices[0][1], vertices[1][0], vertices[1][1]),
                            "stroke":"#ffffff",
                            "stroke-width": str(svg_options.gene_contour_thickness)
                            }
                        if intron_break <= head_start:
                            intron_attribs["stroke-linecap"] = "square"
                        else:
                            intron_attribs["stroke-linecap"] = "round"
                        
                        intron_main.append(etree.Element("path", attrib=intron_attribs))
                        
                        # intron dashed line
                        intron_attribs["stroke"] = "#8c8c8c"
                        intron_attribs["stroke-dasharray"] = "2,3"
                        intron_main.append(etree.Element("path", attrib=intron_attribs))
                        
                        # add the element to the main group
                        intron_elements.append(intron_main)
                        
        
        # Assemble the complete figure
        # intron dashed boxes need to go first
        for exon in exon_elements:
            main_group.append(exon)
        if intron_regions and len(self.cds_regions) > 1:
            for intron in intron_elements:
                main_group.append(intron)
            
            for domain in domain_elements:
                main_group.append(domain)
        # but intron breaks need to be drawn after exons
        else:
            #for exon in exon_elements:
                #main_group.append(exon)
            for domain in domain_elements:
                main_group.append(domain)
            for intron in intron_elements:
                main_group.append(intron)
                

        return main_group


class BGCDomain:
    def __init__(self, protein, ID, env_from, env_to, ali_from, ali_to, hmm_from, hmm_to, score, Evalue):
        self.protein = protein
        self.ID = ID                # domain short name e.g. 'ketoacyl-synt'
        self.env_from = env_from    # Pos. in target seq. at which surr. envelope 
        self.env_to = env_to        #   starts/ends.
        self.ali_from = ali_from    # Position in target sequence at which the 
        self.ali_to = ali_to        #   hit starts/ends
        self.hmm_from = hmm_from    # Position in the hmm at which the hit 
        self.hmm_to = hmm_to        #   starts/ends
        self.score = score          # Only depends on profile HMM and target seq.
        self.Evalue = Evalue        # Based on score and DB size
        
    def get_sequence(self):
        return self.protein.sequence[self.ali_from:self.ali_to+1]
