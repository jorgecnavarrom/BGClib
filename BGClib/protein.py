#!/usr/bin/env python3

"""
BGC library: protein

Protein information from a BGC(locus)
"""

from .annotations.domain import BGCDomain
from .locus import BGCLocus
from .metadata.metabolite import Metabolite
from .data.utilities import *
from .visualization.arrower import ArrowerOpts
from .annotations.hmm import hmmDB
from pathlib import Path
from subprocess import PIPE, Popen
from Bio import SearchIO
from lxml import etree
from operator import itemgetter
from copy import deepcopy
from io import StringIO


class BGCProtein:
    """
    Information about a Protein encoded by a gene which is part of the gene cluster
    """
        
    def __init__(self):
        self.parent_cluster = None  # Parent BGC object
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
        
        # TODO: metabolites should be linked to BGCs
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
            part_two = f"{seq[-remainder:]}\n"
        return f"{part_one}\n{part_two}"
        
        
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
            
        return f">{self.identifier}|{self.protein_id}|{self.gene}|" \
            f"{molecules}\n{self.sequence80(start, end)}"


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
                # Only Python 3.12+ can do this within the f-string
                dom_list = " | ".join(reversed(domain_alias_list))
                return f"< {dom_list}"
            else:
                dom_list = "]-[".join(reversed(domain_alias_list))
                return f"<--[{dom_list}]--|"
        else:
            if simple:
                dom_list = " | ".join(domain_alias_list)
                return f"{dom_list} >"
            else:
                dom_list = "]-[".join(domain_alias_list)
                return f"|--[{dom_list}]-->"
        
        
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
    
    
    # TODO substitute call to ProteinCollection.predict_domains for a direct
    #   prediction using pyhmmer
    def predict_domains(self, hmmdb, domtblout_path="", cpus=1, tc=True, filterdoms=True):
        """
        domtblout_path is the path where the domtableout file will be deposited. 
        If present, a Path-like object is expected
        """
        from .collections.proteins import ProteinCollection

        assert(isinstance(hmmdb, hmmDB))
        if domtblout_path != "":
            assert(isinstance(domtblout_path, Path))
        
        pc = ProteinCollection()
        pc.proteins[self.identifier] = self
        pc.predict_domains(hmmdb, domtblout_path, cpus, tc, filterdoms)
        
        self.classify_sequence(hmmdb)
        
        return
    
    
    def classify_sequence(self, hmmdb = hmmDB()):
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
            
            if self.domain_set & FAS2 & FAS1:
                self.protein_type = "Fatty Acid Synthase"
            elif self.domain_set & FAS1:
                self.protein_type = "Fatty Acid Synthase beta-chain"
            elif self.domain_set & FAS2:
                self.protein_type = "Fatty Acid Synthase alpha-chain"
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
        elif "TIGR03443" in self.domain_set:
            self.protein_type = "alpha-aminoadipate reductase"
            self.role = "other"
            return
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
        base_attribs = {
            "version":"1.1", 
            "baseProfile":"full",
            "width":str(int(self.length/scaling)),
            "height":str(int(2*h + H))
        }

        root = etree.Element(
            "svg", 
            attrib=base_attribs, 
            nsmap={None:'http://www.w3.org/2000/svg'}
        )
        main_group = etree.SubElement(root, "g")
        
        # add title
        title = etree.Element("title")
        title.text = self.identifier
        main_group.append(title)
        
        # genomic locus line
        genomic_locus_attribs = {
            "x1":"0", 
            "y1":str(int(h+H/2)),
            "x2":str(int(self.length/scaling)),
            "y2":str(int(h+H/2)),
            "style":f"stroke:#0a0a0a; stroke-width:{stripe_thickness}"
        }
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
            
            domain_box_attribs = {
                "x":str(int(domain.ali_from/scaling)),
                "y":"0",
                "rx":str(curve_radius),
                "ry":str(curve_radius),
                "width":str(int((domain.ali_to - domain.ali_from)/scaling)),
                "height":str(int(2*H)),
                "fill":f"rgb({color})"
            }
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
                    intron_attribs = {
                        "d":f"M {current_x-offset:d} 0 L {current_x-offset:d} {2*H:d}",
                        "stroke":"#ffffff",
                        "stroke-width":"2"
                    }
                else:
                    # if reverse, mirror position of introns
                    intron_attribs = {
                        "d":f"M {l-current_x+offset:d} 0 L {l-current_x+offset:d} {2*H:d}",
                        "stroke":"#ffffff",
                        "stroke-width":"2"
                    }
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
    

    def arrow_colors(self, mode, hmmdb=hmmDB()):
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
            
            # TODO: remove! this is a hack
            if self.domain_set & FAS1:
                return "#66BF7F"
            elif self.domain_set & FAS2:
                return "#7CC159"
            
            if self.protein_type in {
                "Fatty Acid Synthase beta-chain",
                "Fatty Acid Synthase alpha-chain"
                }:
                if self.protein_type == "Fatty Acid Synthase beta-chain":
                    return "#66BF7F"
                else:
                    return "#7CC159"
            
            if len(self.domain_list) > 0:
                domain = self.domain_list[0]
                try:
                    hex_color = "".join( hex(int(c))[2:] for c in hmmdb.colors[domain.ID] )
                    color = f"#{hex_color}"
                except KeyError:
                    color = "#eff0f1"
            else:
                color = "#eff0f1"
                
            return color
            
        else:
            # "white"
            return "#ffffff"
        

    def arrow_SVG(self, filename, hmmdb=hmmDB(), svg_options=ArrowerOpts(), xoffset=0, yoffset=0):
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
        if svg_options.show_introns:
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
        
    
    def xml_arrow(
            self, 
            hmmdb = hmmDB(), 
            svg_options = ArrowerOpts(), 
            xoffset = 0, 
            yoffset = 0, 
            mirror = False, 
            needs_class_data = True
        ):
        """Creates an gene SVG figure (xml structure)

        Args:
            hmmdb (hmmDB): If BGC's proteins are annotated and need to be \
                painted, this object provides alias and color information
            svg_options (ArrowerOpts): Object that includes all relevant \
                drawing options
            xoffset (int): Offset X (increases when going to the right)
            yoffset (int): Offset Y (increases when going down)
            mirror (bool): True: flip the figure horizontally
            needs_class_data: to be implemented. Highest level call (protein, \
                BGC, BGC list) should define style classes for genes, introns \
                and domains to save even more space

        Returns:
            etree object
        
        Starting from upper left corner (X,Y), the vertices of the arrow are
        (clockwise): A, B, C, D (tip of the arrow), E, F, G
                    C
                    | \\
          (X,Y)A ___|  \\
                |   B   \\ 
                |       / D
                |___F  /
               G    | /
                    |/
                    E

        Where A = (X,Y) = (xoffset, yoffset+h)

        Shorter genes: C, D, E
        Ribbons: A, B, D, F, G
        Shorter Ribbons: A, D, G

        (not considering vertices created by intron regions)
        But note that here the first coordinate to be annotated will be the one
        in the bottom left corner (G or E, depending on whether the arrow is 
        complete or will be a squeezed version)
        
        Domains are similar, except that b and c follow the slope of the arrow
        head. Vertices:
        a, b, c, d (if domain end matches the arrow's), e, f, g
    
        Other relevant numbers:
            L: total arrow length
            l: length of arrow head (45° angle)
            H: height of arrow's body distance(A,G) == distance(B,F)
            h: height of top/bottom arrow edge distance(B,C) == dist.(E,F)
            
        Note that scaling only happens on the x axis
        """
        
        fill_color = self.arrow_colors(svg_options.color_mode, hmmdb)
        
        if mirror:
            flip = self.forward
        else:
            flip = not self.forward and svg_options.original_orientation
        
        L = self.cds_regions[-1][1] - self.cds_regions[0][0]
        
        scaling = svg_options.scaling
        shape = svg_options.shape

        H = svg_options.arrow_height
        H_scaled = H * scaling
        h = 0
        l = 0.5 * H_scaled # needs to be scaled along the x-axis
        if shape == "Arrow":
            h = H/2
            l = H_scaled 
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
        
        # If only one domain and color_mode == 'domains', don't draw domains 
        # even if requested (instead, whole arrow will be filled with domain's 
        # color)
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
        
        # This determines whether we have a normal arrow/ribbon or a squished
        # triangle
        head_start = L - l # In local 'arrow coordinates'
        if head_start <= 0:
            head_start = 0
        center = Y + 0.5*H

        # alpha = (H*(L-start)/l) distance from center to colission with head
        Hl = H/l 
        HL = H/L # for shorter arrows

        #
        # DRAW MAIN GENE AREA
        vertices = list()
        D = [L, center]
        if head_start == 0:
            # short arrow
            if shape == "Ribbon":
                A = [0, Y]
                G = [0, Y + H]
                vertices = [A, D, G]
            else:
                C = [0, 0]
                E = [0, Y + H + h]
                vertices = [C, D, E]
        else:
            # full arrow
            A = [0, Y]
            B = [head_start, Y]
            F = [head_start, Y + H]
            G = [0, Y + H]

            if shape == "Ribbon":
                vertices = [A, B, D, F, G]
            else:
                C = [head_start, 0]
                E = [head_start, Y + H + h]
                vertices = [A, B, C, D, E, F, G]

        # flip if requested and gene is in the reverse strand
        if flip:
            for v in vertices:
                v[0] = L - v[0]
        
        # rescale everything
        rescale = 1 / scaling
        for v in vertices:
            v[0] = v[0]*rescale
        
        # shape properties
        string_vertices = []
        for v in vertices:
            string_vertices.append(
                f"{round(v[0]+xoffset, 2)},{round(v[1]+yoffset, 2)}"
            )
        
        region_attribs = dict()
        region_attribs["class"] = "gene"
        region_attribs["points"] = " ".join(string_vertices)
        region_attribs["fill"] = fill_color

        if needs_class_data:
            region_attribs["stroke-width"] = str(
                svg_options.gene_contour_thickness
            )
            if svg_options.outline:
                region_attribs["stroke"] = "#0a0a0a"
                
        gene = etree.Element("polygon", attrib = region_attribs)
        main_group.append(gene)

        #
        # DRAW INTRONS
        intron_regions = []
        intron_elements = []
        if show_introns and len(self.cds_regions) > 1:
            # calculate intron (inter-CDS) regions.
            # Coordinates are relative to the start of the gene-arrow
            if self.forward:
                cds_start = self.cds_regions[0][0]
                for cds_num in range(len(self.cds_regions) - 1):
                    intron_regions.append(
                        (self.cds_regions[cds_num][1] + 1 - cds_start, 
                        self.cds_regions[cds_num+1][0] - 1 - cds_start)
                    )
            else:
                cds_end = self.cds_regions[-1][1] + 1
                for cds_num in range(len(self.cds_regions)-1, 0, -1):
                    intron_regions.append(
                        (cds_end - self.cds_regions[cds_num][0], 
                        cds_end - self.cds_regions[cds_num-1][1] - 1)
                    )

            for intron_start, intron_end in intron_regions:
                vertices = []

                # draw left-most part of intron
                if intron_start < head_start:
                    g = [intron_start, Y + H]
                    vertices.append(g)
                    a = [intron_start, Y]
                    vertices.append(a)
                else:
                    if head_start == 0:
                        # squeezed triangle
                        alpha = HL*(L-intron_start)
                        if shape == "Ribbon":
                            alpha *= 0.5
                    else:
                        # full arrow/ribbon
                        alpha = L-intron_start
                        if shape == "Arrow":
                            alpha *= Hl
                        else:
                            # alpha is scaled in arrows when dividing by L or
                            # l, but not in ribbons
                            alpha *= rescale
                
                    e_prime = [intron_start, center + alpha]
                    vertices.append(e_prime)
                    c_prime = [intron_start, center - alpha]
                    vertices.append(c_prime)

                # draw right-most part of region
                if intron_end <= head_start:
                    # intron begins and ends before head_start; is a rectangle
                    a_prime = [intron_end, Y]
                    vertices.append(a_prime)
                    g_prime = [intron_end, Y + H]
                    vertices.append(g_prime)
                else:
                    if (intron_start < head_start):
                        b = [head_start, Y]
                        vertices.append(b)
                        c = [head_start, Y-h]
                        vertices.append(c)
                    
                    if head_start == 0:
                        # squeezed triangle
                        alpha = HL*(L-intron_end)
                        if shape == "Ribbon":
                            alpha *= 0.5
                    else:
                        # full arrow/ribbon
                        alpha = L-intron_end
                        if shape == "Arrow":
                            alpha *= Hl
                        else:
                            # alpha is scaled in arrows when dividing by L or
                            # l, but not in ribbons
                            alpha *= rescale
                    c_primeprime = [intron_end, center - alpha]
                    vertices.append(c_primeprime)
                    e_primeprime = [intron_end, center + alpha]
                    vertices.append(e_primeprime)
                        
                    if (intron_start < head_start):
                        e = [head_start, Y + H + h]
                        vertices.append(e)
                        f = [head_start, Y + H]
                        vertices.append(f)

                # flip if requested and gene is in the reverse strand
                if flip:
                    for v in vertices:
                        v[0] = L - v[0]
                
                # rescale everything
                for v in vertices:
                    v[0] = v[0]*rescale
                
                # shape properties
                region_attribs = {
                    "stroke-width": str(svg_options.gene_contour_thickness),
                    "stroke-opacity": "0.85",
                    "class": "intron",
                    "fill": "#e6e6e6",
                    "fill-opacity": "0.5"
                }
                # region_attribs["stroke-dasharray"] = "3,2"
                if svg_options.outline:
                    region_attribs["stroke"] = "#919191"

                string_vertices = []
                for v in vertices:
                    string_vertices.append(
                        f"{round(v[0]+xoffset, 2)},{round(v[1]+yoffset, 2)}"
                    )
                region_attribs["points"] = " ".join(string_vertices)
                intron_element = etree.Element("polygon", attrib=region_attribs)

                intron_elements.append(intron_element)

        # Draw intron elements
        for intron_element in intron_elements:
            main_group.append(intron_element)

        #
        # DOMAINS 
        domain_elements = []
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

            # Convert a domain (start/end) into domain regions (broken 
            # according to gene structure)
            for domain in self.domain_list:
                # get cDNA coordinates
                dstart = domain.ali_from * 3    # inclusive dstart
                dend = (domain.ali_to + 1) * 3  # exclusive dend

                # get genomic coordinates
                offset = 0
                regions = []
                for current_cds, cds in enumerate(forward_regions):
                    # CDS is before domain. Push domain forward
                    if cds[1] < dstart:
                        try:
                            # advance the length of the incoming intron
                            offset = forward_regions[current_cds+1][0] \
                                - cds[1] + 1
                        except IndexError:
                            print(f"Warning, domain {domain.AC} out of range" \
                                  + f" {self.identifier} (incomplete " \
                                  + "GenBank?)")
                            # print(self.identifier)
                            # print(domain.AC)
                            # print(cds, dstart, dend)
                            # print(forward_regions)
                            # print(current_cds)
                            continue
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
                            offset = forward_regions[current_cds+1][0] \
                                - cds[1] + 1
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
                if shape == "Arrow":
                    arrow_collision = L - ( (h-idm)/HL )
                else:
                    # derived from the relationship
                    # beta/idm = L / 0.5H
                    arrow_collision = 2 * idm * L / H
            else:
                if shape == "Arrow":
                    arrow_collision = head_start + (h+idm)/Hl
                else:
                    arrow_collision = head_start + (idm*scaling)
            
            # "linker half height"
            lhh = 0.125 * H
            
            # Work every domain
            for domain in domain_regions:
                ID = domain["ID"]
                DE = domain["DE"]
                alias = domain["alias"]
                AC = domain["AC"]
                regions = domain["regions"]

                vertices = []
                vertices_top = []
                vertices_bottom = []

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
                    elif dbox_start < arrow_collision \
                            and dbox_end > arrow_collision:
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
                        # f = [arrow_collision, center + h - idm]
                        f = [arrow_collision, Y + H - idm]
                        vertices_bottom.append(g)
                        vertices_bottom.append(f)

                        if dbox_end == L:
                            d = [L, center]
                            vertices_top.append(d)
                        else:
                            if head_start == 0:
                                # squeezed triangle
                                alpha = HL*(L-dbox_end)
                                if shape == "Ribbon":
                                    alpha *= 0.5
                            else:
                                # full arrow/ribbon
                                alpha = L-dbox_end
                                if shape == "Arrow":
                                    alpha *= Hl
                                else:
                                    # alpha is scaled in arrows when dividing 
                                    # by L or l, but not in ribbons
                                    alpha *= rescale
                            # alpha = int(alpha)
                        
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
                            # squeezed triangle
                            alpha1 = HL*(L-dbox_start)
                            if shape == "Ribbon":
                                alpha1 *= 0.5
                        else:
                            # full arrow/ribbon
                            alpha1 = L-dbox_start
                            if shape == "Arrow":
                                alpha1 *= Hl
                            else:
                                # alpha is scaled in arrows when dividing 
                                # by L or l, but not in ribbons
                                alpha1 *= rescale
                        # alpha1 = int(alpha1)

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
                                if shape == "Ribbon":
                                    alpha2 *= 0.5
                            else:
                                alpha2 = (L-dbox_end)
                                if shape == "Arrow":
                                    alpha2 *= Hl
                                else:
                                    alpha2 *= rescale

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
                    v[0] = v[0]*rescale

                # Make element
                # General properties of the current domain: color and title
                try:
                    color = hmmdb.colors_hex[ID]
                    color_outline = hmmdb.colors_outline_hex[ID]
                except KeyError:
                    # non-persistent way of having new colors. 
                    # for persistance, need to keep a list of all the new
                    # colors and at some point, append them to the color 
                    # tsv file
                    s_ = (0.45, 0.7)
                    v_ = (0.5, 0.8)
                    color = random_color_tuple((0.0, 1.0), s_, v_)

                    # make darker outline. 
                    # Convert: Hex -> RGB -> HSV, darken -> Hex
                    r = int(f"0x{color[1:3]}",0)
                    g = int(f"0x{color[3:5]}",0)
                    b = int(f"0x{color[5:7]}",0)
                    h_, s_, v_ = rgb_to_hsv(r/255.0, g/255.0, b/255.0)
                    rgb_darker = tuple(int(round(c * 255)) for c in \
                                       hsv_to_rgb(h_, s_, 0.8*v_))
                    color_outline = f"#{rgb_darker[0]:02x}{rgb_darker[1]:02x}{rgb_darker[2]:02x}"
                
                    hmmdb.colors_hex[ID] = color
                    hmmdb.colors_outline_hex[ID] = color_outline

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

        # Draw domains on top
        for domain_element in domain_elements:
            main_group.append(domain_element)

        return main_group
