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
__version__ = "0.3.3"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@westerdijkinstitute.nl"

valid_CBP_types = set({"nrPKS", "rPKS", "NRPS", "t3PKS", "unknown", "other",
        "PKS-NRPS_hybrid", "NRPS-PKS_hybrid", "other_PKS", "unknown_PKS", 
        "no_domains", "NIS"})

# this should all have an alias in CBP_domains.tsv
PKS_domains = set({"SAT", "ketoacyl-synt", "Ketoacyl-synt_C", "KAsynt_C_assoc",
                "Acyl_transf_1", "TIGR04532"})
PKS3_domains = set({"Chal_sti_synt_N", "Chal_sti_synt_C"})
NRPS_domains = set({"Condensation", "AMP-binding", "AMP-binding_C"})
reducing_domains = set({"PKS_ER_names_mod", "KR", "PS-DH"})
NRPS_Independent_Siderophore_domains = set({"IucA_IucC"})
#Terpene_domains = set({"Terpene_synth", "Terpene_synth_C"})
#Squalene_domains = set({"SQHop_cyclase_N", "SQHop_cyclase_C"})


class HMM_DB:
    """
    This class keeps information about HMM databases to be used by the other
    classes
    """
    
    def __init__(self):
        self.db_list = []           # list of paths to hmm databases
        self.alias = {}             # ID to alias
        self.colors = {}            # ID to tuple(r,g,b)
        self.cores = 0              # for hmmer. Remember that it always uses an
                                    # extra core for reading the database
                                    
        self.ID_to_acc = defaultdict(str)   # TODO: load these data...
        self.ID_to_desc = defaultdict(str)
        
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
            with open(alias_file,"r") as f:
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
            with open(colors_file,"r") as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    hmm_ID, colors = line.strip().split("\t")
                    self.colors[hmm_ID] = tuple(colors.split(","))
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
        self.arrow_head_height = 15          # aH: Additional width of arrows' head
                                        #  (i.e. total arrow head = H + 2*aH)
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
        self.valid_color_modes = set({"white", "gray", "random-pastel", "random-dark", "random"})
                                            
        self.outline = True

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
        
        protein_list = []
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
            gbk = Path(_path)
        clusterName = gbk.stem
        
        self.identifier = clusterName
        
        try:
            records = list(SeqIO.parse(str(gbk), "genbank"))
        except ValueError as e:
            print("Error, not able to parse file {}: {}".format(str(e)))
        else:
            self.accession = records[0].id
            self.definition = records[0].description
            self.organism = records[0].annotations["organism"]
            
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
                        
                        self.protein_list.append(protein)
                        locus.protein_list.append(protein)
                        locus.gene_coordinates.append((cds_start,cds_end))
              
                self.loci.append(locus)


    def random_color_tuple(self, s_, v_):
        """
        returns a random tupble of color in RGB according to a range of SV values
        Additional info:
        https://en.wikipedia.org/wiki/HSL_and_HSV
        and http://stackoverflow.com/a/1586291
        """
        
        h = uniform(0, 1) # all possible colors
        s = uniform(s_[0], s_[1])
        v = uniform(v_[0], v_[1])
        
        return tuple(int(c * 255) for c in hsv_to_rgb(h, s, v))
        
        
    def gene_colors(self, mode):
        if mode == "random-pastel":
            s_ = (0.15, 0.4)
            v_ = (0.9, 0.95)
            return self.random_color_tuple(s_, v_)
        elif mode == "random-dark":
            s_ = (0.5, 0.75)
            v_ = (0.4, 0.5)
            return self.random_color_tuple(s_, v_)
        elif mode == "random":
            s_ = (0.6, 0.75)
            v_ = (0.65, 0.85)
            return self.random_color_tuple(s_, v_)
        elif mode == "gray":
            return (180, 180, 180)
        else:
            return (255, 255, 255)
    
    
    def SVG_arrow(self, offsetx, Y, tabs, protein, models=HMM_DB(), svg_options=ArrowerOpts()):
        #tabs, gene, baseX, Y, H, h, l, scaling=30):
        """
        SVG code for an arrow representing a single gene:
            - offset: horizontal position where *locus* starts.
            - Y: Top coordinate of the arrow (in line with the top arrow head)
            - (X,Y) ... upper left (+) or right (-) corner of the arrow
            - L ... arrow length
            - H ... arrow height
            - h ... arrow head edge width
            - l ... arrow head length
            the edges are ABCDEFG starting from (X,Y)     
        """
        
        # remember protein objects' lenght is in aminoacids
        L = (protein.length*3)/svg_options.scaling
        X = offsetx/svg_options.scaling
        h = svg_options.arrow_head_height
        H = svg_options.arrow_height
        l = svg_options.arrow_head_length
        gene_contour_thickness = svg_options.gene_contour_thickness
        internal_domain_margin = svg_options.internal_domain_margin
        domain_contour_thickness = svg_options.domain_contour_thickness
        
        arrow = []
        
        arrow.append("{}<g>".format(tabs*"\t"))
        tabs += 1
        
        if protein.forward:
            head_end = L
            if L < l:
                # squeeze arrow if length shorter than head length
                A = (X,Y-h)
                B = (X+L,Y+H/2)
                C = (X,Y+H+h)
                head_start = 0
                points = [A, B, C]
            else:
                A = (X,Y)
                B = (X+L-l,Y)
                C = (X+L-l,Y-h)
                D = (X+L,Y+H/2)
                E = (X+L-l,Y+H+h)
                F = (X+L-l,Y+H)
                G = (X,Y+H)
                head_start = L - l # relative to the start of the gene, not absolute coords.
                points = [A, B, C, D, E, F, G]
        else:
            head_start = 0
            if L < l:
                # squeeze arrow if length shorter than head length
                A = (X,Y+H/2)
                B = (X+L,Y-h)
                C = (X+L,Y+H+h)
                head_end = L
                points = [A, B, C]
            else:
                A = (X+L,Y)
                B = (X+l,Y)
                C = (X+l,Y-h)
                D = (X,Y+H/2)
                E = (X+l,Y+H+h)
                F = (X+l,Y+H)
                G = (X+L,Y+H)
                head_end = l
                points = [A, B, C, D, E, F, G]
        
        head_length = head_end - head_start
        #rounding error?
        if head_length == 0:
            return []
        
        points_coords = []
        for point in points:
            points_coords.append(str(int(point[0])) + "," + str(int(point[1])))
        
        label = protein.identifier
        if protein.accession != "":
            label = protein.accession
        arrow.append("{}\t<title>{}</title>".format(tabs*"\t", label))
        
        outline = ""
        if svg_options.outline:
            arrow_outline = (10, 10, 10)
            outline = "stroke=\"rgb({})\" ".format(",".join([str(x) for x in arrow_outline]))
        
        arrow.append("{}\t<polygon points=\"{}\" fill=\"rgb({})\" {}stroke-width=\"{:d}\" />".format(tabs*"\t", " ".join(points_coords), ",".join([str(x) for x in self.gene_colors(svg_options.color_mode)]), outline, gene_contour_thickness ))
        
        dH = H - 2*internal_domain_margin
        for domain in protein.domain_list:
            # domain coordinates are also in aminoacids!
            dX = (3*domain.ali_from)/svg_options.scaling
            dL = 3*(domain.ali_to - domain.ali_from)/svg_options.scaling
            
            # Try to get colors from external file
            try:
                domain_color = models.colors[domain.ID]
            except KeyError:
                domain_color = ('255', '255', '255')
                domain_color_outline = ('200', '200', '200')
                print(domain.ID)
            else:
                h_, s, v = rgb_to_hsv(float(domain_color[0])/255.0, 
                        float(domain_color[1])/255.0, float(domain_color[2])/255.0)
                domain_color = tuple(str(int(round(c * 255))) for c in hsv_to_rgb(h_, s*0.7, v*1.1))
                domain_color_outline = tuple(str(int(round(c * 255))) for c in hsv_to_rgb(h_, s, 0.8*v))
            
            arrow.append("{}<g>".format(tabs*"\t"))
            tabs += 1
            
            arrow.append("{}<title>{} (acc)\ndesc</title>".format(tabs*"\t", domain.ID))
            if protein.forward:
                # calculate how far from head_start we (the horizontal guide at 
                #  y=Y+internal_domain_margin) would crash with the slope
                # Using similar triangles (triangle with head_length as base vs
                #  triangle with collision_x as base):
                collision_x = head_length * (h + internal_domain_margin)
                collision_x /= (h + H/2.0)
                collision_x = round(collision_x)
                
                # either option for 'x_margin_offset' works
                #  m = -float(h + H/2)/(head_length) #slope of right line
                #  x_margin_offset = (internal_domain_margin*sqrt(1+m*m))/m
                #  x_margin_offset = -(x_margin_offset)
                #x_margin_offset = internal_domain_margin/sin(pi - atan2(h+H/2.0,-head_length))

                # does the domain fit completely before the arrow head's diagonal?
                if (dX + dL) < head_start + collision_x:
                    d = "{}<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" stroke-linejoin=\"round\" fill=\"rgb({})\" stroke=\"rgb({})\" stroke-width=\"{}\" />".format(tabs*"\t", int(round(X+dX)), 
                        int(round(Y+internal_domain_margin)), dL, dH, ",".join(domain_color), 
                        ",".join(domain_color_outline), domain_contour_thickness)
                    arrow.append(d)
                else:
                    del points[:]
                    
                    # left part of the domain
                    if dX < head_start + collision_x:
                        # add point A
                        points.append((round(X + dX), round(Y + internal_domain_margin)))
                        
                        # add point B
                        points.append((round(X + head_start + collision_x), round(Y + internal_domain_margin)))
                    else:
                        # Domain is a triangle; add point A'
                        start_y_offset = (h + H/2)*(L - dX)
                        start_y_offset /= head_length
                        start_y_offset = int(start_y_offset)
                        points.append((round(X + dX), round(Y + H/2 - start_y_offset)))
                        
                    # handle the rightmost part of the domain
                    if dX + dL >= head_end: # could happen more easily with the scaling
                        # head of the domain's triangle
                        points.append((round(X + head_end), round(Y + H/2)))
                    else:
                        # add points C and D
                        end_y_offset = (2*h + H)*(L - dX - dL)
                        end_y_offset /= 2*head_length
                        end_y_offset = round(end_y_offset)

                        points.append((round(X + dX + dL), round(Y + H/2 - end_y_offset)))
                        points.append((round(X + dX + dL), round(Y + H/2 + end_y_offset)))
                
                    # handle lower part
                    if dX < head_start + collision_x:
                        # add points E and F
                        points.append((round(X + head_start + collision_x), round(Y + H - internal_domain_margin)))
                        points.append((round(X + dX), round(Y + H - internal_domain_margin)))
                    else:
                        # add point F'
                        points.append((round(X + dX), round(Y + H/2 + start_y_offset)))
                
                    del points_coords[:]
                    for point in points:
                        points_coords.append(str(int(point[0])) + "," + str(int(point[1])))
                        
                    d = "{}<polygon points=\"{}\" stroke-linejoin=\"round\" fill=\"rgb({})\" stroke=\"rgb({})\" stroke-width=\"{}\" />".format(tabs*"\t", " ".join(points_coords), ",".join(domain_color), ",".join(domain_color_outline),
                        domain_contour_thickness)
                    arrow.append(d)
            else:
                # Reverse orientation. Note that dX is relative to the _start_
                #  of the protein (now on the right!)
                # calculate how far from head_start we (the horizontal guide at 
                #  y=Y+internal_domain_margin) would crash with the slope
                # Using similar triangles:
                collision_x = head_length * ((H/2) - internal_domain_margin)
                collision_x /= (h + H/2.0)
                collision_x = round(collision_x)
                
                # The whole domain can be represented with a single block
                if L-dX-dL > collision_x:
                    d = "{}<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" stroke-linejoin=\"round\" fill=\"rgb({})\" stroke=\"rgb({})\" stroke-width=\"{}\" />".format(tabs*"\t", int(round(X+L-dX-dL)), 
                        int(round(Y+internal_domain_margin)), dL, dH, ",".join(domain_color), 
                        ",".join(domain_color_outline), domain_contour_thickness)
                    arrow.append(d)
                else:
                    del points[:]
                    
                    print(domain.ID)
                    print(X, collision_x, int(round(X+L-dX-dL)))
                    print()
                
            tabs -= 1
            arrow.append("{}</g>".format(tabs*"\t"))
            
            
        tabs -= 1
        arrow.append("{}</g>".format(tabs*"\t"))
        
        return arrow
    

    def SVG(self, models=HMM_DB(), svg_options=ArrowerOpts(), file_path=Path("./"), extra_label=""):
        """
        Writes an SVG figure for the cluster
        """
        
        h = svg_options.arrow_head_height
        H = svg_options.arrow_height
        scaling = svg_options.scaling
        #hl = svg_options.arrow_head_length
        #gene_contour_thickness = svg_options.gene_contour_thickness
        #idm = svg_options.internal_domain_margin
        #dct = svg_options.domain_contour_thickness
        fontsize = svg_options.fontsize
        margin = svg_options.label_arrow_margin
        stripe_thickness = svg_options.stripe_thickness
        stripe_color = svg_options.stripe_color
        
        svg_data = []

        space_for_labels = 0
        if svg_options.write_id or extra_label != "":
            space_for_labels = fontsize + margin
        
        # TODO: still produce a figure is no loci were annotated but we do have
        # proteins (i.e. consider every protein to be within a "mini-locus")
        #if len(self.loci) == 0 and len(self.protein_list) > 0:
        
        max_width = sum([x.length for x in self.loci])
        # give some space for inter-loci separator
        if len(self.loci) > 1:
            max_width += (len(self.loci)-1)*H
            
        # the extra "margin" for the height is due to the pointy lower parts of 
        # the arrows going further than the coordinate
        svg_data.append("<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" width=\"{:d}\" height=\"{:d}\">".format(int((max_width)/scaling), (2*h + H + space_for_labels + margin)))
        
        tabs = 1
        
        svg_data.append("{}<g>".format(tabs*"\t"))
        tabs += 1
        
        x = 0
        y = 0
        
        svg_data.append("{}<title>{}</title>".format(tabs*"\t",self.identifier))
        if svg_options.write_id or extra_label != "":
            svg_data.append("{}<text x=\"10\" y=\"{}\" font-size=\"{}\">".format(tabs*"\t", fontsize, fontsize))
            y += fontsize + margin
            tabs += 1
            
            label = ""
            if svg_options.write_id:
                label = self.identifier
                if extra_label != "":
                    label += " "
            if extra_label != "":
                label += extra_label
            svg_data.append("{}{}".format(tabs*"\t",label))
            
            tabs -= 1
            svg_data.append("{}</text>".format(tabs*"\t"))
        
        n_locus = 1
        stripe_position = y + h + (H/2)
        for locus in self.loci:
            svg_data.append("{}<g>".format(tabs*"\t"))
            tabs += 1
            
            svg_data.append("{}<line x1=\"{:d}\" y1=\"{:d}\" x2=\"{:2}\" y2=\"{:2}\" style=\"stroke:rgb({}); stroke-width:{:d}\" />".format(tabs*"\t", int(x), int(stripe_position), int(x+locus.length/svg_options.scaling), int(stripe_position), ",".join(map(str,stripe_color)), stripe_thickness))
            
            for p in range(len(locus.protein_list)):
                protein = locus.protein_list[p]
                gene_coordinates = locus.gene_coordinates[p]
                
                svg_data.extend(self.SVG_arrow(x+gene_coordinates[0], y+h, tabs, protein, models, svg_options))
            
            x += locus.length
            # TODO print locus separator mini-fig
            if n_locus < len(self.loci):
                x += H
            
            tabs -= 1
            svg_data.append("{}</g>".format(tabs*"\t"))
            
            n_locus += 1
        
        tabs -= 1
        svg_data.append("{}</g>".format(tabs*"\t"))
        svg_data.append("</svg>")
        
        
        # get this BGC written
        with open(file_path / (self.identifier + ".svg"), "w") as f:
            f.write("{}".format("\n".join(svg_data)))
            
        return
    

class BGCLocus:
    """
    BGC might be divided in more than one locus (because of sequencing, or, as
    is more often the case in fungi, because of genomics)
    
    This class can be used to organize proteins
    """
    def __init__(self):
        self.protein_list = []
        self.gene_coordinates = []  # a list of tuples. end-start might not match
                                    # the corresponding protein lenght due to 
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
        
        self.role = "unknown"       # e.g. [biosynthetic, transporter, tailoring, 
                                    #   resistance, unknown]
        self._CBP_type = "unknown"  # should be one from "valid_CBP_types"
        self.s_cbp = ""             # Sub classes of CBP e.g.: Group V
        self.ss_cbp = ""            # e.g.: Group V2
        
        # e.g.: nrPKS, nrPKS+hrPKS etc.
        # TODO: redifine as "Gene Cluster Architecture"?
        # then GCF would be defined by the list of CBP-signatures
        self._gca = []
        
        
        self.compound_family = ""   # Compound family e.g. "emodin-like"
        self.compound = ""
        self.source = "unknown"     # e.g. MIBiG, curated document etc.
        self.organism = ""
        self.TaxId = ""
        
        self.forward = True
        
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
        return "+".join((self._gca))    
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
            return "|--[" + "]-[".join(domain_alias_list) + "]-->\t{}".format(tag)
        else:
            return "<--["+ "]-[".join(list(reversed(domain_alias_list))) + "]--|\t{}".format(tag)
        
        
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
        
        if len(self.domain_list) < 2:
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
                
            elif len(self.domain_set & set({"TIGR04532","SAT"})) > 0:
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


    def domain_SVG(self, filename, hmmdb):
        """Creates an SVG figure with the domains as boxes
    
        in:
            name: name of the file including path
            hmmdb: contains alias and color information
            
        out:
            an svg file
        """
    
        scaling = 2
        H = 20
        h = H/2
        stripe_thickness = 2
        curve_radius = int(H/10)
        
        data = []
        data.append("<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\">\n".format(int(self.length/scaling), (2*h + H)))
        
        data.append("\t<g>\n")
        data.append("\t<title>{}</title>\n".format(self.identifier))
        data.append("\t<line x1=\"0\" y1=\"{:d}\" x2=\"{:2}\" y2=\"{:2}\" style=\"stroke:rgb({}); stroke-width:{:d}\" />\n".format(int(h+H/2), int(self.length/scaling), int(h+H/2), "10,10,10", stripe_thickness))
        
        for domain in self.domain_list:            
            try:
                color = ",".join([str(x) for x in hmmdb.colors[domain.ID]])
            except KeyError:
                color = "150,150,150"
                
            try:
                title = hmmdb.alias[domain.ID]
            except KeyError:
                title = domain.ID
                
            data.append("\t\t<g>\n")
            data.append("\t\t<title>{}</title>\n".format(title))
            data.append("\t\t\t<rect x=\"{}\" y=\"0\" rx=\"{}\" ry=\"{}\" width=\"{}\" height=\"{}\" fill=\"rgb({})\"/>\n".format(int(domain.ali_from/scaling), curve_radius, curve_radius, int((domain.ali_to - domain.ali_from)/scaling), 2*H, color))
            data.append("\t\t</g>\n")
                
        data.append("\t</g>\n")
        data.append("</svg>\n")
        
        with open(filename, "w") as f:
            for row in data:
                f.write(row)
    

class BGCDomain:
    def __init__(self, protein, ID, env_from, env_to, ali_from, ali_to, hmm_from, hmm_to, score, Evalue):
        self.protein = protein
        self.ID = ID                # e.g. PF00501.27
        self.env_from = env_from    # Pos. in target seq. at which surr. envelope 
        self.env_to = env_to        #   starts/ends.
        self.ali_from = ali_from    # Position in target sequence at which the 
        self.ali_to = ali_to        #   hit starts/ends
        self.hmm_from = hmm_from    # Position in the hmm at which the hit 
        self.hmm_to = hmm_to        #   starts/ends
        self.score = score          # Only depends on profile HMM and target seq.
        self.Evalue = Evalue        # Based on score and DB size
        
    def get_sequence(self):
        return self.protein.sequence[ali_from:ali_to+1]
