#!/usr/bin/env python3

"""
BGC library: BGC

Handles information about a biosynthetic gene cluster
"""

from .protein import BGCProtein
from .collections.proteins import ProteinCollection
from .locus import BGCLocus
from .metadata.organism import Organism
from .metadata.metabolite import Metabolite
from .annotations.hmm import hmmDB
from .visualization.arrower import ArrowerOpts
from .data.utilities import *
from pathlib import Path
from lxml import etree
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import pickle
from operator import itemgetter
from copy import deepcopy

class BGC:
    def __init__(self, gbkfile=None):
        self.identifier = ""        # usually the file name
        
        self.CBPtypes = []          # Core Biosynthetic Protein List: simple
                                    #  5'-3' ordered list of biosynthetic types
        self.CBPtypes_set = set()
        self.CBPcontent = {}        # Core Biosynthetic Protein Content. Every 
                                    #  biosynthetic-type returns a list of 
                                    #  pointers to each Protein object
        
        self.products = set()       # All different 'product' qualifiers  
                                    #  annotated by antiSMASH
        self.contig_edge = False    # antiSMASH v4+ was not able to fully 
                                    #  complete the extension phase. BGC might
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
            self.load(gbkfile)
        
        
    def load(self, gbk) -> None:
        """
        Initializes the object by reading all CDS fields from a GenBank file
        
        Args:
            Input is intended to be a Path() object (but it should work with a 
        string)
        """

        _gbk = Path(gbk)
        assert _gbk.is_file(), f"{gbk} not a valid file"
        self.identifier = _gbk.stem
        
        # # TODO: get lineage? implement the rest of the parsing using gb-io
        # for locus_num, record in enumerate(gb_io.iter(str(_gbk))):
        #     # find organism properties
        #     source_list = list(filter(
        #         lambda feat: feat.kind == 'source', 
        #         record.features
        #     ))
        #     self.organism = Organism()
        #     if len(source_list) == 1:
        #         ft = source_list[0]
        #         print(ft.qualifiers)
        #         # this seems ineficient:
        #         quals = dict((q.key, q.value) for q in ft.qualifiers)
        #         self.organism.taxid = quals['db_xref'].split(":")[-1]
        #         self.organism.fullname = quals['organism']
        #     elif len(source_list) > 1:
        #         print(f"Warning: {_gbk} has more than one 'source' feature")
            
        #     print('info')
        #     print(dir(record))
        #     print(record.accession)
        #     print(record.definition)
        #     print(record.division)
        #     print(record.name)
        #     print("features:\n")
        #     for feature in record.features:
        #         print(feature.kind)
        #         print(feature.location)
        #         print(feature.qualifiers)

        # exit()
        try:
            records = list(SeqIO.parse(str(gbk), "genbank"))
        except ValueError as e:
            e.add_note(f"Error, not able to parse file {str(gbk)}: {e}")
            raise

        self.accession = records[0].id
        self.definition = records[0].description

        if records[0].annotations["organism"] != "":
            self.organism = Organism()
            self.organism.fullname = records[0].annotations["organism"]
            self.organism.lineage = records[0].annotations["taxonomy"]
        
        cds_list = []
        aS5 = True # Whether the BGC prediction comes from antiSMASH 5+

        # traverse all possible records in the file. There's usually only 1
        for locus_num, record in enumerate(records):
            # collect all biosynthetic regions' coordinates (probably only a 
            # single one per gbk, but it's possible to read a whole genome's
            # genbank), then traverse all CDS to see which need to be marked 
            biosynthetic_regions = []
            locus = BGCLocus()
            locus.identifier = f"{self.identifier}~{locus_num}"
            locus.length = len(record)
            cds_num = 0

            for feature in record.features:
                if feature.type == 'source':
                    if "db_xref" in feature.qualifiers:
                        taxon = feature.qualifiers['db_xref'][0]
                        taxon_id = taxon.strip().split(":")[-1]
                        self.organism.taxid = taxon_id

                # antiSMASH <= 4
                if feature.type == "cluster":
                    aS5 = False
                    if "product" in feature.qualifiers:
                        for product in feature.qualifiers["product"]:
                            for p in product.replace(" ","").split("-"):
                                self.products.add(p)
                            
                    if "contig_edge" in feature.qualifiers:
                        if feature.qualifiers["contig_edge"][0] == "True":
                            self.contig_edge = True
                            
                    biosynthetic_regions.append(feature.location)
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

                    biosynthetic_regions.append(feature.location)
                    continue
                            
                if feature.type != "CDS": continue
                
                cds_num += 1
                
                CDS = feature
                
                cds_start = max(0, int(CDS.location.start))
                cds_end = max(0, int(CDS.location.end))

                identifier = f"{self.identifier}~L{locus_num}+CDS{cds_num}"

                product = ""
                if "product" in CDS.qualifiers:
                    product = ", ".join(CDS.qualifiers["product"])
                    
                # NOTE: JGI annotations have a non-standard qualifier:
                # "proteinId" or "proteinID", which is just a number 
                # (sometimes). 
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

                # TODO: test this. It's a bit tricky with antiSMASH 5 now 
                # with regions
                protein_type = ""
                if aS5 and role == "biosynthetic" and "gene_functions" \
                        in CDS.qualifiers:
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
                # TODO: partial compatibility with older antiSMASH results
                elif "sec_met" in CDS.qualifiers:
                    sec_met = {}
                    for item in CDS.qualifiers["sec_met"]:
                        split_item = item.split(": ")
                        if len(split_item) == 2:
                            sec_met[split_item[0]] = split_item[1]
                        else:
                            sec_met[split_item[0]] = ": ".join(
                                split_item[1:]
                                )
                    if "Type" in sec_met:
                        if sec_met["Type"] != "none" \
                                and sec_met["Kind"] == "biosynthetic":
                            protein_type = sec_met["Type"]
                            role = "biosynthetic"

                protein = BGCProtein()
                
                protein.identifier = identifier
                protein.product = product
                protein.protein_id = protein_id
                protein.gene = gene
                protein.role = role
                protein.protein_type = protein_type
                
                if "translation" not in CDS.qualifiers:
                    # print(" Warning. Generating translation for " \
                    #   + f"{identifier} (standard codon table)")
                    # TODO: there should be a way to warn user without 
                    # spamming terminal with prints
                    protein.sequence = CDS.translate(record.seq)
                else:
                    protein.sequence = CDS.qualifiers["translation"][0]
                # TODO: verify sequence doesn't contain stop codons at the 
                # end (anywhere?)

                # TODO: verify that 
                # len(CDS.location.parts) == len(protein.sequence)
                
                # NOTE what happens if no strand info (strand=None)?
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
                protein.cds_regions = tuple(
                    sorted(cds_list, key=itemgetter(0))
                    )

                self.proteins[protein.identifier] = protein
                locus.protein_list.append(protein)
                locus.gene_coordinates.append((cds_start,cds_end))


            # CDSs will probably be ordered, but not guaranteed
            locus.protein_list.sort(key=lambda p: p.cds_regions[0][0])
            self.protein_list.extend(locus.protein_list)

            # Mark CDSs that are inside a biosynthetic region
            biosynthetic_regions.sort(key=lambda r: r.start)
            p_index = 0
            for region in biosynthetic_regions:
                for protein in locus.protein_list[p_index:]:
                    cds_start = protein.cds_regions[0][0]
                    cds_end = protein.cds_regions[-1][1]

                    if cds_end <= region.start:
                        # haven't reach region
                        p_index += 1
                        continue
                    elif cds_start >= region.end:
                        # beyond current region, skip to the next region, but
                        # don't advance p_index
                        break
                    else:
                        p_index += 1
                        protein.in_biosynthetic_region = True

            self.loci.append(locus)
            locus_num += 1

        self.set_CBP_content()

        return
    
    
    # TODO: finish
    def load_fasta(self, fasta: Path) -> None:
        """Creates faux BGCs from a fasta file.
        
        Each entry in the fasta file will be treated as a separate locus.
        Assumes that the fasta file contains a protein sequence.
        """
        return
    
    
    def inter_loci_element(
            self, 
            xoffset: float, 
            yoffset: float, 
            H: float,
            track_thickness: int
        ):
        """Creates an inter-loci figure

        Args:
            xoffset (float): Offset X (increases when going to the right)
            yoffset (float): Offset Y (increases when going down)
            H (float): Arrow main body height
            track_thickness (int): Thickness of the genomic track

        Draws the following SVG figure as xml code to link loci genomic 
        stripes:
        __/ /__
         / /
        
        and returns the corresponding etree element.
        
        xoffset/yoffset correspond to position of point A in an arrow
        """
        
        
        inter_loci_main = etree.Element("g")
        
        inter_loci_attribs = {
            "stroke": "#464646",
            "stroke-width": str(track_thickness)
            }
            
        # left
        inter_loci_attribs["d"] = f"M{int(xoffset)} {int(yoffset + 0.5*H)} " \
            + f"L {int(xoffset+0.375*H)} {int(yoffset + 0.5*H)}"
        inter_loci_main.append(
            etree.Element("path", attrib=inter_loci_attribs)
        )
        
        # right
        inter_loci_attribs["d"] = f"M{int(xoffset+0.625*H)} " \
            + f"{int(yoffset + 0.5*H)} L {int(xoffset+H)} " \
            + f"{int(yoffset + 0.5*H)}"
        inter_loci_main.append(
            etree.Element("path", attrib=inter_loci_attribs)
        )
        
        # slash1
        inter_loci_attribs["d"] = f"M{int(xoffset + 0.25*H)} " \
            + f"{int(yoffset + 0.75*H)} L {int(xoffset + 0.5*H)} " \
            + f"{int(yoffset + 0.25*H)}"
        inter_loci_main.append(
            etree.Element("path", attrib=inter_loci_attribs)
        )
        
        # slash2
        inter_loci_attribs["d"] = f"M{int(xoffset + 0.5*H)} " \
            + f"{int(yoffset + 0.75*H)} L {int(xoffset + 0.75*H)} " \
            + f"{int(yoffset + 0.25*H)}"
        inter_loci_main.append(
            etree.Element("path", attrib=inter_loci_attribs)
        )
        
        return inter_loci_main
    

    # TODO: implement the extra_label parameter. The main issue is how to 
    # calculate the width of the text element to take into account for the SVG 
    # width
    def SVG(
            self, 
            filepath: Path, 
            hmmdb = hmmDB(), 
            svg_options = ArrowerOpts(), 
            extra_label = "", 
            xoffset = 0.0, 
            yoffset = 0.0, 
            mirror = False
        ) -> None:
        """Writes an SVG figure of the cluster

        Args:
            filepath (Path): Name of the output file including path. \
                Pre-existance of subfolders is not checked
            hmmdb (hmmDB): If BGC's proteins are annotated and need to be \
                painted, this object provides alias and color information
            svg_options (ArrowerOpts): Object that includes all relevant \
                drawing options
            extra_label (str): [not used yet]
            xoffset (float): Offset X (increases when going to the right)
            yoffset (float): Offset Y (increases when going down)
            mirror (bool): True: flip the figure horizontally
        """
        
        # Arrows need to be represented in their original orientation. 
        # Be careful as not to modify input svg options though
        bgc_svg_options = deepcopy(svg_options)
        bgc_svg_options.original_orientation = True
        
        H = bgc_svg_options.arrow_height
        margin = bgc_svg_options.topbottom_margin
        
        # length of the SVG figure
        # second term accounts for the inter-loci figure
        L = sum(
            [locus.length/bgc_svg_options.scaling for locus in self.loci]
            ) + H*(len(self.loci)-1)

        # main node
        # Note: we are in dna space (so amino acid sequence must be scaled)
        # width, height, xoffset and yoffset must be corrected or the corners
        # with sharp angles will be chopped off
        thickness = bgc_svg_options.gene_contour_thickness
        Xoffset = xoffset + thickness
        Yoffset = yoffset + thickness
        base_attribs = {
            "version":"1.1", 
            "baseProfile":"full",
            "width":str(int(Xoffset + L + thickness))
        }
        root = etree.Element(
            "svg", 
            attrib=base_attribs, 
            nsmap={None:'http://www.w3.org/2000/svg'}
        )

        fig_height = margin + thickness + H
        if bgc_svg_options.shape == 'Arrow':
            fig_height += H
        root.set('height', str(int(Yoffset + fig_height)))

        inner_tree = self.xml_BGC(
            Xoffset, Yoffset + margin/2, hmmdb, bgc_svg_options, mirror
            )
        root.append(inner_tree)
        
        with open(filepath, "bw") as f:
            f.write(etree.tostring(root, pretty_print=True))


    def xml_BGC(
            self, 
            xoffset = 0.0, 
            yoffset = 0.0, 
            hmmdb = hmmDB(), 
            svg_options = ArrowerOpts(), 
            mirror = False
        ):
        """Creates the internal xml structure of a BGC SVG

        Args:
            hmmdb (hmmDB): If BGC's proteins are annotated and need to be \
                painted, this object provides alias and color information
            svg_options (ArrowerOpts): Object that includes all relevant \
                drawing options
            extra_label (str): [not used yet]
            xoffset (int): Offset X (increases when going to the right)
            yoffset (int): Offset Y (increases when going down)
            mirror (bool): True: flip the figure horizontally
        """
        main_group = etree.Element("g")
        
        bgc_title = etree.Element("title")
        bgc_title.text = self.identifier
        main_group.append(bgc_title)
        
        H = svg_options.arrow_height
        track_thickness = svg_options.stripe_thickness
        h = 0.5*H # body to head vertical tip
        if svg_options.shape == "Ribbon":
            h = 0
        
        current_locus = 1
        Xoffset = xoffset # start of each locus
        loci_list = self.loci
        if mirror:
            loci_list = reversed(self.loci)

        # genomic track vertical position
        y_str = str(int(yoffset + h + 0.5*H))
        for locus in loci_list:
            L = locus.length/svg_options.scaling
            
            line_attribs = {
                "x1": str(int(Xoffset)),
                "y1": y_str,
                "x2": str(int(Xoffset+L)),
                "y2": y_str,
                "stroke": "#464646",
                "stroke-width": str(int(track_thickness))
                }
            line = etree.Element("line", attrib=line_attribs)
            main_group.append(line)
            
            for protein in locus.protein_list:
                gene_start = protein.cds_regions[0][0] 
                gene_end = protein.cds_regions[-1][1]
                
                gene_position = gene_start / svg_options.scaling
                if mirror:
                    gene_length = (
                        gene_end - gene_start - 3
                    )/svg_options.scaling
                   
                    gene_position = L - gene_position - gene_length
                
                inner_tree = protein.xml_arrow(
                    hmmdb, 
                    svg_options, 
                    Xoffset + gene_position, 
                    yoffset, 
                    mirror=mirror
                )
                main_group.append(inner_tree)
                
                
            Xoffset += L
            
            if current_locus < len(self.loci):
                # draw inter-locus fig of length H
                main_group.append(
                    self.inter_loci_element(
                        Xoffset, yoffset + h, H, track_thickness
                    )
                )
            
                Xoffset += H
            current_locus += 1
        
        return main_group


    def classify_proteins(self):
        """Assigns (SM) role and protein_type to every protein in BGC

        Uses predicted domains to classify (Secondary Metabolism) role and 
        protein_type of all its proteins.
        """
        
        for p in self.protein_list:
            p.classify_sequence()
        
        self.set_CBP_content()

        return

    
    def set_CBP_content(self):
        """Annotates internal CBP-related variables

        Role of proteins needs to be set beforehand
        """
        
        del self.CBPtypes[:]        # ordered list
        self.CBPtypes_set.clear()   # set of all CBTs
        self.CBPcontent.clear()     # dict, CBT to list of BGCProtein
        for protein in self.protein_list:
            if protein.role == "biosynthetic":
                self.CBPtypes.append(protein.protein_type)
                try:
                    self.CBPcontent[protein.protein_type].append(protein)
                except KeyError:
                    self.CBPcontent[protein.protein_type] = [protein]
        self.CBPtypes_set = set(self.CBPtypes)


    def predict_domains(
            self, 
            hmmdb: hmmDB, 
            cpus=1, 
            tc=True, 
            filterdoms=True
            ):
        """Annotate domains for this BGC

        Args:
            hmmdb: An hmmDB object pointing to a valid hmm database
            cpus: Number of cpus for hmm detection
            tc: (bool) Use trusted cutoffs (score of the lowest-scoring known\
                true positive)
            filterdoms: (bool) Keep only the highest-scoring domain if over-\
                lapping predictions found  

        Make a a protein collection from the proteins in this BGC and annotate
        their domains using the provided hmm databases
        """
        assert isinstance(hmmdb, hmmDB)
        
        pc = ProteinCollection()
        
        missing_identifier = False

        for protein in self.protein_list:
            assert protein.identifier
            pc.proteins[protein.identifier] = protein
            
        pc.predict_domains(hmmdb, cpus=cpus, tc=tc, filterdoms=filterdoms)
        
        self.attempted_domain_prediction = True
        self.calculate_domain_sets()

        return
    

    def calculate_domain_sets(self) -> None:
        """Annotates the domain set internal variables
        
        Three sets are created:
        - `domain_set`: all domains in BGC
        - `domain_set_core`: only domains found in proteins with 'biosynthetic'
            role
        - `domain_set_complement`: domains outside proteins with 'biosynthetic'
            role
        """
        self.domain_set = set()
        for protein in self.protein_list:
            self.domain_set.update(protein.domain_set)
            if protein.role == "biosynthetic":
                self.domain_set_core.update(protein.domain_set)
            else:
                self.domain_set_complement.update(protein.domain_set)
    

    def clear_domain_predictions(self):
        """Remove all domain information for all the BGC's proteins
        """

        for protein in self.protein_list:
            protein.domain_list = []
            protein.domain_set = set()

        self.domain_set = set()
        self.domain_set_core = set()
        self.domain_set_complement = set()
        
        self.attempted_domain_prediction = False
