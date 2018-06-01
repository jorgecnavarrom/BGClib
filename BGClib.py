#!/usr/bin/env python

"""
BGC library

Classes for data storage and analysis of Biosynthetic Gene Clusters

"""

import sys
from pathlib import Path
from subprocess import PIPE, Popen
import warnings
from io import StringIO
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO
from Bio import SeqIO

__author__ = "Jorge Navarro"
__version__ = "0.2.1"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@westerdijkinstitute.nl"

valid_CBP_types = set({"nrPKS", "rPKS", "NRPS", "t3PKS", "unknown", "other",
        "PKS-NRPS_hybrid", "NRPS-PKS_hybrid", "other_PKS", "unknown_PKS", 
        "no_domains"})

class HMM_DB:
    """
    This class keeps information about HMM databases to be used by the other
    classes
    """
    
    def __init__(self):
        self.db_list = []           # list of paths to hmm databases
        self.alias = {}             # ID to alias
        self.colors = {}            # ID to tuple(r,g,b)
        
        return
    
    def add_database(self, db_path):
        if not db_path.is_file():
            print("Not able to add hmm database (not a file. Wrong path?)")
            return
        elif db_path.suffix.lower() != ".hmm":
            print("Not able to add hmm datase (not a .hmm file)")
            return
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
        return

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



class BGCCollection:
    """
    This class will allow implementarion of collection-wide functions such as 
    single-step prediction of domains and comparisons between collections
    """
    
    def __init__(self):
        pass



class BGC:
    def __init__(self):
        self.identifier = ""    # usually the file name
        self.gca = []       # Gene Cluster Architecture. List of CBP types in the cluster
        self.gcf = []       # should be a list of CBP signatures (numbers?)



class BGCLocus:
    """
    BGC might be divided in more than one locus (because of sequencing, or, as
    is more often the case in fungi, because of genomics)
    
    This class can be used to organize proteins
    """
    def __init__(self):
        proteins = []
        pass



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
        self.length = 0
        
        self.role = "unknown"       # e.g. biosynthetic, transporter, tailoring, 
                                    #   resistance, unknown
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
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.parent_cluster, self.gca, self.accession, self.ref_accession, self.ncbi_id, 
            self.ncbi_ipg_id, self.CBP_type, self.s_cbp, self.ss_cbp, 
            self.compound_family, self.compound, self.source, self.organism, 
            self.TaxId)
    
    
    def fastasizedSeq(self):
        return "\n".join([self._sequence[i:i+80] for i in range(0, self.length, 80)])
        
        
    def fasta(self):
        compound = ""
        if self.compound != "":
            compound = " {}".format(self.compound)
        if self.ref_accession != "":
            return ">{}{}\n{}".format(self.ref_accession, compound, self.fastasizedSeq())
        else:
            return ">{}{}\n{}".format(self.accession, compound, self.fastasizedSeq())


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
            
        # finally, sort by env_from
        self.domain_list.sort(key=lambda x:x.ali_from)
        
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
        self.domain_set = set([d.ID for d in self.domain_list])
        
        return
    
    
    def classify_sequence(self):
        """Classifies a sequence based on its predicted domains according to a 
        set of rules
        """
        
        if len(self.domain_list) == 0:
            if self.attempted_domain_prediction:
                print(" Domain prediction already tried, but no hits found")
                self._CBP_type = "no_domains"
                return
            else:
                self.predict_domains()
                self.attempted_domain_prediction = True
        
        # Basic filtering
        if len(self.domain_set) == 1:
            return "other"
            #if "Trp_DMAT" in domain_set:
                #return "DMAT"
            #else:
                #return "other"
        
        # this should all have an alias in CBP_domains.tsv
        PKS_domains = set({"SAT", "ketoacyl-synt", "Ketoacyl-synt_C", "KAsynt_C_assoc",
                        "Acyl_transf_1", "TIGR04532", "Thioesterase", "Methyltransf_12"})
        PKS3_domains = set({"Chal_sti_synt_N", "Chal_sti_synt_C"})
        NRPS_domains = set({"Condensation", "AMP-binding", "AMP-binding_C"})
        reducing_domains = set({"PKS_ER_names_mod", "KR", "PS-DH"})
        #Terpene_domains = set({"Terpene_synth", "Terpene_synth_C"})
        #Squalene_domains = set({"SQHop_cyclase_N", "SQHop_cyclase_C"})
        
        sequence_type = ""

        # pks/nrps hybrid
        if len(self.domain_set & PKS_domains) > 0 and len(self.domain_set & NRPS_domains) > 0:
            if self.domain_list[0].ID in PKS_domains:
                sequence_type = "PKS-NRPS_hybrid"
            elif self.domain_list[0].ID in NRPS_domains:
                sequence_type = "NRPS-PKS_hybrid"
            else:
                print("PKS/NRPS or NRPS/PKS hybrid?")
                print(self.domain_string({}))
                sequence_type = "NRPS-PKS_hybrid"
        # nrPKS or (hr/prPKS)
        elif len(self.domain_set & PKS_domains) > 0:
            # this is not expected
            #if len(self.domain_set & set({"TIGR04532","SAT"})) > 0 and len(self.domain_set & reducing_domains) > 0:
                #sequence_type = "unknown_PKS"
            
            if len(self.domain_set & set({"TIGR04532","SAT"})) > 0:
                # assume that having a SAT domain is enough for labeling as nrPKS
                # but note that in this category, there seem to be three cases for PT:
                # - PT detected by TIGR04532
                # - DH detected by pfam PF14765 (TIGR doesn't pass --cut_tc
                # - no detection
                # There are either several types of PT domains or not a general model
                # to detect all three+ cases
                sequence_type = "nrPKS"
            elif len(self.domain_set & reducing_domains) > 0:
                sequence_type = "rPKS"
            else:
                sequence_type = "other_PKS"
        elif len(self.domain_set & PKS3_domains) > 0:
            sequence_type = "PKS_III"
        elif len(self.domain_set & NRPS_domains) > 0:
            sequence_type = "NRPS"
        else:
            sequence_type = "other"

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
        
