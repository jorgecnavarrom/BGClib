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
__version__ = "0.1"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@westerdijkinstitute.nl"

valid_CBP_types = set({"nrPKS", "rPKS", "NRPS", "t3PKS", "unknown", "other",
        "PKS-NRPS_hybrid", "NRPS-PKS_hybrid", "other_PKS", "unknown_PKS", 
        "no_domains"})

class HMM_DB:
    """
    This class basically keeps a list of hmm databases to be used later
    """
    
    def __init__(self):
        # list of paths to hmm databases
        self.db_list = []
        
        # TODO should also contain alias for important models and dictionaries
        # to go from AC to ID or AC to DESC
        
        return
    
    def addDatabase(self, db_path):
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


class BGCCollection:
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
        self.parent_cluster = None    # Should point to an object of the BGC class
        
        self.identifier = ""        # Defined by user. should be unique
        
        self.accession = ""         # original NCBI accession if possible
        self.ref_accession = ""     # "RefSeq Selected Product" from NCBI's IPG DB
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
        
        self.hmm_db = None
        self.domain_list = []
        self.domain_set = set()
        self.attempted_domain_prediction = False
        
        return
    
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
        return "\n".join([self._Sequence[i:i+80] for i in range(0, self.length, 80)])
            
    def fasta(self):
        compound = ""
        if self.compound != "":
            compound = " {}".format(self.compound)
        if self.ref_accession != "":
            return ">{}{}\n{}".format(self.ref_accession, compound, self.fastasizedSeq())
        else:
            return ">{}{}\n{}".format(self.accession, compound, self.fastasizedSeq())

    def set_hmmdb(self, hmmdb):
        assert(isinstance(hmmdb, HMM_DB))
        self.hmm_db = hmmdb
        return

    def predict_domains(self, hmmdb):
        assert(isinstance(hmmdb, HMM_DB))
        return
    
    def domain_string(domain_alias, domains, seq_id):
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
        domain_name_list = [x.ID for x in domain_list]
        
        # try to convert names to alias
        domain_alias_list = []
        for d in domain_name_list:
            if d in domain_alias:
                domain_alias_list.append(domain_alias[d])
            else:
                domain_alias_list.append(d)
        
        return "|--[" + "]-[".join(domain_alias_list) + "]-->\t{}".format(seq_id)
    
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
            if self.domain_list[0] in PKS_domains:
                sequence_type = "PKS-NRPS_hybrid"
            elif self.domain_list[0] in NRPS_domains:
                sequence_type = "NRPS-PKS_hybrid"
            else:
                print("PKS/NRPS or NRPS/PKS hybrid?")
                print(" - ".join(self.domain_list))
                sequence_type = "NRPS-PKS_hybrid"
        # nrPKS or (hr/prPKS)
        elif len(domain_set & PKS_domains) > 0:
            # this is not expected
            if len(domain_set & set({"TIGR04532","SAT"})) > 0 and len(domain_set & reducing_domains) > 0:
                sequence_type = "unknown_PKS"
            elif len(domain_set & set({"TIGR04532","SAT"})) > 0:
                # assume that having a SAT domain is enough for labeling as nrPKS
                # but note that in this category, there seem to be three cases for PT:
                # - PT detected by TIGR04532
                # - DH detected by pfam PF14765 (TIGR doesn't pass --cut_tc
                # - no detection
                # There are either several types of PT domains or not a general model
                # to detect all three+ cases
                sequence_type = "nrPKS"
            elif len(domain_set & reducing_domains) > 0:
                sequence_type = "rPKS"
            else:
                sequence_type = "other_PKS"
        elif len(domain_set & PKS3_domains) > 0:
            sequence_type = "PKS_III"
        elif len(domain_set & NRPS_domains) > 0:
            sequence_type = "NRPS"
        else:
            sequence_type = "other"

        self._CBP_type = sequence_type
        
        return


class BGCDomain:
    def __init__(self, protein, AC, env_from, env_to, ali_from, ali_to, hmm_from, hmm_to, score, Evalue):
        self.protein = protein
        self.AC = AC                # e.g. PF00501.27
        self.env_from = env_from    # Pos. in target seq. at which surr. envelope 
        self.env_to = env_to        #   starts/ends.
        self.ali_from = ali_from    # Position in target sequence at which the 
        self.ali_to = ali_to        #   hit starts/ends
        self.hmm_from = hmm_from    # Position in the hmm at which the hit 
        self.hmm_to = hmm_to        #   starts/ends
        self.score = score          # Only depends on profile HMM and target seq.
        self.Evalue = Evalue        # Based on score and DB size
        
