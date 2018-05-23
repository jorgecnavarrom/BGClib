import sys

class BGCProtein:
    """
    Information about a Protein encoded by a gene which is part of the gene cluster
    v0.1
    
    variables:
        - a tuple with all domains' ids

    methods:
        - predict domains. 
            Needs:
                a DomainProfiles object
                a path where to deposit domtable results. .domtable
                a path where to put filtered domtable (non-overlapping) results. .pfd
    
    """
    
    # nrPKS, rPKS, NRPS, t3PKS, PKS/NRPS Hybrids, other_PKS, unknown_PKS
    
    valid_CBP_types = set({"nrPKS", "rPKS", "NRPS", "t3PKS", "unknown", "other",
                                "PKS-NRPS_hybrid", "other_PKS", "unknown_PKS"})
    
    
    def __init__(self):
        # original NCBI accession
        self.Accession = ""
        
        # NCBI reference accession e.g. "RefSeq Selected Product"
        # Reference given by the Identical Protein Groups NCBI database
        self.RefAccession = ""
            
        # NCBI internal id in the Protein database
        self.NCBI_id = ""
        
        # NCBI internal id in the Identical Protein Group database. Should be unique
        self.NCBI_ipg_id = ""
        
        #self.Sequence = None
        self.length = 0
        
        # basic annotation: biosynthetic, transporter, tailoring, resistance, unknown
        self.SM_class = "unknown"
        
        # Should point to an object of the BGC class
        self.parentCluster = ""
        
        # internally, this variable uses an underscore to differentiate from function
        self._CBP_type = "unknown"
        
        # e.g.: nrPKS, nrPKS+hrPKS etc.
        self._GCF = set()
        
        # e.g.: Group V
        self.sGCF = ""
        
        # e.g.: Group V2
        self.ssGCF = ""
        
        # Compound family e.g. "emodin-like"
        self.CompoundFamily = ""
        
        self.Compound = ""
        
        # Data source (e.g. MIBiG, curated document etc.)
        self.Source = "unknown"
        
        # From original file
        self.Organism = ""
        
        self.TaxId = ""
        
        self.forward = True
        
        self.domain_list = []
        self.domain_set = set()
            
        return
    
    @property
    def CBP_type(self):
        return self._CBP_type
    
    @CBP_type.setter
    def CBP_type(self, cbptype):
        try:
            assert cbptype in self.valid_CBP_types
        except AssertionError:
            print("{} is not a valid CBP type".format(cbptype))
        self._CBP_type = cbptype
    
    @property
    def GCF(self):
        return "+".join(sorted(self._GCF))
    
    @GCF.setter
    def GCF(self, gcf):
        if isinstance(gcf, set):
            self._GCF = gcf
        else:
            self._GCF.add(gcf)
    
    @property
    def Sequence(self):
        return self._Sequence
        
    @Sequence.setter
    def Sequence(self, seq):
        # note: does not check if contains valid amino acid characters yet
        self._Sequence = seq.replace("\n", "").strip().replace(" ", "")
        self.length = len(seq)
    
    def annotations(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.parentCluster, self.Accession, self.RefAccession, self.NCBI_id, 
            self.NCBI_ipg_id, self.CBP_type, self.GCF, self.sGCF, self.ssGCF, 
            self.CompoundFamily, self.Compound, self.Source, self.Organism, 
            self.TaxId)
    
    def fastasizedSeq(self):
        return "\n".join([self._Sequence[i:i+80] for i in range(0, self.length, 80)])
            
    def fasta(self):
        if self.RefAccession != "":
            return ">{}\n{}".format(self.RefAccession, self.fastasizedSeq())
        else:
            return ">{}\n{}".format(self.Accession, self.fastasizedSeq())
