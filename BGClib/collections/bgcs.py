#!/usr/bin/env python3

"""
BGC library: BGCCollection

Handles collections of BGC objects
"""

from ..bgc import BGC
from ..collections.proteins import ProteinCollection
from ..annotations.hmm import hmmDB
from multiprocessing import Pool
from pathlib import Path
import pickle

class BGCCollection:
    """Stores a collection of BGCs and implements collection-wide functions
    """
    
    def __init__(self):
        self.bgcs = {}  # key: BGC.identifier
        self.name = ""

    @classmethod
    def from_file(cls, collection_file:Path):
        """Loads a collection from a file
        
        Args:
            collection_file: (Path) path to a pickled collection (usually with
                extension bgccase)

        Returns:
            A BGCCollection object
        """

        assert collection_file.is_file()

        with open(collection_file, "rb") as dc:
            collection_data = pickle.load(dc)

        assert isinstance(collection_data, cls)

        return collection_data


    def __len__(self):
        return len(self.bgcs)
    

    def add_gbk(self, gbk) -> None:
        """Creates a BGC from a gbk file and adds it to the collection

        Args:
            gbk: (str/Path) path to the gbk file
            identifier: (str, optional) internal name for the BGC. By default,
                the stem of the file will be used
        """
        
        bgc = BGC(gbk)
        if bgc.identifier in self.bgcs:
            print(f"Warning: {bgc.identifier} already in BGC Collection. "
                + "Skipping")
        else:
            self.bgcs[bgc.identifier] = bgc

    
    def predict_domains(
        self, 
        hmmdb: hmmDB, 
        cpus = 1, 
        tc = True, 
        filterdoms = True
    ) -> None:
        """Predict domains in the collections' protein sequences

        Args:
            hmmdb: An hmmDB object pointing to a valid hmm database
            cpus: Number of cpus for hmm detection
            tc: (bool) Use trusted cutoffs (score of the lowest-scoring known\
                true positive)
            filterdoms: (bool) Keep only the highest-scoring domain if over-\
                lapping predictions found   


        Compile the protein sequences of the whole collection and predict 
        domains. Assign predicted domains to each protein and filter
        """
        assert isinstance(hmmdb, hmmDB)
        
        # assume all proteins have identifiers
        pc = ProteinCollection()
        for bgc in self.bgcs.values():
            for protein in bgc.protein_list:
                assert protein.identifier
                pc.proteins[protein.identifier] = protein
            
        pc.predict_domains(hmmdb, cpus, tc, filterdoms)
        
        for bgc in self.bgcs.values():
            bgc.attempted_domain_prediction = True
            bgc.calculate_domain_sets()

        return


    def classify_proteins(self, cpus=1) -> None:
        """Use domains to assign type and role for all proteins in collection
        """
        
        with Pool(cpus) as pool:
            for bgc in self.bgcs.values():
                pool.apply_async(bgc.classify_proteins())
            pool.close()
            pool.join()

        return


    def clear_domain_predictions(self, cpus=1):
        """Remove all domain annotations from the proteins in this collection

        This also removes each BGC-level domain information (domain_set, 
            domain_set_core and domain_set_complement)
        """

        with Pool(cpus) as pool:
            for bgc in self.bgcs.values():
                pool.apply_async(bgc.clear_domain_predictions())
            pool.close()
            pool.join()

        return
    

    def clear_protein_roles(self):
        """Removes all protein role annotations from the collection.

        Note that this may affect annotations imported from the original
        GenBank files
        """

        for bgc in self.bgcs.values():
            for protein in bgc.protein_list:
                protein.role = "unknown"
