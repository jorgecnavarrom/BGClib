#!/usr/bin/env python3

"""
BGC library: locus

A BGC attempts to describe the complete pathway towards metabolite synthesis.
This implies that the implicated genes may not be all co-localized in a single
genomic locus. This class attempts to help handling these cases
"""

class BGCLocus:
    """Organizes data that is split in different records in the same file

    BGC might be divided in more than one locus (because of sequencing, or, as
    is more often the case in fungi, because of genomics)
    
    This class can be used to organize proteins
    """
    def __init__(self):
        self.identifier = ""        # e.g. "bgc~L01"
        self.protein_list = []
        self.gene_coordinates = []  # a list of tuples. end-start might not 
                                    # match the corresponding protein length 
                                    # due to splicing, but will be useful for 
                                    # inter-protein region length when making 
                                    # arrower figure.
        self.length = 0

