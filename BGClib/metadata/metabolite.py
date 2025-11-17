#!/usr/bin/env python3

"""
BGC library: metabolite

Metaclass to handle metabolite data (apply to BGC)
"""

class Metabolite:
    def __init__(self):
        self.name = ""                  # e.g "1,3,6,8-tetrahydroxynaphthalene"
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
