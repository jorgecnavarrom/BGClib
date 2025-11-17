#!/usr/bin/env python3

"""
BGC library: organism

Metaclass to handle organism data (apply to BGC)
"""

class Organism:
    def __init__(self):
        self.shortname = ""
        self.fullname = ""
        self.taxid = ""                 # NCBI tax Id
        self.lineage = []
