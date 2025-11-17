"""
BGC library

Classes for data storage and analysis of Biosynthetic Gene Clusters
"""

__author__ = "Jorge Navarro"
__version__ = "0.9.0"
__maintainer__ = "Jorge Navarro"
__email__ = "jorge.c.navarro.munoz@gmail.com"

# Import and expose all public classes
from .annotations.hmm import hmmDB
from .annotations.domain import BGCDomain
from .metadata.organism import Organism
from .metadata.metabolite import Metabolite
from .visualization.arrower import ArrowerOpts
from .locus import BGCLocus
from .protein import BGCProtein
from .bgc import BGC
from .collections.bgcs import BGCCollection
from .collections.proteins import ProteinCollection

# Expose data utilities if needed by users
from .data import utilities

__all__ = [
    'hmmDB',
    'BGCDomain',
    'Organism',
    'Metabolite',
    'ArrowerOpts',
    'BGCLocus',
    'BGCProtein',
    'BGC',
    'BGCCollection',
    'ProteinCollection',
    'utilities',
]