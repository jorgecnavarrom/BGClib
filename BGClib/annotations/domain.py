#!/usr/bin/env python3

"""
BGC library: domain

Handles domain information attached to a protein sequence
"""

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..protein import BGCProtein

class BGCDomain:
    def __init__(self, protein, ID, AC, DE, alias, ali_from, ali_to, \
            hmm_from, hmm_to, hmm_size, score, Evalue, algn_seq):
        self.protein = protein
        self.ID = ID                # domain short name e.g. 'ketoacyl-synt'
        self.AC = AC                # domain accession. e.g. 'PF00109.27'
        self.DE = DE                # domain description. e.g. 
                                   #'Beta-ketoacyl synthase, N-terminal domain'
        self.alias = alias          # domain alias e.g. 'KS'
        self.ali_from = ali_from    # Position in target sequence at which the 
        self.ali_to = ali_to        #   hit starts/ends
        self.hmm_from = hmm_from    # Position in the hmm at which the hit 
        self.hmm_to = hmm_to        #   starts/ends
        self.hmm_size = hmm_size    # Size of the hmm model
        self.score = score          # Only depends on profile HMM & target seq.
        self.Evalue = Evalue        # Depends on score and DB size
        self.algn_seq = algn_seq    # Sequence of hit, aligned to hmm profile
        
    def get_sequence(self):
        return self.protein.sequence[self.ali_from:self.ali_to]

    def get_aligned_sequence(self):
        dashfrom = "-"*self.hmm_from
        dashto = "-"*(self.hmm_size-self.hmm_to)
        return f"{dashfrom}{self.algn_seq}{dashto}"
