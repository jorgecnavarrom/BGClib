#!/usr/bin/env python3

"""
BGC library: utilities

Auxiliary functions and constants
"""

from colorsys import hsv_to_rgb, rgb_to_hsv
from random import uniform

# Don't include Acyl_transf_1 domain as some NRPS have it
PKS_domains = {
    "SAT", "ketoacyl-synt", "Ketoacyl-synt_C", "KAsynt_C_assoc", "TIGR04532"
}
reducing_domains = {"PKS_ER", "KR", "PS-DH"}
PKS3_domains = {"Chal_sti_synt_N", "Chal_sti_synt_C"}
NRPS_domains = {"Condensation", "AMP-binding", "AMP-binding_C"}
NRPS_Independent_Siderophore_domains = {"IucA_IucC"}
Terpene_meroterpenoid_domains = {"mero_tc"}
Terpene_diterpene_domains = {"diterpene_tc"}
Terpene_triterpene_domains = {"SQHop_cyclase_N", "SQHop_cyclase_C"}
Terpene_sesquiterpene_domains = {"TRI5", "Terpene_syn_C_2"}
Terpene_sesquiterpene_bifunc_domains = {"Terpene_syn_C_2", "polyprenyl_synt"}
Terpene_squalene_domains = {"SQS_PSY"}
Terpene_carotenoid_domains = {"TIGR03462"}
Terpene_UbiA_domains = {"UbiA"}
Other_terpene_domains = {"Terpene_synth", "Terpene_synth_C", "Lycopene_cycl",
    "Prenyltrans"}
Terpene_domains =  Terpene_meroterpenoid_domains | Terpene_diterpene_domains \
    | Terpene_triterpene_domains | Terpene_sesquiterpene_domains \
    | Terpene_sesquiterpene_bifunc_domains | Terpene_squalene_domains \
    | Terpene_UbiA_domains | Other_terpene_domains
                   
DMATS_domain = {"Trp_DMAT"}

# Precursors
FAS2 = {"Fas_alpha_ACP" ,"FAS_I_H", "ACPS"}
FAS1 = {"DUF1729", "FAS_meander", "MaoC_dehydrat_N", "MaoC_dehydratas"}
precursor_domains = FAS1 | FAS2

# TODO: add PT domain in this set until we have a better model?
hmmdbs_without_tc = {"FNP_terpene_models"}

# TODO: choose colors for the last ones
valid_CBP_types_fungal = {
    "nrPKS": "#76b7f4", # blue
    "rPKS": "#2c9cdc", # darker blue
    "t3PKS": "#3368c8", #007dfb", # a bit more dark blue
    "NRPS": "#ffc755", # orange
    "NRPS-like": "#ffeb87", # light orange / yellow
    "other_PKS": "#80d8eb", # another blue. Intense-ish 
    "unknown_PKS": "#aaffff", # lighter version of previous
    "PKS-NRPS_hybrid": "#aa007f", # purple
    "PKS-mmNRPS_hybrid": "#9a116f", # darkish purple
    "NRPS-PKS_hybrid": "#a25fe6", # violet
    "NIS": "#c20c28", # red blood
    "Meroterpenoid_synthase": "#b42d2d", # shades of brown/dark orange
    "Diterpene_synthase": "#ba5454", # shades of brown/dark orange
    "Triterpene_synthase": "#e04d37", # shades of brown/dark orange
    "Sesquiterpene_synthase": "#995358", # shades of brown/dark orange
    "Sesquiterpene_bifunctional_synthase": "#f4f4fc", # shades of brown/dark orange
    "Carotenoid_synthase": "#f4815c", # shades of brown/dark orange
    "Squalene_synthase": "#be480a", # shades of brown/dark orange
    "UbiA-type_terpene": "#da746e", # shades of brown/dark orange
    "Terpene_other": "#ca8a73", # shades of brown/dark orange
    "DMATS": "#c5c55d" # "gold"
}

# antiSMASH list retrieved 2020-08-19
valid_CBP_types_antiSMASH = {
    'T1PKS': '#34b4eb',
    'T2PKS': '#999999',
    'T3PKS': '#3368c8',
    'transAT-PKS': '#999999',
    # "Trans-AT PKS fragment, with trans-AT domain not found"
    'transAT-PKS-like': '#999999', 
    'PpyS-KS': '#999999',
    'hglE-KS': '#999999',
    'CDPS': '#cccdac', # "tRNA-dependent cyclodipeptide synthases"
    'PKS-like': '#999999',
    'arylpolyene': '#cdc5c0',
    'resorcinol': '#999999',
    'ladderane': '#cdc69f',
    'PUFA': '#999999',
    'nrps': '#999999',
    'nrps-like': '#999999',
    'thioamide-NRP': '#999999',
    'terpene': '#b44c3a',
    'lanthipeptide': '#5b9950', # Obsolete: split into subclasses in antiSMASH6
    'lipolanthine': '#999999',
    'bacteriocin': '#70b598',
    'betalactone': '#b8cdc5',
    'thiopeptide': '#90b442',
    'linaridin': '#999999',
    'cyanobactin': '#999999',
    'glycocin': '#999999',
    'LAP': '#999999',
    'lassopeptide': '#999999',
    'sactipeptide': '#999999',
    'bottromycin': '#999999',
    'head_to_tail': '#999999',
    'microviridin': '#31b568',
    'proteusin': '#999999',
    'blactam': '#999999',
    'amglyccycl': '#999999',
    'aminocoumarin': '#999999',
    'siderophore': '#c20c28',
    'ectoine': '#999999',
    'butyrolactone': '#b1b0bc',
    'indole': '#999999',
    'nucleoside': '#999999',
    'phosphoglycolipid': '#999999',
    'melanin': '#999999',
    'oligosaccharide': '#999999',
    'furan': '#999999',
    'hserlactone': '#bca8ad',
    'phenazine': '#999999',
    'phosphonate': '#f6f9d4',
    'fused': '#999999',
    'PBDE': '#999999',
    'acyl_amino_acids': '#999999',
    'tropodithietic-acid': '#999999',
    'NAGGN': '#999999',
    'RaS-RiPP': '#999999',
    'fungal-RiPP': '#3cb5a1',
    'TfuA-related': '#999999',
    'other': '#999999',
    'saccharide': '#999999',
    'fatty_acid': '#999999',
    'halogenated': '#999999'
}
valid_CBP_types_antiSMASH_set = set(valid_CBP_types_antiSMASH.keys())


valid_CBP_types_colors = valid_CBP_types_fungal
valid_CBP_types_colors.update(valid_CBP_types_antiSMASH)

valid_CBP_types = valid_CBP_types_fungal.keys() | valid_CBP_types_antiSMASH
# other colors
# "unknown": "#f4f4fc", # super light lilac
# "other": "#fcdcdc", # very light pink
# no_domains": "#ffffff", # white. For stuff like RiPPs (TODO)

role_colors = {
    "biosynthetic":"#f06c6e", # red, rgb(240, 108, 110)
    "tailoring":"#8fc889", # green, rgb(143, 200, 137)
    "transport":"#f0d963", # yellow, rgb(240, 217, 99)
    "regulatory":"#33c1f0", # blue, rgb(51, 193, 240)
    "other":"#eff0f1", # light gray, rgb(239, 240, 241)
    "precursor":"#9797dc", # lilac, rgb(151,151,220)
    "unknown":"#eff0f1", # gray, rgb(220, 220, 220) dcdcdc
    "resistance":"#f0a1ac", # light red, rgb(240, 161, 172) 
    "biosynthetic-additional":"#f0986b" # orange rgb(240, 152, 107)
} 


def random_color_tuple(
    h_:tuple[float, float], 
    s_:tuple[float, float], 
    v_:tuple[float, float]
) -> str:
    """Generates a random hex color given hsv bounds
    
    Args:
        h_: hue. A tuple containing lower and upper bounds
        s_: saturation
        v_: value

    Returns:
        A random hex color

    Each argument is a tuple of lower and upper bounds between 0.0 and 1.0. 
    For hue, these bounds map to a complete round in the cylindrical model 
    (e.g. 0.0 corresponds to red, 0.33 maps to 120° or green, etc.)
    
    Additional info:
    https://en.wikipedia.org/wiki/HSL_and_HSV
    https://stackoverflow.com/a/1586291       
    """
    
    h = uniform(*h_)
    s = uniform(*s_)
    v = uniform(*v_)

    r, g, b = map(int, [c*255 for c in hsv_to_rgb(h, s, v)])

    return f"#{int(r):02x}{int(g):02x}{int(b):02x}"

