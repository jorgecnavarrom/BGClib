#!/usr/bin/env python3

"""
BGC library: hmm

Handles setup of HMM databases
"""

from pathlib import Path
from pyhmmer.hmmer import hmmpress
from pyhmmer.plan7 import HMMFile
from ..data.utilities import rgb_to_hsv, hsv_to_rgb, role_colors

class hmmDB:
    """Handles information about hmm databases
    """
    
    def __init__(self):
        self.db_path_list = []      # list of paths to hmm databases
        self.alias = {}             # ID to alias
        self.colors = {}            # ID to tuple(r,g,b)
        self.colors_hex = {}        # ID to hex color
        self.color_outline = {}     # Precalc. darker version of self.colors
        self.colors_outline_hex = {}
        self.cores = 0              # for hmmer. Remember that it always uses
                                    # an extra core for reading the database
                                    
        self.ID_to_AC = {}          # ID: short name; AC: accession number (w/
        self.ID_to_DE = {}          # version); DE: description
            
        self.ID_to_CL = {}          # Can be found in the Pfam-A.hmm.dat file
        self.ID_to_TP = {}          # Family, Domain, Repeats, Motifs, Coiled-
                                    # coils, Disordered. Also in dat file

        self.clan_to_DE = {}        # Can be found in the Pfam-C file
        
        self.ID_to_role = {}        # manually assigned roles (see role_colors)
        self.domain_signature_to_protein_type = {}  # curated using annotated 
                                                    # proteins from MIBiG
                                                    # Domain Signature = tilde-
                                                    # separated domain_IDs
        
        self.read_domain_colors(
            Path(__file__).parent.parent / "data/domain_color_file_IDh.tsv"
        )
        self.read_domain_roles(
            Path(__file__).parent.parent / "data/SM_domain_roles.tsv"
        )
        self.read_protein_types(
            Path(__file__).parent.parent / "data/protein_types.tsv"
        )
        self.read_domain_alias_file(
            Path(__file__).parent.parent / "data/domains_alias.tsv"
        )
        
        return
    

    def add_included_database(self):
        """Reads hmm profiles included in the library
        """
        
        for hmm in (Path(__file__).parent.parent).glob("data/Domain_models/*.hmm"):
            self.add_database(hmm)
        for hmm in (Path(__file__).parent.parent).glob("data/Domain_models/*.HMM"):
            self.add_database(hmm)
            
        return
    
    
    # TODO: read clan description from Pfam-C if found
    def add_database(self, db_path:Path) -> bool:
        """Adds a database of hmm models

        Args:
            db_path: Path to a hmm file
        
        This method also:

        * Verifies that the database is hmmpressed
        * Stores information about each individual domain (linking its ID to 
        Accession and Description). This information can come from three 
        places. In order of priority:
          - A .dat file (that also links ID to clan accession). Its presence is
        automatically detected (Pfam)
          - A .info file, generated from a previous run
          - The hmm database itself. This method is slow as the whole file needs
        to be read
        """
        
        if not db_path.is_file():
            print(f"Not able to add hmm database (not a file. Wrong path?): "
                + f"{db_path}")
            return False
        elif db_path.suffix.lower() != ".hmm":
            print(f"Not able to add hmm database (not a .hmm file): {db_path}")
            return False
        
        # make sure database is already "pressed" for hmmscan
        try:
            assert Path(db_path.parent / (db_path.name + ".h3i")).is_file()
            assert Path(db_path.parent / (db_path.name + ".h3f")).is_file()
            assert Path(db_path.parent / (db_path.name + ".h3m")).is_file()
            assert Path(db_path.parent / (db_path.name + ".h3p")).is_file()
        except AssertionError:
            with HMMFile(db_path) as hmm:
                hmmpress(hmm, db_path)
            
        # Make sure the files were generated
        try:
            assert Path(db_path.parent / (db_path.name + ".h3i")).is_file()
            assert Path(db_path.parent / (db_path.name + ".h3f")).is_file()
            assert Path(db_path.parent / (db_path.name + ".h3m")).is_file()
            assert Path(db_path.parent / (db_path.name + ".h3p")).is_file()
        except AssertionError:
            exit("Not able to hmmpress the database file")
                    
        # save the domain info file in the same place as the database so it's 
        # accessible to the user
        domain_info_file = db_path.parent / (db_path.name + ".domain_info.tsv")
        db_ID_to_AC = {}
        db_ID_to_DE = {}
        db_ID_to_TP = {}
        db_ID_to_CL = {}
        all_ids = set()

        # first try the .dat file
        dat_file = db_path.parent / (db_path.name + ".dat")
        if dat_file.is_file():
            # print(f"\tFound {dat_file.name} file")
            with open(dat_file) as dat:
                # Storing the whole file in mem. Shouldn't be too big
                records = dat.read().split("//")
                for record in records:
                    if not record.strip():
                        continue

                    # example of record from Pfam:
                    # STOCKHOLM 1.0
                    #=GF ID   14-3-3
                    #=GF AC   PF00244.26
                    #=GF DE   14-3-3 protein
                    #=GF GA   33.2; 33.2;
                    #=GF TP   Repeat
                    #=GF ML   223
                    #=GF CL   CL0020
                    
                    # skip the STOCKHOLM line
                    temp_dict = {}
                    for line in record.strip().splitlines()[1:]:
                        # Skip '#=GF ', then split only once
                        code, value = line.strip()[5:].split(None, 1)
                        temp_dict[code] = value

                    assert temp_dict['ID'], "Error, found model without ID!" \
                        + db_path

                    id = temp_dict['ID']
                    all_ids.add(id)
                    db_ID_to_AC[id] = temp_dict.get('AC', "")
                    db_ID_to_DE[id] = temp_dict.get('DE', '')
                    db_ID_to_TP[id] = temp_dict.get('TP', '')
                    db_ID_to_CL[id] = temp_dict.get('CL', '')

        elif domain_info_file.is_file():
            # print(f"\tFound {domain_info_file.name} file")
            with open(domain_info_file) as dif:
                for line in dif:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    line = line.split('\t') # 'dirty' line with newline chars.
                    if len(line) == 3:
                        id, ac, de = line
                        de = de.strip()
                        tp = ''
                        cl = ''
                    elif len(line) == 5:
                        id, ac, de, tp, cl = line
                        cl = cl.strip()
                    else:
                        exit(f"Error while reading {domain_info_file}: " 
                             + "There should be either 3 (old) or 5 columns")
                    assert id, f"Error: found line in {domain_info_file} " \
                        + "without ID (first column)"
                    
                    all_ids.add(id)
                    db_ID_to_AC[id] = ac
                    db_ID_to_DE[id] = de
                    db_ID_to_TP[id] = tp
                    db_ID_to_CL[id] = cl
                      
        else:
            # Whole hmm file needs to be read, which can take a few moments.
            # Also, no type (TP) or clan (CL) info is expected
            with HMMFile(db_path) as hmm_db:
                for hmm_model in hmm_db:
                    id = hmm_model.name.decode()
                    ac = hmm_model.accession.decode()
                    de = hmm_model.description.decode()

                    all_ids.add(id)
                    db_ID_to_AC[id] = ac
                    db_ID_to_DE[id] = de

        with open(domain_info_file, "w") as dif:
            dif.write("#ID\tAC\tDE\tTP\tCL\n")
            for id in sorted(all_ids):
                ac = db_ID_to_AC.get(id, '')
                de = db_ID_to_DE.get(id, '')
                tp = db_ID_to_TP.get(id, '')
                cl = db_ID_to_CL.get(id, '')
                dif.write(f"{id}\t{ac}\t{de}\t{tp}\t{cl}\n")
        
        self.ID_to_AC.update(db_ID_to_AC)
        self.ID_to_DE.update(db_ID_to_DE)
        self.ID_to_TP.update(db_ID_to_TP)
        self.ID_to_CL.update(db_ID_to_CL)

        self.db_path_list.append(db_path)
        
        return True


    def read_domain_alias_file(self, alias_file: Path) -> None:
        """Loads internal aliases for some key core biosynthetic domains
        
        Args:
            alias_file: path to the internal file
        """
        
        try:
            with open(alias_file) as f:
                for line in f:
                    if line[0] == "#" or line.strip() == "":
                        continue
                    
                    line = line.strip().split("\t")
                    hmm_ID = line[0]
                    hmm_alias = line[1]
                    
                    self.alias[hmm_ID] = hmm_alias
        except FileNotFoundError:
            print(f"Could not open domain alias file ({alias_file})")

        
    def read_domain_colors(self, colors_file:Path, append=True) -> None:
        """Loads domain color information

        Args:
            colors_file: Path to colors tsv file. It should have three columns:
                hmm_ID, RGB, HEX.
        
            append: toggle to false to reset all colors

        Notes: RGB color should be a string with values for each color (0-255)
            separated by a comma, e.g. 160,77,53.
            HEX should be prefixed by a hash symbol, e.g. #a04d35. Case should
            not matter
        """
        if not append:
            self.colors.clear()
            self.color_outline.clear()
        
        try:
            with open(colors_file) as f:
                for line in f:
                    line = line.strip()
                    if line[0] == "#" or line == "":
                        continue
                    
                    line = line.split("\t")
                    hmm_ID = line[0]
                    color_rgb = line[1]
                    
                    r, g, b = map(int, color_rgb.split(","))
                    self.colors[hmm_ID] = r, g, b

                    try:
                        color_hex = line[2]
                    except IndexError:
                        color_hex = f"#{r:02x}{g:02x}{b:02x}"
                    self.colors_hex[hmm_ID] = color_hex
                    
                    # Generate darker outline by going to HSV space
                    h, s, v = rgb_to_hsv(r/255.0, g/255.0, b/255.0)
                    darkRGB = [int(c * 255) for c in hsv_to_rgb(h, s, 0.8*v)]

                    self.color_outline[hmm_ID] = tuple(map(str, darkRGB))
                    r, g, b = darkRGB
                    self.colors_outline_hex[hmm_ID] = f"#{r:02x}{g:02x}{b:02x}"

        except FileNotFoundError:
            print(f"Could not open domain colors file ({colors_file})")
        
        return


    def read_domain_roles(self, domain_roles_file:Path, append=True) -> None:
        """Reads a file linking domains with biosynthetic pathway role

        Args:
            domain_roles_file: A tsv file with two columns: hmm_ID and role
                Where role is one of biosynthetic, tailoring, transport, 
                regulation, other, unknown
            
            append: If False, current information stored in self.the ID_to_role
                will be erased first. Otherwise values will only be overwritten
                for previously assigned IDs
        """
        if not append:
            self.ID_to_role.clear()
        
        try:
            with open(domain_roles_file) as f:
                for line in f:
                    line = line.strip()
                    if line[0] == "#" or line == "":
                        continue
                    
                    # a potential third column may contain a comment
                    line = line.split("\t")
                    ID = line[0]
                    role = line[1]
                    
                    if role not in role_colors:
                        print("Warning. Found unknown role when loading "
                            + f"domain role file:\n{ID}\t{role}")
                    
                    self.ID_to_role[ID] = role
        except FileNotFoundError:
            print(f"Could not open domain role file({domain_roles_file})")
            
        return
    
    
    # TODO is it possible to use and integrate interpro rules?
    def read_protein_types(self, protein_types_file:Path, append=True) -> None:
        """
        input:
            protein_types_file: A tsv file with two columns:
                Domain Signature (tilde separated hmm_IDs) \t protein_type
            
            protein_type comes from Fungal MIBiG BGCs with annotations
        append:
            If False, the domain_signature_to_protein_type dictionary will be 
            erased. Otherwise values will only be overwritten for previously
            assigned values
        """
        if not append:
            self.ID_to_role.clear()
        
        try:
            with open(protein_types_file) as f:
                for line in f:
                    line = line.strip()
                    if line[0] == "#" or line == "":
                        continue
                    
                    stripped_line = line.split("\t")
                    domain_signature = stripped_line[0]
                    protein_type = stripped_line[1]
                    self.domain_signature_to_protein_type[domain_signature] = \
                        protein_type
        except FileNotFoundError:
            print(f"Could not open protein types file({protein_types_file})")
            
        return
