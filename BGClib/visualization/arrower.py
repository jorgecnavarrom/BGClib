#!/usr/bin/env python3

"""
BGC library: arrower

Handles SVG visualization of genetic loci elements
"""

from ..data.utilities import role_colors, random_color_tuple
from pathlib import Path

class ArrowerOpts:
    """Options for Arrower-like figures. 
    
    Colors are linked with domains so see the HMM_DB class for them.
    """
    
    def __init__(self, cfg=None):
        self.shape = 'Arrow'    # Options: 'Arrow' and 'Ribbon'
        self.scaling = 20               # px per bp
        self.topbottom_margin = 2      # upper/bottom space between fig and 
                                        # border

        self.arrow_height = 40          # H: Height of arrows' body
        #self.arrow_head_height = 15    # h: Additional height of arrows' head
                                        #  (i.e. total arrow head = H + 2*aH)
                                        # internally, h = H/2 to preserve 45° 
                                        # angle
                                        
        #self.arrow_head_length = 30    # Distance from start of head till tip.
        
        self.gene_contour_thickness = 2 # note: thickness grows outwards
        
        self.internal_domain_margin = 3
        self.domain_contour_thickness = 1
        
        self.stripe_thickness = 3       # genetic track line
                
        self._color_mode = "white"
        self._valid_color_modes = {
            "white", "gray", "random-pastel", "random-dark", "random", "roles",
            "domains", "none"
        }
        
        self.outline = True
        self.show_domains = True
        self.show_introns = True # Show a gap where introns are

        self.original_orientation = True # If false, all arrows will point forward
        
        if cfg is not None:
            self.load_options(Path(cfg))

        return

    @property
    def color_mode(self):
        return self._color_mode
    @color_mode.setter
    def color_mode(self, cm):
        if cm in self._valid_color_modes:
            self._color_mode = cm
        else:
            print("Color mode not supported; defaulting to 'white'")
            color_modes_txt = "', '".join(self._valid_color_modes)
            print(f"Valid color modes are ['{color_modes_txt}']")
            self._color_mode = "white"


    def load_options(self, cfg:Path) -> None:
        """Loads options for the figures from a config file
        
        Args:
            cfg: a path object to the cfg file

        The file itself consists on option/value pairs (separated by an equal
        symbol). A file with the defaults should be included in the project

        The config file can contain a limited number of options
        """

        try:
            assert cfg.is_file()
        except AssertionError:
            raise FileNotFoundError("Arrower options file: not a valid file")
        
        # collect data
        temp_options = {}
        for line in open(cfg, "r"):
            if line[0] == "#" or line.strip() == "":
                continue
            
            try:
                option, value = line.strip().split("=")
            except ValueError:
                print(f"ArrowerOpts configuration file: ignoring bad line "
                    + f"({line.strip()})")
                continue
            else:
                option = option.strip().lower()
                value = value.strip()

                temp_options[option] = value

        # Check text-based values
        if 'shape' in temp_options:
            shape_value = temp_options['shape']
            try:
                assert shape_value in {'Arrow', 'Ribbon'}
            except AssertionError:
                raise ValueError("Arrower config file: Value for option "
                    + "'shape' should be 'Arrow' or 'Ribbon'")
            else:
                self.shape = shape_value

        if 'color_mode' in temp_options:
            color_mode_value = temp_options['color_mode'].lower()
            if color_mode_value == "grey":
                color_mode_value = "gray"
            try:
                assert color_mode_value in self._valid_color_modes
            except AssertionError:
                vcm_text = ", ".join(self._valid_color_modes)
                raise ValueError("Arrower config file: Value for option "
                    + f"'color_mode' should be one of: {vcm_text}")
            else:
                self._color_mode = color_mode_value

        # Check integer values
        int_opts = {
            'topbottom_margin', 'scaling', 'arrow_height', 
            'gene_contour_thickness', 'internal_domain_margin', 
            'domain_contour_thickness', 'stripe_thickness'
        }
        for option in int_opts:
            if option in temp_options:
                try:
                    value = int(temp_options[option])
                except ValueError as e:
                    e.add_note("Arrower config file: value for option "
                        + f"'{option}' could not be converted to an integer")
                    raise

                try:
                    assert value > 0
                except AssertionError:
                    e.add_note("Arrower config file: value for option "
                        + f"'{option}' should be positive!")
                    raise
                
                setattr(self, option, value)

        # Check boolean values
        bool_opts = {
            'outline', 'show_domains', 'show_introns', 'original_orientation'
        }
        for option in bool_opts:
            if option in temp_options:
                value = temp_options[option].capitalize()
                try:
                    assert value in {'True', 'False'}
                except AssertionError as e:
                    e.add_note("Arrower config file: value for option "
                        + f"'{option}' should be either True or False")
                    raise
                else:
                    setattr(self, option, value == 'True')

        # Finally, warn if weird, unsupported or typo'd options were used
        all_opts = {
            'shape', 'color_mode'
        }
        all_opts.update(int_opts)
        all_opts.update(bool_opts)
        for opt in temp_options:
            if opt not in all_opts:
                print(f"Arrower config file: Warning, unknown option '{opt}'")
                
        return
    
