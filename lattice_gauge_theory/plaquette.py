import numpy as np

class Plaquette(object):
    """ Object representing elementary plaquette in the the lattice. """
    plaquette_number = 0 # counter so that every plaquette is distinguishable
    def __init__(self, start_site):
        """
        Initialize plaquette instance.

        Args:
            start_site (Site): Lattice site in lower-left hand corner of
            plaqeutte.

        Returns:
            None
        """
        Plaquette.plaquette_number += 1
        self.number = Plaquette.plaquette_number
        self._site = start_site
        ######################################
        #   TODO:   FINISH IMPLEMENTATION    #
        ######################################



