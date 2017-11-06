import numpy as np
import random
import itertools
import sys
from collections import Counter


class Lattice(object):
    """ Lattice class. """

    def __init__(self, sites, cell_lengths):
        """
        Initialize a Lattice instance.

        Args:
            sites (List(Site)): List of sites contained in the lattice.
            cell_lengths (np.array(x,y,z)): Vector of cell lengths for
                fundamental cell.

        Returns:
            None
        """
        self.cell_lengths = cell_lengths
        self.sites = sites
        #  self.links = sites.links
        self.number_of_sites = len(self.sites)
        #  self.number_of_links = len(self.links)
        self.site_labels = set([site.label for site in self.sites])
        #  self.site_populations = Counter([site.label for site in self.sites])
        self.enforce_periodic_boundary_conditions()
        self.initialize_site_lookup_table()
        self.nn_energy = False
        #  self.cn_energies = False
        #  self.site_energies = False
        self.jump_lookup_table = False
        for site in self.sites:
            site.p_neighbors = [self.site_with_id(i) for i in site.neighbors]
        self.reset()

    def enforce_periodic_boundary_conditions(self):
        """
        Enforce periodic boundary conditions in each dimension.

        Args:
            None

        Returns:
            None
        """
        for s in self.sites:
            for i in range(3):
                if s.r[i] < 0.0:
                    s.r[i] += self.cell_lengths[i]
                if s.r[i] > self.cell_lengths[i]:
                    s.r[i] -= self.cell_lengths[i]

    def reset(self):
        """
        Reset all time-dependent counters for the lattice and its sites.

        Args:
            None

        Returns:
            None
        """
        self.time = 0.0
        for site in self.sites:
            site.time_occupied = 0.0

    def initialize_site_lookup_table(self):
        """
        Create a lookup table allowing sites in the lattice to be queried using
        'self.site_lookup[n]', where 'n' is the identifying site number.

        Args:
            None

        Returns:
            None
        """
        self.site_lookup = {}
        for site in self.sites:
            self.site_lookup[site.number] = site

    def site_with_id(self, number):
        """
        Select the site with a specific id number.

        Args:
            number (int): The identifying number for a specific site.

        Returns:
            (Site): The site with id number equal to 'number'
        """
        return self.site_lookup[number]

    def update(self, step):
        """
        Update the lattice by accepting a specific step.

        Args:
            step (Step): The step that has been accepted.

        Returns:
            None.
        """
        ######################################################
        # TODO: Implement step to update link variables U_ij #
        ######################################################

    def step(self):
        """
        Randomly select a potential step, then update the lattice state.

        Args:
            None

        Returns:
            None
        """




