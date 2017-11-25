import numpy as np
from lattice_gauge_theory.utils import *


class Lattice(object):
    """
    Lattice class.

    Attributes:
        shape (array-like): Array containing linear extent in each dimension.
        sites (np.arr):

    Methods:
    """

    def __init__(self, shape, link_dtype=int, site_dtype=int,
                 care='links', p_action=True):
        """
        Initializes lattice.

        Args:
            shape (array-like):
                Array containing linear extent in each dimension.
            link_dtype (data-type):
                Data type of information stored on the links of the lattice.
                (DEFAULT: int)
            site_dtype (data-type):
                Data type of information stored on the sites of the lattice.
                (DEFAULT: int)
            p_action (bool):
                Flag specifying if we wish to consider the plaquette (Wilson)
                action.
        """
        self._link_dtype = link_dtype
        self._site_dtype = site_dtype
        self.shape = shape
        self.num_dims = len(shape)
        self.num_sites = np.product(shape)
        self.num_links = self.num_sites * self.num_dims
        if site_dtype is not None:
            self.__dict__['sites'] = np.zeros(shape, dtype=site_dtype)
            self.sites_flat = np.ravel(self.sites)  # flattened for iterating
        if link_dtype is not None:
            self.__dict__['links'] = np.zeros(shape + (len(shape),),
                                              dtype=link_dtype)
            self.links_flat = np.ravel(self.links)

        if p_action:
            # Generate arrays for quickly accessing plaquettes.
            # Potentially memory restrictive.

            # p_sites is an array of link indices for all plaquettes,
            # ordered by sites and direction
            self.p_sites = np.zeros(shape + (self.num_dims,
                                             self.num_dims, 4),
                                    dtype=object)
            for i in multirange(shape):
                for j, k in multirange((self.num_dims, self.num_dims)):
                    if j == k:
                        continue
                    #  import pdb
                    #  pdb.set_trace()
                    self.p_sites[i + (j,k)] = sp_links(i,j,k)
                    for l in range(4):
                        self.p_sites[i+(j,k,l)] = tuple(
                            np.mod(self.p_sites[i+(j,k,l)],
                                   shape+(self.num_dims,))
                        )

            # p_links is an array of link indices for all plaquettes,
            # ordered by link, direction, and sign.
            # So indexing takes the form:
            #   p_links[link, direction, # sign, side of plaquette]
            self.p_links = np.zeros(shape + (self.num_dims,
                                             self.num_dims, 2, 4),
                                            dtype=object)
            for i in multirange(shape):
                for dE, d in multirange((self.num_dims, self.num_dims)):
                    if dE == d:
                        continue
                    for s in {-1, 1}:
                        self.p_links[i+(dE,d,(1-s)//2)] = lp_links(i+(dE,),d,s)
                        for j in range(4):
                            self.p_links[i+(dE,d,(1-s)//2,j)] = tuple(
                                np.mod(self.p_links[i+(dE,d,(1-s)//2,j)],
                                       shape+(self.num_dims,))
                            )
        if care=='gauge':
            self.care = self.links
            self.care_flat = self.links_flat
        elif care=='sites':
            self.care = self.sites
            self.care_flat = self.sites_flat

    def __getitem__(self, idx):
        if not hasattr(idx, '__len__'):
            return self.care_flat[idx]
        elif len(idx) == self.num_dims + 1:
            return self.links[idx]
        elif len(idx) == self.num_dims:
            return self.sites[idx]

    def __setattr__(self, name, val):
        if name == 'links':
            if not val.dtype == self._link_dtype:
                raise ValueError("Value must have type {}".format(
                    self._link_dtype
                ))
            if not val.shape == self.links.shape:
                raise ValueError("Value must have shape {}".format(
                    self.links.shape
                ))
            self.__dict__['links'] = val
            self.__dict__['links_flat'] = np.ravel(self.links)
        elif name == 'sites':
            if not val.dtype == self._site_dtype:
                raise ValueError("Value must have type {}".format(
                    self._site_dtype
                ))
            if not val.shape == self.sites.shape:
                raise ValueError("Value must have shape {}".format(
                    self.sites.shape
                ))
            self.__dict__['sites'] = val
            self.__dict__['sites_flat'] = np.ravel(selv.sites)
        else:
            self.__dict__[name] = val
