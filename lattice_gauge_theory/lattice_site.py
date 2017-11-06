from collections import Counter

class Site(object):
    """ Site class. """

    index = 0
    def __init__(self, number, coordinates, neighbors, energy, label):
        """
        Initialize a lattice Site object.

        Args:
            number (int): An identifying number for this site.
            coordinates (np.array(x,y,z)): The coordinates of this site.
            neighbors (list(int)): A list of the id numbers of neighboring
                sites.
            energy (float): On-site occupation energy.
                (Equal to zero for Abelian gauge theory on links of dual
                lattice)
            label (str): Label for classifying this as a specific site type.

        Returns:
            None.

        Notes:
            There should be a strict 1:1 mapping between sites and site
            numbers.
        """
        self.number = number
        self.index = Site.index
        Site.index += 1
        self.r = coordinates
        self.neighbors = neighbors
        self.p_neighbors = None    # ptr to neighboring sites. init in Lattice
        self.energy = energy
        self.occupation = 0
        #  self.atom = None
        #  self.is_occupied = False
        self.label = label
        #  self.time_occupied = 0.0


