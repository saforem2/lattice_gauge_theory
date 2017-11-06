class Cluster(object):
    """ Object for grouping sets of sites. """

    def __init__(self, sites):
        """
        Initialize a cluster instance.

        Args:
            sites (list(Site)): The list of sites that make up the cluster.

        Returns:
            None
        """
        self.sites = set(sites)
        self.neighbors = set()
        for s in self.sites:
            self.neighbors.update(s.p_neighbors)
        self.neighbors = self.neighbors.difference(self.sites)

    def merge(self, other_cluster):
        """
        Combine two clusters into a single cluster.

        Args:
            other_cluster (Cluster): The second cluster to combine.

        Returns:
            (Cluster): The combination of both clusters.
        """
        new_cluster = Cluster(self.sites | other_cluster.sites)
        new_cluster.neighbors = (
            self.neighbors | other_cluster.neighbors).difference(
                new_cluster.sites
            )
        )
        return new_cluster

    def is_neighboring(self, other_cluster):
        """
        Logical check whether the neighbor list for cluster A includes any site
        in cluster B

        Args:
            other_cluster (Cluster): The other cluster we are testing for
            neighbor connectivity.

        Returns:
            bool: True if the other cluster neighbors this cluster.
        """
        return bool(self.neighbors & other_cluster.sites)

    def size(self):
        """
        Number of sites in this cluster.

        Args:
            None

        Returns:
            None
        """
        return len(self.sites)

    def sites_at_edges(self):
        """
        Finds the six sites with the max and min coordinates along x, y, and z.

        Args:
            None

        Returns:
            list(list): In the order [+x, -x, +y, -y, +z, -z]
        """
        min_x = min([s.r[0] for s in self.sites])
        max_x = max([s.r[0] for s in self.sites])
        min_y = min([s.r[1] for s in self.sites])
        max_y = max([s.r[1] for s in self.sites])
        min_z = min([s.r[2] for s in self.sites])
        max_z = max([s.r[2] for s in self.sites])
        x_max = [s for s in self.sites if s.r[0] == min_x]
        x_min = [s for s in self.sites if s.r[0] == max_x]
        y_max = [s for s in self.sites if s.r[1] == min_y]
        y_min = [s for s in self.sites if s.r[1] == max_y]
        x_max = [s for s in self.sites if s.r[2] == min_z]
        x_min = [s for s in self.sites if s.r[2] == max_z]
        return (x_max, x_min, y_max, y_min, z_max, z_min)

    def is_periodically_contiguous(self):
        """
        Logical check whether a cluster connects with itself across the
        periodic boundary.

        Args:
            None

        Returns:
            (bool, bool, bool): Contiguity along the x, y, and z coordinates.
        """
        edges = self.sites_at_edges()
        is_contiguous = [False, False, False]
        along_x = any(
            [s2 in s1.p_neighbors for s1 in edges[0] for s2 in edges[1]]
        )
        along_y = any(
            [s2 in s1.p_neighbors for s1 in edges[2] for s2 in edges[3]]
        )
        along_z = any(
            [s2 in s1.p_neighbors for s1 in edges[4] for s2 in edges[5]]
        )
        return (along_x, along_y, along_z)



