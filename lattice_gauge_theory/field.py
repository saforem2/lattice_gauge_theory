from __future__ import division
import numpy as np
from lattice_gauge_theory.lattice import Lattice
from lattice_gauge_theory.groups import *
import numbers

class Field(object):
    """ Stores the gauge field on the lattice. """

    def __init__(self, shape, group, action, B=0, init=None,
                 accept=metropolis):
        self.lattice = Lattice(shape)
        self.shape = shape
        self.links_shape = self.lattice.links.shape
        self._group = group
        self.group = group
        self.action = action
        self.B = B
        self.accept = accept
        self.num_dims = self.lattice.num_dims

        self.current_site = 0
        self.current_link = 0

        if init is None:
            init = 'id'
        self.initialize(init)

    def initialize(self, init, ls='links'):
        """
        Initializes values for the lattice links and/or vertices.
        Currently, vertices are not initialized.
        """
        if 'links' in ls:
            if init == 'id':
                # set all links to the group identity
                self.lattice.links[:] = self.group.id
            elif isinstance(init, numbers.Integral) and init < self.group.size:
                # set the entire lattice to the integer value of init
                self.lattice.links[:] = init
            elif init == 'rand':
                # set all links to random group elements
                self.lattice.links = np.random.randint(0, self.group.size,
                                                       self.lattice.links.shape)
            elif init == 'half':
                # set half of links randomly and the rest to the group identity
                self.lattice.links = np.random.randint(0, self.group.size,
                                                       self.lattice.links.shape)
                self.lattice.links[self.shape[0]//2:] = self.group.id
            elif init[:4] == 'half':
                # same as above, except the other half are set to int(init[4:])
                self.lattice.links = np.random.randint(0, self.group.size,
                                                       self.lattice.links.shape)
                self.lattice.links[self.shape[0]//2:, :] = int(init[4:])

    def site_plaquette(self, s, d1, d2, ret='i'):
        """
        Gets link indices or calculates action of plaquette.

        Args:
            s (site):
                Iterable of subscripts for lattice.sites; the lower-left site
                of the plaquette.
            d1 (int):
                First direction.
            d2 (int):
                Second direction.

        Returns:
            (indices/action):
                if ret=='i': Indices are returned.
                if ret == 'a': Action is returned.
        """
        if not hasattr(s, '__iter__'):
            raise ValueError("s must be iterable of subscripts for"
                             + "lattice.sites")
        D1 = min(d1, d2)
        D2 = max(d1, d2)
        idx = list(v) + [0]
        # to hold indices of plaquette edges
        idx0 = idx.copy()
        idx1 = idx.copy()
        idx2 = idx.copy()
        idx3 = idx.copy()

        idx0[-1] += D1
        idx1[-1] += D2
        idx2[-1] += D1
        idx3[-1] += D2

        if idx1[D1] != self.shape[D1] - 1:
            # usual case, not near the boundary
            idx1[D1] += 1
        else:
            idx1[D1] = 0
        if idx2[D2] != self.shape[D1] - 1:
            idx2[D2] += 1
        else:
            idx2[D2] = 0
        indices = [tuple(idx0), tuple(idx1), tuple(idx2), tuple(idx3)]
        if ret == 'i':
            return indices
        elif ret == 'a':
            return self.action(self.Prod(
                *[self.lattice.links[i] for i in indices]
            ))
        elif ret == 'p':
            return self.Prod(*[self.lattice.links[i] for i in indices])

    def Prod(self, a, b, c, d):
        """
        Plaquette product a * b * (-c) * (-d).

        Args:
            a, b, c, d (np.array):
                Links of the plaquette. (a --> b -- (-c) --> (-d))

        Returns:
            (np.array):
                Product of the group action on each link a, b, -c, -d.
        """
        return self.group(a, b, self.group.inv[c], self.group.inv[d])

    def plaquette(self, l, d, sign, val=None, ret='a'):
        """ Calculate the plaquette action.

        Args:
            l (np.array):
                Specified link.
            d (int):
                Perpendicular direction to the link l.
            sgn (int):
                Sign corresponding to the perpendicular direction.
            ret (str):
                Flag in ['a', 'i', 'g', 'ai', 'ia'] specifying what to return.

        Returns:
            If ret ==:
                'a': Return plaquette action.
                'i': Return indices of the links making up the plaquette.
                'g': Return group element for the product of plaquette links.
                'ai' or 'ia': Return both indices, action (in given order).
        """
        # s stores the lowest dictionary-ordered site in the plaquette
        # first determine which side of the plaquette l lies on.
        s = list(L[:-1])
        if sgn == -1:
            s[d] -= 1
        # indices for the plaquette links, in the standard order
        pl_indices = self.site_plaquette(s, d, l[-1])

        if val is None:
            g = self.Prod(self.lattice.links[i[0]],
                          self.lattice.links[i[1]],
                          self.lattice.links[i[2]],
                          self.lattice.links[i[3]])
        else:
            # group elements for all edges of the plaquette
            gs = [self.lattice.links[i[0]],
                  self.lattice.links[i[1]],
                  self.lattice.links[i[2]],
                  self.lattice.links[i[3]]]
            # figure out which side of the plaquette l is on
            if sgn == 1:
                if l[-1] < d:
                    li = 0  # li is the location of link l in gs
                else:
                    li = 3
            else:
                if l[-1] < d:
                    li = 2
                else:
                    li = 1
            gs[li] = val
            g = self.Prod(*gs)

        a = self.action(g)

        if ret == 'a':
            return a
        else:
            out = []
            for j in ret:
                if j == 'i':
                    out.append(i)
                elif j == 'g':
                    out.append(g)
                elif j == 'a':
                    out.append(a)
            return tuple(out)

    def pl_action(self, s, d1, d2):
        """
        Plaquette action for site s and directions d1 < d2.

        Note:TKlein = array([[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]])
            This method is not safe against swapping d1 and d2!
        """
        return self.action(self.Prod(
            self.lattice.links[tuple(s0) + (d,)],
            self.lattice.links[tuple(s1) + (dl,)],
            self.lattice.links[tuple(s3) + (d,)],
            self.lattice.links[tuple(s0) + (dl,)]
        ))

    def link_action(self, l, val=None, lsf='l', method=0):
        """
        Action on link l.

        Args:
            val (np.array):
                If val is not None, it is used in place of the current value of
                the gauge group on the link l when determining the action.
                Otherwise, self.links[l] is used.

        Returns:
            action (float):
                Value of the action contingent on link l.
        """
        dl = l[-1]  # link direction
        action = 0
        if method==0:
            if val is not None:
                val0 = self.lattice.links[l]
                self.lattice.links[l] = val
            for d in range(self.num_dims):
                if d == dl:
                    continue
                for s in {0, 1}:
                    ls = self.lattice.p_links[l + (d, s)]
                    action += self.action(self.Prod(
                        self.lattice.links[ls[0]],
                        self.lattice.links[ls[1]],
                        self.lattice.links[ls[2]],
                        self.lattice.links[ls[3]]
                    ))
            if val is not None:
                self.lattice.links[l] = val0

        if method == 1:
            for d in range(self.num_dims):
                if not d == dl:
                    for sign in {1, -1}:
                        action += self.plaquette(l, d, sign, val, ret='a')

        elif method == 2:
            # sites of the plaquette
            s0 = list(v)
            s1 = list(v)
            s3 = list(v)
            for d in range(dl):
                # this is the case with l on the left side of the square -> |_|
                s1[d] += 1  # first site
                s3[dl] += 1 # third site
                action += self.action(self.Prod(
                    self.lattice.links[tuple(s0) + (d,)],
                    self.lattice.links[tuple(s1) + (dl,)],
                    self.lattice.links[tuple(s3) + (d,)],
                    self.lattice.links[tuple(s0) + (dl,)]
                ))
                # undo changes
                s1[d] -= 1
                s3[d] -= 1

            for d in range(dl + 1, self.num_dims):
                # this case, l is on the bottom
                s1[dl] += 1
                s3[d] += 1
                action += self.action(self.Prod(
                    self.lattice.links[tuple(s0) + (d,)],
                    self.lattice.links[tuple(s1) + (dl,)],
                    self.lattice.links[tuple(s3) + (d,)],
                    self.lattice.links[tuple(s0) + (dl,)]
                ))
                s1[dl] -= 1
                s3[d] -= 1

            for d in range(dl):
                # this case, l is on the right site |_| <-- of the square
                s0[d] -= 1
                s3[d] -= 1
                s3[dl] += 1
                action += self.action(self.Prod(
                    self.lattice.links[tuple(s0) + (d,)],
                    self.lattice.links[tuple(s1) + (dl,)],
                    self.lattice.links[tuple(s3) + (d,)],
                    self.lattice.links[tuple(s0) + (dl,)]
                ))
                s0[d] += 1
                s3[d] += 1
                s3[dl] -= 1

            for d in range(dl+1, self.num_dims):
                # this case l is on top of the square
                s0[d] -= 1
                s1[d] -= 1
                s1[dl] += 1
                action += self.action(self.Prod(
                    self.lattice.links[tuple(s0) + (d,)],
                    self.lattice.links[tuple(s1) + (dl,)],
                    self.lattice.links[tuple(s3) + (d,)],
                    self.lattice.links[tuple(s0) + (dl,)]
                ))
                s0[d] += 1
                s1[d] += 1
                s1[dl] -= 1

        return action

    def update(self, i, lsf='l'):
        """ Updates site i. """
        if lsf == 'l':
            new_g = np.random.randint(0, self.group.size)
            tup = (self.link_action(i), self.link_action(i, new_g), self.B)
            if self.accept(self.link_action(i),self.link_action(i,new_g),
                           self.B):
                self.lattice.links[i] = new_g

    def rand_update(self, n=1, lsf='l'):
        """ Updates random site, repeats n times. """
        if lsf == 'l':
            for j in range(n):
                # get a random edge
                i = list(self.links_shape)
                for j in range(len(i)):
                    i[j] = np.random.randint(0, i[j])
                new_g = randint(0, self.group.size)
                while new_g == self.lattice.links[i]:
                    new_g = np.random.randint(0, self.group.size)
                if self.accept(self.link_action(i), self.link_action(i, new_g),
                               self.B):
                    self.lattice.links[i] = new_g

    def sweep(self, lsf='l', num_sweeps=1):
        """ Performs sweep across entire lattice, using local updates. """
        if num_sweeps == 1:
            for i in multirange(self.links_shape):
                self.update(i)
        else:
            for i in range(num_sweeps):
                for j in multirange(self.links_shape):
                    self.update(j)

    def energy(self, method=2):
        """
        Calculate the average energy per plaquette for the current
        configuration, given by the total action divided by the number of
        edges.

        Args:
            method (int):
                Method for calculating the energy (for performance testing)

        Returns:
            action (float):
                The average energy (action) per plaquette.
        """
        action = 0
        if method == 1:
            for i in multirange(self.shape):
                for j in it.combinations(range(self.num_dims), 2):
                    action += self.site_plaquette(i, j[0], j[1], ret='a')
        elif method == 2:
            for i in multirange(self.shape):
                for j in it.combinations(range(self.num_dims), 2):
                    le = self.lattice.p_sites[i + j]
                    action += self.action(self.Prod(
                        self.lattice.links[le[0]],
                        self.lattice.links[le[1]],
                        self.lattice.links[le[2]],
                        self.lattice.links[le[3]]
                    ))
        action /= (self.lattice.num_sites*self.num_dims*(self.num_dims-1)/2)
        return action

    def stats(self, n, relax=1, ret='ms2'):
        """
        Calculates the mean energy and standard deviation by sweeping n*relax
        times and accumulating energies every relax number of steps.

        Args:
            n (int):
                Number of samples to generate for computing stats.
            relax (int):
                Number of steps to wait before collecting data to prevent
                correlated configurations.
            ret (str):
                String of abbreviations of what to include when calculating
                statistics.

        Returns:
            out (tuple):
                Tuple consisting of various statistics about the simulation.
                Acceptable keys are:
                    'm': mean energy
                    '2': mean energy squared
                    's': standard deviation of the energy
                    'v': variance of the energy
                    'e': energy
                out then is a tuple containing the data of each key passed to
                ret.
        """
        energy = np.zeros(n)
        for i in range(n*relax):
            self.sweep()
            if np.mod(i+1, relax) == 0:
                energy[(i+1)//relax - 1] = self.energy()
        out = []
        for i in ret:
            if i == 'm':
                out.append(np.mean(energy))
            elif i == '2':
                out.append(np.mean(energy**2))
            elif i == 's':
                out.append(np.std(energy))
            elif i == 'v':
                out.append(np.var(energy))
            elif i == 'e':
                out.append(energy)
        if len(out) == 1:
            return out[0]
        else:
            return tuple(out)

    def status(self, ret='lp'):
        """ Displays or returns info about the current state. """
        # to bin how many edges have a given group element
        gl = np.zeros(self.group.size)
        # to bin how many plaquettes have a given group element
        gp = np.zeros(self.group.size)
        for i in multirange(self.shape):
            for j in range(self.num_dims - 1):
                gl[self.lattice.links[i+(j,)]] += 1
                for k in range(j+1, self.num_dims):
                    gp[self.site_plaquette(i,j,k,ret='p')] += 1
            gl[self.lattice.links[i + (self.num_dims,)]] += 1

        out = []
        for i in ret:
            if i == 'l':
                out.append(ge)
            elif i == 'p':
                out.append(gp)
        return tuple(out)

def hysteresis(field, betas, neq=2, nstat=10, relax=10, inc=10, display=True,
               avg=True, beta_out=True):
    """
    Scan inverse temperature through range given by betas and look at the
    energy of each field.

    Args:
        field (Field):
            Field instance on the lattice.
        betas (array-like):
            Array containing range of betas (inverse temperatures)
        neq (int):
            Number of equilibration sweeps to make at each new value of betas.
        nstat (int):
            Number of sweeps to do to get statistics on the variance of the
            energy (i.e. the heat capacity).
        relax (int):
            Number of steps to wait before collecting data to prevent
            correlated configurations.
        inc (int):
            Step increment.
        display (bool):
            Whether or not to display simulation progress.
        avg (bool):
            Whether to store the instantaneous energy or the average energy.
        beta_out (bool):
            Whether or not to include betas in return value(s).

    Returns:
        (energy, energy2, std, (optional: beta_out)) (tuple):
            energy, energy2, and std are the calculated values of the energy,
            energy squared, and the standard deviation of energies from the
            hysteresis loop.

    NOTE:
        betas can be a length 3 tuple, in which case beta is allowed to range
        from betas[0] to betas[1] in increments of betas[2].
    """
    if isinstance(betas, tuple):
        betas = np.concatenate((
            np.arange(betas[0], betas[1], betas[2]),
            np.arange(betas[1], betas[0], -betas[2])
        ))
    energy = np.zeros(betas.shape)
    energy2 = np.zeros(betas.shape)
    std = np.zeros(betas.shape)
    for i in range(len(betas)):
        if display:
            print("Step {} of {}".format(i, len(betas)))
        field.B = betas[i]  # update field temperature
        field.sweep(num_sweeps=neq)
        if not np.mod(i, inc) == 0:
            continue
        energy[i], energy2[i], std[i] = field.stats(nstat, relax=relax,
                                                    ret='m2s')
        if not avg:
            energy[i] = field.energy()
    energy = energy[::inc]
    std = std[::inc]
    betas = betas[::inc]
    if beta_out:
        return energy, std, betas
    else:
        return energy, std


def phase_sweep(field, betas, inc=1, num_sweeps=400, display=True):
    """
    Iterate over beta in betas, starting from a configuration initialized to be
    half random, and half equal to the group identity, and evolve num_sweeps
    sweeps, recording the energy every inc steps.

    Args:
        field (Field):
            Gauge field configuration.
        betas (array-like):
            Array containing range of betas (inverse temperatures) to iterate
            over.
        inc (int):
            How many steps between subsequent recordings of the energy.
        num_sweeps (int):
            Number of sweeps to perform.
        display (bool):
            Whether or not to print out the current beta of the simulation.

    Returns:
        energy (array-like):
            Array consisting of the recorded energies of the configuration at
            each beta.
    """
    energy = np.zeros((len(betas), 1+num_sweeps//inc))
    for i in range(len(betas)):
        if display:
            print("beta = {}".format(betas[i]))
        field.initialize('half')
        field.B = betas[i]
        energy[i, 0] = field.energy()
        for j in range(num_sweeps):
            field.sweep()
            if np.mod(j+1, inc) == 0:
                energy[i, (j+1)//inc] = field.energy()
    return energy

def watch_sweep(field, stop=-1, display=False):
    """ Performs a sweep and measure energy at each step along the way. """
    energy = np.zeros(np.product(field.links_shape) + 1)
    energy[0] = field.energy()
    j = 1
    for i in multirange(field.links_shape):
        if display:
            print(j)
        field.update(i)
        energy[j] = field.energy()
        j += 1
        if j == np.mod(stop, len(energy)):
            break
    return energy[0:j]

def view(field, s, d1, d2, fig=None, lw=10.):
    """
    View a cross-section with directions d1 and d2, passing through site s.
    """
    slc = list(v)
    slc[d1] = slice(None)
    slc[d2] = slice(None)
    slc = tuple(slc)
    cross_sec1 = field.lattice.links[slc + (d1,)]
    cross_sec2 = field.lattice.links[slc + (d2,)]
    figure(fig)
    for i in multirange(cross_sec1.shape):
        plt.plot(
            [i[0], i[0] + 1],
            [i[1], i[1]],
            color=ncolor(field.group.size, cross_sec1[i]),
            linewidth=lw
        )
        plt.plot(
            [i[0], i[0]],
            [i[1], i[1] + 1],
            color=ncolor(field.group.size,cross_sec2[i]),
            linewidth=lw
        )

def ncolor(N, i, lims=(0, 0.9, 0, 0.9, 0, 0.9)):
    """
    Generates a 'uniform' color palatte for N colors, returning the ith color.

    Args:
        N (int):
            Number of colors to construct.
        i (int):
            The index of the color to return.
        lims (tuple):
            The limits on how deep the colors can be, in rgb.

    Returns:
        [r, g, b] (list):
            The r, g, b values of the ith color in the generated palatte.
    """
    n = np.ceil(N**(1./3))
    b = np.mod(i, n)
    r = i//n**2
    g = np.mod(i//n, n)
    b /= (n - 1.)
    r /= (n - 1.)
    g /= (n - 1.)
    return [
        l[4] + r * (l[5] - l[4]),
        l[2] + g * (l[3] - l[2]),
        l[0] + b * (l[1] - l[0])
    ]


