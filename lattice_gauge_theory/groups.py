from __future__ import division
import numpy as np
from lattice_gauge_theory.utils import *

def ZN(N, action=None):
    """
    Generates the group Z_N, and if 'action' is specified, returns the action.
    """
    table = np.mod(np.arange(N) + np.reshape(np.arange(N), (N, 1)), N)
    names = np.array(['(w{})^{}'.format(N, i) for i in range(N)])
    G = IntGroup(table, names=names)
    if action is None:
        return G
    elif action=='delta':
        ZN_action = lambda a: float(a != G.id)
        return G, ZN_action
    elif action == 'U1':
        ZN_action = lambda a: 1 - np.cos(2*np.pi*a/N)
        return G, ZN_action


class IntGroup(object):
    """
    Class for discrete gauge groups.

    Group elements are identified by integers 0, 1, ..., N-1, where N is the
    cardinality of the group.

    Specify a multiplication table to define the group.
    Representations are also interpreted.
    """

    def __init__(self, table, rep=None, names=None):
        """
        Args:
            table (np.array):
                Array of shape [N, N] with integer values from 0 to N-1.
            rep (dict or iterable):
                if dict:
                    dict[i] should be a np.array containing the matrix
                    representation of group element i.
                if iterable:
                    iterable[i] should be a matrix representation of element i.
            names (dict):
                Dictionary for associating names to integers to be used for
                identifying group elements.
        """
        # Error check on table
        if not (table.ndim==2 and table.shape[0] == table.shape[1]):
            raise ValueError("table must be a square array.")
        if not table.dtype == 'int':
            raise ValueError("table must be integer valued.")
        self.table = table
        N = table.shape[0]
        self.N = N # number of group elements
        self.size = N
        # check table vals
        if not (np.amax(table) == N-1 or np.amin(table) == 0):
            raise ValueError("table values must be from 0 to table.shape[0]-1")
        if isinstance(rep, dict):
            # if rep is a dict, convert to array
            r = np.array(N)
            for i in range(N):
                r[i] = dict[i]
            self.rep = r
        else:
            self.rep = np.array(rep)

        # Want self.names[i] to be the names of element i
        # Want self.names_['name'] to be the element with name 'name'
        if isinstance(names, dict):
            self.names_ = np.array(names)
            self.names = np.zeros(N, dtype=object)
            for i in names.keys():
                self.names[names[i]] = i
        elif hasattr(names, '__iter__'):
            self.names = names
            self.names_ = dict()
            for i in range(N):
                self.names_[names[i]] = i

        self.id = np.where(table[0] == 0)[0][0] # group identity element
        # inv[i] stores the group inverse of group element i
        self.inv = np.zeros(N, dtype=int)
        for i in range(N):
            self.inv[i] = np.where(table[i] == self.id)[0][0]
            if not (table[self.inv[i],table[i]]==np.arange(N, dtype=int)).all():
                raise ValueError("Element {} is not invertible".format(i))

        # conjugacy classes
        self.Cl = np.array([
            {self(i, j, self.inv[i]) for i in range(N)} for j in range(N)
        ])

    def __getitem__(self, idx):
        return self.table[idx]

    def __call__(self, *a):
        """ Group product. """
        b = a[0]
        for i in range(1, len(a)):
            b = self.table[b, a[i]]
        return b



########################### Some groups: ##########################
# Klein 4-group:
TKlein = np.array([[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]])
Klein = IntGroup(TKlein)
KDaction = lambda a: float(a!=Klein.id)

# Quaternion group:
# The element identification is 1->0, -1->1, i->2,-i->3,j->4,-j->5,k->6,-k->7
Qnames = ['1','-1','i','-i','j','-j','k','-k']
TQ = np.array([	[0,1,2,3,4,5,6,7],
				[1,0,3,2,5,4,7,6],
				[2,3,1,0,6,7,5,4],
				[3,2,0,1,7,6,4,5],
				[4,5,7,6,1,0,2,3],
				[5,4,6,7,0,1,3,2],
				[6,7,4,5,3,2,1,0],
				[7,6,5,4,2,3,0,1]])
Quaternion = IntGroup(TQ,names=Qnames)
QDaction = lambda a: float(a!=Quaternion.id)



