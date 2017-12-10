from __future__ import division
import numpy as np
import itertools as it

def metropolis(s1, s2, beta):
    """
    Determines acceptance of updated configuration, using the Metropolis rule.

    Args:
        s1 (float):
            Value of the action from the initial configuration.
        s2 (float):
            Value of the action from the updated configuration.
        beta (float):
            Inverse temperature at which both configurations were generated.

    Returns:
        (bool):
            True if the updated configuration is accepted, else False.
    """
    delta_S = s2 - s1
    if delta_S < 0 or np.exp(-beta * delta_S) > np.random.random():
        return True
    return False

def unit_vec(n, d, l=1, dtype=int):
    """ Creates an n-dimensional vector of length l in direction d. """
    vec = np.zeros(n, dtype=dtype)
    vec[d] += l
    return vec

def multirange(t):
    """
    Generates range of indices in multiple dimensions. t should be iterable.
    """
    return it.product(*[range(i) for i in t])

def increment_array(arr, *increments, f=tuple):
    """
    Increments array at indices by amounts specified as pairs contained in
    increments.

    Args:
        arr (array-like):
            Array to be incremented.
        increments (array-like):
            Collection of pairs (idx, val), giving the index of the array to be
            incremented by val.
        f (return-type/method):
            Applied to the incremented array to be returned.
            If f specifies an array-like data type (list, tuple, np.array) the
            returned value is of type f.
            If f is a method, this method is applied before returning the
            incremented array.

    Returns:
        (array-like):
            Incremented array, as a copy of arr, whose data-type is specified
            by return-type of f (default is tuple).

    Examples:
        >>> arr = [1, 2, 3]
        >>> increment_array(arr, (0, 1), (1, 2), (2, 3), f=tuple)
        (2, 4, 6)
        >>> increment_array(arr, (0, 2), f=list)
        [3, 2, 3]
    """
    arr2 = np.array(arr).copy()
    for i in increments:
        arr2[i[0]] += i[1]
    return f(arr2)

def site_plaq_links(site, d1, d2):
    """
    Returns the plaquette links corresponding to the passed site, in the
    directions d1 and d2.
    """
    if d1 < d2:
        return (
            site + (d1,),
            increment_array(site, (d1, 1)) + (d2,),
            increment_array(site, (d2, 1)) + (d1,),
            site + (d2,)
        )
    elif d1 > d2:
        return (
            site + (d2,),
            increment_array(site, (d2, 1)) + (d1,),
            increment_array(site, (d1, 1)) + (d2,),
            site + (d1,)
        )

def link_plaq_links(link, direction, sign):
    """
    Returns the plaquette links corresponding to the passed link, direction and
    sign.
    """
    if sign == 1:
        if direction > link[-1]:
            return (
                link,
                increment_array(link, (link[-1], 1), (-1, direction-link[-1])),
                increment_array(link, (direction, 1)),
                increment_array(link, (-1, direction-link[-1]))
            )
        elif direction < link[-1]:
            return (
                increment_array(link, (-1, direction - link[-1])),
                increment_array(link, (direction,1)),
                increment_array(link, (link[-1], 1), (-1, direction-link[-1])),
                link
            )
    elif sign == -1:
        if direction > link[-1]:
            return (
                increment_array(link, (direction, -1)),
                increment_array(link,
                                (direction, -1),
                                (link[-1], 1),
                                (-1, direction - link[-1])),
                link,
                increment_array(link,
                                (direction, -1),
                                (-1, direction - link[-1]))
            )
        elif direction < link[-1]:
            return (
                increment_array(link,
                                (direction,-1),
                                (-1,direction-link[-1])),
                link,
                increment_array(link,
                                (direction,-1),
                                (link[-1],1),
                                (-1,direction-link[-1])),
                increment_array(link,(direction,-1))
            )
    else:
        raise ValueError("sign must be either +1, or -1")

