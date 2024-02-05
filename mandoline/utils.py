import numpy as np


def expand_array(arr, factor):
    """
    Data reading utility
    ----
    Expand lower resolution 2D array by [factor]
    to broadcast it to a higher level grid.
    This allows broadcasting lower resolution arrays to a higher 
    AMR level grid without altering the data.
    ----
    Example:
    >> expand_array([[1, 2],
                     [3, 4]], factor=2)
    >> [[1, 1, 2, 2],
        [1, 1, 2, 2],
        [3, 3, 4, 4],
        [3, 3, 4, 4]]
    """
    exp = np.repeat(arr, factor).reshape(arr.shape[0], 
                                         arr.shape[1]*factor)
    exp = np.repeat(exp, factor, axis=0).reshape(arr.shape[0]*factor, 
                                                 arr.shape[1]*factor)
    return exp


def sanitize_field_name(fname):
    """
    Remove parentheses from field names
    so that we dont breack file systems
    Y(H) becomes Y_H
    """

    return fname.replace('(','_').replace(')','')
