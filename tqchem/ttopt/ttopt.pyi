from _typeshed import Incomplete
from tqchem.ttopt.maxvol import maxvol as maxvol
from typing import Callable

class TTOpt:
    """Determine most stable conformer for a given input Molecule

    Attributes
    ----------
    n_points_per_variable: list[int]
        Number of points per variable
    function: Callable
        Objective function to be optimized for conformer search
    n_sweeps: int, default=2
        Number of sweeps used in optimization
    rank: int, default=2
        Maximum used rank
    initial_indices: list[list[int]], default=None
        Initial indices for precomputing the initial tensor train
    multiple_indices_function: bool, default=False
        Enable multiple indices evaluation
    """
    n_sweeps: Incomplete
    rank: Incomplete
    seed: Incomplete
    function: Incomplete
    n_torsions: Incomplete
    n_points_per_torsion: Incomplete
    multiple_indices_function: Incomplete
    i_opt: Incomplete
    y_opt_list: Incomplete
    opt_opt_list: Incomplete
    cache: Incomplete
    initial_indices: Incomplete
    def __init__(self, n_points_per_variable: list[int], function: Callable, n_sweeps: int = 2, rank: int = 2, seed: int = 42, initial_indices: list[list[int]] = None, multiple_indices_function: bool = False) -> None: ...
    @property
    def opt_opt(self):
        """Current value of option of the function related to opt-point."""
    @property
    def y_opt(self):
        """Current approximation of optimum of the function of interest."""
    def optimize(self) -> None:
        '''Find the optimum element of the implicitly given multidimensional array.

        This function computes the minimum of the implicitly given
        d-dimensional (d >= 2) array (tensor). The adaptive method based on the
        TT-approximation and the cross-maximum-volume principle are used.

        Returns:
            [np.ndarray, float]: the multi-index that gives the optimum value of the
            tensor (it is 1D np.ndarray of length "d" of int; i.e., "i_opt") and
            the optimum value of the tensor (it is float; i.e., "y_opt") that
            corresponds to the multi-index "i_opt".

        '''
    def comp(self, setsOfIndices):
        """Compute the function for the set of multi-indices (samples).

        Args:
            setsOfIndices (ndarray): 2D array of shape (nSamples, d) containing sets of indices that can be converted into sets of torsion angles.

        Returns:
            energies (ndarray): 1D array (float) of shape(number of new conformers) containing the absolute energies of the new conformers

        Note:
            The set of points (I) should not contain duplicate points. If it
            contains duplicate points (that are not in the cache), then each of
            them will be recalculated without using the cache.

        """
    def comp_opt(self, setsOfIndices, n_sweeps, i_opt: Incomplete | None = None, y_opt: Incomplete | None = None, opt_opt: Incomplete | None = None):
        '''Compute the function for the set of points and save current optimum.

        This helper function (this is wrapper for function "comp") can be
        passed to the optimizer. When making requests, the optimizer must pass
        the grid points of interest (setOfIndices) as arguments, as well as the current
        approximation of the argmin / argmax (i_opt), the corresponding value
        (y_opt) and related option value (opt_opt).

        '''

def ttopt_find(multi_indices, y, opt, i_opt, y_opt, opt_opt):
    """Find the minimum or maximum value on a set of sampled points."""
def ttopt_init(n_torsions, n_points_per_torsion, rank, seed):
    """Build initial approximation for the tensor train
    For simplicity all nodes in the tensor train are represented as
    3-dimensional tensors. The first/last node simply have their
    first/last dimension set to 1.

    edge_dimensions: Dimensions in the tensor train, excluding the leaves
    initial_tensor_train: n_torsions 3-dimensional tensors filled with random numbers
    """
