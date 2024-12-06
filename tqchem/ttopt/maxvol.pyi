from _typeshed import Incomplete

def maxvol(A, e: float = 1.05, k: int = 100):
    '''Compute the maximal-volume submatrix for given tall matrix.

    Args:
        A (np.ndarray): tall matrix of the shape [n, r] (n > r).
        e (float): accuracy parameter (should be >= 1). If the parameter is
            equal to 1, then the maximum number of iterations will be performed
            until true convergence is achieved. If the value is greater than
            one, the algorithm will complete its work faster, but the accuracy
            will be slightly lower (in most cases, the optimal value is within
            the range of 1.01 - 1.1).
        k (int): maximum number of iterations (should be >= 1).

    Returns:
        (np.ndarray, np.ndarray): the row numbers I containing the maximal
        volume submatrix in the form of 1D array of length r and coefficient
        matrix B in the form of 2D array of shape [n, r], such that
        A = B A[I, :] and A (A[I, :])^{-1} = B.

    Note:
        The description of the basic implementation of this algorithm is
        presented in the work: Goreinov S., Oseledets, I., Savostyanov, D.,
        Tyrtyshnikov, E., Zamarashkin, N. "How to find a good submatrix".
        Matrix Methods: Theory, Algorithms And Applications: Dedicated to the Memory of Gene Golub (2010): 247-256.

    '''
def maxvol_rect(A, e: float = 1.1, dr_min: int = 0, dr_max: Incomplete | None = None, e0: float = 1.05, k0: int = 10):
    '''Compute the maximal-volume rectangular submatrix for given tall matrix.

    Within the framework of this function, the original maxvol algorithm is
    first called (see function maxvol). Then additional rows of the matrix
    are greedily added until the maximum allowed number is reached or until
    convergence.

    Args:
        A (np.ndarray): tall matrix of the shape [n, r] (n > r).
        e (float): accuracy parameter.
        dr_min (int): minimum number of added rows (should be >= 0 and <= n-r).
        dr_max (int): maximum number of added rows (should be >= 0). If the
            value is not specified, then the number of added rows will be
            determined by the precision parameter e, while the resulting
            submatrix can even has the same size as the original matrix A.
            If r + dr_max is greater than n, then dr_max will be set such
            that r + dr_max = n.
        e0 (float): accuracy parameter for the original maxvol algorithm
            (should be >= 1). See function "maxvol" for details.
        k0 (int): maximum number of iterations for the original maxvol algorithm
            (should be >= 1). See function "maxvol" for details.

    Returns:
        (np.ndarray, np.ndarray): the row numbers I containing the rectangular
        maximal-volume submatrix in the form of 1D array of length r + dr,
        where dr is a number of additional selected rows (dr >= dr_min and
        dr <= dr_max) and coefficient matrix B in the form of 2D array of shape
        [n, r+dr], such that A = B A[I, :].

    Note:
        The description of the basic implementation of this algorithm is
        presented in the work: Mikhalev A, Oseledets I., "Rectangular
        maximum-volume submatrices and their applications." Linear Algebra and
        its Applications 538 (2018): 187-211.

    '''
