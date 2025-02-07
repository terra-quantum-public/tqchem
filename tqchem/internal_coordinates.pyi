from _typeshed import Incomplete
from tqchem.batframe import bat_component as bat_component, bat_to_xyz as bat_to_xyz, xyz_to_bat as xyz_to_bat
from tqchem.molgraph import MolecularSystem as MolecularSystem, linearized_component_sweep as linearized_component_sweep
from tqchem.puckering import puckering_to_xyz as puckering_to_xyz, xyz_to_puckering as xyz_to_puckering
from typing import Any, Iterable

class InternalCoordinates(dict):
    """Collection of internal coordinates for components of a system.

    Dictionary mapping tuples of atom indices to a representation
    of internal coordinates, i.e. BATFrames or PuckeringFrames.

    For conformer search we do not want to alter all coordinates,
    so Internal coordinates can be filtered for relevant ones (filter_subset)
    and there are coordinates that we don't want to alter, but we need to
    make sure that they are correctly set (constant_coordinates)
    """
    def update(self, coords: dict) -> None:
        """Update internal coordinates with (subset of) internal coordinates"""
    def combine(self, coords: dict) -> None:
        """Combine internal coordinates with (subset of) internal coordinates"""
    def variable_coordinates(self, molecule: MolecularSystem, atom_order: list[int]) -> InternalCoordinates:
        """Get subset of coordinates which are relevant/varied during conformer search

        atom_order is a list of atom ids in the graph and determines the
        order in which the variables will be listed in the subset
        """
    def constant_coordinates(self, molecule: MolecularSystem) -> InternalCoordinates:
        """Get subset of coordinates which need to stay constant during conformer search

        This entails for example improper dihedrals for bonds going into and out
        of a ring, as well as bond lengths and bond angles
        """
    def make_grid(self, **kwargs) -> None:
        """Make a discrete grid for each coordinate of each component."""
    def evaluate_grid(self, indices: dict[tuple[int, ...], list[int]]) -> InternalCoordinates:
        """Evaluate indices on the grids and update the internal coordinates with the new values

        indices: dict[tuple[int, ...], list[int]
            Maps atom indices of the components to a list of indices for each
            coordinate in the component
        """
    def grid_range(self) -> dict[tuple[int, ...], list[range]]:
        """Return grid ranges for each component"""

class MolecularGrid:
    """Grid representation of the conformers of a system

    Attributes
    ----------
    molecule: Molecular System
        Graph representation of the molcule and interface to ase to modify
        the atom positions
    variables: InternalCoordinates
        Subset of the InternalCoordinates of the system
        (mainly Dihedrals and puckering coordinates) which are represented
        on the grid, to perform Conformer Search.
    constant_coordinates: InternalCoordinates
        Subset of InternalCoordinates (Bonds, Angles,
        improper Dihedrals used for the connection of Rings)
        which should not change during the Conformer Search, but need to be set
        to make sure they remain the same.

    Methods
    -------
    __call__:
        Return new molecule based on the gridpoints received on input
    shape:
        Return range of indices for each element of the Molecular Grid
    n_variables:
        Return number of variable elements in the Molecular Grid
    variable_atom_indices:
        Return list of the variable elements in the Molecular Grid
    closest_gridpoint:
        Get grid point closest to the geometry/molecule system provided on input
    """
    molecule: Incomplete
    variables: Incomplete
    constant_coordinates: Incomplete
    def __init__(self, molecule: MolecularSystem = None, shift_by_reference: bool = False, **bond_grids) -> None:
        """Initialize Molecular grid used for conformer search

        Parameters
        ----------
        molecule: MolecularSystem
            Molecule which should be discretized
        shift_by_reference: bool, default=False
            Set dihedrals based on the provided reference structure
        **bond_grids:
            Keywords specifying which grids to use when setting dihedrals
            Arguments are passed to MolecularSystem.set_rotatable_bond_grids
        """
    def __call__(self, superidx: list[int]) -> MolecularSystem:
        """Construct new molecular system based on the incoming indices

        Indices are delinearized to match the structure of the internal coordinates.
        The variables are updated (e.g. with new dihedral values), combined with
        the constants, and converted to a new MolecularSystem with xyz coordinates.

        Parameters
        ----------
        superidx: list[int]
            List of indices, one for each variable coordinate
        """
    def shape(self) -> list[int]: ...
    def n_variables(self): ...
    def closest_gridpoint(self, molecule: MolecularSystem) -> list[int]:
        """Return the gridpoint most similar to the provided molecule"""
    def list_variables(self) -> list[str]:
        """List of variables (Bonds, cycles, inversions) in the molecular grid"""
    def copy(self) -> MolecularGrid:
        """Deepcopy of the molecular grid"""

def linearize(dictionary: dict[Any, Iterable]) -> list[Any]:
    """Converts dictionary values (Iterables) into a one-dimensional list

    Used to convert a dictionary relating components to a list of ranges for
    the component's variables to a list which only contains the ranges

    Example:
    Three components with 2, 1 and 0 dihedrals

    >>> dictionary = {
            (0, 1, 2, 3, 4, 5): [range(0, 12), range(0, 10)],
            (6, 7, 8, 9): [range(0, 8)],
            (10): [],
        }
    >>> linearize(dictionary)
    [range(0, 12), range(0, 10), range(0, 8)]
    """
def indices_to_components(indices: list[int], component_ranges: dict[tuple[int], list[range]]):
    """Assign list of indices to the variables of the components

    component_ranges is a dict mapping each component to a list
    containing as many ranges as there are variables in the component.
    This function returns an equivalent dict, but the ranges are replaced
    with the indices.
    """
def internal_to_xyz(molecule: MolecularSystem, coordinates: InternalCoordinates) -> MolecularSystem:
    """Transform internal coordinates to xyz coordinates.

    Internal Coordinates can consist of multiple components (BAT, puckering, BAT, ...)
    """
def xyz_to_internal(molecule: MolecularSystem) -> InternalCoordinates:
    """Transform xyz coordinates to internal coordinates.

    Internal Coordinates can consist of multiple components (BAT, puckering, BAT, ...)
    """
def molecule_coordinates(molecule: MolecularSystem) -> tuple[InternalCoordinates, InternalCoordinates]:
    """Calculates and returns the internal coordinates for a molecular system"""
