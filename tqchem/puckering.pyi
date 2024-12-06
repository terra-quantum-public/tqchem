import networkx as nx
import numpy as np
from dataclasses import dataclass
from tqchem.batframe import Angle as Angle, BATFrame as BATFrame, Bond as Bond, Dihedral as Dihedral, bat_to_xyz as bat_to_xyz, is_improper_dihedral as is_improper_dihedral, xyz_to_bat as xyz_to_bat
from tqchem.molgraph import MolecularSystem as MolecularSystem, bond_component as bond_component

@dataclass
class PuckeringFrame:
    """Pucker coordinate frame.

    This class represents a frame in puckering coordinates.

    Puckering is implemented according to:
    D. Cremer, J. Pople, https://doi.org/10.1021/ja00839a011

    Attributes
    ----------
    x: np.ndarray
        in-plane x-coordinates of each atom
    y: np.ndarray
        in-plane y-coordinates of each atom
    q: np.ndarray
        (n_atoms - 3) out-of-plane amplitudes
    phi: np.ndarray
        (n_atoms - 3) out-of-plane angles
    batframe: BATFrame
        BATFrame connecting the ring to the rest of the molecule
    cycle: list[int]
        list of atoms in the cycle
    planar_only: bool
        Should only the planar conformation be used?
    """
    x: np.ndarray = ...
    y: np.ndarray = ...
    q: np.ndarray = ...
    phi: np.ndarray = ...
    batframe: BATFrame = ...
    cycle: list = ...
    planar_only: bool = ...
    def update(self, frame) -> None: ...
    def variable_coordinates(self, molecule: MolecularSystem, **kwargs) -> PuckeringFrame:
        """Select coordinates for optimization of conformers.

        This includes the all puckering coordinates and a batframe
        connecting it to possible other parts of the molecule described
        using BAT coordinates.
        """
    def constant_coordinates(self, molecule: MolecularSystem) -> PuckeringFrame:
        """Get subset of coordinates which need to stay constant during conformer search

        Only coordinates of the batframe which are not variable are added
        """
    grid = ...
    def make_grid(self, graph: nx.Graph = None, shift_by_reference: bool = False):
        """Discretize coordinates of puckering frame for optimization."""
    def evaluate_grid(self, indices: list[int]):
        """Evaluate the puckering frame at a given grid index.

        The first parameter in i will be used for the puckering frame
        and the remaining ones for the BAT frame.
        """
    def grid_range(self):
        """Obtain ranges for each of the coordinates."""
    def closest_gridpoint(self, molecule: MolecularSystem) -> list[int]:
        """Return the gridpoint most similar to the provided molecule"""
    def __init__(self, x=..., y=..., q=..., phi=..., batframe=..., cycle=..., planar_only=...) -> None: ...

def reference_coordinate_system(mol: MolecularSystem, cycle: tuple[int, ...]) -> tuple[np.array, np.array, np.array, np.array]:
    """Obtain x0, y0, z0 from the current geometry"""
def puckering_to_xyz(mol: MolecularSystem, frame: PuckeringFrame) -> MolecularSystem:
    """Transform puckering coordinates to xyz coordinates."""
def xyz_to_puckering(mol: MolecularSystem, cycle: tuple[int, ...]) -> PuckeringFrame:
    """Transform xyz coordinates to puckering coordinates."""
def pucker_coordinate_5ring(reference: PuckeringFrame, amplitude: float = 0.8, n_phi: int = 12) -> list[PuckeringFrame]:
    """Generate sets of values representing the most relevant conformers of a 5 ring

    As shown in Figure 1 of https://s3.smu.edu/dedman/catco/ring-puckering.html

    By default, it creates 13 frames that generate all usual conformers,
    where each edge of the ring points up and down once and it includes a planar geometry

    Parameters
    ----------
    reference: PuckeringFrame
        reference frame to build the different ring conformers from
    amplitude: float
        amplitude of the out-of-plane displacement
    n_phi: int
        number of non-planar conformers created
        Specifies the discretization of the pseudorotational cycle shown in figure 1
    """
def pucker_coordinate_6ring_smart(ref: PuckeringFrame, amplitude: float = 0.7, n_phi: float = 12, equatorial: bool = False, tropical: bool = False, planar: bool = True, axial: bool = True) -> list[PuckeringFrame]:
    """Selected generation of puckering points for 6 rings.

    Select which parts of the puckering sphere should be included by setting the flags.
    Please see Fig. 2 of https://s3.smu.edu/dedman/catco/ring-puckering.html for reference.

    axial: include 2 chair conformations
    planar: include planar conformation
    equatorial: include Boat conformation
    tropical: include envelope conformation
    """
def pucker_coordinate_6ring(reference: PuckeringFrame, amplitude: float = 0.7, n_phi: int = 12, n_theta: int = 5) -> list[PuckeringFrame]:
    """Generate sets of values representing the most relevant conformers of a 5 ring

    As shown in Figure 2 of https://s3.smu.edu/dedman/catco/ring-puckering.html

    By default, it creates 39 frames that generate all usual conformers,
    including planar geometry, chair, twist, and boat conformations.
    """
def planar_conformer(reference: PuckeringFrame) -> PuckeringFrame:
    """Returns planar ring conformer"""
def is_planar_ring(cycle: list[int], graph: nx.Graph) -> bool: ...
