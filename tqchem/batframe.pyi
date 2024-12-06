import networkx as nx
from _typeshed import Incomplete
from tqchem.functional_groups import match_atomic_numbers as match_atomic_numbers
from tqchem.molgraph import MolecularSystem as MolecularSystem, bond_component as bond_component, component_type as component_type, get_component_of_atom as get_component_of_atom, get_path as get_path, split_bond_to_subgraphs as split_bond_to_subgraphs

class Bond(tuple):
    """Class representing a bond between two atoms

    based on the tuple class a bond is described by the uid of the
    2 atoms involved

    Attributes
    ----------
    tuple:
        containing uid of atom 0 and atom 1
    """
    def __new__(cls, bond: tuple[int, int], mask: list[int] = None):
        """Initialize Bond as tuple containing 2 atom indices and potentially a mask"""

def bond_mask(bond: Bond, open_graph: nx.Graph) -> list[int]:
    """List of atom indices that should be changed together with the given bond

    The mask should contain the last atom of the bond (the one moved)
    and the connected atoms that rotate with it.
    This mask is passed to set_bond of ase.Atoms
    """

class Angle(tuple):
    """Class representing an angle between three atoms

    based on the tuple class an angle is described by the indices of the
    3 atoms involved

    Attributes
    ----------
    tuple:
        containing uids of atom0 to atom2
    """
    def __new__(cls, angle: tuple[int, int, int], mask: tuple[int, ...] = None):
        """Initialize Angle as tuple containing 3 atom indices and potentially a mask"""

def angle_mask(angle: Angle, open_graph: nx.Graph) -> list[int]:
    """List of atom indices that should change together with the given angle

    The mask should contain the last atom of the angle (the one moved)
    and the connected atoms that should move with it.
    This mask is passed to set_angle of ase.Atoms
    """

class Dihedral(tuple):
    """Class representing an dihedral/torsion angle between four atoms

    based on the tuple class a dihedral is described by the indices of the
    4 atoms involved

    Attributes
    ----------
    tuple:
        containing uids of atom0 to atom3
    """
    def __new__(cls, dihedral: tuple[int, int, int, int], mask: tuple[int, ...] = None):
        """Initialize Dihedral as tuple containing 4 atom indices and potentially a mask"""

class InversionDihedral(Dihedral):
    '''Class representing a dihedral used for Nitrogen inversion

    Dummy class to make Dihedrals used for nitrogen inversion distinguishable
    and to save the 2 dihedral values defining the "normal" and "inverted" state
    '''
    def __new__(cls, dihedral: tuple[int, int, int, int], allowed_values: list[float], mask: tuple[int, ...] = None):
        """Initialize Dihedral as tuple containing 4 atom indices, the dihedral values
        defining the normal and inverted state and optionally a mask"""

def dihedral_mask(dih: Dihedral, open_graph: nx.Graph) -> list[int]:
    """Get a list of atom indices that should rotate together with the given dihedral

    The mask should contain the atom that is rotated around which we rotate
    and the connected atoms that rotate with it.
    This mask is passed to set_dihedral of ase.Atoms
    """
def masked_coordinate(coordinate: Bond | Angle | Dihedral, open_graph: nx.Graph) -> Bond | Angle | Dihedral:
    """Return copy of coordinate containing its appropriate mask"""

class BATFrame(dict):
    """Class/Dictionary storing the bonds, angles and dihedrals and their values

    For a given set of atoms this class contains the information required
    construct a conformer or a frame of an MD trajectory.

    Torsions are implemented as Dihedrals with a mask rotating
    all atoms attached to the end of a given bond.

    Methods
    -------
    angles(self):
        return dictionary relating the angles of the frame to their size
    bonds(self):
        return dictionary relating the bonds of the frame to their lengths
    dihedrals(self):
        return dictionary relating the dihedrals of the frame to their size
    copy(self):
        return a copy of the BATFrame
    """
    def angles(self) -> BATFrame:
        """Return dictionary of all angles and their sizes"""
    def bonds(self) -> BATFrame:
        """Return dictionary of all bonds and their lengths"""
    def dihedrals(self) -> BATFrame:
        """Return dictionary of all dihedrals and their sizes"""
    def copy(self) -> BATFrame:
        """Return copy of the BATFrame"""
    def torsions_from_bonds(self, bonds: list[tuple[int, int]], molecule: MolecularSystem, component: list[int]) -> BATFrame:
        """Find torsions from list of bonds"""
    def variable_coordinates(self, molecule: MolecularSystem, atom_order: list[int], component: list[int]) -> BATFrame:
        """Determine rotatable bonds and collect torsions"""
    def constant_coordinates(self, molecule: MolecularSystem) -> BATFrame:
        """Determine constant coordinates

        These coordinates have to be reset after changing variable bonds.
        This includes for example bond lengths and angles
        """
    grid: Incomplete
    def make_grid(self, graph: nx.Graph, shift_by_reference: bool = False) -> None:
        """Discretize all torsions in the frame"""
    def evaluate_grid(self, indices: list[int]) -> BATFrame:
        """Evaluate the frame at a given grid index"""
    def grid_range(self) -> list[range]:
        """Return the grid range for all torsions in the frame"""
    def closest_gridpoint(self, molecule: MolecularSystem) -> list[int]:
        """Return the gridpoint most similar to the provided molecule"""
    def __add__(self, value: float) -> BATFrame: ...
    def __sub__(self, value: float) -> BATFrame: ...

def is_proper_dihedral(dihedral: Dihedral, graph: nx.Graph) -> bool:
    """Does a dihedral represent a rotation around an actual bond"""
def is_improper_dihedral(dihedral: Dihedral, graph: nx.Graph) -> bool:
    """Does a dihedral represent a rotation around an actual bond"""
def torsion_from_bond(bond: Bond, frame: BATFrame, graph: nx.Graph) -> dict[Dihedral, float]:
    """Get a Dihedral that rotates the subsituents of an entire bond

    The rotation around a bond is achieved by adding all substituents of the
    second atom of the bond to a mask such that they are rotated together

    Parameters
    ----------
    bond: Bond
        Bond for which dihedrals should be found
    frame: BATFrame
        Dictionary mapping Bonds, Angles, Dihedrals to values
    graph: nx.Graph
        Open graph representing the molecule

    Returns
    -------
    dict[Dihedral, float]

    """
def rotatable_bonds(graph: nx.Graph) -> list[tuple[int, int]]:
    """Determine rotatable bonds from graph representation of molecule"""
def xyz_to_bat(mol: MolecularSystem, frame: BATFrame) -> BATFrame:
    """Transform XYZ to BAT coordinates

    Parameters
    ----------
    mol: MolecularSystem
        Molecule as represented in tqchem
    frame: BATFrame
        BATFrame specifying the atoms involved in the BATcoordinates

    Returns
    -------
    frame: BATFrame
        frame with overwritten values for the BATCoordinates
    """
def bat_to_xyz(mol: MolecularSystem, frame: BATFrame) -> MolecularSystem:
    """Transform BAT to XYZ coordinates

    Parameters
    ----------
    mol: MolecularSystem
        Molecule as represented in tqchem
    frame: BATFrame
        Internal coordinates (defined by atom indices) and their values

    Returns
    -------
    new_molecule: MolecularSystem
        Molecule with geometry updated  based on the BATFrame
    """
def improper_dihedral_path_ring(cycle: tuple[int, ...], connecting_atom: int) -> list[int]:
    """Get path for improper dihedral angle in a ring"""
def get_dihedral_path(component: tuple[int, ...], component_sweep: list[tuple[int, ...]], molgraph: MolecularSystem, atom: int, root: int, improper_ring_dihedral: bool = True) -> list[int]:
    """Get path to create dihedral angle for a given atom from MolGraph

    Based on a tree graph, the path to create a dihedral angle for a given atom
    is determined.

    Arguments
    ---------
    component: tuple[int, ...]
        atoms for the component of the whole molecule
    component_sweep: list[tuple[int, ...]
        list of components
    molgraph: MolecularSystem
        Molecular graph and ase.Atoms object
    atom: int
        atom to which the dihedral path leads
    root: int
        root of the path
    improper_ring_dihedral: bool
        if True, the dihedral path for a ring uses improper dihedral

    Returns
    -------
    list[int]:
        path to create dihedral angles for a given atom

    @todo: component and atom only used as component[atom] so we could maybe avoid both arguments?
           is atom an index of the atom in the component?
    """
def bat_component(mol: MolecularSystem, component: list[int] = None) -> BATFrame:
    """Constructor for a BATFrame with values set to None

    mol: MolecularSystem
        Molecular graph and ase.Atoms object
    component: list[int]
        atoms for which bat system will be constructed
    """
def nitrogen_inversion(nitrogen: int, normal_torsion: dict[Dihedral, float], molecule: MolecularSystem) -> None | dict[InversionDihedral, float]:
    '''Return dihedral describing nitrogen inversion

    Assuming we have the following molecule: OCN(C)CC
    The normal torsion is the rotation around this bond: OC-NCC
    We can describe nitrogen inversion by rotating the methyl group at the
    nitrogen atom (OC-N(C)) separately.
    This function returns the Torsion (i.e. {Dihedral: value}) that describes
    this nitrogen inversion.
    Note, that the inversion dihedral needs to be set before the "normal" dihedral
    because otherwise we would overwrite the methyl group set by the normal dihedral
    so that we would rotate both residues of the Nitrogen atom independently.

    None is returned if the nitrogen atom is symmetrically coordinated,
    as for example in -NH2, or only participates in 1 or 2 bonds.
    '''
