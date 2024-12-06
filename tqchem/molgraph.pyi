import ase
import networkx as nx
import numpy as np
from _typeshed import Incomplete
from pathlib import Path
from rdkit import Chem as Chem
from tqchem.chem import ase_from_path as ase_from_path, ase_from_rdkit as ase_from_rdkit, ase_from_smiles as ase_from_smiles, ase_molecule as ase_molecule, determineBonds as determineBonds, element_color as element_color, wiberg_bond_orders as wiberg_bond_orders
from tqchem.functional_groups import amide_bonds_grid as amide_bonds_grid, amino_group_grid as amino_group_grid, find_methyl_groups as find_methyl_groups, iso_propyl_groups_grid as iso_propyl_groups_grid, simple_alkanes_grid as simple_alkanes_grid, single_bonded_chalcogens_grid as single_bonded_chalcogens_grid, t_butyl_groups_grid as t_butyl_groups_grid
from typing import Any

class MolecularSystem:
    """Class representing a system of one or more molecules.

    Serves as the central management class in tqchem, containing an ase.Atoms
    object for energy evaluations and setting dihedral values,
    the molecular graph and an opened version of the graph,
    where either multirings are opened to single rings are opened or
    the entire graph is opened to a tree.
    Additionally it contains the component sweep directing the iteration
    through the components (rings, batframes, bridges) of the graph

    Attributes
    ----------
    graph: nx.Graph
        representation of the molecule by atoms & bonds.
    open_graph: nx.Graph
        graph with removed bonds to remove multirings in case of fragmentation
        or open the entire graph to a tree in case of fragmentation turned off
    component_sweep: list[list[int]]
        List containing the components (rings, bridges, trees) as list of atom indices

    Methods
    -------
    component_type:
        Get type (ring, batframe, bridge) for a given component
    draw:
        Draw molecular graph in initial and/or opened form
    root:
        Get root atom of the graph
    copy:
        Deepcopy of the molecular object
    write:
        Write geometry to file. Wrapper for ase.io.write
    get_positions:
        Get atomic positions. Wrapper for ase.atoms.get_positions
    get_path:
        Get path from the root of the graph to a given atom
    """
    ase: Incomplete
    graph: Incomplete
    open_graph: Incomplete
    component_sweep: Incomplete
    def __init__(self, molecule: MolecularSystem = None) -> None: ...
    def component_type(self, component: list[int]) -> str:
        """For a given component return if it's a batframe, ring or bridge"""
    def draw(self, graph: bool = True, open_graph: bool = False) -> None:
        """Draw the graph and/or the open graph of the system using matplotlib"""
    def root(self) -> int:
        """Return root node of the graph"""
    def copy(self) -> MolecularSystem:
        """Deepcopy of the molecular system"""
    def write(self, filename) -> None:
        """Wrapper for ase.io.write"""
    def get_positions(self):
        """Return atomic positions of the molecule"""
    def get_path(self, atom: int) -> list[int]:
        """Get path from the root node to the given atom/node-id"""
    def add_rotatable_bonds(self, bond_grids=...) -> None:
        """Set bonds to be rotatable, but don't overwrite previously rotatable bonds"""
    def overwrite_rotatable_bonds(self, bond_grids=...) -> None:
        """Unset all rotatable bonds and only set the ones specified in the dict"""
    def set_bond_orders(self) -> None:
        """Determine Bond orders and set them as edge attributes"""
    def set_rotatable_bond_grids(self, include_amide_bonds: bool = False, add_rotatable_bonds: dict[tuple[int, int], np.array] = None, set_rotatable_bonds: dict[tuple[int, int], np.array] = None, amide_grid: np.array = None, amino_grid: np.array = None, alkane_grid: np.array = None, t_butyl_grid: np.array = None, iso_propyl_grid: np.array = None, chalcogen_grid: np.array = None, default_grid: np.array = None, overwrite_grids: np.array = None) -> None:
        """Determine rotatable bonds in the system and assign corresponding bond grids

        Parameters
        ----------
        include_amide_bonds: bool = False,
            Include amide bonds as rotatable bonds
        add_rotatable_bonds: dict[tuple[int, int], np.array] = None,
            Adds rotatable bonds to the default ones in terms of a dict relating bonds
            to an array of allowed dihedral values
        set_rotatable_bonds: dict[tuple[int, int], np.array] = None,
            Overwrites rotatable bonds with a dict relating bonds to an array of
            allowed dihedral values
        amide_grid: np.array, default=[0.0, 180.0]
            Grid of allowed dihedral values for amide bonds
        amino_grid: np.array, default=[0, 60, 120, 180, 240, 300]
            Grid of allowed dihedral values for amino bods
        alkane_grid: np.array, default=[0, 60, 120, 180, 240, 300]
            Grid of allowed dihedral values for alkane bonds
        t_butyl_grid: np.array, default=[0.0, 60.0]
            Grid of allowed dihedral values for t_butyl
        iso_propyl_grid: np.array, default=[0.0, 60.0, 120.0]
            Grid of allowed dihedral values for iso_propyl
        chalcogen_grid: np.array, default=[0, 60, 120, 180, 240, 300]
            Grid of allowed dihedral values for chalcogen bonds
        default_grid: np.array, default=[0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
            Grid of allowed dihedral values for all other bonds
        overwrite_grids: np.array, default=None
            Overwrite all bond specific grids with this grid
        """

def molecule_from_ase(mol_ase: ase.Atoms, single_batframe: bool = False) -> MolecularSystem:
    """Construct MolecularSystem from ase.Atoms object.

    Parameters
    ----------
    mol_ase: ase.Atoms
        The molecular system
    single_batframe: bool = False,
        Use a single batframe for the entire molecule instead of using
        special puckering coordinates for the rings

    Returns
    -------
    MolecularSystem
    """
def molecule_from_rdkit(rd_mol: Chem.Mol, single_batframe: bool = False) -> MolecularSystem:
    """Construct MolecularSystem from rdkit.Mol object.

    Parameters
    ----------
    rd_mol: Chem.Mol
        Rdkit Molecule object
    single_batframe: bool = False,
        Use a single batframe for the entire molecule instead of using
        special puckering coordinates for the rings
    """
def molecule_from_file(file_: str | Path, format_: str = None, single_batframe: bool = False) -> MolecularSystem:
    """Construct MolecularSystem from molecule file.

    Parameters
    ----------
    file_: Path | str,
        Path to molecule file
    format_: str = None,
        File format in case the input is a str or Path and the format
        cannot be inferred from the file suffix.
    single_batframe: bool = False,
        Use a single batframe for the entire molecule instead of using
        special puckering coordinates for the rings
    """
def molecule_from_smiles(smiles: str, single_batframe: bool = False) -> MolecularSystem:
    """Construct MolecularSystem from SMILES string.

    Parameters
    ----------
    smiles: str
        Smiles string
    single_batframe: bool = False,
        Use a single batframe for the entire molecule instead of using
        special puckering coordinates for the rings
    """
def molecularSystem(input_: Path | str | ase.Atoms | Chem.Mol, format_: str = None, single_batframe: bool = False) -> MolecularSystem:
    """Smart factory for MolecularSystem

    Parameters
    ----------
    input_: Path | str | ase.Atoms | Chem.Mol,
        Input for the molecular system
    format_: str = None,
        File format in case the input is a str or Path and the format
        cannot be inferred from the file suffix.
    single_batframe: bool = False,
        Use a single batframe for the entire molecule instead of using
        special puckering coordinates for the rings

    Returns
    -------
    MolecularSystem
    """
def drawGraph(graph: nx.Graph) -> None:
    """Visualize the graph used for the molecule

    Parameters
    ----------
    graph: nx.Graph
        Graph representing the molecule
    """
def tree_from_graph(graph: nx.Graph) -> nx.Graph:
    """Make tree from graph by opening cycles in the graph

    Parameters
    ----------
    graph: nx.Graph
        Graph representing the molecule

    Returns
    -------
    oepn_graph: nx.Graph
        Copy of the Graph with opened cycles

    @todo: a clear greedy algorithm for maintaining the longest path in the tree
           might be preferable. The current implementation just arbitrarily removes
           edges until no cycles are left.
    """
def longest_path_in_tree(tree: nx.Graph) -> list[int]:
    """Find longest path going from one end of the tree to the other

    Parameters
    ----------
    tree: nx.Graph
        Tree representing the molecule

    Returns
    -------
    longest_path: list[int]
        Longest path in the tree as a list of nodes
    """
def ase_to_graph(molecule: ase.Atoms) -> nx.Graph:
    """Convert an ase molecule to a networkx graph

    Parameters
    ----------
    molecule: ase.Atoms
        Atoms representing the molecule
    Returns
    -------
    graph: nx.Graph
        networkx Graph representing the connectivities in the molecule
    """
def remove_cycle_overlap(cycles: list[list[int]]) -> list[list[int]]:
    """Remove overlap between cycles.

    Parameters
    ----------
    cycles: list[list[int]]
        list of lists containing the nodes that each cycle consits of
    Returns
    -------
    orthogonal_cycles: list[list[int]]
        list of lists containing the cycle nodes, but without duplicate nodes
    """
def residual_trees(graph: nx.Graph, cycles_bridges: list[list[int]]) -> list[list[int]]:
    """Find remaining trees in the graph after removing cycles and bridges."""
def cycle_bridge_tree_decomposition(graph: nx.Graph) -> list[list[int]]:
    """Decompose a graph into cycles, bridges, and trees"""
def dfs_sweep(graph: nx.Graph):
    """perform a depth-first search and return the sweep order"""
def in_list_of_lists(element: Any, list_of_lists: list[list[Any]]) -> bool: ...
def index_list_of_lists(element: Any, list_of_lists: list[list[Any]]) -> int: ...
def component_type(component: list[int], graph: nx.Graph) -> str:
    """Determine whether component is a ring, bridge, or tree."""
def sorted_components(components: list[list[int]], sweep: list[int]) -> list[list[int]]:
    """Sort components according to order in sweep"""
def get_root(component_sweep: list[list[int]]) -> int:
    """Get the root of the graph which serves as the beginning of the sweep"""
def sorted_ring(ring: list[int], component_sweep: list[list[int]], molgraph: nx.Graph) -> list[int]:
    """Sort the ring so that the first element connects to previous component."""
def sort_tree_component(component: list[int], sweep: list[int]):
    """Sort the tree component according to the sweep order"""
def open_multiring_graph(graph: nx.Graph) -> nx.Graph:
    """Open multirings (overlapping rings) of a graph to simple rings"""
def sorted_atoms_within_components(components, sweep, graph):
    """Sort atoms within components according to sweep order."""
def get_component_sweep(graph: nx.Graph, single_batframe: bool = False):
    """Determine the component sweep based on the graph representation of the molecule."""
def get_component_of_atom(atom: int, component_sweep: list[list[int]]):
    """Return the component in which the atom is."""
def linearized_component_sweep(component_sweep: list[list[int]]) -> list[int]:
    """Convert component sweep to simple list"""
def get_start_frame(component_sweep: list[list[int]]) -> list[int]:
    """Get first 3 (or less) atoms of the component sweep

    If the component sweep contains 3 or more atoms return these, otherwise
    return as many as possible
    """
def get_path(molgraph: nx.Graph, component_sweep: list[list[int]], root: int, node: int) -> list[int]:
    '''Find path from root to node, or the start frame if the path is shorter than 3.

    This routine provides the path of "root" to "node" used to setup coordinate systems.
    '''
def is_node_in_cycle(graph: nx.Graph, node: int) -> bool: ...
def split_bond_to_components(bond: tuple[int, int], graph: nx.Graph) -> set[int]:
    """Split graph at the specified bond and return the resulting connected components

    Parameters
    ----------
    bond: tuple[int, int]
        Bond where graph should be split

    graph: nx.Graph
        graph that should be split

    Returns
    -------
    components: generator of sets
        A generator of sets of nodes, one for each component of G.
    """
def split_bond_to_subgraphs(bond: tuple[int, int], graph: nx.Graph) -> nx.Graph:
    """Split graph at the specified bond and return the connected subgraphs"""
def bond_component(bond: tuple[int, int], element_in_component: int, graph: nx.Graph) -> set[int]:
    """Split graph at bond and return component that contains the specified element"""
def is_terminal_atom_bond(graph: nx.Graph, edge: tuple[int, int]) -> bool:
    """Is the edge an edge to a single atom at some end of the graph"""
def is_edge_in_cycle(edge: tuple[int, int], cycles: list[list[int]]) -> bool:
    """Are both atoms of an edge inside a single cycle"""
def is_single_bond_from_cycle(edge: tuple[int, int, dict[str, Any]], cycles: list[list[int]]) -> bool:
    """Is exactly one of the atoms in a cycle and the bond order below a threshold

    This is used because single bonds from aromatic cycles can be very close to
    """
def is_single_bond(bond_order: float):
    """Is the bond order below a set threshold"""
def set_rotatable_bond_grids(graph: nx.Graph, bond_type_grids: dict[str, np.ndarray], include_amide_bonds: bool) -> None:
    """Assign grid to rotatable bonds in the graph

    Bonds in cycles are currently excluded because rotations around those
    lead to a rotation of the *entire molecule* around this bond.

    Parameters:
    -----------
    bond_type_grids: dict[str, np.ndarray],
        Specify the grid (allowed dihedral values) for specific bond types:
        default_grid, chalcogen_grid, amino_grid, alkane_grid, t_butyl_grid,
        iso_propyl_grid, amide_grid
    include_amide_bonds: bool
        Include amide bonds as rotatable bonds
    """
