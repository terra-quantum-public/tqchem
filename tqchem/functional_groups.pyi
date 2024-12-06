import networkx as nx
import numpy as np
from collections.abc import Iterable

def match_atomic_numbers(molecule_node: dict[str, int], subgraph_node: dict[str, int | Iterable[int]]):
    """Determine if a node in the molecule matches a node in the subgraph

    Required for checking isomorphisms of graphs.
    We check the atomic numbers and assume that an atomic number of -1
    means an unspecified residue.
    The subgraph node can also have a list of atomic numbers and the graphs
    match the molecule_node's atomic number is in the list.
    """
def match_bond_order_range(molecule_edge: dict[str, int], subgraph_edge: dict[str, int | Iterable[int]]):
    """Determine if the bond order of the molecule is inside a range defined for the subgraph

    Required for checking isomorphisms of graphs.
    """
def bonds_from_isomorphism_mappings(mappings: Iterable[dict[int, int]], subgraph_bonds: list[tuple[int, int]]) -> set[tuple[int, int]]:
    """Converts bonds defined on the subgraph to bonds on the graph

    Parameters
    ----------
    mapping: list[dict[int, int]]
        Mappings for nodes on the graph to nodes on the subgraph
    subgraph_bonds: list[tuple[int, int]]
        list of bonds defined for the subgraph nodes

    Returns
    -------
    bonds: set[tuple[int, int]]
        set of bonds defined on the graph
    """
def find_methyl_groups(graph: nx.Graph) -> set[tuple[int, int]]:
    """Get set of bonds matching an atom being bonded to a methyl group"""
def find_single_bonded_chalcogen(graph: nx.Graph) -> set[tuple[int, int]]:
    """Get set of bonds matching R-O-R or R-S-R bonds"""
def single_bonded_chalcogens_grid(graph: nx.Graph, grid: np.ndarray) -> dict[None]:
    """Returns Dictionary mapping chalcogen single bonds to their grid"""
def find_amide_groups(graph: nx.Graph) -> set[tuple[int, int]]:
    """Get set of bonds matching the C-N bond in amide groups"""
def amide_bonds_grid(graph: nx.Graph, grid: np.ndarray) -> dict[None]:
    """Returns Dictionary mapping amide bonds to their dihedral grid"""
def find_amino_groups(graph: nx.Graph) -> set[tuple[int, int]]:
    """Get set of bonds matching an atom being bonded to a amino group"""
def amino_group_grid(graph: nx.Graph, grid: np.ndarray) -> dict[None]:
    """Returns Dictionary mapping bonds to amino groups to their dihedral grid"""
def find_t_butyl_groups(graph: nx.Graph) -> set[tuple[int, int]]:
    """Get set of bonds matching the R-C bond in t-butyl groups"""
def t_butyl_groups_grid(graph: nx.Graph, grid: np.ndarray) -> dict[None]:
    """Returns Dictionary mapping bonds to t-butyl groups to their dihedral grid"""
def find_iso_propyl_groups(graph: nx.Graph) -> set[dict[int, int]]:
    """Get set of bonds matching the R-C bond in isopropyl groups"""
def iso_propyl_groups_grid(graph: nx.Graph, grid: np.ndarray) -> dict[None]:
    """Returns Dictionary mapping bonds to iso-propyl groups to their dihedral grid"""
def find_alkane_bonds(graph: nx.Graph) -> set[tuple[int, int]]:
    """Get set of bonds matching the simple C-C bond"""
def simple_alkanes_grid(graph: nx.Graph, grid: np.ndarray) -> dict[None]:
    """Returns Dictionary mapping alkane bonds to iso-propyl groups to their dihedral grid"""
