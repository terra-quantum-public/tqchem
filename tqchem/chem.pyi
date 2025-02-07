import ase
import numpy as np
import rdkit.Chem as Chem
from pathlib import Path

def is_smiles_string(string: str) -> bool:
    """Check if string could be a smiles string

    SMILES strings may contain:

        - Letters of the element symbols:
        - b, c, n, o, p, s for aromatic atoms
        - Numbers for rings
        - . - = # $ : / \\ for Bonds
        - [] + - for Atom and charge specification
        - () for branching
        - @ for stereocenters
    """
def is_multiline_string(string: str) -> bool: ...
def is_existing_file(string: str) -> bool:
    """Tries to check if a string could be an existing file

    An OSError is for example raised when the filename is too long.
    As the error codes are OS dependent we can only assume that
    an OSError means that we do not have a valid file.
    """
def ase_from_rdkit(molecule: Chem.Mol, conformer: Chem.Conformer | int = -1) -> ase.Atoms:
    '''Convert rdkit.Mol object to ase.Atoms object

    Note the Molecule needs to be "embedded" or have 3D coordinates'''
def ase_to_rdkit(molecule: ase.Atoms, charge: int = 0) -> Chem.Mol:
    """Convert ase.Atoms object to rdkit.Mol object"""
def ase_from_path(path: Path, format_: str = None) -> ase.Atoms: ...
def ase_from_molecule_file_content(string: str, format_: str = None) -> ase.Atoms:
    """Convert a string containing a molecule file's content to ase.Atoms

    The format could be .xyz, .mol, .sdf, but if none is specified .xyz is assumed
    """
def rdkit_from_smiles(smiles: str) -> Chem.Mol:
    """Convert a smiles string to a rdkit Molecule"""
def ase_from_smiles(smiles: str) -> ase.Atoms:
    """Convert a smiles string to ase.Atoms object using rdkit"""
def ase_to_xyz_content(molecule: ase.Atoms, comment: str = '') -> str:
    """Return string containing the contents of an xyz file"""
def ase_molecule(input_: Path | str | ase.Atoms | Chem.Mol, format_: str = None) -> ase.Atoms:
    """Smart factory for ase.Atoms object

    Parameters
    ----------
    input_: (Path | str | ase.Atoms | Chem.Mol)
        Molecule provided as: Path object, ase.Atoms, rdkit.Chem.Mol, or string
        containing the contents of an xyz file, a Smiles string or a filename
    format_: str, optional
        Optional format string used when a filename or Path or the content
        of a molecule file are provided

    Returns
    -------
    molecule: ase.Atoms
    """
def xyz_contents(input_: Path | str | ase.Atoms | Chem.Mol, comment: str = '') -> str:
    """Given a molecule specification return the contents of the corresponding xyz file"""
def align_molecules(molecules: list[Chem.Mol], reference: Chem.Mol, to_previous: bool = False) -> list[float]:
    """Aligns molecules to a reference molecule.

    Parameters
    ----------
    molecules: list[rdkit.Chem.Mol]
        list of molecules to align
    reference: rdkit.Chem.Mol
        Molecule to align to
    to_previous: bool, optional
        Align to previous molecule
    Returns
    -------
    rmsds: list[float]
        list of root mean square deviations from the reference
    """
def align_xyz_strings(xyzs: list[str], reference: str, to_previous: bool = False) -> tuple[list[str], list[float]]:
    """Aligns molecules to a reference molecule.

    Parameters
    ----------
    xyzs: list[str]
        list of strings containing xyz data
    reference: str
        String containing xyz data for the reference to align against
    to_previous: bool, optional
        Align to previous molecule

    Returns
    -------
    aligned_xyzs: list[str]
        aligned xyz strings
    rmsds: list[float]
        list of root mean square deviations from the reference
    """
def rdkit_molecules_to_xyz(molecules: list[Chem.Mol], path: str | Path) -> None:
    """Write list of rdkit molecules to xyz file

    Parameters
    ----------
    molecules: list[rdkit.Chem.Mol]
        list of molecules to write to file

    path: str | Path
        Path to the xyz file
    """
def adjacencyMatrix(mol: ase.Atoms) -> np.ndarray:
    """Determine adjacency matrix based on atomic distances

    Parameters
    ----------
    mol: ase.atoms.Atoms
        Object representing the molecule

    Returns
    -------
    mat: np.array
        Adjacency matrix for each pair of atoms
    """
def bondMatrix(mol: ase.Atoms, radius_cutoff: float = ...):
    """Determine bond matrix based on atomic distances

    Parameters
    ----------
    mol: ase.atoms.Atoms
        Object representing the molecule
    radius_cutoff: float
        Cutoff for the bond length

    Returns
    -------
    mat: np.array
        matrix which is 1 for bond between atoms and 0 otherwise
    """
def determineBonds(mol: ase.Atoms) -> np.ndarray:
    """Determine the bonds based on the adjacency matrix"""
def element_color(atomic_number: int) -> tuple[float, float, float]:
    """Return element color for a given atomic number"""
def wiberg_bond_orders(molecule: ase.Atoms) -> np.ndarray: ...
def atoms_collided(molecule: ase.Atoms, cutoff: float = 0.5) -> bool:
    """True if two atoms are closer than a set cutoff."""
def adjacency_differs(adjacency1: np.ndarray, adjacency2: np.ndarray) -> bool:
    """Return if 2 adjacency matrices differ"""
def generate_rdkit_conformers(molecule: Chem.Mol, n_conformers: int, threshold: float = 0.125, keep_current: bool = False, threads: int = 4, seed: int = 42) -> list[int]:
    """Generate conformers from rdkit and returns ids of generated ones

    The conformers are stored inside the molecule and a list of conformer
    ids for the conformers generated during this function is returned

    Arguments
    ---------
    molecule: Chem.Mol
        rdkit molecule object
    n_conformers: int,
        Number of conformers to generate. If threshold > 0 we generate
        10 times more to get at least the specified number
    threshold: float, default=0.125
        RMSD threshold below which we consider conformers identical
    keep_current: bool, default=False
        Keep conformers currently stored in molecule
    threads: int, default=4
        Number of threads to use in the conformer generation
    seed: int, default=42
        Random seed

    Returns
    -------
    list of conformer ids for the conformers generated
    """
def calculate_rotational_constants(input_: Path | str | ase.Atoms | Chem.Mol):
    """Calculate rotational constants with respect to the principle axes"""
def calculate_rmsd(input1: Path | str | ase.Atoms | Chem.Mol, input2: Path | str | ase.Atoms | Chem.Mol, align: bool = True, center: bool = True):
    """Calculate RMSD based on Kabsch algorithm for two sets of coordinates"""
