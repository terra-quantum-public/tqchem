import py3Dmol
from ase import Atoms as Atoms
from pathlib import Path
from rdkit import Chem as Chem
from tqchem.chem import xyz_contents as xyz_contents
from tqchem.internal_coordinates import MolecularGrid as MolecularGrid
from tqchem.molgraph import MolecularSystem as MolecularSystem

def plot_3D_molecule(molecule: Path | str | Atoms | Chem.Mol, width: int = 300, height: int = 300) -> py3Dmol.view:
    """Plot molecule using py3Dmol"""
def display_trajectory(xyzs: list[str], time: int = 30, repetitions: int = 1, width: int = 500, height: int = 500) -> py3Dmol.view:
    """Display trajectory using py3Dmol

    Parameters
    ----------
    xyzs: list[str]
        list of xyz contents
    time: int
        time for the animation to complete; larger number means slower animation
    repetitions: int
        animation is repeated as many times; 0 indicated infinite loop
    width: int
        width of the animation window
    height: int
        height of the animation window
    """
def coordinate_screener(mol: MolecularSystem, zoomTo: bool = False, shift_by_reference: bool = True) -> None:
    """Display changes along coordinates using ipython sliders

    ipython widget displaying the molecule and allowing to change
    individual coordinates using sliders

    Parameters
    ----------
    mol: MolecularSystem
        tqchem representation of a molecule
    zoomTo: bool, optional, default False
        zoom to molecule after each update
    """
def coordinate_screener_from_grid(molecule_grid: MolecularGrid, zoomTo: bool = False) -> None:
    """Display changes along coordinates using ipython sliders

    ipython widget displaying the molecule and allowing to change
    individual coordinates using sliders

    Parameters
    ----------
    molecule_grid: MolecularSystem
        Representation of a molecule and discretization of the flexible coordinates
    zoomTo: bool, optional, default False
        zoom to molecule after each update
    """
