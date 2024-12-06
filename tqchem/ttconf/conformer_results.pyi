import pandas as pd
from _typeshed import Incomplete
from pathlib import Path
from tqchem.chem import ase_to_xyz_content as ase_to_xyz_content
from tqchem.internal_coordinates import MolecularGrid as MolecularGrid, MolecularSystem as MolecularSystem

def write_xyz_list(filepath: Path | str, xyzs: list[str]) -> None: ...
def escape_new_lines(string: str) -> str:
    """Escape special characters in string to make it printable as is"""

class ConformerResults:
    """Class to store and process results of a conformer search

    Attributes
    ----------
    dataframe: pd.DataFrame
        Dataframe containing the objective, energy, gridpoints, relaxed & aligned xyz data,
        and rmsd values indicating the fit of the xyz to the minimum energy xyz
    molgrid: MolecularGrid
        MolecularGrid to evaluate the dihedral gridpoints to obtain the unrelaxed
        molecules sampled during the optimization procedure.

    Methods
    -------
    __str__: str
        Return data frame string representation
    minimum_energy_xyz: str, float
        Return xyz string (after gradient optimization) and energy of molecule with lowest energy
    minimum_objective_xyz: str, float
        Return xyz string (after gradient optimization) and objective of molecule with lowest objective
    minimum_energy_unrelaxed_conformer: MolecularSystem, float
        Return Molecular system (before gradient optimization) with lowest energy
    minimum_objective_unrelaxed_conformer: MolecularSystem, float
        Return Molecular system (before gradient optimization) with lowest objective value
    sampled_unrelaxed_conformers: list[MolecularSystem]
        Return list of all sampled molecular systems
    write_sampled_unrelaxed_conformers: None
        Write all sampled conformers to one xyz file with energies in the comment
    minimum_energy_ensemble: list[MolecularSystem]
        Molecular systems within an energy margin around the lowest energy conformer
    write_minimum_energy_ensemble: None
        Write all conformers within an energy margin around the lowest energy conformer
    """
    dataframe: Incomplete
    molgrid: Incomplete
    def __init__(self, dataframe: pd.DataFrame, molecular_grid: MolecularGrid) -> None: ...
    def minimum_energy_xyz(self) -> str:
        """Return xyz content of gradient optimized geometry with lowest energy"""
    def minimum_objective_xyz(self) -> str:
        """Return xyz content of gradient optimized geometry with lowest objective"""
    def minimum_energy_unrelaxed_conformer(self) -> MolecularSystem:
        """Return Molecular system with lowest energy before gradient optimization"""
    def minimum_objective_unrelaxed_conformer(self) -> MolecularSystem:
        """Return Molecular system with lowest objective before gradient optimization"""
    def sampled_unrelaxed_conformers(self) -> list[MolecularSystem]:
        """Return list of all sampled unrelaxed molecular systems"""
    def write_sampled_unrelaxed_conformers(self, filepath: Path | str) -> None:
        """Write all sampled unrelaxed conformers to one xyz file

        Energy and objective (of the relaxed structures) are listed in the comment line"""
    def minimum_energy_ensemble(self, energy_difference: float = 0.2, rmsd_cutoff: float = 0.1) -> list[str]:
        """xyz strings within an energy and rmsd margin around the lowest energy conformer

        The xyz strings contain the gradient optimized structures and their
        energy and objective value are listed in the comment line
        """
    def write_minimum_energy_ensemble(self, filepath: Path | str, energy_difference: float = 0.2, rmsd_cutoff: float = 0.1) -> None:
        """Write xyz within an energy and rmsd margin around the lowest energy conformer

        The xyz strings contain the gradient optimized structures and their
        energy and objective value are listed in the comment line
        """
    def dataframe_to_csv(self, *args):
        """Convert xyz column to raw string and write dataframe to csv

        Wrapper for DataFrame.to_csv which accepts the same keyword arguments
        """
