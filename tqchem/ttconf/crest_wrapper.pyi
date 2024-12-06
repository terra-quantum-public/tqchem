import ase
import numpy as np
from pathlib import Path

def launch_crest(input_structures: list[ase.Atoms], method: str = 'GFN2-xTB', threads: int = 4, optlevel: str = 'normal', charge: int = 0, solvent: str = None) -> tuple[list[ase.Atoms], np.ndarray]:
    """Runs a CREST optimization for a given ensemble of conformers.

    Returns the ensemble and the energies after optimization.
    """
def extract_energy_from_comment(comment):
    """From the comment of the xyz file reads the conformer energy value

    Note: Only extracts the first float and only if it is surrounded
          by whitespace or nothing (see split)
    """
def read_crest_ensemble_xyz(file_: Path) -> tuple[list[ase.Atoms], np.ndarray]:
    """Reads an xyz file with conformer structures and energies as written by CREST"""
