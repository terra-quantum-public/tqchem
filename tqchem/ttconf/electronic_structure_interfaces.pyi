import abc
import ase
import numpy as np
from _typeshed import Incomplete
from abc import ABC
from openmm.app import Simulation
from tqchem.chem import adjacency_differs as adjacency_differs, ase_to_rdkit as ase_to_rdkit, atoms_collided as atoms_collided, bondMatrix as bondMatrix
from tqchem.ttconf.crest_wrapper import launch_crest as launch_crest

DielectricConstant: Incomplete

class ElectronicStructureInterface(ABC, metaclass=abc.ABCMeta):
    """Interface to electronic structure codes

    Needs to implement a single point calculation and a gradient optimization

    Methods
    -------
    __call__:
        Calculate energy for ase.Atoms object (i.e. the molecule)
    _single_point_calculation:
        Calculate energy for the current structure of the ase.Atoms object
    _local_optimization:
        Perform gradient optimization for ase.Atoms object
        and return energy for the optimized structure
    _morse_energy:
        Calculate energy based on a simple morse potential
    """
    def __call__(self, molecule: ase.Atoms) -> float:
        """Calculate energy for molecule

        Parameters
        ----------
        molecule: ase.Atoms
            Molecule object containing atom types and positions

        Returns
        -------
        energy: float
            Single point energy in eV
        """

class TBliteInterface(ElectronicStructureInterface):
    """Interface/Wrapper for ASE/TBLite library to calculate energies

    https://tblite.readthedocs.io/en/latest/users/ase.html

    Attributes
    ----------
    calculator:
        TBLite calculator to perform single point calculations.

    Methods
    -------
    __call__:
        Calculate energy for ase.Atoms object (i.e. the molecule)
    _single_point_calculation:
        Calculate energy for the current structure of the ase.Atoms object
    _local_optimization:
        Perform gradient optimization for ase.Atoms object
        and return energy for the optimized structure
    """
    calculator: Incomplete
    local_optimization: Incomplete
    adjacency: Incomplete
    def __init__(self, method: str, ref_mol: ase.Atoms, charge: int = 0, accuracy: float = 1.0, solvent: str = None, local_optimization: bool = True) -> None:
        """Set up TBlite calculator

        Parameters
        ----------
        method: str
            GFN-xtb method to be used, Options: GFN1-xTB, GFN2-xTB
        ref_mol: ase.Atoms
            Reference molecule for creating a reference adjacency matrix
        charge: int, default=0
            total charge of the system
        accuracy: float, default=1.0
            accuracy used in the TBLite calculator
        solvent: str, default: None
            Add implicit solvation using specified solvent. Not implemented for Tblite
        local_optimization: bool
            perform short relaxation with gradient optimizer
        """

class OpenFFInterface(ElectronicStructureInterface):
    """Interface to OpenFF

    Different Force Fields are supported by changing the OpenMM simulation
    object from which this interface is constructed

    Implementing functions for single point calculations, optimization
    and extracting of coordinates from the openMM system

    Attributes
    ----------
    simulation: Simulation
        OpenMM simulation class for calculating energies and optimized structures
    adjacency: np.array
            Adjacency matrix of the molecular system
    """
    local_optimization: Incomplete
    simulation: Incomplete
    adjacency: Incomplete
    def __init__(self, simulation: Simulation, ref_mol: ase.Atoms, local_optimization: bool = True) -> None:
        """Initialize OpenFF interface

        Parameters
        ----------
        simulation: Simulation
            OpenMM Simulation object
        ref_mol: ase.Atoms
            Reference molecule for creating a reference adjacency matrix
        local_optimization: bool
            perform short relaxation with gradient optimizer
        """

def sageFF_interface(molecule: ase.Atoms, local_optimization: bool = True, charge: int = 0, solvent: str = None) -> OpenFFInterface:
    """Set up OpenMM simulation environment for the Sage OpenFF forcefield

    Parameters
    ----------
    molecule: ase.Atoms
        Molecule to create openMM simulation object from
    local_optimization: bool, default=True
        Enable gradient optimziation
    charge: int, default=0
        Formal charge of the molecule
    solvent: str, default=None
        Add implicit solvation using solvent specified by name
    """
def amberFF_interface(molecule: ase.Atoms, local_optimization: bool = True, charge: int = 0, solvent: str = None) -> OpenFFInterface:
    """Set up OpenMM Simulation environment for Amber Forcefield

    Parameters
    ----------
    molecule: ase.Atoms
        Molecule to create openMM simulation object
    local_optimization: bool, default=True
        Enable gradient optimziation
    charge: int, default=0
        Formal charge of the molecule
    solvent: str, default=None
        Add implicit solvation using solvent specified by name
    """

class CrestEnsembleInterface:
    optlevel: str
    threads: Incomplete
    charge: Incomplete
    solvent: Incomplete
    adjacency: Incomplete
    method: Incomplete
    def __init__(self, method: str, ref_mol: ase.Atoms, charge: int = 0, solvent: str = None, threads: int = 4) -> None:
        """Set up CREST calculator

        Parameters
        ----------
        method: str
            GFN-xtb method to be used, Options: GFN1-xTB, GFN2-xTB
        ref_mol: ase.Atoms
            Reference molecule for creating a reference adjacency matrix
        charge: int, default=0
            Formal charge of the system
        solvent: str, default=None
            Add implicit solvation using solvent specified by name
        threads: int, default=4
            Number of threads used
        """
    def __call__(self, molecules: list[ase.Atoms]) -> np.ndarray:
        """Prepare ensemble of molecules whithout collided atoms and optimize with CREST"""
