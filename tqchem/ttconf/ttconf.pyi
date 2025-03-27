from _typeshed import Incomplete
from abc import ABC
from tqchem.chem import align_xyz_strings as align_xyz_strings, ase_to_xyz_content as ase_to_xyz_content
from tqchem.internal_coordinates import MolecularGrid as MolecularGrid
from tqchem.ttconf.conformer_results import ConformerResults as ConformerResults
from tqchem.ttconf.electronic_structure_interfaces import CrestEnsembleInterface as CrestEnsembleInterface, ElectronicStructureInterface as ElectronicStructureInterface
from tqchem.ttopt.ttopt import TTOpt as TTOpt
from typing import Any, Callable

class EnergyFilter:
    """Class used to filter energies.

    Wrapper storing the reference energy to calculate energy differences
    and apply some filter function (e.g. an exponential for a Boltzmann filter) to it.
    """
    reference_energy: Incomplete
    filter_: Incomplete
    def __init__(self, filter_: Callable[[float], float], reference_energy: float = 0.0) -> None: ...
    def __call__(self, x: float) -> float: ...

class TTConfObjective(ABC):
    """Objective function to be optimized for conformer search"""

class MolecularGridObjective(TTConfObjective):
    """Objective function to be optimized for conformer search

    Update conformer dihedrals, calculate energy and apply filter.
    This object will be directly called in the Optimizer, which
    optimizes based on the filtered value returned.

    Attributes
    ----------
    molgrid: MolecularGrid
        Wrapper for the molecular graph and the grid which is able to
        update internal coordinates and return ase.Atoms objects
    energy_function: ElectronicStructureInterface
        Interface to the electronic structure code, receiving an ase.Atoms
        object and returns the energy
    filter_: Callable[[float], float]
        Function transforming the energy, e.g. making sure it's a positive
        number for example a Boltzmann filter exp(-energy) can be used
    """
    molgrid: Incomplete
    energy_function: Incomplete
    filter: Incomplete
    data: Incomplete
    def __init__(self, molgrid: MolecularGrid, energy_function: ElectronicStructureInterface, filter_: EnergyFilter) -> None:
        """Setup objects"""
    def __call__(self, dihedrals: list[int]) -> float:
        """Update conformer with new dihedrals, calculate and filter energy"""

class EnsembleGridObjective(TTConfObjective):
    """Objective function to be optimized for conformer search

    Update conformer dihedrals, create ensemble of new conformers, calculate its energies and apply filter.
    This object will be directly called in the Optimizer, which
    optimizes based on the filtered value returned.

    Unlike MolecularGridObjective, performs optimization of an ensemble of molecules
    using CREST, which significantly speeds up calculations.

    Attributes
    ----------
    molgrid: MolecularGrid
        Wrapper for the molecular graph and the grid which is able to
        update internal coordinates and return ase.Atoms objects
    energy_function: CrestEnsembleInterface
        Interface to the CREST, receiving an array with ase.Atoms
        objects and returns the energies
    filter_: Callable[[float], float]
        Function transforming the energy, e.g. making sure it's a positive
        number for example a Boltzmann filter exp(-energy) can be used
    """
    energy_function: Incomplete
    molgrid: Incomplete
    filter: Incomplete
    data: Incomplete
    def __init__(self, molgrid: MolecularGrid, energy_function: CrestEnsembleInterface, filter_: EnergyFilter) -> None:
        """Setup objects"""
    def __call__(self, sets_of_dihedrals: list[list[int]]) -> float:
        """Update conformer with new dihedrals, calculate and filter energy"""

class TTopt_Optimizer:
    """Optimizer using TTopt

    The implementation is based on TTconf
    (https://github.com/terra-quantum-io/TTconf)

    Attributes
    ----------
    objective: MolecularGridObjective
        Objective to be optimized, callable returning a float
    ttopt_parameters: dict[str, Any]
        keywordarguments passed to TTopt defining, number of sweeps, rank,
        seed, and number of threads
    """
    objective: Incomplete
    ttopt_parameters: Incomplete
    def __init__(self, objective: TTConfObjective, ttopt_parameters: dict[str, Any] = None) -> None: ...
    def optimize(self) -> ConformerResults:
        """Start the optimization of the objective function using TTopt"""
