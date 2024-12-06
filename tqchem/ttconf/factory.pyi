from tqchem.chem import ase_from_rdkit as ase_from_rdkit, ase_to_rdkit as ase_to_rdkit, atoms_collided as atoms_collided, generate_rdkit_conformers as generate_rdkit_conformers
from tqchem.internal_coordinates import MolecularGrid as MolecularGrid, MolecularSystem as MolecularSystem
from tqchem.ttconf.electronic_structure_interfaces import CrestEnsembleInterface as CrestEnsembleInterface, ElectronicStructureInterface as ElectronicStructureInterface, TBliteInterface as TBliteInterface, amberFF_interface as amberFF_interface, sageFF_interface as sageFF_interface
from tqchem.ttconf.ttconf import EnergyFilter as EnergyFilter, EnsembleGridObjective as EnsembleGridObjective, MolecularGridObjective as MolecularGridObjective, TTopt_Optimizer as TTopt_Optimizer

def indices_from_rdkit_conformers(molgrid: MolecularGrid, n_conformers: int, indices: list[list[int]], charge: int = 0) -> list[list[int]]:
    """Calculates indices closest to the rdkit structures"""
def solve_atom_collisions(molgrid: MolecularGrid, initial_indices: list[int]) -> list[int]:
    """If atoms collide, search for a new set of indices without atom collisions"""
def initial_tt_indices(molecules: list[MolecularSystem], molgrid: MolecularGrid, rank: int, solve_collisions: bool, charge: int = 0) -> list[list[int]]:
    """Return indices to initialize the tensor train in TTOpt

    Creates initial indices from provided molecules, rdkit conformers and random selection.
    """
def molecularGridObjective(molgrid: MolecularGrid, method: str, filter_type: str, charge: int = 0, solvent: str = None, threads: int = 4, reference_energy: float = None, local_optimization: bool = True, ensemble_optimization: bool = False) -> MolecularGridObjective:
    """Construct a MolecularGridObjective"""
def ttconf_optimizer(molecules: MolecularSystem | list[MolecularSystem], method: str = 'gfn2-xtb', n_sweeps: int = 2, rank: int = 2, seed: int = 42, threads: int = 4, charge: int = 0, solvent: str = None, filter_type: str = 'energy', reference_energy: float = None, local_optimization: bool = True, initialize_indices: bool = True, ensemble_optimization: bool = False, solve_collisions: bool = True, **bond_grids) -> TTopt_Optimizer:
    '''Create Optimizer for the parameters provided

    Parameters
    ----------
    molecules: (MolecularSystem | list[MolecularSystem])
        Molecules provided as MolecularSystem
    method: str, default="gfn2-xtb"
        Electronic structure method used
        Options: "GFN1-xTB", "GFN2-xTB", "Sage-FF", "Amber-FF", case insensitive
    n_sweeps: int, default=2
        Number of sweeps through the tensor train
    rank: int, default=2
        Rank used in the tensor train
    seed: int, default=42
        Seed for random number selection
    threads: int, default=4
        Number of threads to use in the ensemble optimization
    charge: int, default=0
        Molecular charge
    solvent: str, default=None
        Solvent name
        Options: "water", "explicit water" and None (i.e. gas phase), case insensitive
    filter_type: str, default="energy"
        Type of function to filter the energy with
        Options:
        "energy": E_i
        "energy difference": E_i - E_0
        "boltzmann energy difference": -exp(-(E_i - E_0))
        where:
        E_0 is the energy of the current molecular system
        E_i is the energy of the system in iteration i
    reference_energy: float, default=None
        Reference energy for filter
    initialize_indices: bool, default=True
        Initialize the initial tensor in TTOpt with the given molecules and rdkit Conformers
    ensemble_optimization: bool, default=False
        Use ensemble optimization of CREST for gradient optimization
    solve_collisions: bool, default=True
        Check for atom collisions for the initial indices and resolves them
    bond_grids:
        Keywords specifying which grids to use when setting dihedrals:

        - include_amide_bonds: bool = False,
            Include amide bonds as rotatable bonds
        - add_rotatable_bonds: dict[tuple[int, int], np.array] = None,
            Adds rotatable bonds to the default ones in terms of a dict relating bonds
            to an array of allowed dihedral values
        - set_rotatable_bonds: dict[tuple[int, int], np.array] = None,
            Overwrites rotatable bonds with a dict relating bonds to an array of
            allowed dihedral values
        - amide_grid: np.array, default=[0.0, 180.0]
            Grid of allowed dihedral values for amide bonds
        - amino_grid: np.array, default=[0, 60, 120, 180, 240, 300]
            Grid of allowed dihedral values for amino bods
        - alkane_grid: np.array, default=[0, 60, 120, 180, 240, 300]
            Grid of allowed dihedral values for alkane bonds
        - t_butyl_grid: np.array, default=[0.0, 60.0]
            Grid of allowed dihedral values for t_butyl
        - iso_propyl_grid: np.array, default=[0.0, 60.0, 120.0]
            Grid of allowed dihedral values for iso_propyl
        - chalcogen_grid: np.array, default=[0, 60, 120, 180, 240, 300]
            Grid of allowed dihedral values for chalcogen bonds
        - default_grid: np.array, default=[0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
            Grid of allowed dihedral values for all other bonds
        - overwrite_grids: np.array, default=None
            Overwrite all bond specific grids with this grid

    Returns
    -------
    TToptOptimizer

    Raises
    ------
    NotImplementedError:
        if provided method or filter_type does not match the allowed ones
    TypeError:
        if molecule input does not match the allowed types
    '''
