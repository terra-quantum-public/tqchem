from argparse import Namespace
from tqchem.chem import ase_from_rdkit as ase_from_rdkit, generate_rdkit_conformers as generate_rdkit_conformers, rdkit_from_smiles as rdkit_from_smiles
from tqchem.functional_groups import match_atomic_numbers as match_atomic_numbers
from tqchem.molgraph import MolecularSystem as MolecularSystem, molecularSystem as molecularSystem
from tqchem.ttconf.conformer_results import ConformerResults as ConformerResults
from tqchem.ttconf.factory import ttconf_optimizer as ttconf_optimizer

def parse_commandline(provided_arguments: list[str] = None) -> Namespace:
    """Parse arguments and provide program description

    Parameters
    ----------
    provided_arguments: list[str], default=None
        If no arguments are provided sys.argv is used. Used for testing.
    """
def generation_parser(subparsers) -> None: ...
def ttconf_parser(subparsers) -> None: ...
def smiles_to_filename(smiles: str) -> str:
    """Replace smile string characters with - or _"""
def molecular_systems(arguments: Namespace) -> list[MolecularSystem]:
    """Create Molecular Systems from the arguments (Smiles and files)"""
def conformer_optimizer(arguments: Namespace) -> ConformerResults:
    """Create optimizer from commandline arguments"""
def write_results(results: ConformerResults, arguments: Namespace) -> None:
    """Write output files from the TTconf run"""
def run_conformer_search(arguments: Namespace):
    """Run conformer search from commandline arguments"""
def generate_conformers(arguments: Namespace):
    """Generate conformers from smiles string or xyz file"""
def main() -> None:
    """Parse arguments and run the function set for the parser"""
