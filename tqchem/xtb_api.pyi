from _typeshed import Incomplete
from tqchem.molgraph import MolecularSystem as MolecularSystem, molecule_from_file as molecule_from_file

XTB_METHODS: Incomplete

def xtb_singlepoint(molecule: MolecularSystem, method: str = 'gfn1-xtb', charge: int = 0):
    '''
    Performs a single-point energy calculation using xTB.

    Parameters
    ----------
    molecule : MolecularSystem
        The molecular structure to perform the calculation on.
    method : str, optional
        The xTB method to use. Defaults to "gfn1-xtb".
    charge : int, optional
        The charge of the molecule. Defaults to 0.

    Returns
    -------
    energy : float or None
        The total energy of the molecule in Eh, or None if the energy could not be extracted.
    '''
def xtb_optimize(molecule: MolecularSystem, method: str = 'gfn1-xtb', charge: int = 0):
    '''
    Optimizes the molecular structure using xTB.

    Parameters
    ----------
    molecule : MolecularSystem
        The molecular structure to perform the optimization on.
    method : str, optional
        The xTB method to use. Defaults to "gfn1-xtb".
    charge : int, optional
        The charge of the molecule. Defaults to 0.

    Returns
    -------
    molecule_opt : MolecularSystem
        The optimized molecular structure.
    energy : float or None
        The total energy of the molecule in Eh, or None if the energy could not be extracted.
    '''
def read_xtb_vibspectrum(filename: str):
    """
    Reads the vibrational spectrum from an xTB output file.

    Parameters
    ----------
    filename : str
        The path to the file containing the vibrational spectrum data.

    Returns
    -------
    frequencies : list
        A list of vibrational frequencies in wave numbers.
    intensities : list
        A list of IR intensities corresponding to the frequencies.
    """
def xtb_ir_spectrum(molecule: MolecularSystem, method: str = 'gfn-xtb', charge: int = 0):
    '''
    Calculates the IR spectrum of a molecule using xTB.

    Parameters
    ----------
    molecule : MolecularSystem
        The molecular structure to calculate the IR spectrum for.
    method : str, optional
        The xTB method to use for the calculation. Defaults to "gfn-xtb".
    charge : int, optional
        The charge of the molecule. Defaults to 0.

    Returns
    -------
    frequencies : list
        A list of vibrational frequencies in wave numbers.
    intensities : list
        A list of IR intensities corresponding to the frequencies.
    '''
