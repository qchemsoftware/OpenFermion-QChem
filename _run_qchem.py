# OpenFermion plugin to interface with Q-CHEM.
# License///

'''Function to initialize Q-CHEM MolecularData by parsing Q-CHEM outputs'''

from openfermion import MolecularData
from openfermionqchem import QchemMolecularData

def run_qchem(molecule,
              run_scf=True,
              run_mp2=False,
              run_ccsd=False,
              verbose=False,
              file_directory =None,
              input_name=None,
              output_name=None):
    """
    This function creates Q-CHEM molecular object and saves Q-CHEM outputs.

    Args:
        * OpenFermion-Q/CHEM v0.0 only needs two arguments
        file_directory: path of Q-CHEM output files.
        output_name: output file name to parse.

        * Upcoming features: Run Q-CHEM through OpenFermion/Q-CHEM. Then,
        molecule: An instance of the MolecularData or QchemMolecularData class.
        run_scf: Optional boolean to run SCF calculation.
        run_mp2: Optional boolean to run MP2 calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        verbose: Boolean whether to print calculation results to screen.
        input_name: provide Q-CHEM input file to run if available.

    Returns:
        molecule: The updated MolecularData object.
    """

    # Return updated MolecularData object
    return QchemMolecularData(qchem_file=output_name, file_directory=file_directory)

if __name__ == "__main__":
    # Set molecule parameters.
    basis        = 'sto-3g'
    multiplicity = 1
    charge       = 0
    bond_length  = 0.8
    # Set calculation parameters.
    run_scf  = 0
    run_mp2  = 0
    run_ccsd = 0
    # need to change openfermion MolecularData for QCHEM
    # QCHEM, no needs to set up initial geometry, basis, charge, ...
    geometry       = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
    molecule       = MolecularData(geometry, basis, multiplicity, charge)
    file_directory = '/Users/yongbin/Desktop/venv/openfermion-qchem/molecules/h2/sto-3g/'+str(bond_length)+'/'
    molecule       = run_qchem(molecule,file_directory=file_directory,output_name='test_qis')
    print(molecule.hf_energy)