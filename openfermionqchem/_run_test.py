"""
Driver to perform unit test.
"""

import numpy as np
from openfermion.chem import MolecularData
from openfermion.linalg import get_sparse_operator, get_ground_state
from openfermion.ops import InteractionOperator
from openfermionqchem import run_qchem
from openfermionqchem import test_files

def print_test(description=None,bond_length=0.0,FCI=0.0,OF=0.0,DD=0.0):

    print("Molecule:      {}".format(description))
    print("Bond distance: {}".format(bond_length))
    print("FCI energy in Hartees for the reference:                              {:10.7f}".format(FCI))
    print("Ground state energy from OpenFermion and Fidelity w.r.t FCI:          {:10.7f}, {:4.2f}".format(OF,OF/FCI))
    print("Ground state energy by Direct Diagonalization and Fidelity w.r.t FCI: {:10.7f}, {:4.2f}".format(DD,DD/FCI))
    print('='*100)

def run_test():
    """
    This function checks code consistencies.

    Prints:
        H2/STO-3G FCI vs Direct Diagonalization energy comparisons.
        H2/3-21G  FCI vs Direct Diagonalization energy comparisons.
        H2/6-311G FCI vs Direct Diagonalization energy comparisons.
    """
    basis        = '3-21g'
    multiplicity = 1
    charge       = 0
    bond_length  = 0.75
    geometry     = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
    molecule     = MolecularData(geometry, basis, multiplicity, charge)
    path         = str(test_files.__path__._path[0])+'/'
    molecule1    = run_qchem(molecule,file_directory=path+'h2-sto3g/',output_name='test_qis')
    molecule2    = run_qchem(molecule,file_directory=path+'h2-321g/',output_name='test_qis')
    molecule3    = run_qchem(molecule,file_directory=path+'h2-6311g/',output_name='test_qis')

    h2_sto3g_fci = -1.1371170673457311
    h2_321g_fci  = -1.1478773770642743
    h2_6311g_fci = -1.1534912078971873

    H1 = InteractionOperator(molecule1.nuclear_repulsion,molecule1.one_body_integrals,0.25*molecule1.two_body_integrals)
    H2 = InteractionOperator(molecule2.nuclear_repulsion,molecule2.one_body_integrals,0.25*molecule2.two_body_integrals)
    H3 = InteractionOperator(molecule3.nuclear_repulsion,molecule3.one_body_integrals,0.25*molecule3.two_body_integrals)

    # Compute ground state energy through OpenFermion
    h2_sto3g_OF, _ = get_ground_state(get_sparse_operator(H1))
    h2_321g_OF, _  = get_ground_state(get_sparse_operator(H2))
    h2_6311g_OF, _ = get_ground_state(get_sparse_operator(H3))

    # Compute ground state energy through Direct Diagonalization of FCI matrix
    h2_sto3g_DD, _ = np.linalg.eigh(molecule1.direct_diagonalization)
    h2_321g_DD, _  = np.linalg.eigh(molecule2.direct_diagonalization)
    h2_6311g_DD, _ = np.linalg.eigh(molecule3.direct_diagonalization)
    h2_sto3g_DD    = min(h2_sto3g_DD) + molecule1.nuclear_repulsion
    h2_321g_DD     = min(h2_321g_DD)  + molecule2.nuclear_repulsion
    h2_6311g_DD    = min(h2_6311g_DD) + molecule3.nuclear_repulsion

    print('='*100)
    print_test(description='H2/STO-3G',bond_length=0.75,FCI=h2_sto3g_fci,OF=h2_sto3g_OF,DD=h2_sto3g_DD)
    print_test(description='H2/3-21G',bond_length=0.75,FCI=h2_321g_fci,OF=h2_321g_OF,DD=h2_321g_DD)
    print_test(description='H2/6-311G',bond_length=0.75,FCI=h2_6311g_fci,OF=h2_6311g_OF,DD=h2_6311g_DD)

if __name__ == "__main__":
    run_test()
