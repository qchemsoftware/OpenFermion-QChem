# Created by Yongbin Kim, Krylov's Group, University of Southern California
# Add License.

"""Class to access Q-CHEM calculation results"""

from functools import reduce
from openfermion.chem import MolecularData
from openfermionqchem._qchem_parser import *
from openfermionqchem._qchem_functiontools import *

class QchemMolecularData(MolecularData):

    """A derived class from openfermion.chem.MolecularData.
    This class is to store electronic structure data of Q-CHEM at a fixed basis set at a fixed geometry
    This class also provides an interface to access the Q-CHEM input or output files (if given).  

    Attributes:
        as same as openfermion.chem.MolecularData class.
    """

    # v0.0
    def __init__(self, geometry=None, basis=None, multiplicity=None,
                charge=0, qchem_file="", file_directory=""):
        self._qchem_file = qchem_file
        self._file_directory = file_directory
        # Upcoming version will combine functions to read Q-CHEM output files.
        self.read_out()
        self.read_fchk()
        self.read_extra()
        self.qchem_basic()

    def read_out(self):
        """ Save molecular energies e.g) HF, MP2, and CCSD """
        with open(self._file_directory+self._qchem_file+'.inp.out','rb') as ofile:
            self._qchem_output = ofile.read().decode()
    def read_fchk(self):
        """ Save general molecular data e.g) charge, multiplicity, orbital coefficients etc. """
        with open(self._file_directory+self._qchem_file+'.inp.fchk','rb') as ffile:
            self._qchem_fchk = ffile.read().decode()
    def read_extra(self):
        """
        Two electron integrals 
        Molecular orbital energies
        CCSD single and double amplitudes
        compatibility: Q-CHEM 5.4+
        """
        with open(self._file_directory+'two_body_int_for_qis.dat','rb') as eri:
            self._qchem_eri = eri.read().decode()
        with open(self._file_directory+'mo_ene_for_qis.dat','rb') as orb_map:
            self._qchem_omap = orb_map.read().decode()
        with open(self._file_directory+'cc_t1_for_qis.dat','rb') as cc_t1:
            self._qchem_cc_t1 = cc_t1.read().decode()
        with open(self._file_directory+'cc_t2_for_qis.dat','rb') as cc_t2:
            self._qchem_cc_t2 = cc_t2.read().decode()

    def qchem_basic(self):
        self._qchem_data = get_qchem_basic_data(self._qchem_fchk)

    @property
    def geometry(self):
        """
        A list of tuples giving the coordinates of each atom in angstrom.
        e.g) [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
        """
        self._geometry = self._qchem_data['geometry']
        return self._geometry

    @property
    def basis(self):
        """ A string giving the basis set. e.g) sto-3g. """
        self._basis = self._qchem_data['basis']
        return self._basis

    @property
    def charge(self):
        """ An integer giving total molecular charge. """
        self._charge = self._qchem_data['charge']
        return self._charge

    @property
    def multiplicity(self):
        """ An integer giving the spin multiplicity. """
        self._multiplicity = self._qchem_data['multiplicity']
        return self._multiplicity

    @property
    def n_atoms(self):
        """ An integer giving total number of atoms in the molecule. """
        self._n_atoms = self._qchem_data['n_atoms']
        return self._n_atoms

    @property
    def n_electrons(self):
        """ An integer giving total number of electrons in the molecule. """
        self._n_electrons = self._qchem_data['n_electrons']
        return self._n_electrons

    @property
    def atoms(self):
        """ A list of the atoms in the molecule. """
        self._atoms = self._qchem_data['atoms']
        return self._atoms

    @property
    def protons(self):
        """ A list of the nuclear charges in the molecule. """
        self._protons = self._qchem_data['protons']
        return self._protons

    @property
    def n_orbitals(self):
        """ An integer giving total number of spatial orbitals. """
        self._n_orbitals = self._qchem_data['n_orbitals']

        return self._n_orbitals
    @property
    def n_qubits(self):
        """ An integer giving the number of qubits that would be needed."""
        self._n_qubits = self._qchem_data['n_orbitals'] * 2
        return self._n_qubits

    @property
    def mo_maps(self):
        """ A dictionary giving the map between the orbital energy and corresponding index. """
        self._mo_maps = get_correct_orb_order(self._qchem_omap)
        return self._mo_maps

    @property
    def hf_energy(self):
        """ Hartree-Fock energy. """
        self._hf_energy = get_hf_energy(self._qchem_output)
        return self._hf_energy

    @property
    def nuclear_repulsion(self):
        """ Nuclear repulsion energy. """
        self._nuclear_repulsion = get_nuclear_repulsion(self._qchem_output)
        return self._nuclear_repulsion

    @property
    def canonical_orbitals(self):
        """ Numpy array giving the Hartree-Fock canonical orbital coefficients (represented on AO basis). """
        self._canonical_orbitals = get_mo_coeff(self._qchem_fchk, self._qchem_data['n_orbitals'])
        return self._canonical_orbitals

    @property
    def orbital_energies(self):
        """ Numpy array givig the canonical orbital energies. """
        self._orbital_energies = get_mo_energy(self._qchem_fchk)
        return self._orbital_energies

    @property
    def overlap_integrals(self):
        """ Numpy array giving the overlap integrals (represented on AO basis). """
        self._overlap_integrals = get_overlap(self._qchem_fchk, self._qchem_data['n_orbitals'])
        return self._overlap_integrals

    @property
    def h_core(self):
        """ Numpy array giving the core Hamiltonian """
        self._h_core = get_h_core(self._qchem_fchk, self._qchem_data['n_orbitals'])
        return self._h_core

    @property
    def one_body_integrals(self):
        """ 
        Numpy array giving the one-electron integrals (spin orbital representation). 
        hpq = MO_dagger * H_core * MO
        """
        mo      = self.canonical_orbitals
        h_core  = self.h_core
        spatial = reduce(np.dot, (mo.T, h_core, mo))
        self._one_body_integrals = get_one_body_integrals(spatial,2*self._qchem_data['n_orbitals'])
        return self._one_body_integrals

    @property
    def two_body_integrals(self):
        """ Numpy array giving the anti-symmetrized two-electron integrals (spin orbital representation). """
        self._two_body_integrals = get_two_body_integrals(self._qchem_eri,self.mo_maps)
        return self._two_body_integrals

    @property
    def mp2_energy(self):
        """ MP2 energy. """
        self._mp2_energy = get_mp2_energy(self._qchem_output)
        return self._mp2_energy

    @property
    def ccsd_energy(self):
        """ CCSD energy. """
        self._ccsd_energy = get_ccsd_energy(self._qchem_output)
        return self._ccsd_energy

    @property
    def ccsd_single_amps(self):
        """ Numpy array giving CCSD singles. T1_ia """
        self._ccsd_single_amps = get_ccsd_t1_amps(self._qchem_cc_t1,self.mo_maps)
        return self._ccsd_single_amps

    @property
    def ccsd_double_amps(self):
        """ Numpy array giving CCSD doubles. T2_ijab """
        self._ccsd_double_amps = get_ccsd_t2_amps(self._qchem_cc_t2,self.mo_maps)
        return self._ccsd_double_amps

    @property
    def direct_diagonalization(self):
        """ Numpy array giving FCI matrix. """
        n_eles,n_spin_orbs = self.n_electrons, 2*self.n_orbitals
        one_body, two_body = self.one_body_integrals, self.two_body_integrals
        self._fci_matrix   = build_fci_matrix(n_eles,n_spin_orbs,one_body,two_body)
        return self._fci_matrix