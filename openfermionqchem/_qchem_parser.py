"""
Driver to parse Q-CHEM output, checkpoint, and integral files
compatibility: Q-CHEM 5.4+
Add License.
"""

import re
import numpy as np
import copy
np.set_printoptions(precision=10)

def periodic_table(proton):
    # The Periodic Table as a python list and dictionary.
    periodic_table = [  #
        '?', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
        'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
        'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
        'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
        'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
        'Lr']

    periodic_hash_table = {}
    for atomic_number, atom in enumerate(periodic_table):
        periodic_hash_table[atomic_number] = atom

    return periodic_hash_table.get(proton)

def parse(output, from_where=0, letter_type=None):
    """
    Parse Q-CHEM output file

    Args:
        output: Q-CHEM output or checkpoint or integrals.
        from_where: starting position of specific letter. Default=0
        letter_type: A string of letter to find

    Returns:
        An integer giving the position of given letter
    """
    output = output[from_where:]
    prev   = 0
    for letter in re.finditer(letter_type, output):
        if letter.start() > prev + 1:
            return letter.start() + from_where
        prev = letter.end()
    print('raise error')

def fill_full_matrix(lower_triangle, nbas):
    """
    Convert lower triangular matrix into full nbas x nbas matrix

    Args:
        lower_triangle: Numpy array of lower triangular matrix.
        nbas: An integer of the number of basis functions.

    Returns:
        Numpy array of nbas x nbas matrix.
    """
    mat = np.zeros((nbas,nbas))
    k   = 0
    for i in range(nbas):
        for j in range(i+1):
            mat[i,j] = float(lower_triangle[j+k])
            mat[j,i] = float(lower_triangle[j+k])
        k += i+1
    return mat

def get_qchem_basic_data(fchk):
    """
    Parse Q-CHEM checkpoint file and save molecular information into Python dictionary

    Args:
        fchk: Q-CHEM checkpoint file.
    
    Returns:
        Python dictionary containing parsed data.
    """
    data = {}
    bohr_to_angs = 0.529177249

    num = fchk.find('Number of atoms')
    data['basis'] = str(fchk[:num].split()[-1]).lower()
    data['n_atoms'] = int(fchk[num:num+100].split()[4])

    num = fchk.find('Charge')
    data['charge'] = int(fchk[num:num+100].split()[2])

    num = fchk.find('Multiplicity')
    data['multiplicity'] = int(fchk[num:num+100].split()[2])

    num = fchk.find('Number of alpha electrons')
    data['n_alpha'] = int(fchk[num:num+100].split()[5])
    num = fchk.find('Number of beta electrons')
    data['n_beta']  = int(fchk[num:num+100].split()[5])

    num = fchk.find('Number of electrons')
    data['n_electrons'] = int(fchk[num:num+100].split()[4])

    num1 = fchk.find('Current cartesian coordinates')
    num2 = fchk.find('Nuclear charges')
    num3 = fchk.find('Number of basis functions')    
    coord, atoms, temp = [], [], []
    geometry, protons  = [], []
    for xyz in range(6, len(fchk[num1:num2].split())):
        temp.append(float(fchk[num1:num2].split()[xyz]))
        if xyz % 3 == 2:
            coord.append(temp)
            temp = []
    for atm in range(5, len(fchk[num2:num3].split())):
        protons.append(int(float(fchk[num2:num3].split()[atm])))
        atom = periodic_table(int(float(fchk[num2:num3].split()[atm])))
        atoms.append(atom)
        # decimal point needed to be modified
        geometry.append((atom,tuple(round(float(xyz)*bohr_to_angs,5) for xyz in coord[atm-5])))
    data['geometry'] = geometry
    data['atoms'] = atoms
    data['protons'] = protons
    data['n_orbitals'] = int(fchk[num3:num3+100].split()[5])
    return data

def correct_orbital_order(mo_order,orbs,omin):
    """
    Rearrange the orbital indices in orbital energy ascending order.
    """
    temp = copy.deepcopy(orbs)
    for i in range(1,len(orbs)):
        for j in range(0,len(orbs)-1):
            if float(temp[j]) > float(temp[j+1]):
                temp[j+1],temp[j] = temp[j],temp[j+1]
    for i in range(len(orbs)):
        mo_order[i+omin] = temp.index(orbs[i])+omin # alpha
        mo_order[i+omin+len(orbs)] = temp.index(orbs[i])+omin+len(orbs) # beta
    return mo_order

def get_correct_orb_order(orb_map):
    """
    Map between molecular orbital energies and corresponding orbital indices.
    Q-CHEM CCMAN2 orbital indices are not in ascending order.

    Args:
        orb_map: Q-CHEM molecular orbital energy file.

    returns:
        Python dictionary giving the map between orbital energy and corresponding index.
    """
    onum1     = orb_map.find('O')
    onum2     = orb_map.find('V')
    occ_order = orb_map[onum1:onum2].split()[3:]
    vir_order = orb_map[onum2:].split()[3:]
    n_occ     = int(len(occ_order))
    n_vir     = int(len(vir_order))
    mo_map    = {}
    mo_map    = correct_orbital_order(mo_map,occ_order[:n_occ//2],0)
    mo_map    = correct_orbital_order(mo_map,vir_order[:n_vir//2],n_occ)
    return mo_map

def get_hf_energy(output):
    """ Parse Q-CHEM output and return HF energy. """
    num = output.find('Total energy in the final basis set')
    return float(output[num:num+100].split()[8])

def get_nuclear_repulsion(output):
    """ Parse Q-CHEM output and return nuclear repulsion energy. """
    num = output.find('Nuclear Repulsion Energy')
    return float(output[num:num+100].split()[4])

def get_mo_coeff(fchk, nbas):
    """ 
    Parse Q-CHEM checkpoint file and return molecular orbital coefficients.
    
    Args:
        fchk: Q-CHEM checkpoint file.
        nbas: The number of basis functions.

    Returns:
        Numpy array of molecular orbital coefficients.
    """
    num1      = fchk.find('Alpha MO coefficients')
    num2      = parse(fchk, from_where=num1, letter_type='Alpha Orbital Energies')
    aorbs     = [float(coeff) for coeff in fchk[num1:num2].split()[6:]]
    return np.array(aorbs).reshape(nbas, nbas).transpose()

def get_mo_energy(fchk):
    """ Parse Q-CHEM checkpoint file and return molecular orbital energies. """
    num1      = fchk.find('Alpha Orbital Energies')
    num2      = parse(fchk, from_where=num1, letter_type='Total SCF Density')
    aorbs_ene = [float(ene) for ene in fchk[num1:num2].split()[6:]]
    return np.array(aorbs_ene)

def get_overlap(fchk, nbas):
    """ Parse Q-CHEM checkpoint file and return overlap matrix. """
    num1 = fchk.find('Overlap Matrix')
    num2 = parse(fchk, from_where=num1, letter_type='Core Hamiltonian Matrix')
    return fill_full_matrix(fchk[num1:num2].split()[5:], nbas)

def get_h_core(fchk, nbas):
    """ Parse Q-CHEM checkpoint file and return core Hamiltonian. """
    num   = fchk.find('Core Hamiltonian Matrix')
    ncomp = int(fchk[num:num+100].split()[5])
    hcore = [float(coeff) for coeff in fchk[num:].split()[6:6+ncomp]]
    return fill_full_matrix(hcore,nbas)

def sign(p,q,r,s):
    if p%2 != r%2 or q%2 != s%2:
        return -1.0
    return +1.0

def get_one_body_integrals(spatial,norbs):
    """
    One electron integrals from spatial to spin orbital representation
    """
    one_body = np.zeros((norbs, norbs))
    for p in range(norbs//2):
        for q in range(norbs//2):
            # Populate 1-body coefficients. Require p and q have same spin.
            if abs(spatial[p, q]) > 10**-15:
                one_body[2 * p, 2 * q] = spatial[p, q]
                one_body[2 * p + 1, 2 * q + 1] = spatial[p, q]
    return one_body

def swap_block(spin,n_max):
    if spin < n_max/2: # alpha block in bare QCHEM two-electron integrals
        return int(spin*2) # alpha becomes even index
    # beta block in bare QCHEM two-electron integrals
    return int(2*(spin-n_max/2)+1) # beta becomes odd index

def qchem_two_body(mo_maps,pmax,qmax,rmax,smax,pmin,qmin,rmin,smin,raw,two_body,permute=None):
    """
    Transformations:

        i) Q-CHEM two electron integrals storage: 1D array and antisymmetrized <pq||rs>
            -> 1D to 4D

        ii) orbital indices: O: occ, V: vir, A: alpha, B: beta
        [OA1, OA2, ..., OAn, OB1, OB2, ..., OBn, VA1, VA2, ..., VAm, VB1, VB2, ..., VBm]
            -> [OA1, OB1, OA2, OB2, ..., OAn, OBn, VA1, VB1, VA2, VB2, ..., VAm, VBm]

        iii) permutations:
        OOOO and VVVV no need permutations
        OOOV = <pq||rs> -> -<pq||sr>, OOVO -> +<rs||pq>, OVOO -> -<sr||pq>, VOOO
        OOVV = <pq||rs> -> +<rs||pq>, VVOO
        OVOV = <pq||rs> -> -<pq||sr>, OVVO -> +<qp||sr>, VOVO -> -<qp||rs>, VOOV
        OVVV = <pq||rs> -> -<qp||rs>, VOVV -> +<rs||pq>, VVOV -> -<rs||qp>, VVVO
    """

    for p in range(pmax):
        for q in range(qmax):
            for r in range(rmax):
                for s in range(smax):
                    # correct orbital order
                    i,j,k,l = mo_maps[p+pmin],mo_maps[q+qmin],mo_maps[r+rmin],mo_maps[s+smin]
                    # swap spin block; [aa..abbb..b] -> [abab..ab]
                    i,j,k,l = swap_block(i-pmin,pmax),swap_block(j-qmin,qmax),swap_block(k-rmin,rmax),swap_block(l-smin,smax)
                    idx = (p*qmax*rmax*smax) + (q*rmax*smax) + (r*smax) + s
                    two_body[i+pmin,j+qmin,k+rmin,l+smin] += float(raw[idx])
                    if permute == 'ooov':
                        two_body[i+pmin,j+qmin,l+smin,k+rmin] += -1.0*float(raw[idx])
                        two_body[k+rmin,l+smin,i+pmin,j+qmin] += +1.0*float(raw[idx])
                        two_body[l+smin,k+rmin,i+pmin,j+qmin] += -1.0*float(raw[idx])
                    elif permute == 'oovv':
                        two_body[k+rmin,l+smin,i+pmin,j+qmin] += +1.0*float(raw[idx])
                    elif permute == 'ovov':
                        two_body[i+pmin,j+qmin,l+smin,k+rmin] += -1.0*float(raw[idx])
                        two_body[j+qmin,i+pmin,l+smin,k+rmin] += +1.0*float(raw[idx])
                        two_body[j+qmin,i+pmin,k+rmin,l+smin] += -1.0*float(raw[idx])
                    elif permute == 'ovvv':
                        two_body[j+qmin,i+pmin,k+rmin,l+smin] += -1.0*float(raw[idx])
                        two_body[k+rmin,l+smin,i+pmin,j+qmin] += +1.0*float(raw[idx])
                        two_body[k+rmin,l+smin,j+qmin,i+pmin] += -1.0*float(raw[idx])
    return two_body

def make_openfermion_format(qchem_two_body):
    """
    Convert Q-CHEM two electron integral to OpenFermion format.
    <pq||rs> -> <pq||sr>
    """
    n           = len(qchem_two_body)
    of_two_body = np.zeros((n,n,n,n))
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    of_two_body[p,q,r,s] = qchem_two_body[p,q,s,r]
    return of_two_body

def get_two_body_integrals(vpqrs,mo_maps):
    """
    Parse Q-CHEM two electron integral file and return Numpy array of two electron integrals.

    Args:
        vqprs: Q-CHEM two electron integral file.
        mo_maps: Mapping between orbital index and orbital energy.

    Returns:
        Numpy array of two electron integrals, Vpqsr.
    """
    # read Q-CHEM two-electron integrals
    # O: occupied V: virtual
    num1 = vpqrs.find('OOOO')
    num2 = vpqrs.find('OOOV')
    num3 = vpqrs.find('OOVV')
    num4 = vpqrs.find('OVOV')
    num5 = vpqrs.find('OVVV')
    num6 = vpqrs.find('VVVV')
    oooo = vpqrs[num1:num2].split()[6:]
    ooov = vpqrs[num2:num3].split()[6:]
    oovv = vpqrs[num3:num4].split()[6:]
    ovov = vpqrs[num4:num5].split()[6:]
    ovvv = vpqrs[num5:num6].split()[6:]
    vvvv = vpqrs[num6:].split()[6:]

    n_occ    = int(vpqrs[num1:num1+100].split()[2])
    n_vir    = int(vpqrs[num3:num3+100].split()[4])
    n_orb    = int(n_occ+n_vir)
    two_body = np.zeros((n_orb,n_orb,n_orb,n_orb))

    two_body = qchem_two_body(mo_maps,n_occ,n_occ,n_occ,n_occ,0,0,0,0,oooo,two_body)
    two_body = qchem_two_body(mo_maps,n_occ,n_occ,n_occ,n_vir,0,0,0,n_occ,ooov,two_body,permute='ooov')
    two_body = qchem_two_body(mo_maps,n_occ,n_occ,n_vir,n_vir,0,0,n_occ,n_occ,oovv,two_body,permute='oovv')
    two_body = qchem_two_body(mo_maps,n_occ,n_vir,n_occ,n_vir,0,n_occ,0,n_occ,ovov,two_body,permute='ovov')
    two_body = qchem_two_body(mo_maps,n_occ,n_vir,n_vir,n_vir,0,n_occ,n_occ,n_occ,ovvv,two_body,permute='ovvv')
    two_body = qchem_two_body(mo_maps,n_vir,n_vir,n_vir,n_vir,n_occ,n_occ,n_occ,n_occ,vvvv,two_body)

    # <pq||rs> -> <pq||sr>
    return make_openfermion_format(two_body)

def get_mp2_energy(output):
    """ Parse Q-CHEM output and return MP2 energy. """
    if 'MP2         total energy' in output:
        num = output.find('MP2         total energy')
        return float(output[num:num+100].split()[4])
    elif 'MP2 energy' in output:
        num = output.find('MP2 energy')
        return float(output[num:num+100].split()[3])

def get_ccsd_energy(output):
    """ Parse Q-CHEM output and return CCSD energy. """
    num = output.find('CCSD total energy')
    return float(output[num:num+100].split()[4])

def get_ccsd_t1_amps(cc_t1,mo_maps):
    """ Parse and return Q-CHEM CCSD T1 amplitude"""
    # t1_ia -> stored as an array i*amax+a array
    raw     = cc_t1[:].split()[5:]
    n_occ   = int(cc_t1[:].split()[3])
    n_vir   = int(cc_t1[:].split()[4])
    n       = n_occ + n_vir
    t1_amps = np.zeros(((n,n)))

    pmax,qmax = n_occ,n_vir
    pmin,qmin = 0,n_occ
    for p in range(pmax):
        for q in range(qmax):
            i,a = mo_maps[p+pmin],mo_maps[q+qmin]
            i,a = swap_block(i-pmin,pmax),swap_block(a-qmin,qmax)
            idx = (p * qmax) + q
            t1_amps[a+qmin,i+pmin] = float(raw[idx])
    return t1_amps

def get_ccsd_t2_amps(cc_t2,mo_maps):
    """ Parse and return Q-CHEM T2 amplitude """
    # t2_ijab -> stored as an array i*(jmax amax bmax) + j *amax bmax + a*bmax +b
    raw     = cc_t2[:].split()[7:]
    n_occ   = int(cc_t2[:].split()[3])
    n_vir   = int(cc_t2[:].split()[5])
    n       = n_occ + n_vir
    t2_amps = np.zeros(((n,n,n,n)))

    # p,q = occ & r,s = vir
    pmax,qmax,rmax,smax = n_occ,n_occ,n_vir,n_vir
    pmin,qmin,rmin,smin = 0,0,n_occ,n_occ
    for p in range(pmax):
        for q in range(qmax):
            for r in range(rmax):
                for s in range(smax):
                    # correct orbital order
                    i,j,a,b = mo_maps[p+pmin],mo_maps[q+qmin],mo_maps[r+rmin],mo_maps[s+smin]
                    # swap spin block; [aa..abbb..b] -> [abab..ab]
                    i,j,a,b = swap_block(i-pmin,pmax),swap_block(j-qmin,qmax),swap_block(a-rmin,rmax),swap_block(b-smin,smax)
                    idx = (p*qmax*rmax*smax) + (q*rmax*smax) + (r*smax) + s
                    t2_amps[a+rmin,i+pmin,b+smin,j+qmin] += float(raw[idx])
    return 0.5*t2_amps