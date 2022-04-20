import numpy as np
import itertools
from openfermionqchem._qchem_parser import make_openfermion_format

def fact(n):
    """ Compute n! """
    if n > 1:
        return n * fact(n-1)
    else:
        return 1

def bubble_sort(list1,list2):
    """ Simple bubble sort function """
    n_eles = len(list1)
    m,n    = 0,0 # number of permutations
    for i in range(1,n_eles):
        for j in range(0, n_eles-1):
            if list1[j] > list1[j+1]:
                m += 1
                list1[j+1],list1[j] = list1[j],list1[j+1]
            if list2[j] > list2[j+1]:
                n += 1
                list2[j+1],list2[j] = list2[j],list2[j+1]
    return m,n

def determinants(n_spin_orbs, n_eles, spin_flip=False):
    # n_determinants = n_spin_orbs choose n_eles
    n_dets     = int(fact(n_spin_orbs) / (fact(n_eles) * fact(n_spin_orbs - n_eles)))
    if spin_flip:
        # single-SF -> n_alpha * n_beta_orbs choose one
        n_dets = n_eles * (n_spin_orbs / 2)
    # itertools.product(), computes the cartesian product of input iterables.
    # e.g) print list(product([1,2,3],repeat = 2))
    # [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
    all_permute = list(itertools.product(range(n_spin_orbs),repeat=n_eles))
    dets        = []
    for conf in all_permute:
        if set(conf) in dets: # remove repetition
            continue
        elif len(set(conf)) < n_eles: # pauli-exclusion principle
            continue
        else:
            if spin_flip:
                n_beta = len([spin for spin in conf if spin % 2 == 1])
                if n_beta == 0 or n_beta > 1:
                # if sum([spin % 2 for spin in conf]) == 0 or sum([spin % 2 for spin in conf]) == n_eles:
                    continue
                else:
                    dets.append(set(conf))
            else:
                dets.append(set(conf))

    if len(dets) != n_dets:
        print('check determinants')
        exit()

    return dets

def compare_dets(bra,ket):
    if len(set(bra).difference(set(ket))) > 2:
        return 3, 1, [], -1, -1, -1, -1

    bra,ket = list(bra),list(ket)
    if len(set(bra).difference(set(ket))) == 0:
        m,n         = bubble_sort(bra,ket)
        p1,p2,a1,a2 = -1,-1,-1,-1
        differ_by   = 0

    elif len(set(bra).difference(set(ket))) == 1:
        p2,a2  = -1,-1
        p1     = list(set(bra)-set(ket))[0] # reference
        a1     = list(set(ket)-set(bra))[0] # substituted
        idx_p1 = bra.index(p1)
        idx_a1 = ket.index(a1)
        # move the diff spin orbital to the right most position
        bra.remove(p1); ket.remove(a1)
        bra.append(p1); ket.append(a1)
        # only permute common indices
        m,n = bubble_sort(bra[:-1],ket[:-1])
        # extra permutation due to the diff spin orbital
        m += (len(bra)-idx_p1-1)
        n += (len(ket)-idx_a1-1)
        differ_by = 1

    elif len(set(bra).difference(set(ket))) == 2:
        bra_complement = [x for x in bra if x not in ket] # maintain orders
        ket_complement = [x for x in ket if x not in bra]
        p1,p2          = bra_complement[0],bra_complement[1]
        a1,a2          = ket_complement[0],ket_complement[1]
        # move first
        idx_p1 = bra.index(p1)
        idx_a1 = ket.index(a1)
        bra.remove(p1); ket.remove(a1)
        bra.append(p1); ket.append(a1)
        e_m = (len(bra)-idx_p1-1)
        e_n = (len(ket)-idx_a1-1)
        # move second
        idx_p2 = bra.index(p2)
        idx_a2 = ket.index(a2)
        bra.remove(p2); ket.remove(a2)
        bra.append(p2); ket.append(a2)
        e_m += (len(bra)-idx_p2-1)
        e_n += (len(ket)-idx_a2-1)
        m,n  = bubble_sort(bra[:-2],ket[:-2])
        # extra permutation due to the diff spin orbital
        m += e_m
        n += e_n
        differ_by = 2

    sign   = pow(-1,m) * pow(-1,n)
    common = list(set(bra)&set(ket))    # common index
    return differ_by, sign, common, p1, p2, a1, a2

def build_fci_matrix(n_eles,n_spin_orbs,one_body,two_body):
    # dimension of FCI matrix = n_determinants x n_determinants
    # n_determinants = n_spin_orbs choose n_eles
    dets       = determinants(n_spin_orbs, n_eles)
    fci_matrix = np.zeros((len(dets),len(dets)))
    # <pq||sr> -> <pq||rs>
    two_body   = make_openfermion_format(two_body)
    """
    Slater Rules
    H = Z + V
    1. <P1P2P3...Pn|H|P1P2P3...Pn>, two determinants are identical
        one_body: <Pi|z|Pi>, i from 1 to n
        two_body: 0.5 * <PiPj||PiPj>, i and j from 1 to n
    2. <P1P2P3...Pn|H|A1P2P3...Pn>, two determinants differ by one; P1->A1
        one_body: <P1|z|A1>
        two_body: <P1Pj||A1Pj>, j from 1 to n
    3. <P1P2P3...Pn|H|A1A2P3...Pn>, two determinants differ by two; P1->A1 and P2->A2
        one_body: no contributions
        two_body: <P1P2||A1A2>
    """
    for bra in range(len(dets)):
        for ket in range(len(dets)):
            differ_by,sign,common,p1,p2,a1,a2 = compare_dets(dets[bra],dets[ket])
            if differ_by == 0:
                for i in common:
                    fci_matrix[bra,ket] += sign * one_body[i,i]
                    for j in [ele for ele in common if ele != i]:
                    # for j in common:
                        fci_matrix[bra,ket] += 0.5 * sign * two_body[i,j,i,j]
                        # fci_matrix[bra][ket] += 0.25 * sign * two_body[i][j][i][j]
                        # fci_matrix[bra][ket] -= 0.25 * sign * two_body[i][j][j][i]

            elif differ_by == 1:
                fci_matrix[bra,ket] += sign * one_body[p1,a1]
                for j in common:
                    fci_matrix[bra,ket] += sign * two_body[p1,j,a1,j]
                    # fci_matrix[bra][ket] += 0.5 * sign * two_body[p1][j][a1][j]
                    # fci_matrix[bra][ket] -= 0.5 * sign * two_body[p1][j][j][a1]

            elif differ_by == 2:
                fci_matrix[bra,ket] += sign * two_body[p1,p2,a1,a2]
                # fci_matrix[bra][ket] += 0.5 * sign * two_body[p1][p2][a1][a2]
                # fci_matrix[bra][ket] -= 0.5 * sign * two_body[p1][p2][a2][a1]

            else:
                continue

    return fci_matrix