import numpy as np

import parser
import parameters
import contactmap


def parser():
    file = '../pdb/1fd3.pdb'
    res_list = par.get_residues(file, 0, 0, 5)
    print(res_list)


def parameters_dist():
    file = '../pdb/1r6j.pdb'
    dist = parameters.get_residue_dists(file, altloc='A', first_to_remove=5)
    dist_alt = parameters.get_residue_dists(file, altloc='B', first_to_remove=5)
    result = []
    for d1, d2 in zip(dist,dist_alt):
        result.append(np.allclose(d1,d2, atol=1e-3))
    print(result)


def parameters_disulfide():
    file = '../pdb/1fd3.pdb'
    dist = parameters.get_disulfide_bonds(file)
    print(dist)


def parameters_angles():
    file = '../pdb/1fd3.pdb'
    dist = parameters.get_residues_angles(file)
    print(dist)


def parameters_dih():
    file = '../pdb/1ucs.pdb'
    dist = parameters.get_residues_dih(file)
    print(dist)


def cm():
    file = '../pdb/1r6j.pdb'
    cc = contactmap.get_contact_map(file, 4.5)
    contactmap.visualize_contact_map(cc[::-1,:], f'Contact map: {file[-8:-4]}')


def save_para():
    file = '../pdb/1ucs.pdb'
    thr = 4.5
    altloc = 'A'
    first_to_remove = 0
    parameters.save_parameters(file, threshold=thr, altloc=altloc,
                                first_to_remove=first_to_remove)
