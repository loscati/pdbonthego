import math
from typing import List

import numpy as np

from .parser import get_residues

def get_residue_dists(
    file: str, 
    altloc: str ='A', 
    mod: int=0, 
    ch: int=0, 
    first_to_remove:int =0
    ) -> List[List[float]]:
    '''Distances between C_alpha atoms of subsequent residues
    This corresponds to harmonic bond terms in the Go-like Hamiltonian

    Args:
        file (str): absolute/relative path for pdb file
        altloc (str) Specifies the protein configuration
        mod (int): selects the wanted model (must be => 0)
        ch (int): selects the wanted chain (must be => 0)
        first_to_remove (int): number of residues to remove from the beginning
        of the chain (e.g. because they are added artificially to make the
        protein crystallize)

    Returns:
        list: list of the distances between C_alpha atoms
    '''

    res_list = get_residues(file, mod, ch, first_to_remove)
    calpha_dists = []
    for i in range(1, len(res_list)):
        Ca0 = res_list[i-1]['CA']
        Ca1 = res_list[i]['CA']

        # Taking the altloc picked by the user
        if Ca0.is_disordered():
            Ca0.disordered_select(altloc)
        if Ca1.is_disordered():
            Ca1.disordered_select(altloc)

        # minus operator has been overwritten to return distances
        calpha_d = Ca1 - Ca0
        calpha_dists.append(calpha_d)

    return calpha_dists


def get_disulfide_bonds(
    file: str, 
    altloc: str ='A', 
    mod: int=0, 
    ch: int=0, 
    first_to_remove:int =0
    ) -> List[float]:
    '''Disulfide bridges bond lengths
    It search for pair of Cysteines if present in the chain

    Args:
        file (str): absolute/relative path for pdb file
        altloc (str) Specifies the protein configuration
        mod (int): selects the wanted model (must be => 0)
        ch (int): selects the wanted chain (must be => 0)
        first_to_remove (int): number of residues to remove from the beginning
        of the chain (e.g. because they are added artificially to make the
        protein crystallize)

    Returns:
        list: list of disulfide bond lengths
    '''

    res_list = get_residues(file, mod, ch, first_to_remove)
    cys_list = []
    # Checking for pairs of Cysteins
    for res in res_list:
        if res.get_resname() == 'CYS':
            cys_list.append(res)

    n_cys = len(cys_list)
    if len(cys_list) <= 1:
        # no bonds with only one CYS
        return []

    dis_bridges = []
    for i in range(n_cys):
        for j in range(i+1, n_cys):
            SG1 = cys_list[i]['SG']
            SG2 = cys_list[j]['SG']
            bond_len = SG2 - SG1
            # Saving only actual bonds, which have length of about 2.05 AA
            if bond_len < 2.15:
                dis_bridges.append(bond_len)

    return dis_bridges


def get_residues_angles(
    file: str, 
    altloc: str ='A', 
    mod: int=0, 
    ch: int=0, 
    first_to_remove:int =0,
    degrees: bool =True
    ) -> List[float]:
    '''Angles between successive residues, hence C_alpha

    Angle definition involves three successive C_alphas and makes use of the
    dot product:
    vec{dr} = vec{C_alpha^{i}} - vec{C_alpha^{i-1}}
    vec{dp} = vec{C_alpha^{i+1}} - vec{C_alpha^{i}}
    theta_i^0 = pi - acos( (dr_{i-1,i} \dot dr_{i,i+1})/modules )

    pi - acos() is used because we do not need the angle as defined by the
    dot product, instead we seek the angle form by the 'bending' of two
    subsequent residues (described by their C_alpha).

    Args:
        file (str): absolute/relative path for pdb file
        altloc (str) Specifies the protein configuration
        mod (int): selects the wanted model (must be => 0)
        ch (int): selects the wanted chain (must be => 0)
        first_to_remove (int): number of residues to remove from the beginning
        of the chain (e.g. because they are added artificially to make the
        protein crystallize)
        degrees (bool): choose between radians and degrees, if set to False
        the result will be in radians

    Returns:
        list: list of angles in radians or degrees
    '''

    res_list = get_residues(file, mod, ch, first_to_remove)
    angles = []
    for i in range(1, len(res_list) - 1):
        Ca0 = res_list[i-1]['CA']
        Ca1 = res_list[i]['CA']
        Ca2 = res_list[i+1]['CA']

        # Taking the altloc picked by the user
        if Ca0.is_disordered():
            Ca0.disordered_select(altloc)
        if Ca1.is_disordered():
            Ca1.disordered_select(altloc)
        if Ca2.is_disordered():
            Ca2.disordered_select(altloc)

        dr = Ca1.get_coord() - Ca0.get_coord()
        dp = Ca2.get_coord() - Ca1.get_coord()
        dr_mag = np.sqrt(dr.dot(dr)) # mag = magnitude
        dp_mag = np.sqrt(dp.dot(dp))

        # math.acos(): [-1,1] -> [0, pi]
        # math.acos(-1) = pi; math.acos(1) = 0
        theta = math.pi - math.acos(dr.dot(dp)/(dr_mag*dp_mag)) # radians

        angles.append(theta)

    if degrees:
        angles = [rad*180/math.pi for rad in angles]

    return angles


def get_residues_dih(
    file: str, 
    altloc: str ='A', 
    mod: int=0, 
    ch: int=0, 
    first_to_remove:int =0,
    degrees: bool =True
    ) -> List[float]:
    '''Dihedral angles between successive residues

    For details about the computation see:
    https://en.wikipedia.org/wiki/Dihedral_angle#In_polymer_physics

    Args:
        file (str): absolute/relative path for pdb file
        altloc (str) Specifies the protein configuration
        mod (int): selects the wanted model (must be => 0)
        ch (int): selects the wanted chain (must be => 0)
        first_to_remove (int): number of residues to remove from the beginning
        of the chain (e.g. because they are added artificially to make the
        protein crystallize)
        degrees (bool): choose between radians and degrees, if set to False
        the result will be in radians

    Returns:
        list: list of dihedral angles in radians or degrees
    '''

    res_list = get_residues(file, mod, ch, first_to_remove)
    dihedrals = []
    for i in range(2, len(res_list) - 1):
        Ca0 = res_list[i-2]['CA']
        Ca1 = res_list[i-1]['CA']
        Ca2 = res_list[i]['CA']
        Ca3 = res_list[i+1]['CA']

        # Taking the altloc picked by the user
        if Ca0.is_disordered():
            Ca0.disordered_select(altloc)
        if Ca1.is_disordered():
            Ca1.disordered_select(altloc)
        if Ca2.is_disordered():
            Ca2.disordered_select(altloc)
        if Ca3.is_disordered():
            Ca3.disordered_select(altloc)

        u1 = Ca1.get_coord() - Ca0.get_coord()
        u2 = Ca2.get_coord() - Ca1.get_coord()
        u3 = Ca3.get_coord() - Ca2.get_coord()
        u2_mag = np.sqrt(u2.dot(u2))

        n1 = np.cross(u1,u2)
        n2 = np.cross(u2,u3)

        # math.atan2: R --> [-pi, pi]
        dih = math.atan2(u2_mag*u1.dot(n2), n1.dot(n2))
        dihedrals.append(dih)

    if degrees:
        dihedrals = [rad*180/math.pi for rad in dihedrals]

    return dihedrals


def save_parameters(
    file: str, 
    altloc: str ='A', 
    mod: int=0, 
    ch: int=0, 
    first_to_remove:int =0, 
    degrees: bool =True
    ) -> None:
    '''Save on a .dat file all parameters for the Go-like hamiltonian.
    For example is the pdb file name is "1ucs.pdb" and the configuration
    selecteds is 'A', then all parameters are saved in the "1ucs-A.dat" file

    Args:
        file (str): absolute/relative path for pdb file
        threshold (float): max distance for which two heavy atoms are considered
        'in contact'
        altloc (str) Specifies the protein configuration
        mod (int): selects the wanted model (must be => 0)
        ch (int): selects the wanted chain (must be => 0)
        first_to_remove (int): number of residues to remove from the beginning
        of the chain (e.g. because they are added artificially to make the
        protein crystallize)
        degrees (bool): choose between radians and degrees, if set to False
        the result will be in radians

    Returns:
        None
    '''

    res_list = get_residues(file, mod, ch, first_to_remove)
    dist = get_residue_dists(file, altloc=altloc,
                            first_to_remove=first_to_remove)
    disulfide = get_disulfide_bonds(file, altloc=altloc,
                                    first_to_remove=first_to_remove)
    angles = get_residues_angles(file, altloc=altloc,
                                first_to_remove=first_to_remove,
                                degrees=degrees)
    dih = get_residues_dih(file, altloc=altloc,
                            first_to_remove=first_to_remove,
                            degrees=degrees)
    # contact_map = get_contact_map(file, altloc=altloc,
    #                               threshold=threshold,
    #                               first_to_remove=first_to_remove)

    name_protein = file[-8:-4] # get unique protein ID of 4 characters
    
    with open(f'{name_protein}-{altloc}.dat', 'w') as fout:
        # Info about the file
        fout.write(f'INFO\nProtein: {name_protein},\nAltloc: {altloc},\nNumber of residues removed in the beginning: {first_to_remove},\nModel number: {mod},\nChain number: {ch},\n')
        fout.write('Sequence:\n')
        for res in res_list:
            fout.write(f'{res.get_resname()} ')
        fout.write('\n\n----------------------------------------\n')
        fout.write('PARAMETERS\n\n')

        fout.write('C_alpha distances (Angstroms):\n')
        fout.write('[')
        for d in dist:
            fout.write(f'{d:.3f}, ')
        fout.write(']\n\n')

        fout.write('Disulfide bond lengths (Angstroms):\n')
        fout.write('[')
        for d in disulfide:
            fout.write(f'{d:.3f}, ')
        fout.write(']\n\n')

        fout.write('Angles between residues (degrees):\n')
        fout.write('[')
        for a in angles:
            fout.write(f'{a:.3f}, ')
        fout.write(']\n\n')

        fout.write('Dihedral angles (degrees):\n')
        fout.write('[')
        for d in dih:
            fout.write(f'{d:.3f}, ')
        fout.write(']\n\n')

    return None


if __name__ == '__main__':
    pass
