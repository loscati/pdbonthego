import math

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import parser
import parameters
import tester

def get_contact_map(file, threshold, altloc='A', mod=0, ch=0,
                    first_to_remove=0) -> '2D numpy array':
    '''2D map of distances between C_alphas which have, in their residue, at
    least a pair of heavy atoms closer than threshold

    From the contact map, disulfide bridges are excluded

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

    Retuns:
        2D numpy array: a 2 dimentional numpy array containg C_alpha distances
        of residues 'in contact'. All dimentions are equal to the number of
        residues in the protein. Following the convention, the diagonal is
        placed from bottom left to top right; where the first residue pair is
        placed at the bottom left. Moreover, the matrix is symmetric w.r.t.
        the anti-diagonal
    '''
    res_list = parser.get_residues(file, mod, ch, first_to_remove)
    n = len(res_list)
    cm = np.zeros((n,n))

    # Contacts are computed between residues which are separated by at least
    # three other residues in the polypeptide chain [noel2012]
    for row in range(0, n-4):
        for col in range(row+4, n):

            res1 = res_list[row]
            res2 = res_list[col]
            stop_res1 = False
            for atom1 in res1:
                # Contact map takes into account only heavy atoms
                if atom1.get_name() == 'H':
                    continue
                if atom1.is_disordered():
                    atom1.disordered_select(altloc)
                for atom2 in res2:
                    if atom2.get_name() == 'H':
                        continue
                    if atom2.is_disordered():
                        atom2.disordered_select(altloc)

                    distance = atom1 - atom2

                    # TODO:
                    # Mettere opzione per considerare ponti di disofuro
                    # come contatti o meno
                    
                    # Check if the pair of atoms are involved
                    # in a disulfide bond
                    # if (atom1.get_name() == 'SG'
                    #     and atom2.get_name() == 'SG'):
                    #     print(distance)
                    #     print('PONTE')
                    #     continue

                    if distance < threshold:

                        ####### DEBUG
                        # if row == 4 or row == 6 or row == 10:
                        #     print(f'Residues: ({row+1},{col+1})')
                        #     print(f'Contact atom from res {row+1}: {atom1.get_name()}')
                        #     print(f'Contact atom from res {col+1}: {atom2.get_name()}')
                        #     print(f'Distance (not between CA): {distance:.4f}\n\n')
                        #######

                        Ca1 = res1['CA']
                        Ca2 = res2['CA']
                        if Ca1.is_disordered():
                            Ca1.disordered_select(altloc)
                        if Ca2.is_disordered():
                            Ca2.disordered_select(altloc)

                        cm[row,col] = Ca2 - Ca1
                        stop_res1 = True
                        break

                if stop_res1:
                    break

    # Obtain standard representation of contact maps, e.g. [noel2016]:
    # no subtraction of the diag because cm has always zeros in the diag
    cm = cm + cm.T
    return cm


def visualize_contact_map(cm, title) -> None:
    '''Visualizing contact map cm using an heat map

    Args:
        cm (2D numpy array): contact map
        title (str): title of the plot

    Returns:
        None
    '''
    
    matplotlib.rcParams.update({'font.size': 16})

    fig, ax = plt.subplots()
    ax.set_title(title, pad=15)

    # fix y ticks and labels: display tick and label every 10 residues
    res_n = cm.shape[0]
    fix_yloc = [i for i in range(res_n%5, res_n - res_n%5, 5)]
    fix_ylabel = [i for i in range(res_n - res_n%5, res_n%5, -5)]
    ax.set_yticks(fix_yloc)
    ax.set_yticklabels(fix_ylabel)

    # x ticks and labels
    ax.set_xticks([i for i in range(4, res_n, 5)])
    ax.set_xticklabels([i for i in range(5, res_n, 5)])

    ax.grid()
    ax.set_xlabel('Residues')
    ax.set_ylabel('Residues')
    
    # Obtain standard representation of contact maps, e.g. [noel2016]:
    im = ax.imshow(cm[::-1,:], cmap = 'binary')
    # Lateral legend bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, ax=ax, cax=cax)

    plt.show()


if __name__ == '__main__':
    tester.cm()
