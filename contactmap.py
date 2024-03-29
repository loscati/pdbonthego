"""Compute and draw the interaction matrix, or contact map from a PDB structure"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pdb_parser import get_residues

def get_contact_map(file, model_id, chain_id, threshold, altloc="A", to_include=None, to_ignore=None):
    """
    Matrix with the interaction network, or contact map, extracted
    from a PDb structure. The contact definition is purely geometric.
    Two residues are said to be in contact if at least one pair of
    non-hydrogen atoms is spatially closer than the threshold.

    The matrix entry is the distance between the alpha carbons of the
    respective residues.

    Disulfide bridges are not considered.

    Parameters
    --------
        See pdb_parser.get_residues for more details
        threshold : float
            Distance below which two non-hydrogen atoms are considered in contact
        altloc : str
            Alternative location for disordered atoms. Default: 'A'

    Returns
    -------
        contact_map : array_like
            Two dimensional numpy array with entries the alpha carbon distances
            of contacting residues. All dimensions are equal to the number of
            residues in the protein
    """
    res_list = get_residues(file, model_id, chain_id, to_include, to_ignore)
    n_residues = len(res_list)
    contact_map = np.zeros((n_residues,n_residues))

    # Contacts are computed between residues which are separated by at least
    # three other residues in the polypeptide chain [noel2012]
    for row in range(0, n_residues-4):
        for col in range(row+4, n_residues):

            res1 = res_list[row]
            res2 = res_list[col]
            stop_res1 = False
            for atom1 in res1:
                # Contact map takes into account only heavy atoms.
                # In PDBs, Hydrogen are highlighted by a string with
                # initial char equal to H. But there are atoms such as
                # Oxygen in ARG that are identified as OH
                if atom1.get_name().startswith("H"):
                    continue
                if atom1.is_disordered():
                    atom1.disordered_select(altloc)
                for atom2 in res2:
                    if atom2.get_name().startswith("H"):
                        continue
                    if atom2.is_disordered():
                        atom2.disordered_select(altloc)

                    # the minus operator between Bio.PDB.Atom.Atom
                    # objects is overloaded to return the distance
                    # between them
                    distance = atom1 - atom2
                    if distance < threshold:
                        Ca1 = res1['CA']
                        Ca2 = res2['CA']
                        if Ca1.is_disordered():
                            Ca1.disordered_select(altloc)
                        if Ca2.is_disordered():
                            Ca2.disordered_select(altloc)

                        contact_map[row,col] = Ca2 - Ca1
                        stop_res1 = True
                        break
                if stop_res1:
                    break
    return contact_map


def visualize_contact_map(
    cm: np.ndarray, 
    cmap: str, 
    title: str, 
    reflected: bool =True
    ) -> matplotlib.figure.Figure:
    '''Visualizing contact map cm using an heat map

    Args:
        cm (2D numpy array): contact map
        cmpa (str): cmap selection passed to imshow
        title (str): title of the plot

    Returns:
        Matplotlib Figure object
    '''
    matplotlib.rcParams.update({
        'font.size': 30,
        'text.usetex': True
    })

    fig, ax = plt.subplots(figsize=(8,7))
    ax.set_title(title, pad=10, fontsize=30)

    if reflected:
        # fix y ticks and labels: display tick and label every 10 residues
        res_n = cm.shape[0]
        fix_yloc = [i for i in range(res_n%5, res_n - res_n%5, 5)]
        fix_ylabel = [i for i in range(res_n - res_n%5, res_n%5, -5)]
        ax.set_yticks(fix_yloc)
        ax.set_yticklabels(fix_ylabel, fontsize=20)

        # x ticks and labels
        ax.set_xticks([i for i in range(4, res_n, 5)])
        ax.set_xticklabels([i for i in range(5, res_n, 5)], fontsize=20)

        # Obtain standard representation of contact maps, e.g. [noel2016]:
        im = ax.imshow(cm[::-1,:], cmap=cmap)
    else:
        im = ax.imshow(cm, cmap=cmap)   


    ax.grid(alpha=.5)
    ax.set_xlabel('Residue index', labelpad=10)
    ax.set_ylabel('Residue index', labelpad=10)

    # Lateral legend bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    c = fig.colorbar(im, ax=ax, cax=cax)
    c.set_label('\AA', rotation=0, labelpad=10)

    plt.show()

    return fig


if __name__ == '__main__':
    pass
