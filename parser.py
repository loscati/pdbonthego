from typing import List

import Bio.PDB as pdb
from Bio.PDB.Selection import unfold_entities

# TODO: better handle warning by PDBParser (see github on PDB warning)

def get_residues(
    file: str, 
    mod: int, 
    ch: int, 
    first_to_remove: int
    ) -> List[pdb.Residue.Residue]:
    '''Residues from pdb file
    It select only proteinogenic residues excluding all water molecules and
    ions.

    Args:
        file (str): absolute/relative path for pdb file
        mod (int): selects the wanted model (must be => 0)
        ch (int): selects the wanted chain (must be => 0)
        first_to_remove (int): number of residues to remove from the beginning
        of the chain (e.g. because they are added artificially to make the
        protein crystallize)

    Returns:
        list: list of proteinogenic residues excluding all water molecules
        and ions
    '''

    parser = pdb.PDBParser()
    name_protein = file[-8:-4] # get unique protein ID of 4 characters
    # This assumes that the file name is the protein ID
    structure = parser.get_structure(name_protein, file)

    # Info to print while calling the parser
    # print(f'Parsing: {name_protein}')
    # print('Models: ', len(list(structure.get_models())))
    # print('Chains: ', len(list(structure.get_chains())))

    # Unpacking the selected chain
    models = unfold_entities(structure, 'M')
    chains = unfold_entities(models[mod], 'C')
    res_list_full = unfold_entities(chains[ch], 'R')
    # filtering out all but proteinogeneic residues
    proteinogenic_res = \
    ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET',
    'ASN','PYL','PRO','GLN','ARG','SER','THR','SEC','VAL','TRP','TYR']
    res_list = []
    for res in res_list_full:
        if res.get_resname() in proteinogenic_res:
            res_list.append(res)

    return res_list[first_to_remove:]


if __name__ == '__main__':
    pass
