# PDBontheGō
Set of functions for extracting parameters from a PDB file to construct a coarse-grained, Gō-model, of the analyzed protein. They allow you to compute: bond lengths, dihedral angles, contact map (gemetrical definition, more details in `contactmap.py`), etc.

The porject is divided into:
+ `parser` which contains a `Bio.PDB.PDBParser()` wrapper from the [Biopython project](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ) and it is used to extract residues objects from the PDB file
+ `parameters` contains functions to extract bond lengths, bond angles, dihedral angles and disulfite bonds lengths
+ `contactmap` helps you to compute the native contact map as define [here](https://www.nature.com/articles/s41598-019-44928-3) and plot it.

Each function has a detailed docstring explaning what it does, arguments and its result.

### Installation
```bash
git clone https://github.com/LeonardoSalicari/pdbonthego.git
```

### Usage example
Clone this repository into the project directory
```Python
import parser
import parameters
import contactmap 

prot_name = '1ucs' # PDB name example
path_pdb = f'/path/to/pdbfile/{prot_name}.pdb'

# Parameters to select model-chain-altloc (see Biopython wiki and pdb101.rcsb.org/)
altloc = 'A' # alternative location
model = 0 # model index
chain = 0 # chain index
residues_to_remove = 0 # from the beginning of the chain (N-terminus)
threshold = 4.5 # threshold for the contact map, same lenght units of the PDB file (usually Angstrom)

res = parser.get_residues(path_pdb, altloc, model, chain, residues_to_remove)
bonds = parameters.get_residue_dists(path_pdb, altloc, model, chain, residues_to_remove)
contact_map = contactmap.get_contact_map(path_pdb, threshold, altloc, model, chain, residues_to_remove)
```

In order to select the right model-chain-altloc you have to know the structure of the PDB file, [this site](https://pdb101.rcsb.org/) contain a complete introduction to this datastructure.
More info about how `parser`, and all other function which uses it, handle these arguments can be found in the *The Structure Object* section of the [Bio.PDB package](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ).

### Compatibility
The package is tested with Python 3.9.5 and `biopython` version 1.79

### TODO
- [ ] create a `setup.py` file 
