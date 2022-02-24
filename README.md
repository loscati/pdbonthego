# PDBontheGō
Set of functions used to extract parameters from a PDB file in order to construct a coarse-grained Gō-model of the analyzed protein, such as bond lengths, dihedral angles, contact map etc.

Functions are grouped into four files characterized as follow:
+ `parser` contains a wrapper to `Bio.PDB.PDBParser()` of the [Biopython project](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ) used to extract only residues from the PDB file
+ `parameters` has functions to extract bond length, bond angles, dihedral angles and disulfite bonds
+ `contactmap` computes the native contact map as define [here](https://www.nature.com/articles/s41598-019-44928-3) and plots it.

Each function has a detailed docstring explaning what it does, arguments and its result.

## Installation
```bash
git clone https://github.com/LeonardoSalicari/pdbonthego.git
```

## Usage example
```Python
import sys
sys.path.insert(0, "/path/to/pdbonthego/directory")
import parser
import parameters
import contactmap as cm

prot_name = '1ucs' # PDB name example
path_pdb = f'/path/to/pdbfile/{prot_name}.pdb'

# Parameters to select the wanted model-chain-altloc
altloc = 'A' # alternative location
model = 0 # model index
chain = 0 # chain index
residues_to_remove = 0 # from the beginning of the chain (N-terminus)
threshold = 4.5 # threshold for the contact map, same lenght units of the PDB file (usually Angstrom)

res = parser.get_residues(path_pdb, altloc, model, chain, residues_to_remove)
bonds = parameters.get_residue_dists(path_pdb, altloc, model, chain, residues_to_remove)
contact_map = cm.get_contact_map(path_pdb, threshold, altloc, model, chain, residues_to_remove)
```

In order to select the right model-chain-altloc you have to know the structure of a PDB file, [this site](https://pdb101.rcsb.org/) contain a complete introduction to this type of files.
More info about how `parser`, and all other function which uses it, handle these arguments can be found in the *The Structure Object* section of the [Bio.PDB package](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ).

## Tests
The package is tested with Python 3.8 and `biopython` version 1.78

## TODO
- [ ] create a `setup.py` file 
