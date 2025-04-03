from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors

smiles_list = [
    "NC(=O)[C@@H]1CC[C@@H]2CN1C(=O)N2OS(=O)(=O)O",
    "C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(SCCNC=N)C[C@H]12",
    "C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(S[C@@H]3CN[C@H](C(=O)N(C)C)C3)[C@H](C)[C@H]12"
]

# Konwersja SMILES do RDKit Mol
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles)]

# Generowanie fragmentów BRICS
fragments = set()
for mol in mols:
    for frag in BRICS.BRICSDecompose(mol):
        fragments.add(frag)

# *Tworzenie nowych związków*
new_mols = []
for frag1 in fragments:
    for frag2 in fragments:
        for frag3 in fragments:
            mol1 = Chem.MolFromSmiles(frag1)
            mol2 = Chem.MolFromSmiles(frag2)
            mol3 = Chem.MolFromSmiles(frag3)

            if mol1 and mol2 and mol3:
                # *Łączenie trzech fragmentów*
                new_mol = Chem.CombineMols(mol1, mol2)
                new_mol = Chem.CombineMols(new_mol, mol3)

                # *Filtracja (np. reguła Lipińskiego)*
                mw = Descriptors.MolWt(new_mol)
                logp = Descriptors.MolLogP(new_mol)
                h_donors = Descriptors.NumHDonors(new_mol)
                h_acceptors = Descriptors.NumHAcceptors(new_mol)

                if mw < 500 and logp < 5 and h_donors <= 5 and h_acceptors <= 10:
                    new_mols.append(new_mol)

# *Zapisywanie do SDF*
w = Chem.SDWriter("generated_ligands.sdf")
for mol in new_mols:
    if mol:
        w.write(mol)
w.close()

print(f"Wygenerowano {len(new_mols)} nowych związków po filtracji.")