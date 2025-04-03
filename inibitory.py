from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, rdFMCS

# Wczytanie pliku SDF z ligandami (np. inhibitorem KPC-2)
suppl = Chem.SDMolSupplier("polaczone_bialka_inh.sdf")
mols = [mol for mol in suppl if mol is not None]

# Generowanie fragmentów za pomocą BRICS (strategiczne cięcie molekuł)
fragments = set()
for mol in mols:
    for frag in BRICS.BRICSDecompose(mol):
        fragments.add(frag)

# Łączenie fragmentów w nowe cząsteczki
new_mols = []
for frag1 in fragments:
    for frag2 in fragments:
        new_mol = Chem.CombineMols(Chem.MolFromSmiles(frag1), Chem.MolFromSmiles(frag2))
        new_mols.append(new_mol)

# Zapisywanie wyników do pliku SDF
w = Chem.SDWriter("generated_ligands.sdf")
for mol in new_mols:
    if mol is not None:
        w.write(mol)
w.close()

print(f"Wygenerowano {len(new_mols)} nowych związków.")