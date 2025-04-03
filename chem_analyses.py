from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator
import numpy as np

# Wczytaj bibliotekę białek (SDF)
protein_sdf = "polaczone_bialka.sdf"  # <- podmień na swoją nazwę pliku
proteins = Chem.SDMolSupplier(protein_sdf)


# Wczytaj inhibitory (SMILES lub SDF)
inhibitor_smiles = [
    "NC(=O)[C@@H]1CC[C@@H]2CN1C(=O)N2OS(=O)(=O)O",
    "C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(SCCNC=N)C[C@H]12",
    "C[C@@H](O)[C@H]1C(=O)N2C(C(=O)O)=C(S[C@@H]3CN[C@H](C(=O)N(C)C)C3)[C@H](C)[C@H]12"
]
# rdFingerprintGenerator.GetMAC # użyć tego
inhibitors = [Chem.MolFromSmiles(smiles) for smiles in inhibitor_smiles]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
inhibitor_fps = [gen.GetFingerprint(mol) for mol in inhibitors]

# Wczytaj fingerprinty białek
protein_fps = [gen.GetFingerprint(mol) for mol in proteins if mol is not None]


# Oblicz średnie podobieństwo dla każdego inhibitora
similarity_scores = []

for i, inhibitor_fp in enumerate(inhibitor_fps):
    similarities = [DataStructs.TanimotoSimilarity(inhibitor_fp, protein_fp) for protein_fp in protein_fps]
    avg_similarity = np.mean(similarities)  # Średnie podobieństwo do wszystkich białek
    similarity_scores.append((inhibitor_smiles[i], avg_similarity))

# Posortuj według średniego podobieństwa
similarity_scores.sort(key=lambda x: x[1], reverse=True)

# Wypisz wynik
print("Najbardziej podobny inhibitor do całej biblioteki białek:")
for i, (inhibitor, avg_sim) in enumerate(similarity_scores):
    print(f"{i+1}. Inhibitor: {inhibitor}, Średnie podobieństwo: {avg_sim:.3f}")

# Wybierz najlepszy inhibitor
best_inhibitor = similarity_scores[0]
print(f"\nNajlepiej dopasowany inhibitor: {best_inhibitor[0]} (Podobieństwo: {best_inhibitor[1]:.3f})")