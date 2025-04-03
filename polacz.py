from rdkit import Chem
from rdkit.Chem import SDWriter, rdFingerprintGenerator

# Generator fingerprintów Morgana
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)

# Wczytaj pliki SDF
sdf_files = ["plik1.sdf", "plik2.sdf", "plik3.sdf"]
all_molecules = []
seen_fps = set()

# Wczytaj cząsteczki z plików SDF
for sdf_file in sdf_files:
    supplier = Chem.SDMolSupplier(sdf_file)
    for mol in supplier:
        if mol is not None:
            fp = gen.GetFingerprint(mol)  # Generowanie fingerprintu Morgana
            fp_bytes = fp.ToBinary()  # Konwersja fingerprintu do bajtów (porównywanie)

            # Sprawdzanie, czy fingerprint już występuje
            if fp_bytes not in seen_fps:
                seen_fps.add(fp_bytes)
                all_molecules.append(mol)

# Zapisywanie unikalnych cząsteczek do nowego pliku SDF
output_writer = SDWriter("polaczone_bialka.sdf")
for mol in all_molecules:
    output_writer.write(mol)
output_writer.close()

print(f"Połączono i zapisano {len(all_molecules)} unikalnych cząsteczek do pliku 'polaczone_bialka.sdf'")
