# 🔬 Inhibitory KPC-2 – analiza wstępna

W projekcie skupiłyśmy się na enzymie **KPC-2** 🧪, odpowiedzialnym za rozkład β-laktamów i rozwój oporności na antybiotyki. Jego działanie stanowi poważne wyzwanie kliniczne, dlatego naszym celem było znalezienie potencjalnych **inhibitorów** tego enzymu.

## ✅ Etapy projektu:

- Stworzenie biblioteki potencjalnych inhibitorów (PubChem, ChEMBL) 
- Dokowanie molekularne wybranych związków z enzymem KPC-2  
- Analiza właściwości związków za pomocą narzędzi ADMET AI  

🔁 Projekt ten stanowi punkt wyjścia do dalszych analiz, takich jak symulacje dynamiki molekularnej.

## ⚙️ Uruchamianie
Aby uruchomić program, użyj pliku **`KPC-2-project.ipynb`**.

## 📦 Wymagane biblioteki (Python)

- `rdkit`  
- `biopython`  
- `numpy`  
- `matplotlib`  
- `py3Dmol`  
- `openbabel`  
- `json`, `os`, `glob`, `csv` *(standardowe moduły Pythona)*

**Uwaga**: W repozytorium znajdują się również pliki wynikowe, takie jak `merged_affinities.txt`, które umożliwiają odtworzenie kluczowych kroków bez potrzeby generowania tysięcy plików ponownie.
