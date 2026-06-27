# Mini Project – Struktur Data Bioinformatika (BIF1223)

**Integrasi Struktur Data untuk Pipeline Analisis Sederhana**  
IPB University · Semester Genap 2025/2026

---

## Dataset

| Atribut | Detail |
|---|---|
| Organisme | 10 spesies vertebrata (manusia, tikus, sapi, ayam, zebrafish, dll.) |
| Tipe data | Mitokondria lengkap (*complete mitochondrial genome*) |
| Sumber | NCBI RefSeq Nucleotide |
| Format | FASTA |
| Ukuran | ~16–20 kb per sekuens |

**Alasan pemilihan:** Mitokondria vertebrata memiliki ukuran kecil-sedang, tersedia di RefSeq, dan menunjukkan variasi GC Content yang bermakna antar spesies sehingga analisis perbandingan menjadi informatif.

---

## Struktur Folder

```
project/
├── download_ncbi.py   # Unduh data dari NCBI via Biopython Entrez
├── analysis.py        # Pipeline utama: baca, analisis, visualisasi, CSV
├── data/
│   └── mitochondria_vertebrates.fasta
├── output/
│   ├── gc_analysis_results.csv
│   └── gc_content_plot.png
├── requirements.txt
└── README.md
```

---

## Cara Menjalankan

```bash
# Langkah 1 – Unduh data dari NCBI
python download_ncbi.py

# Langkah 2 – Jalankan pipeline analisis
python analysis.py
```

Program akan menghasilkan:
- `output/gc_analysis_results.csv` – hasil analisis tabular
- `output/gc_content_plot.png`    – grafik GC Content

---

## Dependensi

| Library | Kegunaan |
|---|---|
| `biopython` | Entrez (unduh NCBI) + SeqIO (parsing FASTA) |
| `pandas` | DataFrame + ekspor CSV |
| `matplotlib` | Visualisasi bar chart |

---

## Kontak

Nama  : Akmal Fauzan Nurantho 

NIM   : G0401241060

Prodi : Bioinformatika, IPB University
