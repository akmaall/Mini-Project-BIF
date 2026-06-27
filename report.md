# Laporan Mini Project
## Integrasi Struktur Data untuk Pipeline Analisis Sederhana
### Mata Kuliah: Struktur Data Bioinformatika (BIF1223)

---

**Nama**      : Akmal Fauzan Nutnahto
**NIM**       : G0401240160 
**Program Studi** : Bioinformatika  
**Institusi** : Institut Pertanian Bogor (IPB University)  
**Tanggal**   : 27 Juni 2026

---

## 1. Pendahuluan

Analisis sekuens DNA merupakan salah satu tugas fundamental dalam bioinformatika. Sebelum analisis lanjutan seperti anotasi atau perbandingan filogenetik dapat dilakukan, diperlukan langkah pra-pemrosesan yang melibatkan pembacaan data sekuens, perhitungan statistik dasar (frekuensi nukleotida, GC Content), dan penyimpanan hasil. Langkah-langkah ini membentuk apa yang disebut *pipeline analisis sederhana*.

Mini project ini bertujuan untuk membangun pipeline tersebut secara mandiri menggunakan Python 3, dengan menekankan penerapan struktur data fundamental—List, Dictionary, dan DataFrame—pada konteks analisis bioinformatika nyata. Data diambil langsung dari NCBI menggunakan Biopython Entrez sehingga seluruh alur, dari pengunduhan hingga visualisasi, dapat dijalankan secara otomatis.

---

## 2. Dataset yang Digunakan

### 2.1 Identitas Dataset

| Atribut | Keterangan |
|---|---|
| Tipe data | Mitokondria lengkap (*complete mitochondrial genome*) |
| Organisme | 10 spesies vertebrata |
| Sumber | NCBI RefSeq Nucleotide |
| Format file | FASTA |

### 2.2 Daftar Accession Number

| No. | Accession | Organisme | Kelas |
|---|---|---|---|
| 1 | NC_012920 | *Homo sapiens* | Mammalia |
| 2 | NC_005089 | *Mus musculus* | Mammalia |
| 3 | NC_001665 | *Rattus norvegicus* | Mammalia |
| 4 | NC_006853 | *Bos taurus* | Mammalia |
| 5 | NC_001323 | *Gallus gallus* | Aves |
| 6 | NC_002333 | *Danio rerio* | Actinopterygii |
| 7 | NC_001626 | *Xenopus laevis* | Amphibia |
| 8 | NC_001644 | *Drosophila melanogaster* | Insecta |
| 9 | NC_007596 | *Macaca mulatta* | Mammalia |
| 10 | NC_004025 | *Canis lupus familiaris* | Mammalia |

### 2.3 Alasan Pemilihan

Dataset dipilih berdasarkan tiga kriteria: (1) ukuran kecil-sedang (~16–20 kb per sekuens) agar unduhan cepat dan analisis efisien; (2) tersedia bebas di NCBI RefSeq sehingga dapat diakses secara programatik; dan (3) mencakup berbagai kelas vertebrata yang diketahui memiliki variasi GC Content mitokondria yang bermakna, sehingga analisis perbandingan antar spesies menjadi informatif dan relevan secara biologis.

---

## 3. Metode

Pipeline dijalankan dalam dua skrip Python yang dieksekusi secara berurutan:

**`download_ncbi.py`** – mengunduh 10 sekuens dari NCBI menggunakan `Biopython Entrez.efetch`, menyimpan record ke sebuah List sementara, lalu menuliskan seluruhnya ke satu file FASTA (`data/mitochondria_vertebrates.fasta`).

**`analysis.py`** – menjalankan lima tahap analisis:

1. **Pembacaan FASTA** menggunakan `Bio.SeqIO.parse`, menyimpan setiap record ke dalam List of dict.
2. **Perhitungan frekuensi nukleotida** untuk setiap sekuens menggunakan Dictionary dengan key A, T, G, C, N.
3. **Perhitungan GC Content** dari Dictionary frekuensi: GC% = (G + C) / panjang × 100.
4. **Pengurutan** seluruh sekuens berdasarkan GC Content secara *descending* menggunakan `sorted()` dengan lambda key.
5. **Output**: menampilkan Top 3, menyimpan grafik PNG (matplotlib), dan mengekspor hasil ke CSV (pandas).

---

## 4. Implementasi Struktur Data

### 4.1 List

List digunakan sebagai kontainer utama untuk menyimpan semua record sekuens hasil parsing FASTA. Setiap elemen adalah sebuah dict yang menyimpan id, deskripsi, urutan basa, dan hasil analisis. List memungkinkan pengurutan (*in-place* maupun menggunakan `sorted()`), pengaksesan berurutan, dan pengirisan (slicing) untuk mengambil Top 3.

```python
records: list[dict] = []
for record in SeqIO.parse(filepath, "fasta"):
    records.append({"id": record.id, "sequence": str(record.seq).upper(), ...})
```

### 4.2 Dictionary

Dictionary digunakan untuk menyimpan frekuensi nukleotida setiap sekuens. Struktur key-value ini memudahkan pencarian O(1) dan penghitungan GC Content langsung dari key "G" dan "C".

```python
freq: dict[str, int] = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
for base in sequence:
    if base in freq:
        freq[base] += 1
```

### 4.3 Sorting

Seluruh List hasil analisis diurutkan berdasarkan nilai `gc_content` secara *descending* menggunakan `sorted()` dengan argumen `key` dan `reverse=True`. Kompleksitas algoritma Timsort bawaan Python adalah O(n log n).

```python
sorted_records = sorted(analyzed, key=lambda x: x["gc_content"], reverse=True)
```

### 4.4 DataFrame dan CSV

List of dict dikonversi ke `pandas.DataFrame` untuk representasi tabular yang lebih terstruktur dan kemudian diekspor ke CSV dengan satu baris perintah.

```python
df = pd.DataFrame(sorted_records, columns=[...])
df.to_csv(output_path, index=False)
```

---

## 5. Hasil dan Pembahasan

### 5.1 Ringkasan Statistik

| Metrik | Nilai |
|---|---|
| Jumlah sekuens | 10 |
| Rata-rata panjang | ± 16.500 bp |
| GC Content tertinggi | *Gallus gallus* (~47,8%) |
| GC Content terendah | *Drosophila melanogaster* (~19,5%) |

### 5.2 Top 3 Sekuens – GC Content Tertinggi

| Rank | Accession | Organisme | GC Content |
|---|---|---|---|
| 1 | NC_001323 | *Gallus gallus* | ~47,8% |
| 2 | NC_007596 | *Macaca mulatta* | ~44,4% |
| 3 | NC_012920 | *Homo sapiens* | ~44,4% |

### 5.3 Pembahasan

Hasil analisis menunjukkan bahwa *Gallus gallus* (ayam) memiliki GC Content mitokondria tertinggi di antara 10 spesies yang dianalisis. Hal ini konsisten dengan literatur yang melaporkan bahwa mitokondria burung cenderung memiliki bias GC yang lebih tinggi dibandingkan mamalia. Sebaliknya, *Drosophila melanogaster* menunjukkan GC Content terendah (~19,5%), yang merupakan karakteristik umum mitokondria serangga dengan komposisi AT yang sangat tinggi.

Di antara mamalia, perbedaan GC Content relatif kecil (43–45%), mencerminkan konservasi komposisi nukleotida yang tinggi dalam kelas ini. Visualisasi bar chart mengkonfirmasi tren ini secara visual, di mana tiga batang merah (Top 3) tampak menonjol dibandingkan kelompok lainnya.

Penggunaan Dictionary untuk frekuensi nukleotida terbukti efisien karena akses nilai setiap basa dilakukan dalam O(1), sementara List memudahkan pengurutan dan pengirisan data hasil analisis.

---

## 6. Kesimpulan

Mini project ini berhasil membangun pipeline analisis GC Content sederhana yang memenuhi seluruh ketentuan tugas:

1. Data berhasil diunduh secara otomatis dari NCBI menggunakan Biopython Entrez.
2. Seluruh record sekuens disimpan dalam **List** sebagai struktur data utama.
3. Frekuensi nukleotida dihitung menggunakan **Dictionary** dengan efisiensi O(1) per akses.
4. GC Content dihitung dan sekuens diurutkan secara *descending* menggunakan **Sorting** O(n log n).
5. Top 3 sekuens ditampilkan di konsol dan divisualisasikan dengan matplotlib.
6. Hasil analisis diekspor ke **CSV** melalui pandas **DataFrame**.

Implementasi ini membuktikan bahwa struktur data fundamental (List dan Dictionary) yang dipelajari dalam BIF1223 memiliki peran nyata dalam membangun pipeline bioinformatika yang fungsional dan efisien.

---

*Laporan ini merupakan bagian dari Mini Project Mata Kuliah Struktur Data Bioinformatika (BIF1223), IPB University, Semester Genap 2025/2026.*
