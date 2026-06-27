"""
analysis.py
===========
Pipeline analisis sederhana bioinformatika

Langkah:
  1. Membaca file FASTA lalu menyimpan record ke dalam List
  2. Menghitung frekuensi nukleotida menggunakan Dictionary
  3. Menghitung GC Content setiap sekuens
  4. Mengurutkan sekuens berdasarkan GC Content (descending)
  5. Menampilkan 3 sekuens dengan GC Content tertinggi
  6. Visualisasi GC Content dengan matplotlib
  7. Menyimpan hasil analisis ke file CSV
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

matplotlib.use("Agg")

# Konfigurasi path
FASTA_FILE   = Path("data")   / "mitochondria_vertebrates.fasta"
OUTPUT_CSV   = Path("output") / "gc_analysis_results.csv"
OUTPUT_PLOT  = Path("output") / "gc_content_plot.png"

Path("output").mkdir(parents=True, exist_ok=True)


# STEP 1 – Membaca FASTA dan menyimpan ke dalam List
def load_fasta(filepath: Path) -> list[dict]:
    """
    Membaca file FASTA dan mengembalikan List of dict.
    Setiap elemen berisi: id, description, dan sequence (string).

    Struktur Data: LIST befungsi untuk menampung seluruh record sekuens
    sehingga bisa diiterasi, diurutkan, dan diiris (slice).
    """
    records: list[dict] = []   # ← LIST utama

    for record in SeqIO.parse(filepath, "fasta"):
        records.append({
            "id"         : record.id,
            "description": record.description,
            "sequence"   : str(record.seq).upper(),
        })

    return records


# STEP 2 – Menghitung frekuensi nukleotida dengan Dictionary
def nucleotide_frequency(sequence: str) -> dict[str, int]:
    """
    Menghitung frekuensi A, T, G, C dari sebuah sekuens.

    Struktur Data: DICTIONARY berfungsi sebagai key = simbol nukleotida, value = jumlah kemunculan.
    Cocok untuk pencarian O(1) dan representasi mapping.
    """
    freq: dict[str, int] = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
    for base in sequence:
        if base in freq:
            freq[base] += 1
        else:
            freq["N"] += 1   # karakter ambiguous dihitung sebagai N
    return freq


# STEP 3 – Menghitung GC Content
def gc_content(freq: dict[str, int], length: int) -> float:
    """Menghitung persentase GC Content dari dictionary frekuensi."""
    if length == 0:
        return 0.0
    gc = freq["G"] + freq["C"]
    return round((gc / length) * 100, 4)


# STEP 4 – Mengurutkan sekuens berdasarkan GC Content
def sort_by_gc(analyzed: list[dict]) -> list[dict]:
    """
    Mengurutkan list hasil analisis berdasarkan gc_content secara descending.

    Struktur Data: SORTING (built-in sort dengan key function) untuk Kompleksitas O(n log n) – Python Timsort.
    """
    return sorted(analyzed, key=lambda x: x["gc_content"], reverse=True)


# STEP 5 – Pipeline utama
def run_pipeline(fasta_path: Path) -> list[dict]:
    """Menjalankan seluruh pipeline analisis dan mengembalikan hasil terurut."""

    # 1. Load ke List
    print("=" * 65)
    print("  PIPELINE ANALISIS GC CONTENT – BIF1223 Mini Project")
    print("=" * 65)
    records: list[dict] = load_fasta(fasta_path)
    print(f"\n[1] {len(records)} sekuens berhasil dibaca dari '{fasta_path}'\n")

    # 2 & 3. Hitung frekuensi + GC Content; hasilkan List of dict lengkap
    analyzed: list[dict] = []
    for rec in records:
        seq   = rec["sequence"]
        freq  = nucleotide_frequency(seq)        # Dictionary frekuensi
        gc    = gc_content(freq, len(seq))
        label = rec["description"].split(",")[0]  # ambil bagian pertama deskripsi

        analyzed.append({
            "id"             : rec["id"],
            "label"          : label,
            "length"         : len(seq),
            "A"              : freq["A"],
            "T"              : freq["T"],
            "G"              : freq["G"],
            "C"              : freq["C"],
            "N"              : freq["N"],
            "gc_content"     : gc,
            "sequence_preview": seq[:50] + "..." if len(seq) > 50 else seq,
        })

    print("[2] Frekuensi nukleotida dihitung menggunakan Dictionary")
    print("[3] GC Content dihitung untuk setiap sekuens\n")

    # 4. Urutkan berdasarkan GC Content descending
    sorted_records: list[dict] = sort_by_gc(analyzed)
    for rank, rec in enumerate(sorted_records, start=1):
        rec["rank"] = rank

    print("[4] Sekuens diurutkan berdasarkan GC Content (tinggi → rendah)")
    print(f"    {'Rank':<5} {'ID':<18} {'GC Content':>11}  {'Panjang':>9}")
    print("    " + "-" * 48)
    for rec in sorted_records:
        print(f"    {rec['rank']:<5} {rec['id']:<18} {rec['gc_content']:>10.4f}%  {rec['length']:>8,} bp")

    return sorted_records


# STEP 6 – Menampilkan Top 3
def print_top3(sorted_records: list[dict]) -> None:
    """Menampilkan 3 sekuens dengan GC Content tertinggi."""
    top3: list[dict] = sorted_records[:3]   # slice List → top 3

    print("\n" + "=" * 65)
    print("  [5] TOP 3 SEKUENS – GC CONTENT TERTINGGI")
    print("=" * 65)
    for rec in top3:
        print(f"\n  Rank #{rec['rank']}")
        print(f"    ID          : {rec['id']}")
        print(f"    Deskripsi   : {rec['label']}")
        print(f"    Panjang     : {rec['length']:,} bp")
        print(f"    GC Content  : {rec['gc_content']:.4f}%")
        print(f"    Frekuensi   : A={rec['A']}  T={rec['T']}  G={rec['G']}  C={rec['C']}  N={rec['N']}")
        print(f"    Pratinjau   : {rec['sequence_preview']}")


# STEP 7 – Visualisasi dengan matplotlib
def plot_gc_content(sorted_records: list[dict], output_path: Path) -> None:
    """
    Membuat dua grafik:
      (a) Bar chart GC Content semua sekuens (terurut)
      (b) Bar chart komposisi nukleotida Top 3
    """
    ids    = [r["id"] for r in sorted_records]
    gc_vals = [r["gc_content"] for r in sorted_records]
    colors  = ["#e74c3c" if r["rank"] <= 3 else "#3498db" for r in sorted_records]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.patch.set_facecolor("#1a1a2e")

    # Panel kiri: GC Content seluruh sekuens
    ax1 = axes[0]
    ax1.set_facecolor("#16213e")
    bars = ax1.bar(ids, gc_vals, color=colors, edgecolor="white", linewidth=0.5)

    for bar, val in zip(bars, gc_vals):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.2,
            f"{val:.2f}%",
            ha="center", va="bottom", fontsize=7.5, color="white"
        )

    ax1.set_title("GC Content per Sekuens (Descending)", color="white", fontsize=13, pad=12)
    ax1.set_xlabel("Accession ID", color="white", fontsize=10)
    ax1.set_ylabel("GC Content (%)", color="white", fontsize=10)
    ax1.tick_params(colors="white", rotation=35)
    ax1.spines[:].set_color("#444")
    ax1.set_ylim(0, max(gc_vals) + 5)

    # legenda warna
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor="#e74c3c", label="Top 3"),
                       Patch(facecolor="#3498db", label="Others")]
    ax1.legend(handles=legend_elements, facecolor="#1a1a2e", labelcolor="white", fontsize=9)

    # Panel kanan: Komposisi nukleotida Top 3
    ax2 = axes[1]
    ax2.set_facecolor("#16213e")

    top3       = sorted_records[:3]
    top3_ids   = [r["id"] for r in top3]
    bases      = ["A", "T", "G", "C"]
    base_colors = ["#2ecc71", "#e67e22", "#e74c3c", "#9b59b6"]

    x      = range(len(top3))
    width  = 0.18

    for i, (base, col) in enumerate(zip(bases, base_colors)):
        freqs = [r[base] / r["length"] * 100 for r in top3]
        offset = (i - 1.5) * width
        rects = ax2.bar([xi + offset for xi in x], freqs, width,
                        label=base, color=col, edgecolor="white", linewidth=0.5)
        for rect, fv in zip(rects, freqs):
            ax2.text(rect.get_x() + rect.get_width() / 2,
                     rect.get_height() + 0.2,
                     f"{fv:.1f}%", ha="center", va="bottom",
                     fontsize=7, color="white")

    ax2.set_title("Komposisi Nukleotida – Top 3 Sekuens", color="white", fontsize=13, pad=12)
    ax2.set_xlabel("Accession ID", color="white", fontsize=10)
    ax2.set_ylabel("Frekuensi Relatif (%)", color="white", fontsize=10)
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(top3_ids, color="white")
    ax2.tick_params(colors="white")
    ax2.spines[:].set_color("#444")
    ax2.legend(facecolor="#1a1a2e", labelcolor="white", fontsize=9)

    plt.tight_layout(pad=2.5)
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="#1a1a2e")
    plt.close()
    print(f"\n[6] Visualisasi disimpan ke '{output_path}'")


# STEP 8 – Menyimpan ke CSV menggunakan pandas DataFrame
def save_csv(sorted_records: list[dict], output_path: Path) -> None:
    """
    Mengonversi List of dict ke pandas DataFrame lalu menyimpan ke CSV.

    Struktur Data: DATAFRAME untuk Representasi tabular dua dimensi, setiap kolom merupakan
        Series yang efisien untuk operasi vektorisasi.
    """
    columns = ["rank", "id", "label", "length",
               "A", "T", "G", "C", "N",
               "gc_content", "sequence_preview"]

    df = pd.DataFrame(sorted_records, columns=columns)
    df.to_csv(output_path, index=False, encoding="utf-8")
    print(f"[7] Hasil analisis disimpan ke '{output_path}'")
    print(f"\n    Preview tabel:\n")
    print(df[["rank", "id", "length", "gc_content"]].to_string(index=False))


# Entry point
if __name__ == "__main__":
    if not FASTA_FILE.exists():
        raise FileNotFoundError(
            f"File FASTA tidak ditemukan: '{FASTA_FILE}'\n"
            "Jalankan 'python download_ncbi.py' terlebih dahulu."
        )

    sorted_records = run_pipeline(FASTA_FILE)
    print_top3(sorted_records)
    plot_gc_content(sorted_records, OUTPUT_PLOT)
    save_csv(sorted_records, OUTPUT_CSV)

    print("\n" + "=" * 65)
    print("  Pipeline selesai. Semua output tersimpan di folder output/")
    print("=" * 65)
