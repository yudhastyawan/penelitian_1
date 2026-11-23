# Perbandingan forward aktif v6 dengan COMPOSTI

COMPOSTI dijalankan untuk model acuan sintetis Benchmark A dan B dengan API
`reflectivityCPP.compute_displ_Levin`. API itu juga dipanggil oleh skrip
reproduksi [Mreyen dkk. (2026)](https://academic.oup.com/gji/article/244/2/ggaf498/8373849)
dalam [arsip Zenodo 17807710](https://zenodo.org/records/17807710). Uji ini
menggunakan model, geometri, dan pencitraan proyek MFJ–SMVI. Kesamaan dengan
studi Mreyen ada pada solver dan API sintesis.

## Kondisi pengujian

- Gaya titik vertikal di permukaan, perpindahan vertikal, gelombang langsung dan
  refleksi ganda, tanpa gangguan sintetik. Spektrum sumber yang sama
  `(f/25 Hz)^2 exp(−(f/25 Hz)^2)` dikalikan satu kali. Parameter
  Vp, Vs, densitas, ketebalan, Qp, dan Qs disalin dari model acuan A/B.
- 56 offset yang sama (4–114 m) dan 31 frekuensi native 8–55 Hz untuk
  citra MFJ. Jejak waktu dibentuk dari grid 0,5 Hz dalam pita yang sama
  dengan rekaman dua detik. Kurva pada gambar dibagi puncak absolut jejak
  COMPOSTI di offset masing-masing. Metrik jejak memakai data panjang penuh
  dan skala perpindahan yang sama.
- Spektrum COMPOSTI dikonjugasi agar sesuai dengan konvensi waktu v6
  sebelum pencitraan. Keluaran v6 dibagi 2π untuk normalisasi gaya titik.
  Tidak ada penyelarasan amplitudo atau fase yang diestimasi dari data.
- Metrik adalah L2 relatif terhadap COMPOSTI. Kedua sisi memakai prosedur
  pencitraan MFJ aktif yang sama sebelum fusi aktif–pasif.

## Hasil

| Pengaturan v6 | Spektrum kompleks A | Spektrum kompleks B | Citra MFJ A | Citra MFJ B |
| --- | ---: | ---: | ---: | ---: |
| Bawaan | 6,404% | 8,939% | 15,813% | 22,794% |
| Pergeseran frekuensi kompleks dihilangkan, grid bawaan | 0,510% | 0,571% | 0,026% | 0,055% |
| Pergeseran dihilangkan, 8.192 titik log dan batas 4 kali | 0,514% | 0,578% | 0,024% | 0,055% |
| Pergeseran dihilangkan, rentang slowness dan taper COMPOSTI | 0,090% | 0,102% | 0,023% | 0,033% |

| Kondisi jejak | Benchmark | 4 m | 20 m | 60 m | 114 m |
| --- | --- | ---: | ---: | ---: | ---: |
| v6 bawaan | A | 1,71% | 6,76% | 15,88% | 24,75% |
| v6 bawaan | B | 2,71% | 11,15% | 25,03% | 35,05% |
| Pergeseran dihilangkan, grid bawaan | A | 0,757% | 0,019% | 0,014% | 0,039% |
| Pergeseran dihilangkan, grid bawaan | B | 0,683% | 0,025% | 0,034% | 0,074% |

Penghilangan pergeseran numerik `η = 0,005ω` pada grid bawaan yang sama
menjelaskan sebagian besar perbedaan dengan COMPOSTI untuk kanal aktif pada
dua model ini. Redaman material melalui Qp dan Qs tetap dipakai dalam kontrol.
Menyamakan rentang slowness dan taper lebih lanjut mendekatkan spektrum kompleks.
Ketika jumlah titik pada kontrol tersebut digandakan dari 4.096 menjadi 8.192
per segmen, perubahan L2 spektrum berbobot sumber berada di bawah 0,005% pada
kedua model.

Citra COMPOSTI masih memperlihatkan energi pada kecepatan fase tinggi.
Perbandingan mendukung keberadaan fitur tersebut pada dua solver, sementara
identitas moda dan kegunaannya untuk estimasi Vp memerlukan pengujian lanjutan.
Uji PyProp8 sebelumnya memakai pergeseran frekuensi yang sama di kedua solver,
sehingga uji itu mengisolasi kualitas solver ketika pergeseran disamakan.
Uji COMPOSTI di sini memakai pergeseran bawaan masing-masing sebagai pembanding
awal, lalu mengubahnya dalam kontrol. Hasil berlaku untuk model acuan aktif A
dan B. Kanal pasif, stabilitas di seluruh ruang prior, dan inversi ulang
memerlukan pemeriksaan tersendiri.

## Reproduksi

COMPOSTI menggunakan commit
`58ccb64749b7621d761e2fd108ff505be7ad3562` dari
[repositori penulis](https://github.com/mniskanen/composti). Commit tersebut
adalah versi upstream yang digunakan pada audit ini. Arsip reproduksi Mreyen
mengonfirmasi API dan menunjuk ke COMPOSTI tanpa pin commit.

Ekstensi upstream dikompilasi secara serial dengan Apple Clang di macOS ARM64.
Pada checkout terpisah `tmp/composti_external/upstream`, perubahan portabilitas
adalah mengganti panggilan Bessel orde 0 dan 1 dari `std::cyl_bessel_j` ke
`::jn`, serta menghapus opsi kompilasi dan penautan `-fopenmp`.
Kompilasi memakai Python 3.10, NumPy 1.26.4, Cython 3.0.12, dan eigency 5.0.1.1
dalam `tmp/composti_external/venv`. Saya menjalankan
`CC=/usr/bin/clang CXX=/usr/bin/clang++ python setup.py build_ext --inplace`
dari direktori `src/reflectivityCPP` dengan interpreter lingkungan tersebut.

Skrip demonstrasi bawaan COMPOSTI juga selesai dalam mode gambar noninteraktif.
Demonstrasi itu menampilkan pembanding SPECFEM3D secara visual tanpa ambang numerik.

Dari akar proyek, jalankan ulang perbandingan dengan:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python diagnostics_review_v6/benchmark_composti/compare.py
```

`summary.json` memuat metrik lengkap. Berkas `benchmark_*_comparison.npz`
memuat jejak, spektrum kompleks, citra dan grid. Panel berlabel “v6
diselaraskan” memakai kontrol dengan pergeseran dihilangkan dan rentang
slowness serta taper yang disamakan.
