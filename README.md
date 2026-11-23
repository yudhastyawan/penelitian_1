# MFJ–SMVI: implementasi perbaikan lanjutan v6

Kode ini merupakan kerangka eksperimen sintetik. Hasil inversi dan intervalnya
belum layak disebut hasil penelitian tervalidasi. Evaluasi, alasan perubahan,
dan bukti pemeriksaan dicatat di `EVALUASI_DAN_PERBAIKAN.md`.

## Naskah dan artefak yang tersimpan

[Naskah lengkap edisi terintegrasi (PDF)](output/pdf/Naskah_Lengkap_Penelitian_MFJ_SMVI_v6_Edisi_Terintegrasi.pdf)
memuat hasil v6, perbandingan COMPOSTI, diagram, dan pseudocode. Sumber
[Markdown](output/pdf/Naskah_Lengkap_Penelitian_MFJ_SMVI_v6_Edisi_Terintegrasi.md)
dan aset gambar tersedia untuk membangun ulang dokumen dari akar proyek.

- `results_v6/` menyimpan observasi dan sampel sintetis A/B, gambar,
  manifest, serta snapshot kode untuk menelusuri hasil.
- `results_v6_validation/` dan `diagnostics_review_v6/` memuat bukti
  verifikasi numerik, pemeriksaan restart dan noise baru, serta
  [pembandingan COMPOSTI](diagnostics_review_v6/benchmark_composti/RINGKASAN.md).
- Arsip v1–v5, log mesin, dan PDF revisi lama tersimpan di workspace lokal.
  Rujukan ke arsip dalam naskah menggunakan jalur workspace tersebut.
  Manifest hash mempertahankan catatan run asli, termasuk log lokal.
  Dua tes historis otomatis dilewati bila artefak v4 belum tersedia.

## Menjalankan

Gunakan `.venv`; lingkungan `env/` lama mewarisi SciPy native yang bermasalah.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q
OPENBLAS_NUM_THREADS=1 .venv/bin/python validate_reassessment.py --output-dir results_v6_validation
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u run_benchmark.py --quick --output-dir results_smoke_baru
```

Folder output harus baru/kosong. Driver menolak menimpa eksperimen lama.
`--quick` memeriksa integrasi; tetap memakai frekuensi native, tetapi anggaran
inversinya kecil. Menghilangkan `--quick` memperbesar anggaran dan resolusi;
keberhasilan eksekusi tetap tidak otomatis membuktikan konvergensi.

Pemeriksaan ulang dua benchmark yang disimpan di `results_v6/` memakai:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u run_benchmark.py --quick --iterations 60 --components 4 --mc-samples 4 --refinement-steps 30 --benchmarks-only --parallel-benchmarks --output-dir results_v6
```

Perintah tersebut mencakup dua inversi, evaluasi geometri/kernel, gelombang,
spektrum, sampel, dan diagnostik kuadratur. `--benchmarks-only` menyatakan secara
eksplisit bahwa inversi ablation, MCMC, dan stress tambahan tidak dijalankan.
Untuk menjalankan seluruh rangkaian, hilangkan opsi tersebut dan gunakan folder
output baru. Hasil tidak boleh dibandingkan sebagai klaim keunggulan metode
ketika pembanding belum konvergen atau tidak memiliki anggaran yang setara.

## Laporan dan verifikasi

```sh
.venv/bin/python generate_pdf_manuscript.py --results-dir results_v6 --diagnostic-report
```

Ini hanya menulis Markdown diagnostik di folder hasil. Penerbitan naskah/PDF
penelitian memerlukan pemeriksaan `research_quality.py` yang lolos; uji cepat
atau eksperimen tanpa validasi posterior akan ditolak. Generator juga memeriksa
manifest dan hash kode/dependensi. Pandoc/pdflatex hanya diperlukan untuk PDF.

`validate_reassessment.py` membandingkan respons compliance P–SV dengan
propagator independen PyProp8 1.1.5 pada kedua model. Pemeriksaan ini memiliki
cakupan frekuensi dan bilangan gelombang terbatas; bukan validasi gelombang
lapangan. `disba` menyediakan referensi dispersi Rayleigh yang terpisah.

Backend diferensiasi otomatis menghitung model numerik yang sama dengan Numba.
Tes memeriksa spektrum, nilai target, dan gradien finite difference, termasuk
kanal aktif/pasif serta mask kecepatan. Tidak ada cabang moda yang dipaksakan,
pemindahan posterior pembanding, atau inisialisasi dengan parameter sejati.

## Perbaikan lanjutan yang diterapkan

- NaN/Inf, massa transport negatif/tanpa total positif, observasi kompleks dan grid
  yang tidak valid ditolak. Pemeriksaan terjadi sebelum loop transport/forward.
- Skala Wasserstein dinyatakan dalam m/s. Default likelihood ditentukan satu kali
  dari grid parent penuh; pemotongan Rayleigh-window mempertahankan skala yang
  sama. Untuk membandingkan operator dengan grid parent berbeda, set
  `SMVIConfig.wasserstein_scale_mps` secara eksplisit ke nilai yang sama.
- Bobot trapesium dihitung pada setiap pita frekuensi kanal yang dipertahankan,
  termasuk endpoint-nya. Kanal pasif 1–12 Hz mengintegrasikan 11 Hz.
- Prediksi menggunakan draw yang seimbang per komponen, lalu memakai bobot
  campuran untuk menggabungkannya. Default minimum 128 draw, maksimum 512.
  Jumlah draw, bootstrap MC, standard error dan ketercapaian target presisi
  disimpan. Batas anggaran tidak otomatis dianggap sebagai konvergensi.
- Pembanding prior dan peringatan hasil tersimpan dalam JSON eksperimen dan
  selalu dirender ulang. Simpan catatan manual di
  `DIAGNOSTIC_EDITORIAL_NOTES.md` dalam folder hasil; generator membacanya
  tanpa mengubah file catatan tersebut.

Perbandingan terhadap prior adalah pemeriksaan diagnostik satu titik; bukan
referensi posterior yang terkalibrasi. Galat Monte Carlo adalah ketelitian
estimasi dari distribusi campuran yang dipasang, berbeda dari spread respons
model, ketidakpastian noise, dan validasi posterior.

## Batas interpretasi

- Noise pasif adalah perturbasi EGF fenomenologis, bukan simulasi interferometri.
- Densitas dan Q dianggap diketahui. Sigma likelihood komposit belum dikalibrasi.
- SMVI memakai bobot sama dan kernel seluruh parameter komponen. Uji analitik
  dua moda tidak membuktikan konvergensi posterior pada masalah seismik.
- Grid frekuensi native menggantikan interpolasi spektrum yang terlalu jarang.
  Kedua kanal memakai node yang sama dalam overlap untuk mencegah alternasi bobot.
  Kepadatan grid gambar tidak memperbaiki aliasing akuisisi 48 kanal, dx 2 m.
- Spektrum rata-rata prediksi disimpan terpisah dari spektrum pada model mean.
- Referensi pembanding adalah versi tanpa noise dari transformasi/geometri
  masing-masing; metrik ini mengukur ketahanan noise, bukan peringkat metode.

## Referensi fungsi dan persamaan

- [Indeks fungsi](docs/REFERENSI_FUNGSI.md) menghubungkan peran setiap fungsi
  dengan definisi kode, persamaan yang diterapkan, dan rujukannya.
- [Persamaan sesuai implementasi](docs/PERSAMAAN_KODE.md) memuat formula untuk
  sintesis P–SV, MFJ, fusi, target inversi, SMVI, dan verifikasi.

Ringkasan formulasi tersedia di [MATHEMATICAL_FORMULATION.md](MATHEMATICAL_FORMULATION.md).
Tujuan penelitian dan rencana validasi ada di [proposal](PHD_RESEARCH_PROPOSAL.md).
