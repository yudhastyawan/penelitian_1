# Proposal Penelitian: Inversi Spektrum Aktif-Pasif MFJ-SMVI

Status revisi: koreksi metodologi setelah audit implementasi, 21 November 2026.
Dokumen ini menyatakan sasaran penelitian; hasil validasi harus dirujuk dari
`results_v6/benchmark_verification_summary.json` beserta status konvergensinya.
Hasil v6 masih diagnostik; kelayakan publikasi diperiksa terpisah. Angka dan klaim
versi sebelumnya disimpan dalam folder arsip v3 dan v4.

## Masalah dan tujuan

Penelitian menguji apakah pencocokan spektrum frekuensi-kecepatan fase tanpa
picking dapat mengurangi kesalahan identifikasi moda pada lapisan lunak terjebak,
serta menambah informasi tentang Vp dan kedalaman batas jenuh melalui bagian
spektrum berkecepatan tinggi. Integrasi sumber aktif dan pasif diuji pada domain
pengolahan yang sama dengan ketidakpastian parameter yang dinyatakan eksplisit.

## Hipotesis yang harus diuji

1. Penggabungan pita aktif-pasif menambah resolusi parameter dibandingkan masing-
   masing kanal, jika forward operator dan noise model konsisten.
2. Informasi spektral di luar jendela Rayleigh dapat mempersempit ketidakpastian Vp.
3. Campuran Gaussian berkorelasi dapat mengaproksimasi posterior yang tidak cukup
   diwakili satu Gaussian, dengan biaya dan akurasi yang diukur secara adil.
4. Geometri nested 2D dapat menawarkan trade-off cakupan radial dan azimut yang
   berguna. Tidak ada konfigurasi yang dinyatakan terbaik sebelum hasil pengujian.

## Metode dan eksperimen

Model utama menggunakan rekursi compliance P-SV dan kuadratur Sommerfeld dengan
pengurangan asimtot statik. Validasi dimulai dari solusi half-space analitik,
kesetaraan pembagian lapisan homogen, sensitivitas terhadap compliance, serta
konvergensi kuadratur. Referensi Rayleigh independen diperoleh dari disba.

Benchmark A memiliki LVL pada 4-10 m dengan Vs 150 m/s. Benchmark B memiliki
batas muka air tanah pada 3.5 m dengan Vp 550 menjadi 1650 m/s. Geologi berlapis
1D dan densitas/Q yang diketahui adalah asumsi eksperimen, bukan kemampuan
estimasi yang telah dibuktikan.

Empat desain memakai 48 kanal aktif dan 25 stasiun pasif: linear, L-array,
cincin konsentris, dan nested spiral. Setiap desain menggunakan koordinatnya
sendiri dan jarak pasangan aktual. Metode imaging konvensional menjalankan
kernel sebenarnya pada gather, bukan citra yang diberi artefak buatan.

Ablasi wajib menjalankan active-only, passive-only, joint, Rayleigh-window,
Gaussian tunggal, picked-curve, dan MCMC secara langsung. Sampel hasil tidak
boleh digeser untuk menunjukkan mekanisme yang diharapkan. Anggaran komputasi,
seed, model prior, chain, dan sampel disimpan.

Stress-test melemahkan sumber pada domain frekuensi-offset sebelum pencitraan,
menambahkan noise dengan SNR terukur, dan menghitung ulang bobot fusi sebenarnya.
Label SNR, kopling, serta keterbatasan bandwidth tidak boleh ditentukan hanya
melalui judul gambar.

## Tahapan pembuktian

1. Regresi numerik dan eksperimen cepat untuk memeriksa seluruh alur.
2. Studi konvergensi kuadratur aktif dan pasif, sensitivitas hyperparameter,
   multi-start, multi-seed, serta kalibrasi posterior pada dataset berulang.
3. MCMC dengan diagnostik dan eksplorasi moda yang memadai; bandingkan biaya per
   akurasi, bukan jumlah iterasi semata. Tidak ada target speedup yang dianggap hasil.
4. Validasi data lapangan terhadap logging dan piezometer, dengan pengujian efek
   heterogenitas lateral dan ketidakpastian parameter nuisance.
5. Publikasi metode dan perangkat lunak setelah bukti tersebut tersedia.

## Pedoman akuisisi sementara

Geometri nested saat ini memiliki 9 sensor inti (spasi inline 10.5 m dengan
stagger melintang) serta 16 sensor spiral radius 52-320 m. Enam tembakan aktif
berada 4, 10, dan 20 m di luar masing-masing ujung bentangan 94 m. Pasangan
pasif tidak dibuang ke bin radial kasar secara default.

Pemilihan instrumen, sumber, durasi, dan respons frekuensi harus didasarkan
pada noise floor, hasil uji lapangan, serta spesifikasi terverifikasi. Klaim
spasi inti menjamin seluruh pita 8-14 Hz atau pelat tumbuk tertentu menjamin
energi frekuensi tertentu belum dianggap terbukti oleh kode sintetik.
