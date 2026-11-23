# Evaluasi ulang menyeluruh dan perbaikan v5

Tanggal: 21 November 2026. Dokumen ini memisahkan kerusakan implementasi,
pemeriksaan yang sudah dilakukan, dan pertanyaan ilmiah yang belum terjawab.

## Kesalahan pekerjaan sebelumnya

Saya terlalu cepat menganggap keluaran uji singkat sebagai perbaikan yang sudah
memadai. Lolos unit test dan berhasil membuat PDF tidak menjamin spektrum,
inversi, atau interval ketidakpastian dapat dipercaya. Penggantian besar pada
v4 juga tidak disertai bukti independen dan pemeriksaan hasil yang cukup.
Karena itu naskah/PDF v4 ditarik dari keluaran aktif, bukan diperhalus lagi lalu
dinyatakan selesai. Salinan kode, dokumen, tes, dan hasilnya dipertahankan di
`archive_v4_before_full_reassessment/`.

## Temuan dan tindakan

| Bagian | Masalah yang ditemukan | Perbaikan yang diterapkan |
|:--|:--|:--|
| Spektrum | Sintesis aktif hanya 9 frekuensi dan pasif 8 frekuensi pada uji v4; interpolasi ke grid lain dapat mengubah atau menghilangkan ridge | Fusi memakai irisan native; pada bagian dalam overlap kedua kanal harus tersedia di frekuensi yang sama, sehingga tidak ada alternasi kanal akibat grid berbeda. Grid benchmark dan diagnostik diperapat secara eksplisit |
| Bobot frekuensi | Setiap baris memiliki bobot sama meskipun jarak frekuensinya berbeda | Rata-rata misfit menggunakan bobot trapesium dalam Hz; mask ablation mengalikan bobot tersebut |
| Gelombang | Seleksi setiap indeks ke-k pada jarak pasangan tidak teratur mengabaikan ujung apertur; halfwidth pasif 18.5 m membuat jejak bertumpuk | Pemilihan mendekati radius merata, kedua ujung disertakan, halfwidth 0.35 kali celah minimum, jendela waktu diperpanjang |
| Konsistensi sumber | Bentuk sumber pasif untuk FFT dan sintesis frekuensi berbeda | Keduanya memakai bentuk sumber pasif yang sama; realisasi noise time gather tetap ditandai sebagai ilustrasi terpisah |
| Akuisisi | Grid gambar halus mudah disalahartikan sebagai bebas aliasing | Garis Nyquist spasial sederhana ditampilkan pada spektrum aktif; 48 kanal dan dx 2 m tetap, tanpa perubahan instrumen tersembunyi |
| SMVI | Kernel hanya pada mean; attraction dinormalisasi per kolom; bobot komponen memiliki dinamika tambahan | Kernel pada seluruh parameter bebas komponen, faktor 1/K, repulsi pada ruang yang sama, bobot komponen sama |
| Gradien | Finite difference dengan langkah relatif besar mahal untuk grid rapat | Backend AD CPU float64 menghitung rekursi dan kuadratur yang sama; spektrum, target dan gradien dibandingkan dengan Numba/finite difference |
| Inisialisasi | Uji singkat dimulai terlalu jauh dari distribusi yang memberi fit baik | Opsi refinemen L-BFGS dari kandidat Latin hypercube berdasarkan data/prior; tidak memakai parameter sejati |
| Prediksi | Spektrum pada mean parameter mudah disebut sebagai mean prediksi | Keduanya dihitung dan disimpan terpisah; gambar memakai mean prediksi ketika tersedia, dengan label yang sesuai |
| Pembanding transformasi | Plane-wave dan Bessel dinilai terhadap MFJ bersih, sehingga referensi menguntungkan MFJ | Masing-masing dinilai terhadap transformasi bersihnya sendiri pada geometri yang sama. Metrik menjadi konsistensi noise, bukan peringkat metode |
| Solver eksternal | Adapter PyProp8 salah mengurai nilai return lalu menyembunyikan exception | Adapter diperbaiki dengan unit SI, return shape benar dan kegagalan eksplisit. Adapter elastic tidak disamakan dengan solver ber-Q tanpa mencocokkan moduli |
| Integritas hasil | Output lama berpotensi bercampur ketika memakai folder yang sama | Driver menolak output yang sudah berisi; hash kode, tes, dependensi dan konfigurasi disimpan; perubahan kode saat run ditolak |
| Laporan | Status “completed” saja cukup untuk membuat naskah penelitian | Gerbang penerbitan terpisah memeriksa mode, validasi resolusi, posterior independen, kalibrasi dan fisika. Laporan diagnostik memiliki jalur terpisah |

## Bukti yang sudah diperoleh

Tes otomatis: **31 pemeriksaan lolos**, termasuk:

- Adapter waveform PyProp8 dijalankan langsung (output 64 × 2, seluruhnya finite).
  Uji waktu tempuh homogen juga memeriksa arah propagasi/konvensi FFT.

- Respons compliance P–SV berlapis dibandingkan dengan propagator PyProp8 1.1.5
  pada dua model, tiga frekuensi dan enam bilangan gelombang per model.
  Galat relatif maksimum A: `6.0685e-7`; B: `3.7046e-6`.
  Pemeriksaan memakai moduli kompleks yang sama dan sumber sangat dekat permukaan.
  Cakupannya compliance di titik uji, bukan seluruh waveform atau data lapangan.
- Spektrum dan log posterior backend AD cocok dengan Numba; gradien diverifikasi
  melalui finite difference kecil. Kanal fused, aktif, pasif, dan mask kecepatan
  diperiksa, bukan hanya satu kasus.
- SMVI empat komponen diuji pada target dua moda yang diketahui. Pemeriksaan
  meliputi mean, simpangan baku, fraksi sampel kedua moda dan bobot sama.
- Kernel diuji pada komponen yang mean-nya sama tetapi kovariansnya berbeda,
  sehingga tes gagal jika pembaruan kembali hanya memakai mean.
- Fusi tidak mengarang frekuensi/irisan baru, wiggle mencakup kedua ujung apertur
  tanpa overlap, serta laporan uji cepat ditolak oleh gerbang publikasi.
- Pemeriksaan terdahulu tentang W2 diskrit, transformasi prior/Jacobian,
  konsistensi kovarians dan sampling, MCMC tanpa clipping, SNR pada gather,
  pembagian lapisan homogen dan kuadratur tetap dijalankan.

Bukti fisika numerik tersimpan di
`results_v5_validation/independent_physics_validation.json`.
Spektrum diagnostik rapat memakai frekuensi native 0.25 Hz dan grid kecepatan
2 m/s di `diagnostics_review_v5/`. Gambar tersebut berasal dari sintesis langsung,
bukan pembesaran atau penghalusan citra lama. Manifest menyimpan hash dan
menandai bahwa ini diagnosis forward, bukan hasil inversi.

## Eksperimen ulang

Konfigurasi run dua benchmark: 4 komponen, 60 iterasi, 4 sampel Monte Carlo per
komponen per iterasi, refinemen awal maksimum 30 langkah, seed 2026/2027.
Grid aktif native berselang 2 Hz di luar overlap, 0.5 Hz pada overlap, serta
node akhir 55 Hz; grid pasif 0.5 Hz, dan kecepatan 4 m/s.
Kedua benchmark dijalankan pada proses numerik terpisah. Mode tetap ditandai diagnostik. Data, sampel dan manifest berada di `results_v5/`. Snapshot kode dan
dependensi sesuai hash run tersedia di `results_v5/source_snapshot/`.

Rangkaian ini menjalankan dua inversi, pemeriksaan geometri/kernel, figure dan
kuadratur. Inversi ablation, MCMC, dan stress tambahan tidak termasuk dalam run
ini dan tidak akan dinyatakan sudah dilakukan. Pengujian implementasinya tetap
ada, tetapi itu berbeda dari bukti perbandingan ilmiah yang memadai.

Metrik v4 dan v5 tidak boleh dibandingkan sebagai persentase peningkatan:
grid, bobot frekuensi dan realisasi data berubah. Kecocokan hasil v5 harus dibaca
terhadap input v5 yang disimpan, bukan terhadap gambar atau angka lama.

Run awal v5 dihentikan ketika smoke check tambahan mendeteksi bentuk vektor
sumber adapter eksternal yang salah. Bentuk input diperbaiki, tes langsung
ditambahkan. Pemeriksaan input berikutnya mendeteksi risiko alternasi kanal
pada grid overlap tidak cocok; kebijakan fusi diperbaiki, tes regresi
ditambahkan, dan run diulang dari awal. Run terhenti dipertahankan di
`results_v5_interrupted_adapter_check/` dan
`results_v5_interrupted_overlap_check/` dengan status `interrupted`; keduanya
tidak memiliki laporan hasil selesai.

## Batas yang masih nyata

Perubahan ini belum membuktikan bahwa masalah seismik sudah terpecahkan atau
posterior sudah konvergen. Densitas/Q diasumsikan diketahui. Noise EGF pasif
masih fenomenologis. Sigma likelihood komposit belum dikalibrasi. Referensi
MCMC independen yang memadai, sensitivitas prior, studi resolusi seluruh domain,
pengulangan dataset dan validasi lapangan masih diperlukan sebelum klaim ilmiah.

Area berkecepatan tinggi tidak boleh langsung diberi label Guided-P hanya dari
penampakan spektrum. Identifikasi moda memerlukan bukti tambahan. Pada dx 2 m,
komponen 150 m/s memiliki batas Nyquist sederhana 37.5 Hz; plot lebih rapat tidak
menghapus persoalan tersebut. Gelombang band-limited juga dapat memiliki ringing.
Pola demikian tidak boleh dihilangkan hanya agar gambar tampak lebih rapi.

Naskah baru ditahan sampai pemeriksaan ilmiah lolos. Perintah untuk memperoleh
laporan yang sesuai status bukti tersedia di `README.md`.

## Hasil run dan evaluasi terhadap prior

Kedua benchmark selesai dalam 841.7 detik waktu total. Data input, prediksi tersimpan, angka laporan, jumlah iterasi dan snapshot kode telah diperiksa ulang.

| Benchmark | J pada model tengah prior | J pada model mean hasil | J rata-rata prediksi sampel | Status |
|:--|--:|--:|--:|:--|
| benchmark_a | 0.36943 | 0.63742 | 0.30536 | belum tervalidasi |
| benchmark_b | 0.38388 | 0.37867 | 0.46983 | belum tervalidasi |

J lebih kecil berarti kecocokan komposit lebih baik. Model tengah prior adalah satu titik pembanding, bukan distribusi posterior yang sudah diketahui. Prior pada benchmark ini informatif dan mengelilingi model sintetik; perbandingan ini tidak membuktikan kalibrasi posterior.

**Temuan yang tidak boleh disembunyikan:** pada A, spektrum dari rata-rata parameter memberi misfit lebih besar daripada titik tengah prior, walaupun rata-rata prediksi sampel lebih baik. Pada B, rata-rata prediksi sampel lebih buruk daripada titik tengah prior; RMSE profil Vs juga belum membaik terhadap pembanding tersebut. Karena itu hasil ini tidak diterbitkan sebagai bukti keberhasilan inversi.

Penerbitan penelitian tetap diblokir. Pemeriksaan resolusi seluruh domain, restart independen, referensi posterior yang memadai dan kalibrasi data berulang belum selesai. Semua array dipertahankan agar kekurangan tersebut bisa diperiksa dan ditindaklanjuti.

Bukti rinci: `results_v5_validation/saved_result_verification.json`.

## Audit lanjutan: celah yang belum tertangkap 31 tes

Pemeriksaan berikut mempertahankan kode numerik, array hasil v5, dan hash run.
Temuan baru tidak berarti perbaikan di bagian ini sudah diterapkan.

1. **Skala Wasserstein berubah pada ablation jendela kecepatan.** Dengan dua
   spektrum yang sama dan W2 sekitar 50.04 m/s, kontribusi penalti W2 menjadi
   sekitar 24.57 kali lebih besar ketika batas kecepatan diubah dari 2450 menjadi
   598 m/s. Penyebabnya penyebut memakai span grid yang sudah dipotong. Perbandingan
   full-window dan Rayleigh-window belum mengisolasi efek informasi spektrum.
   Inversi ablation ini belum termasuk run v5 yang selesai.
2. **Bobot frekuensi ablation meluas ke luar pita.** Pada konfigurasi quick,
   pemaskeran bobot global memberi bobot 1.25 Hz pada titik pasif 12 Hz; trapesium
   khusus pita pasif memberi 0.25 Hz. Jumlah bobot menjadi 12 Hz, padahal rentang
   1–12 Hz hanya 11 Hz. Bobot perlu dihitung ulang pada support kanal terpilih.
3. **NaN dapat menggantung transport.** Subprocess telah menyelesaikan pemanggilan
   W2 valid dan masuk kasus NaN, tetapi tidak kembali sampai dihentikan timeout
   5 detik. Loop tidak memajukan indeks pada massa NaN. Observasi belum menjalani
   pemeriksaan finite di constructor likelihood. Data v5 yang diperiksa tetap finite.
4. **Estimasi prediktif 16 sampel belum stabil.** Komponen A terwakili [4,4,6,2]
   dan B [2,4,5,5], walaupun bobot campuran masing-masing 0.25. Evaluasi ulang 128
   sampel seimbang memberi J A 0.28951 (sebelumnya 0.30536), B 0.44736 (sebelumnya
   0.46983). Simpangan baku J pada resampling 16-draw dari pool tersebut sekitar
   0.036 dan 0.043. Ini sensitivitas Monte Carlo dari pool empiris, bukan interval
   kredibel yang terkalibrasi. Arah hasil B masih kurang baik terhadap prior pada
   estimasi 128-draw; besarnya selisih tidak boleh dibaca tanpa galat Monte Carlo.
5. **Catatan kritis hilang saat laporan dibuat ulang.** Reproduksi pada salinan
   temporer menunjukkan `--diagnostic-report` berhasil tetapi menghapus bagian
   perbandingan terhadap prior serta peringatan hasil yang ditambahkan manual.
   Laporan asli tetap dipertahankan. Temuan itu harus berasal dari data terstruktur
   agar regenerasi laporan tidak menghilangkannya.

Prioritas sebelum penggunaan lanjutan: validasi input finite, konsistensi skala
likelihood antar-ablation, bobot pita kanal, kestabilan prediksi, dan regenerasi
laporan lengkap. Bukti terstruktur: `audit_followup/additional_audit_findings.json`.

## Implementasi lanjutan v6

Kelima temuan audit lanjutan telah diterapkan dalam kode v6: validasi input
transport, skala Wasserstein yang tetap sepanjang ablation, kuadratur pada pita
kanal, prediksi campuran seimbang dengan galat Monte Carlo dan budget adaptif,
serta regenerasi laporan dari pembanding/peringatan terstruktur. Catatan editorial
manual kini mempunyai file terpisah yang dipertahankan oleh generator.

42 tes lolos, termasuk regresi kasus yang ditemukan. Uji integrasi dua model
(2 iterasi, dua komponen) selesai; hasilnya berada di `results_v6_smoke/` dan
tetap merupakan smoke test. Run diagnostik 60 iterasi per model menggunakan
empat komponen serta empat sampel optimasi per komponen selesai di
`results_v6/`. Hasil v5 dan snapshot kode sebelumnya tetap dipertahankan.

### Hasil pelaksanaan v6

Kedua benchmark selesai dalam 1994.2 detik. Run tetap diagnostik, dengan 60 iterasi pada masing-masing model. Inversi ablation, MCMC dan stress tambahan tidak termasuk run ini. Isolasi perubahan skala dan bobot ablation diperiksa melalui kasus terkontrol serta tes gradien/backend.

| Benchmark | Draw prediktif | Draw per komponen | J prediktif | MC standard error J | J prior midpoint | Target presisi |
|:--|--:|:--|--:|--:|--:|:--|
| benchmark_a | 256 | [64, 64, 64, 64] | 0.297636 | 0.006747 | 0.369432 | terpenuhi |
| benchmark_b | 512 | [128, 128, 128, 128] | 0.463700 | 0.007018 | 0.383878 | terpenuhi |

Semua draw prediktif yang disimpan telah dievaluasi ulang. Input, mean spektrum, spread respons, MC standard error, riwayat bootstrap, metrik, dan peringatan pembanding berhasil direproduksi. Target presisi ini menilai ketelitian estimasi numerik dari q hasil optimasi, bukan validasi q terhadap posterior sebenarnya.

**Temuan B tetap terlihat:** J prediktif lebih besar daripada prior midpoint; RMSE profil Vs juga masih lebih besar daripada pembanding tersebut. Peringatan tersebut berasal dari JSON eksperimen dan tidak dihapus saat regenerasi. Hasil ini tidak dinyatakan sebagai keberhasilan ilmiah atau kalibrasi posterior.

Kasus transport terkontrol yang sebelumnya berbeda sekitar 24.57 kali sekarang memiliki rasio kontribusi W2 sekitar 1.00000000009, pada skala tetap yang sama. Bobot pita pasif berjumlah 11 Hz dengan bobot endpoint 12 Hz sebesar 0.25 Hz. NaN ditolak melalui ValueError.

Bukti tersimpan di `results_v6_validation/`: log 42 tes, kontrol deterministik, validasi compliance independen, replay setiap benchmark, dan pemeriksaan regenerasi laporan/penerbitan. Snapshot source yang cocok dengan hash run berada di `results_v6/source_snapshot/`. Laporan terbaru: `results_v6/DIAGNOSTIC_VERIFICATION_REPORT.md`.
