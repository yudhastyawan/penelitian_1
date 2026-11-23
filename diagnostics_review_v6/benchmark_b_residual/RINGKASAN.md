# Diagnostik residual Benchmark B v6

Seluruh 512 sampel prediktif (128 per komponen) diputar ulang. Hasil fusi cocok dengan spektrum prediksi tersimpan. Nilai misfit prior juga cocok dengan laporan v6. Perhitungan memakai bobot frekuensi, mask kecepatan, serta skala Wasserstein 348 m/s dari Benchmark B.

| Spektrum | J model tengah prior | J prediksi rata-rata | Perubahan J |
|:--|--:|--:|--:|
| Gabungan | 0,38388 | 0,46370 | +0,07982 |
| Aktif | 0,40471 | 0,51392 | +0,10920 |
| Pasif | 0,22725 | 0,14485 | −0,08241 |

Setiap kanal dinilai terhadap observasinya sendiri pada dukungan frekuensi dan mask v6. Angka aktif dan pasif mengevaluasi model dari inversi gabungan. Ablasi dengan refit per kanal belum dilakukan.

| Pita pada spektrum gabungan | Perubahan penalti NCC | Perubahan penalti Wasserstein | Kontribusi ke perubahan J |
|:--|--:|--:|--:|
| 1–<8 Hz | −0,00524 | −0,00073 | −0,00597 |
| 8–12 Hz | −0,00362 | +0,00041 | −0,00321 |
| >12–30 Hz | +0,03537 | −0,03000 | +0,00537 |
| >30–42 Hz | +0,03160 | +0,05555 | +0,08715 |
| >42–55 Hz | +0,00122 | −0,00473 | −0,00351 |

Kontribusi berbobot terbesar berada di >30–42 Hz. Perbaikannya pada beberapa pita lain mengimbangi sebagian kenaikan tersebut. Penurunan NCC gabungan terutama terkonsentrasi pada kecepatan ≤680 m/s, dengan kontribusi +0,05886 dari total perubahan penalti NCC +0,05932. Penguraian kecepatan menggunakan selisih kuadrat dari spektrum yang dinormalisasi dengan norma L2 per irisan, sehingga jumlahnya sama dengan penalti NCC. Transport Wasserstein dapat memindahkan massa melintasi batas jendela kecepatan dan di sini hanya diuraikan menurut frekuensi.

| Komponen | J pusat | J spektrum rata-rata komponen | J campuran tanpa komponen ini |
|:--|--:|--:|--:|
| 1 | 0,46906 | 0,42237 | 0,43914 |
| 2 | 0,62679 | 0,61183 | 0,39973 |
| 3 | 0,41891 | 0,38548 | 0,47201 |
| 4 | 0,34157 | 0,33392 | 0,50593 |

Komponen 2 memberi kecocokan terburuk. Mengeluarkannya dan menormalisasi ulang bobot memperbaiki misfit campuran dari 0,46370 ke 0,39973. Angka tersebut tetap di atas prior 0,38388. Spektrum rata-rata komponen 4 secara tersendiri lebih cocok daripada prior. Perbandingan ini mengukur kecocokan pada data yang dipakai untuk inversi dan tidak menetapkan bobot baru. Empat pusat komponen memperkirakan Vp lapisan ketiga sekitar 1.398–1.606 m/s, di bawah acuan 1.780 m/s.

Nilai rinci ada di `residual_summary.json`. Spektrum per kanal, pusat dan rata-rata komponen ada di `spectra.npz`. Skala yang sama dipakai untuk seluruh perbandingan. Analisis ini melokalisasi selisih kecocokan, sedangkan penyebab fisik dan kestabilan komponen memerlukan pengujian lanjutan.
