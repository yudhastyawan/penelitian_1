# Restart Benchmark B pada observasi v6

Tiga restart memakai seed 2028, 2029, dan 2030. Seed 2027 adalah run asli. Setiap restart menggunakan konfigurasi SMVI, 60 iterasi, geometri, dan observasi asli yang sama. Seed mengubah inisialisasi, sampling optimasi, dan sampling prediksi.

Misfit prior spektrum gabungan **0.38388**. Δ positif berarti misfit prediksi di atas prior. Pita aktif memakai $30<f\leq42$ Hz dengan bobot v6.

| Seed | J gabungan | Galat baku MC | Draw | Target presisi | Δ aktif | Kontribusi pita aktif | Δ pasif |
|--:|--:|--:|--:|:--|--:|--:|--:|
| 2027 | 0.46370 | 0.00702 | 512 | terpenuhi | +0.10920 | +0.10013 | -0.08241 |
| 2028 | 0.34061 | 0.00956 | 128 | terpenuhi | -0.03545 | -0.02339 | -0.07629 |
| 2029 | 0.45133 | 0.00591 | 256 | terpenuhi | +0.09602 | +0.08792 | -0.09505 |
| 2030 | 0.37075 | 0.00679 | 256 | terpenuhi | +0.00102 | +0.01669 | -0.08443 |

## Parameter target

| Seed | Parameter | Rata-rata | Kuantil 95% campuran | Kuantil 95% prior | Acuan |
|--:|:--|--:|--:|--:|--:|
| 2027 | Vp_2 | 1703.19 | [1350.49, 2024.64] | [1221.25, 2028.75] | 1650.00 |
| 2027 | Vp_3 | 1507.64 | [1369.93, 1748.11] | [1371.25, 2178.75] | 1780.00 |
| 2027 | h_1 | 4.14 | [3.19, 5.22] | [2.10, 5.71] | 3.50 |
| 2028 | Vp_2 | 1421.31 | [1248.10, 1701.35] | [1221.25, 2028.75] | 1650.00 |
| 2028 | Vp_3 | 1522.58 | [1399.53, 1726.18] | [1371.25, 2178.75] | 1780.00 |
| 2028 | h_1 | 4.33 | [3.24, 5.55] | [2.10, 5.71] | 3.50 |
| 2029 | Vp_2 | 1444.56 | [1244.96, 1832.72] | [1221.25, 2028.75] | 1650.00 |
| 2029 | Vp_3 | 1541.66 | [1392.80, 1764.93] | [1371.25, 2178.75] | 1780.00 |
| 2029 | h_1 | 4.79 | [3.59, 5.49] | [2.10, 5.71] | 3.50 |
| 2030 | Vp_2 | 1327.09 | [1226.55, 1528.16] | [1221.25, 2028.75] | 1650.00 |
| 2030 | Vp_3 | 1584.97 | [1379.36, 1927.87] | [1371.25, 2178.75] | 1780.00 |
| 2030 | h_1 | 4.25 | [2.95, 5.39] | [2.10, 5.71] | 3.50 |

Satuan Vp dalam m/s dan h dalam m. Prior ketiga parameter berbentuk seragam pada batas yang disimpan. Kuantil campuran mengukur sebaran parameter bersyarat pada satu observasi dan satu hasil pengoptimasian.

## Pemeriksaan operator pada pita 30–42 Hz

Uji grid memakai spektrum bersih dari model yang tetap dan data acuan sintetik. Pembandingan gelombang memakai model berlapis Benchmark B, empat offset, serta pita 30–42 Hz. Lihat hasil JSON masing-masing untuk grid, frekuensi, dan galat.

## Cakupan

Perbandingan terhadap prior memakai observasi v6 yang sama. Ablasi dengan refit kanal terpisah dan pengukuran cakupan kuantil pada model acuan berulang merupakan eksperimen tersendiri.

## Evaluasi pada realisasi noise baru

Keempat hasil fitting dibekukan dan dinilai pada delapan realisasi noise baru dari model sintetik B yang sama. ΔJ adalah selisih berpasangan terhadap model tengah prior untuk setiap realisasi.

| Seed fitting | Prediksi mengungguli prior | Median ΔJ | Rentang ΔJ |
|--:|--:|--:|--:|
| 2027 | 6/8 | -0.06913 | [-0.14792, +0.11381] |
| 2028 | 8/8 | -0.10005 | [-0.15540, -0.00041] |
| 2029 | 6/8 | -0.07042 | [-0.16179, +0.09358] |
| 2030 | 7/8 | -0.10718 | [-0.15559, +0.03183] |

## Interpretasi

Pada observasi fitting, dua dari tiga restart mengungguli prior dan satu lebih buruk. Kontribusi misfit pita aktif >30–42 Hz berubah tanda antar-seed, sementara prediksi pasif membaik pada semua seed. Galat Monte Carlo prediksi lebih kecil daripada rentang J antar-restart. Kondisi untuk menganggap pola residual aktif stabil belum terpenuhi, sehingga ablasi belum dijalankan.

Pada realisasi noise baru, setiap hasil fitting mengungguli prior pada sedikitnya enam dari delapan realisasi. Hasil tersebut berlaku pada model sintetik dan mekanisme noise yang sama. Pada empat hasil fitting, rata-rata Vp lapisan 3 semuanya di bawah acuan, dan rentang kuantil 95% memuat acuan pada satu hasil. Keempat hasil memakai observasi yang sama, jadi frekuensi ini tidak mengukur cakupan berulang.

Pada spektrum aktif bersih 30–42 Hz, perapatan grid dari 2 Hz/4 m/s menjadi 0,25 Hz/2 m/s mengubah besarnya J. Peringkat model tengah prior terhadap lima model hasil v6 tetap sama. Uji gelombang band-limited terhadap PyProp8 pada empat offset memberi selisih relatif L2 maksimum sekitar 0,38%.
