---
title: "Inversi Spektrum Seismik Aktif dan Pasif dengan MFJ–SMVI"
subtitle: "Naskah lengkap penelitian dan lampiran | v6 dengan verifikasi COMPOSTI dan cakupan v3"
lang: id
papersize: a4
geometry: [left=20mm, right=20mm, top=18mm, bottom=18mm]
fontsize: 11pt
colorlinks: true
header-includes:
  - \usepackage{pdflscape}
  - \usepackage{graphicx}
  - |
    \usepackage{xcolor}
    \definecolor{alurbg}{HTML}{E8F1F4}
    \definecolor{alurfg}{HTML}{163D50}
    \newcommand{\kotakalur}[1]{\colorbox{alurbg}{\parbox[c][2.7em][c]{.193\linewidth}{\centering\small\sffamily\bfseries\color{alurfg}#1}}}
    \newcommand{\panahalur}{\hspace{.3em}{\color{alurfg}$\rightarrow$}\hspace{.3em}}
    \newcommand{\alur}[4]{\par\medskip\noindent\kotakalur{#1}\panahalur\kotakalur{#2}\panahalur\kotakalur{#3}\panahalur\kotakalur{#4}\par\medskip}
---

**Status ilmiah:** verifikasi sintetik. Konvergensi posterior dan validasi lapangan masih terbuka.

\clearpage

# Abstrak

Penelitian ini menguji inversi spektrum aktif–pasif tanpa picking dengan MFJ–SMVI pada dua model P–SV berlapis. Misfit prediksi A dan B adalah 0,29764 dan 0,46370, dibanding model tengah prior 0,36943 dan 0,38388. Pada delapan realisasi noise baru dari model B, tiap prediksi mengungguli prior pada enam hingga delapan realisasi. Pada uji aktif COMPOSTI, menghilangkan pergeseran frekuensi kompleks menurunkan selisih citra MFJ dari 15,81–22,79% menjadi 0,026–0,055%. Estimasi kecepatan kompresi lapisan ketiga masih sensitif terhadap optimasi.

**Kata kunci:** gelombang permukaan, spektrum dispersi, MFJ, SMVI.

# Abstract

This study investigates picking-free inversion of active and passive seismic dispersion images with a phase-only Modified Frequency–Bessel transform and Stein Mixture Variational Inference. A layered P–SV forward solver is tested on a trapped low-shear-velocity layer (A) and a compressional-velocity contrast at a shallow water table (B). The mean predictive misfits are 0.29764 and 0.46370, against prior-midpoint values of 0.36943 and 0.38388. On eight new noise realizations of the same B model, each fixed prediction outperforms the prior midpoint in six to eight cases. Removing the complex-frequency shift reduces the active MFJ discrepancy relative to COMPOSTI from 15.81–22.79% to 0.026–0.055%. This test covers the two active reference models. Recovery of the third-layer compressional velocity remains sensitive to optimization.

**Keywords.** Surface waves, dispersion images, picking-free inversion, P–SV, MFJ, SMVI.

# Kesenjangan dan pertanyaan penelitian

1. **Ambiguitas moda pada lapisan lambat.** Cabang Rayleigh dapat mendekat dan dominasi energi dapat beralih antar-moda. *Picking* maksimum per frekuensi berpotensi membentuk lintasan yang menyambungkan cabang berbeda. Seberapa jauh inversi citra utuh memperbaiki pemulihan $V_{s,2}$ dan ketebalan lapisan lambat terhadap inversi kurva dengan asumsi moda fundamental? [Forbriger (2003)](https://doi.org/10.1046/j.1365-246X.2003.01929.x).
2. **Sensitivitas terhadap kompresi.** Respons Rayleigh dan fitur spektral berkecepatan tinggi memiliki sensitivitas parameter yang dapat berbeda. Apakah jendela kecepatan yang lebih luas memperbaiki pemulihan $V_p$ dan kedalaman kontras, dan dapatkah fitur tersebut diidentifikasi sebagai moda terpandu? [Liu dkk. (2024)](https://doi.org/10.11720/wtyht.2024.1387).
3. **Konsistensi citra dan akuisisi.** Pencitraan terpengaruh apertur, sampling spasial, dan pilihan kernel. Bagaimana geometri memengaruhi cakupan radial dan azimut, dan apakah fusi aktif–pasif meningkatkan pemulihan parameter dibanding inversi tiap kanal? [Foti dkk. (2018)](https://doi.org/10.1007/s10518-017-0206-7).
4. **Ketidakpastian dan biaya komputasi.** Tukar imbangan antara ketebalan dan kecepatan dapat menghasilkan beberapa daerah dengan kecocokan yang serupa. Kapan campuran Gaussian berkorelasi memberi aproksimasi yang lebih baik daripada satu Gaussian, jika dinilai terhadap pembanding MCMC yang cukup mengeksplorasi distribusi target?

Tujuan penelitian adalah membangun dan memverifikasi operator maju, menguji empat rancangan geometri, mengukur pemulihan parameter pada A/B, serta menguji kontribusi tiap komponen melalui ablasi. [Lampiran C](#lampiran-c) memuat landasan, protokol lapangan, rancangan eksperimen, dan keterlacakan bagian v3.

\clearpage

# Latar belakang dan tujuan

Penelitian menguji apakah pencocokan spektrum frekuensi-kecepatan fase tanpa
picking dapat mengurangi kesalahan identifikasi moda pada lapisan lunak terjebak,
serta menambah informasi tentang Vp dan kedalaman batas jenuh melalui bagian
spektrum berkecepatan tinggi. Integrasi sumber aktif dan pasif diuji pada domain
pengolahan yang sama dengan ketidakpastian parameter yang dinyatakan eksplisit.

## Hipotesis penelitian

1. Penggabungan pita aktif-pasif menambah resolusi parameter dibandingkan masing-
   masing kanal, jika forward operator dan noise model konsisten.
2. Informasi spektral di luar jendela Rayleigh dapat mempersempit ketidakpastian Vp.
3. Campuran Gaussian berkorelasi dapat mengaproksimasi posterior yang tidak cukup
   diwakili satu Gaussian, dengan biaya dan akurasi yang diukur secara adil.
4. Geometri nested 2D dapat menawarkan trade-off cakupan radial dan azimut yang
   berguna. Tidak ada konfigurasi yang dinyatakan terbaik sebelum hasil pengujian.

# Desain eksperimen dan akuisisi

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

Rencana ablasi mencakup inversi aktif, pasif, gabungan, jendela Rayleigh, Gaussian tunggal, dan pembanding kurva hasil picking. Rencana pembanding posterior mencakup MCMC dengan rantai independen. Dokumentasi tiap eksperimen memuat anggaran komputasi, seed, prior, dan keluaran inferensi.

Uji stres yang direncanakan melemahkan sumber pada domain frekuensi-offset sebelum pencitraan, menambahkan noise dengan SNR terukur, dan menghitung ulang bobot fusi.

## Geometri rancangan

Geometri nested saat ini memiliki 9 sensor inti (spasi inline 10.5 m dengan
stagger melintang) serta 16 sensor spiral radius 52-320 m. Enam tembakan aktif
berada 4, 10, dan 20 m di luar masing-masing ujung bentangan 94 m. Pasangan
pasif tidak dibuang ke bin radial kasar secara default.

Pemilihan instrumen, sumber, durasi, dan respons frekuensi harus didasarkan
pada noise floor, hasil uji lapangan, serta spesifikasi terverifikasi. Klaim
spasi inti menjamin seluruh pita 8-14 Hz atau pelat tumbuk tertentu menjamin
energi frekuensi tertentu belum dianggap terbukti oleh kode sintetik.

`\clearpage`{=latex}

## Pembentukan observasi sintetik dan model gangguan

Untuk kanal $j$ dan offset $r$, gather bersih di domain frekuensi adalah $U_j(f,r)=S_j(f)G_{zz}(f,r)$. Kanal aktif memakai $S_a(f)=x e^{-x}$ dengan $x=[f/(25\,\mathrm{Hz})]^2$. Kanal pasif memakai $S_p(f)=[1+(f/(6\,\mathrm{Hz}))^{1.5}]^{-1}$. Kanal pasif mewakili respons Green vertikal dengan gangguan fenomenologis. Simulasi sumber ambien terdistribusi dan korelasi silang belum termasuk dalam benchmark.

**Kanal aktif.** Setiap offset menerima faktor amplitudo $A_r$ dengan $\log A_r\sim\mathcal N(0,0.35^2)$, dan pergeseran fase $\phi_r\sim\mathcal N(0,0.08^2)$ rad. Satu pasangan faktor per offset berlaku di semua frekuensi. Gather setelah modulasi adalah $\tilde U_a=A_r e^{i\phi_r}U_a$.

**Kanal pasif.** Gather setelah modulasi adalah $\tilde U_p=e^{i\phi_p(f,r)}U_p$ dengan $\phi_p=0.045\,a_r 2\pi f r/(260\,\mathrm{m/s})$. Faktor $a_r$ adalah rata-rata $\cos[2(\theta-25^\circ)]$ untuk pasangan yang masuk kelompok jarak $r$. Jumlah pasangan dalam kelompok itu dinyatakan sebagai $N_r$.

**Noise aditif.** Pada kedua kanal, $U_{\mathrm{obs},j}=\tilde U_j+b_j(r)\,R_j(f)(X+iY)/\sqrt{2}$, dengan $X,Y\sim\mathcal N(0,1)$ independen untuk tiap sel frekuensi-offset. $R_j(f)$ adalah RMS amplitudo $\tilde U_j$ di seluruh offset pada frekuensi itu. Faktor $b_a=0.08$ dan $b_p(r)=0.045/\sqrt{\max(N_r,1)}$ mengatur skala relatif gangguan kompleks.

**Pengacakan.** Pembuatan observasi menggunakan satu generator NumPy per benchmark, dengan urutan kanal aktif lalu pasif. Seed observasi A dan B berturut-turut 2026 dan 2027. Optimasi awal memakai generator terpisah dengan seed masing-masing 2026 dan 2027. Restart B memakai seed optimasi 2028–2030 dengan observasi B tetap. Evaluasi tambahan memakai seed gangguan 3101–3108 pada model acuan B. Prediksi maju kandidat memakai operator yang sama dalam keadaan tanpa gangguan.

`\clearpage`{=latex}

# Diagram alir penelitian

Dua bagian diagram mengikuti data dari perancangan benchmark sampai pelaporan. Pada setiap kotak, input dan output merujuk pada objek yang diteruskan atau disimpan. Evaluasi kandidat pada tahap 08 menggunakan kembali operator tahap 02, 04, dan 05. Tahap 11 menjalankan ulang tahap 07–10 ketika menguji restart.

```{=latex}
\noindent\includegraphics[width=\linewidth,height=.76\textheight,keepaspectratio]{output/pdf/figures/diagram_alir_mfj_smvi_1.pdf}
\clearpage
\noindent\includegraphics[width=\linewidth,height=.88\textheight,keepaspectratio]{output/pdf/figures/diagram_alir_mfj_smvi_2.pdf}
\clearpage
```

# Pseudocode alur kerja

[Lampiran B](#lampiran-metode) merinci diagram alir, input, output, dan pseudocode untuk setiap metode.

Pseudocode berikut merangkum jalur yang dipakai untuk benchmark. Nama prosedur bersifat deskriptif. Hubungan input dan output mengikuti diagram alir.

## Pembuatan observasi dan evaluasi kandidat

```text
OBSERVASI(model_acuan, properti_tetap, geometri, grid, seed_data)
    rng <- inisialisasi_generator_acak(seed_data)
    untuk kanal dalam [aktif, pasif]:
        G <- Green_PSV_Sommerfeld(model_acuan, geometri[kanal], grid)
        U <- kalikan_spektrum_sumber(G, kanal)
        U_obs <- tambahkan_gangguan(U, kanal, rng)
        S_obs[kanal] <- citrakan_MFJ(U_obs, grid_kecepatan, offset)
    D_obs, bobot_fusi <- normalisasi_dan_fusi_native(S_obs)
    simpan gather, spektrum kanal, bobot fusi, D_obs dan konfigurasi
    kembalikan D_obs

TARGET(z, D_obs, prior, geometri, grid, pengaturan)
    m, log_jacobian <- transformasi_ke_domain_fisik(z, prior)
    untuk kanal dalam [aktif, pasif]:
        U_syn <- sintesis_PSV_dan_sumber(m, properti_tetap, kanal)
        S_syn[kanal] <- citrakan_MFJ(U_syn, grid_kecepatan, offset)
    D_syn <- normalisasi_dan_fusi_native(S_syn)
    J, NCC, W2 <- misfit_berbobot(D_obs, D_syn, mask, skala_W2)
    log_target <- -J / (2 * sigma_D^2) + log_jacobian
    kembalikan log_target, J, NCC, W2, D_syn
```

## Optimasi campuran Gaussian

```text
INFERENSI(D_obs, ruang_parameter, pengaturan, seed_fit)
    kandidat <- transformasikan_Latin_hypercube_ke_laten(seed_fit)
    pilih K pusat dengan skor target tinggi dan jarak awal
    perhalus setiap pusat dengan L-BFGS-B hingga batas langkah
    inisialisasi Cholesky dan bobot setiap komponen = 1 / K
    untuk iterasi = 1 sampai batas_iterasi:
        tarik sampel normal yang dipakai bersama semua komponen
        untuk setiap komponen k dan sampel normal eps:
            z <- pusat[k] + Cholesky[k] * eps
            g <- gradien TARGET(z) - gradien log_campuran(z)
            akumulasi g untuk pusat dan outer(g, eps) untuk Cholesky
        rata-ratakan gradien dari sampel tiap komponen
        ubah gradien diagonal ke koordinat log diagonal
        arah <- transport_Stein(pusat, Cholesky, gradien)
        perbarui pusat dan Cholesky dengan Adam
        catat ELBO, misfit dan pusat komponen
    ambil sampel parameter laten dari campuran akhir
    kembalikan komponen, sampel fisik dan riwayat optimasi
```

## Prediksi dan verifikasi

```text
PREDIKSI(komponen, prior, operator, aturan_presisi)
    tetapkan jumlah sampel awal yang sama per komponen
    ulangi:
        tambah sampel yang diperlukan per komponen
        teruskan tiap sampel fisik melalui operator kandidat
        D_pred <- jumlah berbobot dari mean spektrum tiap komponen
        hitung galat baku RMS mean citra dan galat baku bootstrap J
        bandingkan J dengan J dari setengah sampel secara nested
        jika kriteria presisi terpenuhi atau anggaran habis: berhenti
        gandakan anggaran sampel hingga batas maksimum
    evaluasi D_pred, spektrum model rata-rata, pusat terbaik dan prior
    kembalikan spektrum prediksi, galat baku, metrik dan status presisi

EVALUASI_LANJUTAN(observasi_awal, konfigurasi, prediksi, acuan)
    ulangi INFERENSI dan PREDIKSI dengan seed baru pada observasi awal
    urai residual per kanal, pita frekuensi dan wilayah kecepatan
    hasilkan observasi baru dari acuan dan seed noise independen
    nilai prediksi tetap pada observasi baru
    bandingkan estimasi parameter terhadap prior dan acuan
    uji pengaruh resolusi grid dan bandingkan jejak dengan propagator lain
    laporkan bukti, variasi antar-run dan eksperimen yang masih terbuka
    simpan naskah, ringkasan metrik, sampel dan gambar
```

# Formulasi dan metode komputasi

## Model dan fungsi Green P-SV {#green-psv}

Medium terdiri dari tiga lapisan homogen isotropik di atas half-space. Parameter
inversi adalah empat $V_s$, empat $V_p$, dan tiga ketebalan $h$. Densitas dan
faktor kualitas ditetapkan sebagai nuisance parameters yang diketahui dalam
benchmark sintetik; keduanya tidak diestimasi dari data.

Dengan konvensi waktu $e^{-i\omega t}$, kecepatan kompleks adalah
$\tilde V=V(1-i/(2Q))$. Frekuensi evaluasi menggunakan
$\tilde\omega=\omega+i\eta$, $\eta=0.005\omega$. Akar vertikal
$\nu_\alpha=\sqrt{k^2-\tilde\omega^2/\tilde V_\alpha^2}$ dipilih dengan bagian
riil tidak negatif untuk memenuhi peluruhan/radiasi keluar. Pengaruh pergeseran numerik $\eta$ pada hasil aktif diperiksa melalui [pembandingan COMPOSTI](#uji-composti).

Rekursi impedansi P-SV dihitung dari half-space menuju permukaan menggunakan
partisi perpindahan dan traksi serta propagator peluruhan
$E_j=\operatorname{diag}(e^{-\nu_{p,j}h_j},e^{-\nu_{s,j}h_j})$:

$$R_j=-(T_{u,j}-K_{j+1}U_{u,j})^{-1}(T_{d,j}-K_{j+1}U_{d,j}),$$
$$K_j=(T_{d,j}+T_{u,j}E_jR_jE_j)(U_{d,j}+U_{u,j}E_jR_jE_j)^{-1},$$
$$G_{zz}(\omega,k)=-[K_1^{-1}]_{22}.$$

Respons Green pada offset $r$ dihitung dari compliance lengkap dalam domain
bilangan gelombang. Untuk setiap $\omega=2\pi f$, respons sebelum pembobotan
sumber adalah

$$G_{zz}(\omega,r)=\int_0^{k_{\max}}[kG_{zz}(\omega,k)-A]J_0(kr)\,dk
+\frac{A}{r},\qquad A=\frac{1-\tilde\nu_1}{\tilde\mu_1}.$$

$\tilde\nu_1$ dan $\tilde\mu_1$ masing-masing adalah rasio Poisson dan modulus
geser kompleks lapisan permukaan. Pengurangan $A$ meredam ekor integran dan
$A/r$ memulihkan kontribusi statik analitik. Gather bersih kanal $j$ diperoleh
dengan mengalikan respons terhadap spektrum sumber satu kali,
$U_j(f,r)=S_j(f)G_{zz}(\omega,r)$.

\newpage

**Kuadratur v6.** Untuk setiap frekuensi real $f$, dengan $r_{\min}$ offset
minimum pada kanal yang dihitung dan $c_{\min}$ batas bawah grid kecepatan fase,
batas bilangan gelombang ditetapkan sebagai

$$v_{\mathrm{num}}=\max(c_{\min}/2,40\,\mathrm{m/s}),\qquad
k_{\mathrm{split}}=\frac{\omega}{5000\,\mathrm{m/s}},\qquad
k_{\max}=\max\!\left(\frac{3\omega}{v_{\mathrm{num}}},\frac{8}{r_{\min}}\right).$$

Segmen linier berawal dari $10^{-10}\,\mathrm{m}^{-1}$ sebagai pendekatan
batas nol. Segmen ini dibentuk dengan 33 titik hingga $k_{\mathrm{split}}$,
lalu titik terakhirnya dibuang karena sudah menjadi titik pertama segmen
logaritmik hingga $k_{\max}$. Jumlah titik logaritmik adalah
$\max(64,\lfloor1024d\rfloor)$, dengan $d$ parameter densitas kuadratur.
Run utama memakai $d=2$, sehingga grid gabungan memuat 2.080 titik
(32 linier dan 2.048 logaritmik). Bobot trapesium titik interior adalah
$w_i=(k_{i+1}-k_{i-1})/2$. Pada kedua ujung, bobotnya setengah selisih ke titik
terdekat. Integral dihampiri dengan jumlah
$\sum_i w_i[k_iG_{zz}(\omega,k_i)-A]J_0(k_i r)$.
Grid dan bobot ini dipakai bersama oleh operator maju dan jalur autodiff.

## Transformasi MFJ dan pembanding

Xi dkk. (2021) memperkenalkan pencitraan menggunakan kernel Hankel untuk analisis fungsi Green empiris. Implementasi v6 memakai kernel satu arah dengan normalisasi fase, whitening, taper, dan bobot radial berikut.

Whitening per trace diterapkan sebagai
$\hat U_m=U_m/(|U_m|+0.02\max_l|U_l|+\epsilon)$.
Dengan taper ruang $w_m$ dan bobot trapesium radial $\Delta r_m$:

$$I(f,c)=\left|\sum_m w_m\hat U_m\,
\frac{H_0^{(2)}(2\pi f r_m/c)}{|H_0^{(2)}(2\pi f r_m/c)|}
\sqrt{r_m}\,\Delta r_m\right|.$$

Pembanding plane-wave dan standing-Bessel menggunakan kernel
$e^{-ikr}$ serta $J_0(kr)$ pada gather yang sama. Tidak ada artefak buatan,
penggeseran punggung moda, atau pengaburan citra untuk menentukan pemenang.
Kurva Rayleigh dari `disba` (Luu, v0.7.0) dihitung dari model acuan untuk pembanding dispersi. Identifikasi moda Guided-P pada zona kecepatan tinggi masih perlu diuji.

## Fusi aktif-pasif

Spektrum dinormalisasi min-max per frekuensi, lalu dipangkatkan dengan $p$.
Inversi utama A/B, termasuk restart B, memakai $p=1{,}35$. Uji perbandingan
geometri dan kernel memakai $p=1{,}40$. Kedua kanal
disalin hanya pada frekuensi native yang benar-benar disintesis dalam pita 1-55 Hz.
Grid fusi memakai frekuensi pasif di bawah batas overlap dan frekuensi aktif
di atas batas overlap. Di bagian dalam overlap, hanya frekuensi native yang
tersedia pada kedua kanal yang digabungkan. Kebijakan ini mencegah alternasi
kanal akibat perbedaan grid. Driver mensintesis kedua kanal pada node overlap
yang sama (0.5 Hz). Tidak ada interpolasi yang menciptakan irisan
data baru. Pada overlap 8-12 Hz, koefisien
apertur adalah $\alpha=\sin^2[\pi(f-8)/8]$. Di luar overlap,
$\alpha=0$ di bawah 8 Hz dan $\alpha=1$ di atas 12 Hz.

$$v=\max[(\operatorname{mean}D/\max D)^2,10^{-4}],\qquad
w_a=\frac{\alpha/v_a}{\alpha/v_a+(1-\alpha)/v_p},\quad w_p=1-w_a.$$

Besaran $v$ adalah proxy bentuk spektral, bukan estimator noise variance yang
terkalibrasi. Fusi tidak dapat menggunakan kanal pasif di atas batas rekamannya
12 Hz. Klaim penyelamatan otomatis pita 12-14 Hz tidak dibuat.

## Likelihood komposit dengan Wasserstein yang benar

Misfit NCC adalah satu dikurangi cosine similarity pada setiap irisan frekuensi.
Rata-rata frekuensi memakai bobot trapesium dalam Hz agar perubahan kepadatan
grid tidak otomatis mengubah bobot relatif pita aktif dan pasif.
Distribusi diskrit untuk transport adalah $p_q\propto\max(D_q,0)^2+10^{-12}$.
Transport monoton menghitung biaya kuadratik:

$$W_2^2(p,q)=\int_0^1|F_p^{-1}(u)-F_q^{-1}(u)|^2\,du.$$

Ini berbeda dari integral kuadrat selisih CDF. Untuk dua point masses dengan
jarak $d$, implementasi harus memberikan $W_2=|d|$.

$$J=(1-\lambda)\overline{(1-\mathrm{NCC})}
+\lambda\frac{\overline{W_2^2}}{s_W^2},
\qquad \log L=-J/(2\sigma_D^2).$$

$s_W$ merupakan skala transport tetap dalam m/s, yang disimpan pada konfigurasi
benchmark. Jika tidak diberikan, likelihood menetapkannya satu kali sebagai
0.15 kali span grid parent penuh sebelum mask kecepatan diterapkan. Full-window
dan Rayleigh-window pada parent yang sama wajib memakai nilai $s_W$ yang sama.
Jika grid parent berbeda, eksperimen perbandingan harus menetapkan skala tersebut
secara eksplisit. Helper misfit mandiri mempunyai default absolut 200 m/s.

Bobot frekuensi trapesium dihitung ulang pada support kanal/mask terpilih;
setiap pita kontigu memakai endpoint-nya sendiri. Bobot global yang sekadar
dikalikan mask tidak dipakai untuk ablation karena dapat memasukkan lebar
frekuensi dari luar pita yang dipilih.

Pada irisan dengan $f\le12$ Hz, misfit v6 memakai kecepatan sampai 680 m/s. Aturan ini juga mencakup bagian overlap dan mengikuti nilai frekuensi aktual.
Pada eksperimen Rayleigh-only, kedua spektrum dibatasi dengan mask kecepatan yang
sama. Nilai $\sigma_D$ adalah hyperparameter generalized likelihood; interval
hasil inferensi belum boleh dianggap terkalibrasi secara statistik.

## Prior dan transformasi parameter

Prior seragam didefinisikan pada domain fisik yang memenuhi batas kotak dan
$V_p\ge1.45V_s$. Kecepatan $V_s$ dan ketebalan memakai transformasi sigmoid.
Untuk setiap $V_p$, batas bawah kondisional adalah
$b_p(V_s)=\max(V_{p,\min},1.45V_s)$, sehingga
$V_p=b_p+(V_{p,\max}-b_p)\sigma(z_p)$.

Jacobian bersifat triangular. Determinannya dihitung dari hasil kali diagonal,
termasuk lebar interval $V_p$ kondisional. Tidak ada clipping model setelah
sampling dan tidak ada regularisasi laten tersembunyi. Inisialisasi memakai Latin hypercube, seleksi menurut target laten, dan refinemen L-BFGS-B sebelum pembaruan campuran. Langkah ini menggunakan data dan prior.

## Inferensi campuran berkorelasi

$$q(z)=\frac{1}{K}\sum_{k=1}^K\mathcal N(z;\mu_k,L_kL_k^T).$$

$L_k$ adalah matriks triangular bawah dengan diagonal positif. Sampling,
densitas, skor, dan optimasi memakai matriks yang sama. Parameter bebas
komponen adalah $\Theta_k=(\mu_k,\log\operatorname{diag}L_k,
\operatorname{strictLower}L_k)$. Kernel dihitung atas seluruh $\Theta$,
bukan hanya mean. Bobot campuran sama dan tidak diubah melalui aturan ad hoc.

Untuk gradien ELBO per komponen $g_j$, pembaruan Stein adalah

$$\Phi_i=\frac{1}{K}\sum_j
[k(\Theta_j,\Theta_i)g_j+\alpha\nabla_{\Theta_j}k(\Theta_j,\Theta_i)].$$

Kernel RBF menggunakan $k=\exp(-\|\Theta_j-\Theta_i\|^2/h)$ dengan
$h=\operatorname{median}_{i<j}\|\Theta_i-\Theta_j\|^2/\log(K+1)$.
Tidak ada normalisasi per kolom yang mengubah persamaan tersebut. Parameter
$\alpha$ adalah kekuatan regularisasi entropi pada distribusi komponen.
Koordinat log-diagonal menjaga Cholesky positif; ini merupakan parameterisasi
implementasi, bukan klaim identik dengan setiap detail implementasi GeoPVI.

Sampel normal standar dibagikan antar-komponen pada setiap iterasi. Gradien
pathwise menggunakan skor target dikurangi skor campuran, diikuti Adam.
Gradien target dapat dihitung melalui diferensiasi otomatis CPU float64 dari
rekursi dan kuadratur yang sama; spektrum, nilai target, serta gradien diperiksa
terhadap backend Numba dan finite difference. Backend AD bukan solver fisika
independen dan tidak memakai surrogate moda.

Inisialisasi memakai Latin hypercube dan, jika diaktifkan, refinemen L-BFGS
berdasarkan data serta prior yang sama. Nilai model sejati tidak diberikan kepada
optimizer. Densitas/Q tetap merupakan asumsi diketahui dalam eksperimen sintetik.
Pembanding Gaussian tunggal memakai $K=1$, tanpa repulsi atau pemindahan sampel.

Rujukan metodologi:
[Zhao dan Curtis (2026), bagian 2.4](https://arxiv.org/html/2610.23603v1).
Uji target dua moda memeriksa mean, simpangan baku, fraksi kedua moda dan bobot
campuran. Keberhasilan uji tersebut tidak membuktikan konvergensi posterior seismik.

## Pembanding MCMC yang direncanakan

MCMC stretch-move mengevaluasi proposal tanpa clipping. Evaluasi di luar domain
valid ditolak melalui log posterior, bukan diproyeksikan ke batas. Awal MCMC
independen dari pusat SMVI. Chain lengkap, acceptance rate, estimasi waktu
autokorelasi, dan ESS disimpan. Panjang chain di bawah 50 kali estimasi waktu
autokorelasi ditandai tidak memadai; pemeriksaan tersebut sendiri belum membuktikan
konvergensi atau eksplorasi seluruh moda.

Metrik spektrum dihitung ulang pada posterior mean dan komponen terbaik secara
terpisah. Komponen terbaik dari sekumpulan pusat bukan global MAP. RMSE relatif
profil adalah $100\sqrt{\operatorname{mean}[((\hat V_s-V_s)/V_s)^2]}$;
MAPE dilaporkan terpisah. Coverage 11 parameter pada satu model tidak sama dengan
kalibrasi CI dari pengulangan dataset.

## Visualisasi dan status verifikasi

Wiggle dipilih mendekati jarak radial merata dengan kedua ujung apertur disertakan.
Halfwidth adalah 0.35 kali celah minimum jejak terpilih sehingga jejak tidak saling
menutupi. Normalisasi per trace untuk tampilan tidak boleh dibaca sebagai informasi
amplitudo absolut. Source spectrum pasif pada FFT dan sintesis frekuensi sama;
noise time gather adalah ilustrasi terpisah dari realisasi frekuensi yang di-fit.

Pada geometri aktif reguler, garis $c=2\Delta x f$ menunjukkan batas Nyquist
spasial sederhana. Memperhalus grid gambar tidak mengatasi aliasing akuisisi.
Geometri lapangan tetap 48 kanal dengan jarak 2 m; perbaikan ini tidak diam-diam
mengganti jumlah alat atau jarak receiver.

Rata-rata spektrum dari sampel campuran dihitung terpisah dari spektrum pada
rata-rata parameter fisik. Keduanya disimpan dan diberi label berbeda. Rata-rata
prediktif dengan jumlah draw terbatas adalah diagnostik Monte Carlo.

`research_quality.py` memisahkan eksekusi selesai dari kelayakan publikasi.
Uji cepat, konvergensi yang belum diperiksa dari restart independen, referensi
posterior yang belum tervalidasi, uji resolusi yang hilang, dan kalibrasi dataset
berulang yang belum dilakukan menahan penetapan hasil sebagai penelitian tervalidasi. Draf ini mendokumentasikan hasil diagnostik dengan status tersebut.

## Ketelitian estimasi prediksi (v6)

Perhitungan prediksi tidak lagi memilih 16 titik dari kumpulan sampel campuran.
Untuk setiap komponen k diambil jumlah draw yang sama, n_k, dari Gaussian laten
menggunakan Cholesky yang tersimpan, kemudian dipetakan ke model fisik. Estimator
mean adalah jumlah w_k kali mean spektrum per komponen. Standard error mean
spektrum dihitung dari akar jumlah w_k^2 Var(D_k)/n_k. Spread respons model
campuran dihitung terpisah; keduanya tidak boleh dipertukarkan.

Budget bertambah adaptif dari minimum yang dikonfigurasi sampai maksimum.
Diagnostik mencakup standard error RMS citra, bootstrap standard error J, serta
perbedaan J terhadap estimator dari setengah draw secara nested. Target default
0.01 adalah target ketelitian numerik, bukan ambang keberhasilan inversi.
Apabila target belum tercapai pada batas budget, flag precision_target_met tetap
false dan ditampilkan pada laporan serta gerbang publikasi.

Bootstrap dilakukan di dalam masing-masing komponen dan digabung dengan bobot
campuran, sehingga tidak menambahkan fluktuasi komposisi komponen yang tidak
perlu. Galat ini bersifat kondisional pada q hasil optimasi; tidak membuktikan
kalibrasi q terhadap posterior sebenarnya atau noise pengukuran.


# Parameter inversi dan konfigurasi benchmark

## Ruang parameter fisik

Vektor inversi berisi 11 parameter dengan urutan $(V_{s,1:4},V_{p,1:4},h_{1:3})$. Indeks mengikuti urutan dari permukaan. Lapisan 4 adalah half-space dan tidak memiliki parameter ketebalan. Kedalaman antarmuka lapisan ke-$j$ dihitung sebagai jumlah ketebalan hingga $h_j$. Pada Benchmark B, $h_1$ merepresentasikan kedalaman kontras sintetik muka air tanah.

| Parameter | Satuan | Prior A | Acuan A | Prior B | Acuan B |
| :-- | --: | --: | --: | --: | --: |
| $V_{s,1}$ | m/s | 170–350 | 260 | 115–240 | 175 |
| $V_{s,2}$ | m/s | 95–240 | 150 | 160–320 | 235 |
| $V_{s,3}$ | m/s | 270–510 | 390 | 250–460 | 350 |
| $V_{s,4}$ | m/s | 500–860 | 680 | 440–750 | 590 |
| $V_{p,1}$ | m/s | 380–760 | 540 | 360–820 | 550 |
| $V_{p,2}$ | m/s | 320–720 | 480 | 1200–2050 | 1650 |
| $V_{p,3}$ | m/s | 620–1200 | 880 | 1350–2200 | 1780 |
| $V_{p,4}$ | m/s | 1100–1900 | 1480 | 1550–2450 | 1980 |
| $h_{1}$ | m | 2.2–6.5 | 4 | 2–5.8 | 3.5 |
| $h_{2}$ | m | 3.5–9.5 | 6 | 3.8–9.8 | 6.5 |
| $h_{3}$ | m | 9–20 | 14 | 9–20 | 14 |

Prior memakai densitas seragam pada domain fisik yang dibatasi interval tabel dan syarat $V_p\ge1.45V_s$ untuk tiap lapisan. Nilai acuan dipakai saat sintesis dan evaluasi pemulihan parameter. Optimizer menerima observasi, prior, serta properti yang ditetapkan.

## Properti medium yang ditetapkan

| Benchmark / lapisan | Densitas (kg/m$^3$) | Qs | Qp |
| :-- | --: | --: | --: |
| A / 1 | 1780 | 25 | 50 |
| A / 2 | 1620 | 18 | 40 |
| A / 3 | 1960 | 35 | 75 |
| A / 4 | 2180 | 55 | 110 |
| B / 1 | 1660 | 22 | 45 |
| B / 2 | 1980 | 28 | 85 |
| B / 3 | 2060 | 38 | 95 |
| B / 4 | 2200 | 55 | 120 |

## Pengaturan inferensi yang digunakan

| Pengaturan | A | B |
| :-- | --: | --: |
| Seed run awal | 2026 | 2027 |
| Rentang grid kecepatan (m/s) | 120–1248 | 130–2450 |
| Langkah grid kecepatan (m/s) | 4 | 4 |
| Densitas kuadratur $d$ | 2 | 2 |
| Pangkat normalisasi fusi $p$ | 1,35 | 1,35 |
| Skala Wasserstein $s_W$ (m/s) | 169.2 | 348.0 |
| Komponen × iterasi | 4 × 60 | 4 × 60 |
| Sampel gradien per komponen/iterasi | 4 | 4 |
| Kandidat awal / batas refinemen per pusat | 32 / 30 | 32 / 30 |
| Adam untuk pusat / Cholesky | 0,04 / 0,015 | 0,04 / 0,015 |
| Bobot Wasserstein / skala likelihood | 0,35 / 0,08 | 0,35 / 0,08 |
| Koefisien tolakan Stein | 0,2 | 0,2 |
| Sampel parameter | 1200 | 1200 |
| Sampel prediksi min–maks | 128–512 | 128–512 |
| Bootstrap prediksi / toleransi | 32 / 0,01 | 32 / 0,01 |

Gradien memakai autodiff. Pita fusi 1–55 Hz dengan overlap 8–12 Hz. Mask 680 m/s berlaku hingga 12 Hz. Batas iterasi pada tabel adalah anggaran komputasi, sedangkan toleransi prediksi mengatur galat Monte Carlo.

Run menggunakan 4 komponen, 60 iterasi, dan 4 sampel Monte Carlo per komponen pada setiap iterasi optimasi. Kanal aktif menyintesis 31 frekuensi dari 8 sampai 55 Hz. Kanal pasif menyintesis 24 frekuensi dari 0.5 sampai 12 Hz. Waktu eksekusi total 1994.2 detik. Run v6 mencakup dua benchmark utama serta perbandingan geometri dan kernel. Inversi ablasi, MCMC, dan uji stres tambahan belum dijalankan dalam run ini.

# Hasil dan pembahasan

## Perbandingan geometri dan kernel pencitraan

Keempat rancangan memakai 25 stasiun pasif dan menghasilkan 300 pasangan per rancangan. Cakupan radial dihitung dari proporsi 29 bin berukuran 20 m dalam rentang 5–585 m yang berisi setidaknya satu pasangan. Momen azimut orde dua, $|\langle e^{2i\theta}\rangle|$, merangkum konsentrasi arah pasangan. Angka mendekati nol menunjukkan pembatalan momen orde dua pada sampel arah tersebut.

| Geometri | Jarak min. (m) | Jarak maks. (m) | Bin terisi (%) | Momen azimut |
| :-- | --: | --: | --: | --: |
| Linear | 25,00 | 600,00 | 79,3 | 1,000 |
| L-array | 25,00 | 424,26 | 69,0 | 0,357 |
| Cincin | 19,13 | 560,00 | 65,5 | 0,000 |
| Nested 2D | 10,92 | 562,63 | 89,7 | 0,104 |

Nested 2D menyediakan jarak minimum terpendek dan proporsi bin radial terisi tertinggi. Bentangan linear menjangkau jarak maksimum terbesar dan mempunyai arah pasangan yang paling terkonsentrasi. Geometri cincin memiliki momen azimut orde dua paling kecil. Ketiga metrik menerangkan aspek rancangan yang berbeda.

Perbandingan citra dispersi menggunakan model sintetik A, dengan seed gangguan 2029 untuk aktif dan 2030 untuk pasif. NCC pada overlap 8–12 Hz dihitung sebagai kesamaan kosinus seluruh jendela. MAE spektral adalah rata-rata selisih absolut terhadap citra bersih pada grid gabungan. Setiap kombinasi geometri dan kernel dibandingkan dengan citra bersihnya sendiri. Perbandingan ini memakai pangkat normalisasi fusi $p=1{,}40$ untuk citra bising dan bersih.

| Geometri / kernel | NCC overlap | MAE spektral |
| :-- | --: | --: |
| Linear / Plane-wave | 0,9591 | 0,0588 |
| Linear / Bessel berdiri | 0,9454 | 0,0772 |
| Linear / MFJ | 0,9593 | 0,0588 |
| L-array / MFJ | 0,9488 | 0,0547 |
| Cincin / MFJ | 0,9789 | 0,0369 |
| Nested 2D / MFJ | 0,7115 | 0,0809 |

Dalam realisasi uji ini MFJ pada geometri cincin memberi kesamaan paling tinggi terhadap referensi bersihnya sendiri. Nested 2D mengisi lebih banyak bin radial dan memiliki kesamaan citra yang lebih rendah pada overlap. Penilaian akurasi dispersi antar-kernel memerlukan referensi fisik bersama. Penilaian kinerja akuisisi memerlukan beberapa realisasi noise dan rentang model yang lebih luas.

## Benchmark A

| Evaluasi spektrum | Misfit | NCC | W2 (m/s) |
| :-- | --: | --: | --: |
| Model tengah prior | 0.36943 | 0.79136 | 138.29 |
| Model sintetik acuan | 0.17391 | 0.94367 | 105.97 |
| Spektrum pada model rata-rata | 0.63742 | 0.78033 | 201.14 |
| Spektrum komponen dengan misfit terendah | 0.23352 | 0.90679 | 118.93 |
| Spektrum prediksi rata-rata | 0.29764 | 0.86808 | 131.65 |

Prediksi memakai 256 sampel dengan 64 sampel per komponen. Galat baku Monte Carlo untuk misfit sebesar 0.00675, dan galat baku RMS spektrum sebesar 0.00836. Target presisi 0.01 terpenuhi. Pilihan komponen terbaik didasarkan pada misfit di pusat komponen.

| Parameter | Acuan | Rata-rata | Kuantil 2,5% | Kuantil 97,5% |
| :-- | --: | --: | --: | --: |
| Vs 1 | 260.00 | 271.14 | 232.87 | 318.67 |
| Vs 2 | 150.00 | 148.93 | 121.78 | 177.98 |
| Vs 3 | 390.00 | 412.79 | 339.88 | 482.89 |
| Vs 4 | 680.00 | 669.07 | 520.45 | 811.12 |
| Vp 1 | 540.00 | 563.14 | 461.70 | 671.73 |
| Vp 2 | 480.00 | 542.59 | 429.65 | 692.69 |
| Vp 3 | 880.00 | 938.26 | 785.04 | 1097.22 |
| Vp 4 | 1480.00 | 1550.56 | 1203.12 | 1868.08 |
| h 1 | 4.00 | 4.43 | 3.56 | 5.40 |
| h 2 | 6.00 | 6.51 | 4.63 | 8.91 |
| h 3 | 14.00 | 15.25 | 11.05 | 18.99 |

Kuantil di atas berasal dari campuran hasil estimasi dan belum dikalibrasi. Kecepatan dalam m/s dan ketebalan dalam m. RMSE relatif profil kecepatan geser sebesar 14.09%, dibanding 15.43% untuk pembanding prior.

## Benchmark B

| Evaluasi spektrum | Misfit | NCC | W2 (m/s) |
| :-- | --: | --: | --: |
| Model tengah prior | 0.38388 | 0.90193 | 332.82 |
| Model sintetik acuan | 0.33542 | 0.91342 | 310.79 |
| Spektrum pada model rata-rata | 0.37867 | 0.80480 | 295.17 |
| Spektrum komponen dengan misfit terendah | 0.34157 | 0.77306 | 259.12 |
| Spektrum prediksi rata-rata | 0.46370 | 0.81067 | 343.31 |

Prediksi memakai 512 sampel dengan 128 sampel per komponen. Galat baku Monte Carlo untuk misfit sebesar 0.00702, dan galat baku RMS spektrum sebesar 0.00464. Target presisi 0.01 terpenuhi. Pilihan komponen terbaik didasarkan pada misfit di pusat komponen.

| Parameter | Acuan | Rata-rata | Kuantil 2,5% | Kuantil 97,5% |
| :-- | --: | --: | --: | --: |
| Vs 1 | 175.00 | 186.35 | 132.73 | 234.49 |
| Vs 2 | 235.00 | 239.26 | 187.01 | 290.95 |
| Vs 3 | 350.00 | 354.84 | 290.69 | 418.74 |
| Vs 4 | 590.00 | 603.54 | 477.57 | 717.56 |
| Vp 1 | 550.00 | 627.69 | 517.91 | 764.10 |
| Vp 2 | 1650.00 | 1703.19 | 1350.49 | 2024.64 |
| Vp 3 | 1780.00 | 1507.64 | 1369.93 | 1748.11 |
| Vp 4 | 1980.00 | 2153.95 | 1803.08 | 2377.34 |
| h 1 | 3.50 | 4.14 | 3.19 | 5.22 |
| h 2 | 6.50 | 6.94 | 4.35 | 9.58 |
| h 3 | 14.00 | 15.23 | 10.66 | 19.60 |

Kuantil di atas berasal dari campuran hasil estimasi dan belum dikalibrasi. Kecepatan dalam m/s dan ketebalan dalam m. RMSE relatif profil kecepatan geser sebesar 10.52%, dibanding 9.82% untuk pembanding prior.

## Interpretasi terhadap prior

Pada A, rata-rata spektrum prediksi memperbaiki misfit terhadap prior, sementara spektrum yang dihitung pada rata-rata parameter masih kurang cocok. Pada B, rata-rata spektrum prediksi dan RMSE relatif profil kecepatan geser masih lebih buruk daripada pembanding prior. Model tengah prior merupakan satu titik dalam ruang parameter yang dibatasi. Evaluasi model sintetik acuan memakai data bising yang sama dan berfungsi sebagai pembanding diagnostik.

## Kestabilan Benchmark B dan prediksi pada noise baru

Tiga restart memakai seed 2028–2030 dan observasi asli yang sama. Seed memengaruhi inisialisasi dan sampling optimasi. Misfit prior gabungan adalah 0,38388. Kolom perubahan kanal dan kontribusi pita menyatakan selisih terhadap prior. Nilai positif berarti kecocokan menurun.

| Seed | J gabungan | Galat baku MC | $\Delta$ aktif | $\Delta$ pita aktif | $\Delta$ pasif |
| :-- | --: | --: | --: | --: | --: |
| 2027 | 0,46370 | 0,00702 | 0,10920 | 0,10013 | -0,08241 |
| 2028 | 0,34061 | 0,00956 | -0,03545 | -0,02339 | -0,07629 |
| 2029 | 0,45133 | 0,00591 | 0,09602 | 0,08792 | -0,09505 |
| 2030 | 0,37075 | 0,00679 | 0,00102 | 0,01669 | -0,08443 |

Kontribusi pita aktif di atas 30 hingga 42 Hz berubah arah antar-seed. Dua
dari tiga restart memberi misfit prediksi gabungan lebih kecil dari prior
pada observasi fitting. Prediksi pasif membaik pada semua seed.

### Evaluasi pada noise baru

Untuk menilai kinerja di luar observasi fitting, prediksi dari setiap run
dipertahankan tetap dan dinilai pada delapan realisasi noise baru (seed
3101–3108) dari model acuan dan generator gangguan Benchmark B yang sama.
Untuk setiap realisasi, $\Delta J=J_{\mathrm{pred}}-J_{\mathrm{prior}}$
dihitung terhadap model tengah prior pada observasi yang sama.
Nilai negatif menunjukkan perbaikan.

| Seed optimasi | $\Delta J<0$ | Median $\Delta J$ | Min. $\Delta J$ | Maks. $\Delta J$ |
| :-- | --: | --: | --: | --: |
| 2027 | 6/8 | -0,06913 | -0,14792 | +0,11381 |
| 2028 | 8/8 | -0,10005 | -0,15540 | -0,00041 |
| 2029 | 6/8 | -0,07042 | -0,16179 | +0,09358 |
| 2030 | 7/8 | -0,10718 | -0,15559 | +0,03183 |

Perbaikan terkecil pada seed optimasi 2028 sebesar 0,00041. Median dan
rentang di atas merangkum variasi akibat noise pada satu model acuan B.
Generalisasi ke model geologi lain dan data lapangan memerlukan pengujian
lanjutan.

| Seed | Vp lapisan 3 (m/s) | Kuantil 95% (m/s) |
| :-- | --: | --: |
| 2027 | 1507,64 | 1369,93–1748,11 |
| 2028 | 1522,58 | 1399,53–1726,18 |
| 2029 | 1541,66 | 1392,80–1764,93 |
| 2030 | 1584,97 | 1379,36–1927,87 |

Nilai acuan $V_{p,3}$ adalah 1780 m/s. Semua rata-rata berada di bawah acuan. Satu dari empat rentang kuantil memuat nilai acuan pada observasi yang sama. Analisis ini mendukung kehati-hatian dalam menafsirkan pemulihan parameter dan belum mengukur cakupan pada dataset berulang.

## Pemeriksaan khusus pita aktif 30–42 Hz

Pada spektrum sintetik bersih dan model tetap, grid dari 2 Hz/4 m/s hingga 0,25 Hz/2 m/s mengubah nilai misfit. Model tengah prior tetap lebih cocok daripada model parameter rata-rata dan keempat pusat komponen v6 pada pita uji itu. Pembandingan jejak band-limited memakai propagator PyProp8 (Valentine dan Sambridge, 2022) dan kuadratur bilangan gelombang terpisah. Selisih relatif L2 maksimum sekitar 0,38% pada empat offset yang diuji (4, 20, 60, 114 m). Perbandingan PyProp8 ini memakai pergeseran frekuensi kompleks yang sama pada kedua solver.


\clearpage

## Verifikasi forward aktif dengan COMPOSTI {#uji-composti}

**Protokol.** COMPOSTI (Niskanen dan Lähivaara, 2023) dijalankan melalui
fungsi Levin yang dipanggil dalam kode reproduksi Mreyen dkk. (2026).
Kedua solver memakai model acuan A/B, gaya titik vertikal di permukaan,
spektrum sumber 25 Hz, dan 56 offset identik (4–114 m), tanpa noise.
Spektrum kompleks dan citra MFJ aktif dievaluasi pada 31 frekuensi native
di rentang 8–55 Hz. Jejak waktu dibentuk dari grid seragam 0,5 Hz di pita
yang sama. Keluaran COMPOSTI dikonjugasi untuk konvensi waktu v6, dan
spektrum v6 dibagi $2\pi$ untuk normalisasi gaya titik satuan. Tidak
dilakukan fitting skala, fase, atau pergeseran waktu.

**Selisih relatif L2 terhadap COMPOSTI pada seluruh grid aktif.** Citra
dibentuk dengan operator MFJ yang sama sebelum fusi.

| Pengaturan v6 | Spektrum A | Spektrum B | MFJ A | MFJ B |
| :-- | --: | --: | --: | --: |
| Bawaan | 6,404% | 8,939% | 15,813% | 22,794% |
| $\eta=0$, grid bawaan | 0,510% | 0,571% | 0,026% | 0,055% |
| $\eta=0$, grid diperhalus | 0,514% | 0,578% | 0,024% | 0,055% |
| $\eta=0$, rentang dan taper COMPOSTI | 0,090% | 0,102% | 0,023% | 0,033% |

**Selisih relatif L2 jejak vertikal pada empat offset.** Metrik memakai
seluruh jejak dua detik. Gambar 9–10 pada Lampiran A menampilkan jendela
waktu 0–0,95 s dan memakai satu skala per offset untuk ketiga kurva.

| Benchmark dan pengaturan | 4 m | 20 m | 60 m | 114 m |
| :-- | --: | --: | --: | --: |
| A, bawaan | 1,708% | 6,756% | 15,877% | 24,748% |
| A, $\eta=0$, grid bawaan | 0,757% | 0,019% | 0,014% | 0,038% |
| B, bawaan | 2,712% | 11,148% | 25,028% | 35,049% |
| B, $\eta=0$, grid bawaan | 0,683% | 0,025% | 0,034% | 0,074% |

Penghilangan pergeseran numerik $\eta=0{,}005\omega$ pada grid bawaan
menjelaskan sebagian besar selisih. Atenuasi material melalui $Q_p$ dan
$Q_s$ tetap dipakai. Kontrol grid memperbanyak titik logaritmik menjadi
8.192 dan melipatempatkan batas integrasi. Kontrol tambahan menyamakan
rentang slowness dan taper dengan COMPOSTI. Kesepakatan pada A/B mendukung
solver P–SV untuk kanal aktif dalam kondisi tersebut. Energi berkecepatan
fase tinggi juga tampak pada COMPOSTI untuk B. Atribusi *guided P* masih
memerlukan uji moda. Kanal pasif, kestabilan pada ruang prior, dan dampak
pada hasil inversi berada di luar cakupan pembandingan ini.

# Verifikasi numerik dan keterbatasan

Terhadap PyProp8, galat relatif maksimum compliance adalah $6{,}07\times10^{-7}$ pada A dan $3{,}70\times10^{-6}$ pada B di titik yang diuji. Kuadratur lolos pemeriksaan pada frekuensi dan offset terpilih. Sebanyak 42 tes perangkat lunak v6 lulus, dan prediksi serta galat Monte Carlo direproduksi dari sampel tersimpan.

Model sintetik mengasumsikan medium 1D dengan densitas dan faktor kualitas diketahui. Restart pada observasi yang sama menguji sensitivitas terhadap seed. Evaluasi resolusi dan kalibrasi memerlukan beberapa model acuan, refit per kanal, dan pembanding posterior yang memadai.

\clearpage

# Kesimpulan dan tahapan berikutnya

Untuk Benchmark B, prediksi gabungan dari setiap run mengungguli model tengah prior pada enam hingga delapan dari delapan realisasi noise baru pada model acuan yang sama. Pada observasi fitting, kecocokan aktif tetap sensitif terhadap seed. Pemulihan parameter memberi hasil berbeda: rata-rata $V_{p,3}$ berada di bawah acuan pada keempat run, dan hanya satu rentang kuantil 95% yang mencakup acuan. Pengujian manfaat fusi bagi resolusi $V_p$ dan kedalaman muka air tanah memerlukan refit per kanal dan pengulangan pada model acuan yang berbeda.

## Status hipotesis dan bukti

| Hipotesis | Bukti yang tersedia | Uji yang masih diperlukan |
| :-- | :-- | :-- |
| **H1.** Fusi menambah resolusi parameter. | Prediksi gabungan B membaik terhadap model tengah prior pada sebagian besar realisasi noise baru. Efek fusi pada resolusi parameter belum diukur. | Jalankan ulang inversi aktif, pasif, dan gabungan dengan prior, skala misfit, dan desain yang sama. Nilai prediksi kedua kanal pada observasi baru dan ukur galat parameter terhadap model acuan. |
| **H2.** Spektrum di luar jendela Rayleigh mempersempit ketidakpastian $V_p$. | Keempat rata-rata $V_{p,3}$ di bawah nilai acuan. Satu dari empat rentang kuantil 95% memuatnya pada observasi fitting yang sama. | Uji spektrum penuh dan jendela Rayleigh dengan skala Wasserstein yang sama. Verifikasi atribusi moda dan ukur galat, lebar, serta cakupan interval pada beberapa model acuan. |
| **H3.** Campuran Gaussian memberi pendekatan posterior yang lebih baik. | Empat run campuran menguji sensitivitas optimasi terhadap seed. Keunggulan atas satu Gaussian belum diukur. | Bandingkan campuran dan satu Gaussian dengan anggaran komputasi setara. Ukur akurasi terhadap pembanding MCMC dengan diagnostik eksplorasi moda. |
| **H4.** Nested 2D memberi trade-off cakupan radial dan azimut. | Nested mengisi 89,7% bin radial, dengan momen azimut 0,104 (linear 1,000, cincin 0,000). NCC terhadap citra bersih masing-masing adalah 0,7115 untuk nested dan 0,9789 untuk cincin. | Ulangi pada beberapa model dan noise, gunakan referensi fisik bersama, dan ukur pemulihan parameter dari setiap geometri. |

## Tahapan berikutnya

1. Uji sensitivitas $\eta$ dan konvergensi kuadratur pada kanal pasif dan model dalam ruang prior. Ulangi inversi B untuk menilai dampaknya pada $V_p$ dan kedalaman batas lapisan.
2. Uji sensitivitas hyperparameter, multi-start, multi-seed, serta kalibrasi posterior pada dataset berulang.
3. MCMC dengan diagnostik dan eksplorasi moda yang memadai; bandingkan biaya per
   akurasi, bukan jumlah iterasi semata. Tidak ada target speedup yang dianggap hasil.
4. Validasi data lapangan terhadap logging dan piezometer, dengan pengujian efek
   heterogenitas lateral dan ketidakpastian parameter nuisance.
5. Publikasi metode dan perangkat lunak setelah bukti tersebut tersedia.

\clearpage

# Sumber dan keterlacakan

Jalur di bawah relatif terhadap folder proyek. Folder hasil utama adalah `results_v6/`. Tabel konfigurasi, ringkasan metrik, dan gambar berasal dari artefak run yang sudah selesai.

| Isi naskah | Berkas sumber dan penghasil |
|:--|:--|
| Parameter dan observasi A/B | Berkas `benchmark_*_observations.npz` di hasil utama. Pembuatnya `models.py` dan `run_benchmark.py`. |
| Konfigurasi, geometri, kernel dan metrik utama | Berkas `benchmark_verification_summary.json` dan `run_manifest.json` di hasil utama. Pembuatnya `run_benchmark.py`, `acquisition_design.py`, dan `posterior_prediction.py`. |
| Sampel inferensi dan gambar utama | Berkas `benchmark_*_inference.npz` dan gambar `fig0*` serta `benchmark_*` di hasil utama. Pembuat gambar `visualize.py`. |
| Panel B dengan judul revisi | Berkas `benchmark_b_*_v6.png` di `output/pdf/figures/`. Skrip `render_b_labels.py` di folder yang sama. |
| Restart B dan evaluasi noise baru | Folder `diagnostics_review_v6/benchmark_b_restarts/`. Lihat `RINGKASAN.md`, `seed_*/summary.json` dan `new_noise_evaluation.json`. Skrip `run.py`, `summarize.py`, dan `check_new_noise.py` di folder yang sama. |
| Resolusi dan gelombang | Folder restart B. Lihat `resolution_30_42Hz.json` dan `independent_waveform_30_42Hz.json`. Skrip `check_resolution.py` dan `check_waveform.py`. |
| Verifikasi aktif COMPOSTI | Folder `diagnostics_review_v6/benchmark_composti/`. Lihat `summary.json`, `benchmark_*_comparison.npz`, dan `RINGKASAN.md`. Skrip `compare.py` mereproduksi hasil dan mencatat commit upstream. |
| Compliance, replay dan tes | Folder `results_v6_validation/`. Lihat `independent_physics_validation.json`, `benchmark_*_saved_verification.json` dan `pytest_results.txt`. Pemeriksa `validate_reassessment.py`, `verify_saved_results.py` dalam `diagnostics_review_v6/`, dan `tests/`. |


Salinan kode yang digunakan untuk run utama tersimpan di `results_v6/source_snapshot/`. Manifest menyimpan hash sumber dan `results_v6_validation/verified_artifact_sha256.json` menyimpan hash artefak yang diperiksa.

\clearpage

# Daftar pustaka

- Xi, C., Xia, J., Mi, B., Dai, T., Liu, Y., dan Ning, L. (2021). Modified frequency–Bessel transform method for dispersion imaging of Rayleigh waves from ambient seismic noise. *Geophysical Journal International*, 225(2), 1271–1280. [Artikel](https://academic.oup.com/gji/article/225/2/1271/6081096).
- Zhao, X., dan Curtis, A. (2026). Bayesian geophysical inversion using Stein mixture variational inference for multimodal or complex solutions. Pracetak arXiv:2610.23603, versi 1. [Pracetak](https://arxiv.org/abs/2610.23603v1).
- Valentine, A. P., dan Sambridge, M. (2022). pyprop8: A lightweight code to simulate seismic observables in a layered half-space. *Journal of Open Source Software*, 7(76), 4217. [Artikel](https://joss.theoj.org/papers/10.21105/joss.04217).
- Mreyen, A.-S., Michel, H., dan Nguyen, F. (2026). Full wavefield surface wave analysis with Bayesian Evidential Learning. *Geophysical Journal International*, 244(2), ggaf498. [Artikel](https://doi.org/10.1093/gji/ggaf498) · [Kode dan data reproduksi](https://zenodo.org/records/17807710).
- Niskanen, M., dan Lähivaara, T. (2023). COMPOSTI: A Python-based program for seismic trans-dimensional inversion. *SoftwareX*, 21, 101298. [Artikel](https://doi.org/10.1016/j.softx.2022.101298).
- Luu, K. *disba*, versi 0.7.0. Perangkat lunak untuk perhitungan dispersi gelombang permukaan. [Repositori pengembang](https://github.com/keurfonluu/disba).

- Forbriger, T. (2003). Inversion of shallow-seismic wavefields: I. Wavefield transformation. *Geophysical Journal International*, 153(3), 719–734. [Artikel](https://doi.org/10.1046/j.1365-246X.2003.01929.x).
- Foti, S., dkk. (2018). Guidelines for the good practice of surface wave analysis: a product of the InterPACIFIC project. *Bulletin of Earthquake Engineering*, 16, 2367–2420. [Artikel](https://doi.org/10.1007/s10518-017-0206-7).
- Liu, T., Sun, C.-Y., dan Cai, R.-Q. (2024). Forward modeling and wave field analysis of seismic surface waves and guided P-waves. *Geophysical and Geochemical Exploration*, 48(4), 986–995. [Artikel](https://doi.org/10.11720/wtyht.2024.1387).
- Goodman, J., dan Weare, J. (2010). Ensemble samplers with affine invariance. *Communications in Applied Mathematics and Computational Science*, 5(1), 65–80. [Artikel](https://doi.org/10.2140/camcos.2010.5.65).

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\phantomsection
\addcontentsline{toc}{section}{Lampiran A. Gambar hasil}
\section*{Lampiran A. Gambar 1. Geometri akuisisi}
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{results_v6/fig0_field_acquisition_design_comparison.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 2. Perbandingan citra dispersi}
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{results_v6/fig0b_dispersion_image_comparison.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 3. Gelombang dan fusi spektrum - Benchmark A}
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{results_v6/benchmark_a_fig1_waveforms_and_fusion.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 4. Spektrum dan riwayat optimasi - Benchmark A}
\noindent\small Garis putus-putus menunjukkan kurva dispersi Rayleigh model acuan dari \texttt{disba}.\par
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{results_v6/benchmark_a_fig2_fc_spectrum_comparison.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 5. Profil dan distribusi parameter estimasi - Benchmark A}
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{results_v6/benchmark_a_fig3_posterior_profiles_and_marginals.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 6. Gelombang dan fusi spektrum - Benchmark B}
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{output/pdf/figures/benchmark_b_gelombang_v6.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 7. Spektrum dan riwayat optimasi - Benchmark B}
\noindent\small Garis putus-putus menunjukkan kurva dispersi Rayleigh model acuan dari \texttt{disba}.\par
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{output/pdf/figures/benchmark_b_spektrum_v6.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 8. Profil dan distribusi parameter estimasi - Benchmark B}
\begin{center}
\includegraphics[width=\linewidth,height=.80\textheight,keepaspectratio]{\detokenize{output/pdf/figures/benchmark_b_profil_v6.png}}
\end{center}
\end{landscape}
```

`\clearpage`{=latex}

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 9. Perbandingan jejak vertikal dan citra MFJ dengan COMPOSTI - Benchmark A}
\noindent\small Baris atas menggunakan skala yang sama untuk semua kurva pada setiap offset. Panel v6 diselaraskan memakai $\eta=0$ dengan rentang slowness dan taper COMPOSTI.\par
\begin{center}
\includegraphics[width=\linewidth,height=.78\textheight,keepaspectratio]{\detokenize{diagnostics_review_v6/benchmark_composti/benchmark_a_composti.png}}
\end{center}
\end{landscape}
```

```{=latex}
\begin{landscape}
\thispagestyle{empty}
\section*{Gambar 10. Perbandingan jejak vertikal dan citra MFJ dengan COMPOSTI - Benchmark B}
\noindent\small Baris atas menggunakan skala yang sama untuk semua kurva pada setiap offset. Panel v6 diselaraskan memakai $\eta=0$ dengan rentang slowness dan taper COMPOSTI.\par
\begin{center}
\includegraphics[width=\linewidth,height=.78\textheight,keepaspectratio]{\detokenize{diagnostics_review_v6/benchmark_composti/benchmark_b_composti.png}}
\end{center}
\end{landscape}
```

# Lampiran B. Diagram alir dan pseudocode per metode {#lampiran-metode}

Lampiran ini mengikuti eksperimen sintetik v6, mulai dari geometri hingga evaluasi. Setiap halaman metode memuat alur, input, output, dan pseudocode. Nama prosedur ditulis untuk membantu pembacaan. Nama berkas sumber menunjukkan lokasi implementasinya.

```{=latex}
\alur{Geometri dan\newline model acuan}{Observasi sintetis\newline dan citra MFJ}{Fusi dan inversi\newline campuran SMVI}{Prediksi dan\newline evaluasi}
```

**Jalur observasi.** Model acuan dan geometri menghasilkan gather aktif dan pasif. Gangguan sintetik ditambahkan sebelum pencitraan. Fusi kedua citra menghasilkan spektrum observasi $D_{\mathrm{obs}}$.

**Jalur kandidat.** SMVI mengusulkan parameter laten $z$. Setiap kandidat diterjemahkan menjadi model fisik $m$, lalu melewati operator sintesis, MFJ, dan fusi untuk menghasilkan $D_{\mathrm{syn}}$. Skor kecocokan terhadap $D_{\mathrm{obs}}$ memberi arah pembaruan komponen.

**Jalur evaluasi.** Campuran akhir menghasilkan sampel parameter dan prediksi spektrum. Evaluasi B memakai restart pada observasi yang sama, prediksi tetap pada noise baru, serta uji terarah resolusi dan gelombang. Verifikasi aktif A/B mencakup pembanding COMPOSTI.

| Bagian | Metode | Keluaran utama |
|--:|:--|:--|
| [B.1](#metode-b1) | Geometri akuisisi | Offset, azimut, dan jumlah pasangan |
| [B.2](#metode-b2) | Sintesis P–SV dan gangguan | Gather kompleks $U(f,r)$ |
| [B.3](#metode-b3) | Pencitraan MFJ dan kernel pembanding | Citra $I(f,c)$ |
| [B.4](#metode-b4) | Normalisasi dan fusi | Spektrum $D(f,c)$ dan bobot kanal |
| [B.5](#metode-b5) | Prior, NCC dan Wasserstein | Target dan gradien di ruang laten |
| [B.6](#metode-b6) | SMVI | Campuran Gaussian dan riwayat optimasi |
| [B.7](#metode-b7) | Prediksi dan ringkasan parameter | Spektrum prediksi, kuantil, dan galat MC |
| [B.8](#metode-b8) | Evaluasi dan verifikasi | Pembanding, sensitivitas dan batas kesimpulan |

**Notasi.** $f$ adalah frekuensi, $r$ offset, $c$ kecepatan fase, dan $k$ bilangan gelombang. $m$ berisi empat $V_s$, empat $V_p$, dan tiga ketebalan. $\rho$, $Q_s$, dan $Q_p$ diberikan untuk setiap lapisan sintetik.

`\clearpage`{=latex}

## B.1. Geometri akuisisi {#metode-b1}

**Tujuan.** Membangun dua kumpulan offset untuk operator maju dan mendeskripsikan cakupan setiap rancangan pasif.

```{=latex}
\alur{Koordinat stasiun\newline dan tembakan}{Hitung jarak\newline serta azimut}{Kelompokkan offset\newline yang berimpit}{Offset, jumlah\newline pasangan, azimut}
```

**Input.** Posisi tembakan dan penerima aktif, koordinat stasiun pasif, serta pilihan pengelompokan jarak. Run utama memakai 48 kanal aktif dan rancangan pasif nested dengan 25 stasiun.

**Output.** Offset aktif $r_a$, offset pasif $r_p$, jumlah pasangan $N_r$, dan faktor arah $a_r$. Metrik rancangan meliputi jarak minimum dan maksimum, bin radial terisi, serta momen azimut orde dua.

```text
GEOMETRI(tembakan, penerima_aktif, stasiun_pasif)
    r_a <- urutkan offset absolut unik dari seluruh tembakan

    untuk setiap pasangan stasiun pasif i < j:
        r <- jarak Euclidean stasiun i dan j
        theta <- sudut vektor dari stasiun i ke j
        simpan (r, theta) jika r > 1 m

    kelompok <- gabungkan jarak yang berimpit menurut presisi v6
    untuk setiap kelompok jarak r:
        N_r <- jumlah pasangan dalam kelompok
        a_r <- rata-rata cos(2 * (theta - 25 derajat))

    kembalikan r_a, jarak unik pasif, N_r, a_r

METRIK_RANCANGAN(daftar_jarak, daftar_azimut)
    hitung jarak minimum dan maksimum
    hitung proporsi 29 bin radial 5 sampai 585 m yang terisi
    momen_orde_2 <- abs(rata-rata exp(2i * azimut))
    kembalikan metrik geometri
```

Pada run utama, opsi bin radial kasar dinonaktifkan. Bin histogram pada metrik rancangan hanya dipakai untuk merangkum cakupan jarak. Offset dan faktor arah masuk ke metode 2, sedangkan offset juga menentukan bobot pencitraan pada metode 3.

**Implementasi.** `forward_mfj.py` dan `acquisition_design.py`.

`\clearpage`{=latex}

## B.2. Sintesis P–SV dan gangguan {#metode-b2}

**Tujuan.** Menghitung gather kompleks dari model berlapis untuk setiap kanal. Jalur observasi menambahkan gangguan, sedangkan jalur kandidat memakai gather bersih.

```{=latex}
\alur{Model berlapis\newline dan grid}{Compliance P--SV\newline dan integral}{Spektrum sumber\newline serta gangguan}{Gather kompleks\newline tiap kanal}
```

**Input.** $m$, densitas dan $Q$, offset, frekuensi, spektrum sumber, serta seed jika membuat observasi. Sumber aktif memakai skala frekuensi 25 Hz, dan kanal pasif memakai bentuk spektral dengan skala 6 Hz.

**Output.** $U_a(f,r)$ dan $U_p(f,r)$ pada grid frekuensi masing-masing.

```text
SINTESIS(model, kanal, f, r, gangguan, rng)
    untuk setiap frekuensi f_i:
        bentuk grid k linier–logaritmik dan bobot trapesium
        untuk setiap k:
            mulai dari impedansi half-space
            lakukan rekursi ke atas melalui lapisan berhingga
            ambil compliance vertikal Gzz(f_i, k)
        kurangi ekor statik A dari k * Gzz(f_i, k)
        integralkan dengan kernel J0(k * r), pulihkan suku A / r
        kalikan dengan spektrum sumber S_kanal(f_i)

    jika gangguan dipakai dan kanal aktif:
        kalikan amplitudo log-normal dan fase acak per offset
    jika gangguan dipakai dan kanal pasif:
        modulasi fase menurut jarak, frekuensi, dan faktor arah a_r
    jika gangguan digunakan:
        tambah noise Gaussian kompleks berskala RMS per frekuensi
        gunakan 0,08 untuk aktif dan 0,045 / sqrt(max(N_r,1)) untuk pasif

    kembalikan gather kompleks
```

Batas grid dan bobot dijabarkan pada [formulasi P–SV](#green-psv). Observasi A dan B memakai seed 2026 dan 2027. Gangguan pasif bersifat fenomenologis. Jejak waktu dalam lampiran dihitung pada grid FFT tersendiri untuk ilustrasi gelombang. Inversi memakai gather kompleks pada grid frekuensi native.

**Implementasi.** `forward_mfj.py`. Persamaan sumber dan gangguan lengkap tercantum dalam bagian pembentukan observasi pada naskah utama.

`\clearpage`{=latex}

## B.3. Pencitraan MFJ dan kernel pembanding {#metode-b3}

**Tujuan.** Mengubah getaran terhadap offset menjadi citra frekuensi–kecepatan fase dengan menguji penjumlahan koheren pada setiap kecepatan kandidat.

```{=latex}
\alur{Gather kompleks\newline dan offset}{Whitening, taper,\newline bobot radial}{Proyeksi kernel\newline pada grid c}{Citra per kanal\newline pada grid f--c}
```

**Input.** Gather $U(f,r)$, offset $r$, grid $c$, taper ruang, dan kernel. $J_0$ dan $H_0^{(2)}$ diterapkan pada argumen $2\pi f r/c$.

**Output.** Citra $I_a(f,c)$ dan $I_p(f,c)$ yang dinormalisasi per frekuensi.

```text
CITRA(U, f, r, c, jenis_kernel)
    Delta_r <- bobot trapesium offset
    taper <- jendela ruang pada offset

    untuk setiap frekuensi f_i:
        eps <- 0,02 * max_offset(abs(U[f_i])) + 1e-30
        u <- U[f_i] / (abs(U[f_i]) + eps)
        u <- taper * u * sqrt(r) * Delta_r
        untuk setiap kecepatan c_j:
            x <- 2 * pi * f_i * r / c_j
            jika MFJ: K <- Hankel_2_orde_0(x) / abs(Hankel_2_orde_0(x))
            jika plane-wave: K <- exp(-i * x)
            jika Bessel berdiri: K <- J0(x)
            I[f_i, c_j] <- abs(jumlah_offset(u * K))
        bagi I[f_i] dengan maksimum irisan jika maksimum positif

    kembalikan I
```

Pada uji perbandingan, ketiga kernel dipakai pada geometri linear. MFJ juga dipakai pada rancangan L-array, cincin, dan nested. Pemrosesan sebelum proyeksi sama; perbedaan kernel mengubah cara fase dan amplitudo antar-offset digabungkan. Citra bersih pembanding dibuat dengan kombinasi geometri dan kernel yang sama. Kurva Rayleigh dari `disba` digunakan untuk overlay model acuan pada gambar benchmark.

**Implementasi.** `forward_mfj.py`, `acquisition_design.py`, dan `visualize.py`.

`\clearpage`{=latex}

## B.4. Normalisasi dan fusi aktif–pasif {#metode-b4}

**Tujuan.** Membentuk spektrum gabungan dengan dukungan pasif di frekuensi rendah, aktif di frekuensi tinggi, serta transisi pada 8–12 Hz.

```{=latex}
\alur{Dua citra MFJ\newline pada grid native}{Normalisasi\newline per frekuensi}{Pilih node cocok\newline dan hitung bobot}{Citra gabungan\newline dan bobot kanal}
```

**Input.** Citra per kanal, frekuensi native, pangkat normalisasi $p$, pita fusi 1–55 Hz, dan overlap 8–12 Hz.

**Output.** $D(f,c)$, citra aktif dan pasif pada grid gabungan, bobot $w_a(f),w_p(f)$, serta proksi bentuk tiap kanal.

```text
FUSI(I_a, I_p, f_a, f_p, p, overlap=[8,12])
    untuk masing-masing kanal:
        normalisasi min-max tiap irisan dan pangkatkan dengan p

    f_bb <- node pasif <= 8
            + node yang sama pada kedua kanal untuk 8 < f < 12
            + node aktif >= 12
    urutkan, ambil node unik, batasi ke 1 sampai 55 Hz
    salin hanya irisan yang memiliki node native yang cocok
    normalisasi min-max kembali pada irisan yang tersedia

    untuk setiap f pada f_bb:
        alpha <- transisi sin^2 dari pasif menuju aktif
        v_j <- max((mean(D_j) / (max(D_j) + eps))^2, 1e-4)
        jika satu kanal tersedia: gunakan bobot penuh pada kanal itu
        jika keduanya tersedia:
            normalkan alpha/v_a dan (1-alpha)/v_p agar berjumlah 1
        D[f] <- minmax(w_a * D_a[f] + w_p * D_p[f])

    kembalikan D, D_a, D_p, w_a, w_p
```

Gunakan $p=1{,}35$ untuk inversi A/B dan $p=1{,}40$ untuk perbandingan geometri dan kernel. Irisan konstan dipetakan ke nol oleh normalisasi. Proksi $v_j$ merangkum bentuk spektral; satuan dan interpretasinya berbeda dari varians noise pengukuran. Metode yang sama menghasilkan $D_{\mathrm{obs}}$ dan $D_{\mathrm{syn}}$.

**Implementasi.** `fusion.py`.

`\clearpage`{=latex}

## B.5. Prior, NCC dan Wasserstein {#metode-b5}

**Tujuan.** Mengukur kecocokan citra kandidat dan membentuk target untuk optimasi di ruang laten.

```{=latex}
\alur{Parameter laten\newline dan batas prior}{Transformasi dan\newline prediksi spektrum}{Mask, NCC dan\newline transport $W_2$}{Misfit, skor target,\newline serta gradien}
```

**Input.** $z$, prior fisik, $D_{\mathrm{obs}}$, geometri dan grid, serta pengaturan mask dan misfit.

**Output.** Misfit $J$, NCC, $W_2$ agregat, log target laten dan gradiennya. Transformasi sigmoid menerapkan batas fisik dan syarat $V_p\geq1.45V_s$ per lapisan.

```text
TARGET(z, prior, D_obs, grid, mask)
    m, log_det <- transformasi_sigmoid_dan_Jacobian(z, prior)
    D_syn <- FUSI(CITRA(SINTESIS(m, gangguan=tidak)))
    hitung bobot trapesium pada dukungan frekuensi terpilih

    untuk setiap frekuensi f dengan bobot positif:
        terapkan mask kecepatan yang sama pada observasi dan prediksi
        untuk f <= 12 Hz, pertahankan c <= 680 m/s
        a, b <- irisan observasi dan prediksi yang dipertahankan
        L_NCC <- 1 - cosine_similarity(a, b)
        p, q <- normalkan massa max(a,0)^2 + eps dan max(b,0)^2 + eps
        W2_sq <- biaya kuadratik transport monoton p ke q pada grid c

    J <- (1 - lambda) * mean_berbobot(L_NCC)
         + lambda * mean_berbobot(W2_sq) / skala_W2^2
    skor <- -J / (2 * sigma_D^2) + log_det
    gradien <- autodiff(skor terhadap z)
    kembalikan J, NCC, W2, skor, gradien
```

NCC menilai kesamaan bentuk irisan. Wasserstein memberi penalti perpindahan massa sepanjang sumbu kecepatan. V6 memakai $\lambda=0.35$, $\sigma_D=0.08$, dan skala transport 169,2 m/s untuk A serta 348 m/s untuk B. Saat membatasi jendela pada model yang sama, pertahankan skala dari grid induk. Nilai $\sigma_D$ mengatur kekuatan target komposit.

**Implementasi.** `models.py`, `smvi_inversion.py` dan `autodiff_likelihood.py`.

`\clearpage`{=latex}

## B.6. Optimasi SMVI {#metode-b6}

**Tujuan.** Mengestimasi campuran Gaussian berkorelasi di ruang laten. Setiap komponen mempunyai pusat dan faktor Cholesky, dengan bobot sama.

```{=latex}
\alur{Target laten\newline dan seed}{Inisialisasi dan\newline refinemen pusat}{Skor, interaksi Stein\newline dan langkah Adam}{Campuran akhir\newline dan riwayat}
```

**Input.** Target metode 5, seed, jumlah komponen $K$, anggaran iterasi, dan konfigurasi Adam. V6 memakai $K=4$, 60 iterasi, serta empat sampel gradien per komponen per iterasi.

**Output.** Pusat $\mu_k$, faktor Cholesky $L_k$, bobot $1/K$, dan riwayat ELBO serta misfit.

```text
SMVI(target, seed, K, anggaran)
    kandidat <- Latin_hypercube_dalam_koordinat_sigmoid(seed)
    petakan kandidat ke laten dan hitung skor target
    pilih K pusat berskor tinggi dengan aturan jarak awal
    perhalus tiap pusat dengan L-BFGS-B
    inisialisasi Cholesky dan tetapkan bobot 1/K

    untuk setiap iterasi sampai anggaran habis:
        tarik sampel normal yang dipakai bersama antar-komponen
        untuk setiap komponen k dan setiap sampel eps:
            z <- mu[k] + L[k] * eps
            g <- grad_log_target(z) - grad_log_campuran(z)
            akumulasi g untuk pusat dan outer(g, eps) untuk Cholesky
        rata-ratakan gradien; terapkan aturan rantai log-diagonal
        v <- arah Stein dengan kernel RBF pada parameter komponen
        perbarui pusat dan Cholesky menggunakan Adam dan arah v
        batasi log-diagonal Cholesky, catat ELBO dan misfit

    kembalikan mu, L, bobot, riwayat
```

Kernel Stein memakai pusat, log-diagonal, dan bagian bawah Cholesky di luar diagonal. Inisialisasi v6 memakai 32 kandidat, dengan batas 30 iterasi refinemen per pusat. Batas 60 iterasi mengatur anggaran optimasi. Kestabilan prosedur diperiksa melalui restart pada metode 8.

**Implementasi.** `SteinMixtureVariationalInference` di `smvi_inversion.py`.

`\clearpage`{=latex}

## B.7. Ringkasan parameter dan prediksi spektrum {#metode-b7}

**Tujuan.** Merangkum campuran yang sudah diestimasi dan mengukur presisi Monte Carlo dari spektrum prediksi rata-rata.

```{=latex}
\alur{Campuran akhir\newline di ruang laten}{Sampel fisik dan\newline prediksi per sampel}{Gabungkan, bootstrap,\newline periksa presisi}{Profil, spektrum,\newline dan galat baku}
```

**Input.** Campuran akhir, transformasi fisik, operator prediksi, anggaran 128–512 draw, dan toleransi 0,01.

**Output.** Sampel parameter, kuantil, spektrum prediksi rata-rata, sebaran spektrum, galat baku, dan status presisi.

```text
RINGKAS_DAN_PREDIKSI(campuran, operator, konfigurasi)
    tarik 1200 sampel parameter dan transformasikan ke fisik
    hitung rata-rata parameter dan kuantil per parameter
    hitung spektrum pada parameter rata-rata
    nilai semua pusat dan pilih pusat dengan J terkecil

    mulai dengan jumlah draw yang seimbang per komponen
    ulangi:
        tambahkan draw hingga anggaran tahap ini terpenuhi
        untuk setiap draw: hitung D_syn dari model fisik
        D_pred <- jumlah berbobot mean spektrum per komponen
        hitung RMS galat baku citra dan galat baku J (bootstrap)
        hitung selisih J antara setengah sampel dan semua sampel
        jika kedua galat baku <= toleransi
          dan selisih J <= 2*toleransi: catat target tercapai, berhenti
        jika anggaran maksimum tercapai: catat status, berhenti
        gandakan anggaran hingga batas maksimum

    kembalikan ringkasan parameter dan diagnostik prediksi
```

$\operatorname{mean}[D(m)]$ dan $D(\operatorname{mean}[m])$ dapat berbeda karena operator maju nonlinier. Pusat yang dipilih berdasarkan $J$ adalah pembanding di antara pusat campuran. Ketelitian Monte Carlo pada spektrum mengukur stabilitas estimasi prediksi bersyarat pada campuran yang sudah dioptimasi.

**Implementasi.** `posterior_prediction.py` dan `smvi_inversion.py`.

`\clearpage`{=latex}

## B.8. Evaluasi dan verifikasi v6 {#metode-b8}

**Tujuan.** Menguji kecocokan, pemulihan parameter, sensitivitas terhadap titik awal dan noise, serta pemeriksaan numerik pada pita aktif yang bermasalah.

```{=latex}
\alur{Observasi, acuan,\newline prior dan prediksi}{Bandingkan misfit\newline dan profil}{Uji seed, noise,\newline grid dan propagator}{Ringkasan bukti\newline dan batas validasi}
```

**Input.** Observasi v6, model acuan, prediksi tersimpan, model tengah prior, seed tambahan, grid pembanding, propagator PyProp8, dan solver COMPOSTI.

**Output.** Metrik per kanal dan pita, galat parameter, variasi antar-seed, serta batas interpretasi hasil.

```text
VERIFIKASI(observasi, acuan, prior, hasil_v6)
    nilai prediksi, model tengah prior, dan acuan pada observasi
    uraikan selisih misfit B menurut pita dan kanal
    evaluasi pusat, prediksi per komponen dan campuran tanpa satu komponen

    untuk seed_fit dalam [2028, 2029, 2030]:
        latih ulang SMVI pada observasi B yang sama
        hitung prediksi dan bandingkan dengan prior serta acuan

    untuk seed_noise dari 3101 sampai 3108:
        buat observasi baru pada model acuan B yang sama
        nilai semua prediksi tetap dan model tengah prior

    bandingkan misfit model tetap pada grid bertingkat (30-42 Hz)
    bandingkan jejak berlapis dengan PyProp8 di empat offset
    untuk model acuan A dan B:
        sintesis jejak aktif COMPOSTI dan v6 dengan sumber identik
        hitung selisih spektrum, jejak, dan citra MFJ relatif L2
        ulangi v6 dengan eta = 0, kemudian samakan rentang dan taper
    simpan metrik, konfigurasi, sampel, dan referensi artefak
```

**Cakupan.** Evaluasi noise baru menguji prediksi dari hasil fitting yang sudah ada. Uji pemulihan pada beberapa model acuan, refit aktif/pasif/gabungan, dan pembanding posterior MCMC masih merupakan agenda penelitian. Perbandingan kurva Rayleigh memeriksa fitur dispersi, sedangkan PyProp8 dan COMPOSTI memeriksa respons berlapis pada kondisi uji masing-masing.

**Implementasi.** `validate_reassessment.py` dan dua subfolder di `diagnostics_review_v6/`, yaitu `benchmark_b_residual/` dan `benchmark_b_restarts/`. Skrip COMPOSTI tersedia di `diagnostics_review_v6/benchmark_composti/compare.py`.

**Rujukan dalam naskah.** Bagian utama memuat persamaan, pengaturan eksperimen, dan hasil numerik.


\clearpage

# Lampiran C. Landasan, protokol, dan keterlacakan v3 {#lampiran-c}

## C.1. Hubungan fisika gelombang dengan implementasi

**Sumber, batas, dan respons P–SV.** Untuk gaya titik vertikal pada medium isotropik berlapis, simetri silindris mereduksi masalah menjadi kopling P dan SV untuk setiap bilangan gelombang horizontal $k$. Kondisi batas pada antarmuka adalah kontinuitas perpindahan dan traksi. Half-space memakai solusi yang meluruh atau meradiasikan energi ke bawah. Di permukaan, sumber memberi traksi normal dengan traksi geser nol. Jika $K$ memetakan perpindahan ke traksi, compliance vertikal diperoleh dari elemen vertikal invers $K$. Pada [formulasi P–SV](#green-psv), tandanya mengikuti konvensi traksi yang dipakai kode.

Partisi yang dipakai dalam rekursi v6 dapat ditulis dengan $\gamma=2k^2-\tilde\omega^2/\tilde V_s^2$. Tanda atas berlaku untuk gelombang turun, dan tanda bawah untuk gelombang naik.

$$U_{d/u}=\begin{pmatrix}k & \pm\nu_s\\ \pm\nu_p & k\end{pmatrix},
\qquad
T_{d/u}=-\tilde\mu\begin{pmatrix}\pm2k\nu_p & \gamma\\ \gamma & \pm2k\nu_s\end{pmatrix}.$$

Pada half-space, $K_N=T_dU_d^{-1}$ menjadi kondisi awal rekursi.

**Integral dan moda.** V6 mengintegrasikan $kG_{zz}$ setelah pengurangan asimtot statik Boussinesq dan mengembalikan suku statik secara analitik. Untuk normalisasi perpindahan per satuan gaya titik dengan konvensi Fourier 2D standar, keluaran sintesis dibagi $2\pi$, seperti pada [pembandingan COMPOSTI](#uji-composti). Moda normal dapat dipelajari melalui akar kondisi batas homogen. Moda bocor memerlukan pemeriksaan radiasi dan struktur spektral yang sesuai. Citra MFJ yang terang pada kecepatan tinggi memberi hipotesis untuk diuji dengan pelacakan cabang atau respons terhadap perubahan ketebalan, $V_p$, dan offset.

**Fase dan waktu tiba.** Untuk konvensi maju $e^{-i\omega t}$, fase gelombang keluar berperilaku seperti $e^{ikr}$. Evaluasi dengan konvensi invers NumPy memakai spektrum positif yang dikonjugasi, lalu diberi fase tunda $e^{-i\omega t_0}$. Kode memakai $t_0=0{,}045$ s pada ilustrasi aktif dan $0{,}12$ s pada ilustrasi pasif. Untuk cabang terisolasi dengan fase $k(\omega)r$, pendekatan fase stasioner memberi waktu amplop

$$t_{\mathrm{amplop}}\approx t_0 + r\,\frac{dk}{d\omega}=t_0+\frac r{U_g},\qquad
U_g=\frac{c}{1-\frac{f}{c}\frac{dc}{df}},\quad c=\frac{\omega}{k}.$$

Citra $f$–$c$ menampilkan kecepatan fase, sementara kemiringan amplop jejak berkaitan dengan kecepatan grup. Interferensi beberapa cabang dapat menyulitkan pembacaan amplop.

**Tiga jenis sampling.** Jejak ilustratif disintesis langsung pada grid FFT dengan $\Delta f=1/(N_t\Delta t)$. Durasi periodik FFT perlu mencakup kedatangan dan gelombang susulan. Interpolasi spektrum kompleks pada grid jarang berisiko salah mengikuti fase ketika perubahan antar-node melebihi $\pi$. Resolusi kuadratur $k$ mengatur akurasi integral maju. Sampling receiver $\Delta x$ pada array reguler memberi batas Nyquist spasial $c\ge 2\Delta x f$. Ketiganya perlu diperiksa secara tersendiri.

\clearpage

## C.2. Arsitektur kode dan evolusi eksperimen

| Tahap | Modul dan fungsi | Input → output |
| :-- | :-- | :-- |
| Pengaturan eksperimen | `run_benchmark.py` | Konfigurasi, seed dan model A/B → observasi, hasil inversi, gambar dan manifest |
| Model dan geometri | `models.py`, `acquisition_design.py` | Parameter lapisan, prior dan koordinat → model dan ukuran cakupan |
| Sintesis dan pencitraan | `forward_mfj.py` | Model, sumber dan offset → gather kompleks, citra MFJ aktif dan pasif |
| Fusi | `fusion.py` | Citra dan frekuensi native → citra gabungan dan bobot per kanal |
| Inferensi | `smvi_inversion.py`, `autodiff_likelihood.py` | Model kandidat, prior dan observasi → target, gradien, dan campuran Gaussian |
| Prediksi, visualisasi dan pemeriksaan | `posterior_prediction.py`, `visualize.py`, `validate_reassessment.py` | Campuran, observasi dan referensi → prediksi, gambar, galat numerik |
| Uji tambahan | `ablation_experiments.py` dan `diagnostics_review_v6/` | Kanal terpisah, mask, pembanding atau noise tambahan → uji terkontrol dan diagnostik |

**Pustaka dan peran.** Numba digunakan pada jalur sintesis dan pencitraan. PyTorch menyediakan diferensiasi otomatis untuk optimasi. `disba` menghasilkan kurva pembanding Rayleigh model acuan. PyProp8 dan COMPOSTI dipakai untuk verifikasi independen dengan lingkup yang dilaporkan. Adapter dan pustaka yang tersimpan dalam arsip perlu ditelusuri ke manifest sebelum dikaitkan dengan hasil tertentu.

| Generasi | Perubahan metodologis dan pelajaran |
| :-- | :-- |
| v1 | Prototipe awal dengan diskretisasi kasar dan struktur kovariansi yang lebih sederhana. |
| v2–v3 | Penambahan sintesis modal, varian Hankel, kovariansi Cholesky, geometri nested, dan demonstrasi ablasi. Cabang pada sintesis modal lama ditentukan secara pendekatan. |
| v4–v5 | Penilaian ulang sintesis dan tampilan, pemeriksaan independen serta konsistensi grid native. Status dan keterbatasan run terdahulu dicatat dalam `EVALUASI_DAN_PERBAIKAN.md`. |
| v6 dan pemeriksaan lanjutan | Sintesis compliance P–SV dengan kuadratur Sommerfeld, transport kuadratik, prediksi dari sampel campuran, restart, evaluasi noise baru, dan perbandingan COMPOSTI terbatas pada forward aktif. |

Parameterisasi, data dan definisi metrik berubah antargenerasi. Penilaian manfaat perubahan memerlukan pembandingan berpasangan pada input serta definisi metrik yang sama. Arsip v3 menyediakan jejak desain eksperimen dan alasan metodologis. Bukti kuantitatif saat ini terlacak ke run dan diagnostik v6.

\clearpage

## C.3. Protokol pilot akuisisi lapangan

Geometri pilot mengikuti opsi nested 2D yang dikaji pada benchmark. Rancangan referensi memakai 48 kanal aktif dengan spasi 2 m dan enam lokasi sumber pada offset 4, 10, dan 20 m dari tiap ujung. Array pasif memiliki 9 sensor inti dan 16 sensor di bagian luar. Rancangan linear, L-array dan cincin tetap menjadi pembanding. Keputusan lapangan mengikuti keterukuran sinyal, respons instrumen, akses lokasi dan hasil pilot. [Panduan Foti dkk. (2018)](https://doi.org/10.1007/s10518-017-0206-7).

| Komponen | Parameter yang diuji pada pilot | Keluaran pemeriksaan |
| :-- | :-- | :-- |
| Sensor aktif dan digitizer | Respons pada pita sasaran, kopling, waktu trigger, laju sampel, filter dan batas saturasi | Rekaman uji, respons instrumen, pita tanpa clipping |
| Sensor pasif | Respons terkalibrasi, noise instrumen, sinkronisasi, durasi, orientasi sensor jika 3C | Spektrum noise, statistik clock, stabilitas korelasi terhadap durasi |
| Sumber aktif dan pelat | Uji sumber, pelat dan ulangan tembakan dengan geometri tetap | Spektrum sumber dan SNR per pita serta offset yang terukur |
| Posisi array | Koordinat sumber dan receiver, bukaan, distribusi pasangan radial dan azimut | Peta as-built, daftar offset, cakupan bin, jejak perubahan geometri |

1. **Survei pendahuluan.** Dokumentasikan kondisi lokasi, respons sensor, dan kebisingan. Tetapkan pita serta offset yang ingin diukur. Jadwalkan pilot tanpa mengunci merek alat, waktu rekam, atau ambang kualitas sebelum data uji tersedia.
2. **Instalasi dan pemetaan.** Pasang bentangan dan array, ukur posisi aktual, dokumentasikan kopling, periksa sinkronisasi dan uji jalur rekam. Hitung kembali cakupan pasangan dari koordinat aktual.
3. **Rekam uji aktif dan pasif.** Bandingkan konfigurasi sumber dan ulangan pada offset yang sama. Untuk pasif, bandingkan stabilitas korelasi antarsegmen dan menurut arah pasangan. Catat waktu tembakan agar segmen bising ambien dapat dipilih secara konsisten.
4. **QC sebelum pembongkaran.** Periksa clipping, kesalahan waktu, kelengkapan kanal, kurva sensitivitas instrumen, rasio sinyal terhadap gangguan, dan konsistensi citra antarulangan. Laporkan frekuensi overlap yang benar-benar didukung kedua kanal. Simpan data mentah dan metadata sebelum mengubah tata letak.

Nilai sensor, durasi, tebal pelat, dan ambang NCC pada v3 merupakan opsi desain historis. Keandalan konfigurasi pada lokasi penelitian bergantung pada pengukuran pilot. Pada setup sintetik v6, dukungan kanal pasif berakhir pada 12 Hz.

\clearpage

## C.4. Ablasi dan uji ketahanan yang perlu dijalankan

| Pembandingan | Protokol pengujian | Metrik |
| :-- | :-- | :-- |
| Aktif, pasif dan fusi | Jalankan inversi terpisah pada observasi tiap kanal dan gabungan dengan prior serta anggaran yang sama. Skor kedua kanal pada data baru. | Galat profil, galat ketebalan dan kecocokan prediksi |
| Picking dan citra utuh | Terapkan picking pada citra yang sama. Cocokkan kurva yang dipilih dengan pembanding Rayleigh dan dokumentasikan bagian yang ambigu. Bandingkan dengan refit citra penuh. | Keberhasilan pemulihan LVL, galat $V_{s,2}$ dan $h_2$ |
| Jendela Rayleigh dan pita kecepatan luas | Tetapkan mask sebelum fitting, gunakan skala Wasserstein yang sama, ulangi pada beberapa model B dan periksa atribusi moda. | Galat $V_p$, lebar dan cakupan interval, resolusi kedalaman |
| Satu Gaussian, campuran dan MCMC | Optimalkan $K=1$ dan $K>1$ dengan anggaran setara. Jalankan referensi ensemble independen dan periksa ESS serta eksplorasi moda. | Akurasi, statistik antarparameter, biaya per evaluasi dan waktu |

**Uji ketahanan.** Rancang faktor kekuatan sumber, noise, variasi kopling, dan geometri. Ubah spektrum sumber aktif di domain frekuensi–offset sebelum pencitraan. Tambahkan gangguan dengan SNR terukur, jalankan ulang MFJ serta fusi, lalu fit kandidat dengan prosedur yang sama. Gunakan model dan seed berulang untuk menilai perbedaan. Laporkan SNR aktual, bobot tiap kanal, galat pemulihan, dan performa per pita. Pisahkan kasus yang menguji noise saja dari kasus dengan perubahan geometri.

**Status dan sumber selisih dengan v3.** Run v6 yang dilaporkan belum menjalankan rangkaian ilmiah lengkap di atas. Kode arsip v3 menggeser sampel aktif dan pasif setelah inversi, membentuk sampel pembanding Rayleigh dengan generator acak dan membangun pembanding Gaussian dari komponen SMVI. MCMC lama dimulai dari pusat SMVI dan memakai chain singkat. Stress-test v3 memodifikasi citra setelah fusi serta menetapkan bobot transisi melalui kurva yang sudah ditentukan. Detail tersebut berada di `archive_v3_before_scientific_fixes/ablation_experiments.py`. Karena itu angka ablasi, pengurangan ketidakpastian, dan penyelamatan pita dari v3 memiliki status historis. Implementasi saat ini tersedia di `ablation_experiments.py` dan pengujian ilmiahnya perlu dipisahkan dari uji perangkat lunak. [Goodman dan Weare (2010)](https://doi.org/10.2140/camcos.2010.5.65).

\clearpage

## C.5. Peta cakupan naskah lama dan naskah ini

| Bagian v3 | Lokasi atau perlakuan dalam edisi ini |
| :-- | :-- |
| Abstrak Indonesia dan Inggris | Keduanya tersedia di bagian awal, dengan angka dan pembatasan dari hasil v6. |
| Bab I, §1.1–1.4 | Bagian kesenjangan, pertanyaan, tujuan, hipotesis, serta matriks bukti. |
| Bab II, §2.1–2.3b | Formulasi P–SV dan C.1. Rekursi serta asimtot dipertahankan dengan formulasi v6. Sintesis modal heuristik dari arsip dijelaskan sebagai tahap historis. Teori fase, grup dan sampling dikembalikan. |
| Bab II, §2.4–2.7 | Metode komputasi dan Lampiran B. MFJ ditulis sebagai varian berbasis fase. Proksi fusi dan transport kuadratik mengikuti implementasi v6. |
| Bab III, §3.1–3.4 | Model dan properti A/B, observasi, geometri dan kode tercakup dalam metode, konfigurasi benchmark, diagram, Lampiran B dan C.2. |
| Bab III, §3.5–3.6 | C.3 memulihkan topik instrumen dan SOP dalam bentuk rancangan pilot dengan parameter yang harus diukur. |
| Bab IV, §4.1–4.2 | Hasil geometri dan kernel dengan referensi bersih tiap kombinasi serta riwayat desain di C.2. |
| Bab IV, §4.3–4.4 | Benchmark A/B dengan gambar dan estimasi v6. Uji restart, noise baru, PyProp8 dan COMPOSTI memperjelas batas interpretasi. |
| Bab IV, §4.5–4.6 | C.4 memulihkan empat pembandingan ilmiah dan protokol stres. Angka lama terlacak ke arsip, dengan keterbatasan metodologis dijelaskan pada C.4. |
| Bab V dan daftar pustaka | Kesimpulan berbasis hasil v6, tahap lanjutan, serta rujukan teori, instrumen, metode dan verifikasi yang dipakai dalam edisi ini. |

**Sumber historis.** Folder dan berkas berikut menyimpan naskah v3 dengan rincian, angka, dan redaksi asli.

- Folder `archive_v3_before_scientific_fixes/`.
- Berkas `NASKAH_LENGKAP_PENELITIAN_PHD_MFJ_SMVI.md`.
