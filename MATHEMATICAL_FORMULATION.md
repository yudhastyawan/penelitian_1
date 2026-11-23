# Formulasi Matematis MFJ-SMVI: Implementasi Numerik Terkoreksi

Status: eksperimen sintetik; validasi lapangan dan kalibrasi posterior belum selesai.

## 1. Model dan fungsi Green P-SV

Medium terdiri dari tiga lapisan homogen isotropik di atas half-space. Parameter
inversi adalah empat $V_s$, empat $V_p$, dan tiga ketebalan $h$. Densitas dan
faktor kualitas ditetapkan sebagai nuisance parameters yang diketahui dalam
benchmark sintetik; keduanya tidak diestimasi dari data.

Dengan konvensi waktu $e^{-i\omega t}$, kecepatan kompleks adalah
$\tilde V=V(1-i/(2Q))$. Frekuensi evaluasi menggunakan
$\tilde\omega=\omega+i\eta$, $\eta=0.005\omega$. Akar vertikal
$\nu_\alpha=\sqrt{k^2-\tilde\omega^2/\tilde V_\alpha^2}$ dipilih dengan bagian
riil tidak negatif untuk memenuhi peluruhan/radiasi keluar.

Rekursi impedansi P-SV dihitung dari half-space menuju permukaan menggunakan
partisi perpindahan dan traksi serta propagator peluruhan
$E_j=\operatorname{diag}(e^{-\nu_{p,j}h_j},e^{-\nu_{s,j}h_j})$:

$$R_j=-(T_{u,j}-K_{j+1}U_{u,j})^{-1}(T_{d,j}-K_{j+1}U_{d,j}),$$
$$K_j=(T_{d,j}+T_{u,j}E_jR_jE_j)(U_{d,j}+U_{u,j}E_jR_jE_j)^{-1},$$
$$G_{zz}(\omega,k)=-[K_1^{-1}]_{22}.$$

Sintesis aktual mengintegrasikan compliance tersebut, bukan cabang modal yang
kecepatan dan amplitudonya ditentukan melalui aturan pendekatan:

$$U(r,\omega)=S(\omega)\left[\int_0^{k_{\max}}
(kG_{zz}(\omega,k)-A)J_0(kr)\,dk+\frac{A}{r}\right],
\qquad A=\frac{1-\tilde\nu_1}{\tilde\mu_1}.$$

Pengurangan asimtot statik mengurangi kontribusi ekor integral. Suku $A/r$
mengembalikan integral statik analitik. Kuadratur memakai segmen linier pada
bilangan gelombang sangat rendah dan segmen logaritmik sampai batas terhingga.
Jumlah titik dan batas integrasi merupakan parameter numerik yang wajib diuji
konvergensinya; implementasi ini tidak mengklaim bebas aliasing secara mutlak.

## 2. Transformasi MFJ dan pembanding

Whitening per trace diterapkan sebagai
$\hat U_m=U_m/(|U_m|+0.02\max_l|U_l|+\epsilon)$.
Dengan taper ruang $w_m$ dan bobot trapesium radial $\Delta r_m$:

$$I(f,c)=\left|\sum_m w_m\hat U_m\,
\frac{H_0^{(2)}(2\pi f r_m/c)}{|H_0^{(2)}(2\pi f r_m/c)|}
\sqrt{r_m}\,\Delta r_m\right|.$$

Pembanding plane-wave dan standing-Bessel menggunakan kernel
$e^{-ikr}$ serta $J_0(kr)$ pada gather yang sama. Tidak ada artefak buatan,
penggeseran punggung moda, atau pengaburan citra untuk menentukan pemenang.
Referensi dispersi Rayleigh dari `disba` bersifat independen dari integrator.
Belum tersedia ekstraksi pole Guided-P yang tervalidasi; overlay pendekatan
Guided-P telah dihapus.

## 3. Fusi aktif-pasif

Spektrum dinormalisasi min-max per frekuensi dan dipangkatkan 1.35. Kedua kanal
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

## 4. Likelihood komposit dengan Wasserstein yang benar

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

Mask kecepatan pasif menggunakan nilai frekuensi $f\le12$ Hz, bukan indeks baris.
Pada eksperimen Rayleigh-only, kedua spektrum dibatasi dengan mask kecepatan yang
sama. Nilai $\sigma_D$ adalah hyperparameter generalized likelihood; interval
hasil inferensi belum boleh dianggap terkalibrasi secara statistik.

## 5. Prior dan transformasi parameter

Prior seragam didefinisikan pada domain fisik yang memenuhi batas kotak dan
$V_p\ge1.45V_s$. Kecepatan $V_s$ dan ketebalan memakai transformasi sigmoid.
Untuk setiap $V_p$, batas bawah kondisional adalah
$b_p(V_s)=\max(V_{p,\min},1.45V_s)$, sehingga
$V_p=b_p+(V_{p,\max}-b_p)\sigma(z_p)$.

Jacobian bersifat triangular. Determinannya dihitung dari hasil kali diagonal,
termasuk lebar interval $V_p$ kondisional. Tidak ada clipping model setelah
sampling dan tidak ada regularisasi laten tersembunyi. Inisialisasi hanya
menggunakan Latin hypercube, tanpa mengakses parameter sebenarnya.

## 6. Inferensi campuran berkorelasi

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

## 7. MCMC dan pelaporan

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

## 8. Visualisasi dan gerbang penerbitan

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
berulang yang belum dilakukan memblokir naskah penelitian. Laporan diagnostik
boleh dibuat dengan status tersebut terlihat jelas.

## 9. Ketelitian estimasi prediksi (v6)

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
