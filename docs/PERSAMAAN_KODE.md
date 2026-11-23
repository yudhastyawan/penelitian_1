# Persamaan yang digunakan dalam kode MFJ–SMVI v6

[Indeks fungsi dan tautan kode](REFERENSI_FUNGSI.md)

Notasi $`f`$ dalam Hz, kecepatan dalam m/s, jarak dalam m, dan frekuensi sudut $`\omega=2\pi f`$ dalam rad/s. Persamaan di bawah mengikuti kode yang digunakan pada eksperimen v6. Rujukan mendasari metode dan definisi matematis. Ambang, mask, gain, dan aturan penggabungan yang ditetapkan oleh proyek dicatat sebagai pilihan implementasi.

Keterangan simbol berlaku per bagian. Tanda tilde mengikuti besaran yang dijelaskan, misalnya modulus kompleks atau target tak ternormalisasi. Huruf $`i`$ dalam eksponensial kompleks adalah satuan imajiner, sedangkan subskrip $`i,j,k,m`$ menunjukkan indeks. Tanda garis atas menyatakan rerata jika disertai indeks ensemble atau frekuensi dan konjugasi pada spektrum FFT. Tanda topi menandai hasil estimasi, koreksi, atau normalisasi sesuai konteks. Simbol $`\pi`$ adalah konstanta lingkaran. Huruf $`e`$ sebagai basis eksponensial adalah bilangan Euler.

## Model berlapis dan Jacobian

Untuk $`N`$ lapisan termasuk half-space, vektor parameter memiliki $`3N-1`$ entri. Antarmuka ke-$`j`$ berada di $`z_j=\sum_{i=1}^j h_i`$, $`1\le j<N`$. Untuk medium isotropik, $`\mu=\rho V_s^2`$ dan $`\lambda=\rho(V_p^2-2V_s^2)`$. Dengan $`r=V_p/V_s`$, rasio Poisson adalah $`\nu=(r^2/2-1)/(r^2-1)`$. Besaran profil pada kedalaman $`z`$ diambil dari lapisan yang memuat $`z`$. Titik tepat pada antarmuka masuk ke lapisan di bawahnya. Rata-rata dan kuantil dihitung pada ensemble di setiap kedalaman. [Aki dan Richards, *Quantitative Seismology*](https://www.ldeo.columbia.edu/~richards/Aki_Richards.html) · [Kode model](../models.py#L56).

Jika $`a_i`$ dan $`b_i`$ adalah batas fisik, transformasi dari ruang laten adalah

```math
s_i=\sigma(z_i)=\frac1{1+e^{-z_i}},\qquad m_i=a_i+(b_i-a_i)s_i.
```

Pada koordinat $`V_{p,i}`$, batas bawah menjadi $`a_{p,i}=\max(V_{p,i}^{\min},1.45V_{s,i})`$. Transformasi invers menggunakan $`\log u_i-\log(1-u_i)`$ dengan $`u_i=(m_i-a_i)/(b_i-a_i)`$. Karena $`V_s`$ dipetakan lebih dahulu, Jacobian triangular dan diagonalnya $`(b_i-a_i)\sigma(z_i)\sigma(-z_i)`$. Implementasi log-determinan memakai

```math
\log|\det Dm(z)|=\sum_i[\log(b_i-a_i)-\mathrm{softplus}(-z_i)-\mathrm{softplus}(z_i)].
```

Batas rasio 1,45, ordering parameter, dan densitas serta $`Q`$ tetap adalah konfigurasi proyek. Faktor Jacobian mengikuti perubahan variabel untuk prior seragam di domain fisik. [Stan Reference Manual, *constraint transforms*](https://mc-stan.org/docs/reference-manual/transforms.html) · [Transformasi v6](../models.py#L197).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`N,i,j`$ | Jumlah lapisan termasuk half-space, indeks lapisan atau koordinat, dan indeks antarmuka. | Tanpa satuan |
| $`h_i,z_j,z`$ | Tebal lapisan berhingga ke-$`i`$, kedalaman kumulatif antarmuka ke-$`j`$, dan kedalaman evaluasi profil. | m |
| $`V_s,V_p,V_{s,i},V_{p,i}`$ | Kecepatan gelombang S dan P, dengan subskrip $`i`$ untuk lapisan ke-$`i`$. | m/s |
| $`\rho,\mu,\lambda`$ | Densitas massa, modulus geser, dan parameter Lamé pertama. | kg/m³ untuk $`\rho`$, Pa untuk $`\mu,\lambda`$ |
| $`r,\nu`$ | Rasio $`V_p/V_s`$ dan rasio Poisson elastik. Simbol $`r`$ pada bagian berikutnya menyatakan jarak. | Tanpa satuan |
| $`z_i,z,m_i,m(z)`$ | Koordinat laten ke-$`i`$, vektor laten, parameter fisik ke-$`i`$, dan pemetaan vektor laten ke vektor fisik. $`z`$ di sini adalah koordinat laten saat menjadi argumen transformasi. | Laten tanpa satuan. Fisik mengikuti jenis parameter |
| $`a_i,b_i,a_{p,i},V_{p,i}^{\min}`$ | Batas bawah dan atas koordinat fisik, batas bawah kondisional $`V_{p,i}`$, dan batas bawah nominal $`V_{p,i}`$ pada prior. | Sesuai parameter. Kedua batas $`V_p`$ dalam m/s |
| $`\sigma(z_i),s_i,u_i`$ | Fungsi sigmoid, keluaran sigmoid koordinat ke-$`i`$, dan posisi relatif parameter fisik pada interval terbuka (0,1) untuk transformasi balik. | Tanpa satuan |
| $`Dm(z),\det Dm(z)`$ | Matriks turunan $`m`$ terhadap $`z`$ dan determinan yang dipakai untuk perubahan variabel prior. | Mengikuti skala koordinat fisik yang dipakai kode |
| $`Q,Q_p,Q_s`$ | Faktor kualitas gelombang. Subskrip $`p,s`$ merujuk ke gelombang P dan S. | Tanpa satuan |
| $`\log,\mathrm{softplus}(x)`$ | Logaritma natural dan fungsi $`\log(1+e^x)`$, yang dievaluasi dengan bentuk stabil dalam kode. | Tanpa satuan |

## Geometri dan sampling spasial

Untuk koordinat $`(x_i,y_i)`$, $`r_{ij}=\sqrt{(x_j-x_i)^2+(y_j-y_i)^2}`$ dan $`\theta_{ij}=\mathrm{atan2}(y_j-y_i,x_j-x_i)`$. Ada $`N(N-1)/2`$ pasangan unik. Untuk stasiun luar nested, $`\theta_k=k\pi(3-\sqrt5)+0.35`$ dan radius berjarak geometris. Faktor gangguan arah per bin adalah rerata $`\cos[2(\theta_{ij}-25^\circ)]`$. Metrik geometri mencakup persentase bin radial terisi dan $`|\mathrm{mean}\exp(2i\theta_{ij})|`$. Angka Nyquist berdasarkan jarak minimum pada array irregular adalah proksi. Untuk spasi reguler, $`k\le\pi/\Delta x`$ memberi $`c\ge2\Delta x f`$. [Kode geometri](../forward_mfj.py#L198) · [Pembandingan rancangan](../acquisition_design.py#L125) · [Foti dkk., panduan akuisisi gelombang permukaan](https://doi.org/10.1007/s10518-017-0206-7).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`(x_i,y_i),(x_j,y_j),r_{ij}`$ | Koordinat dua stasiun dan jarak horizontal di antara keduanya. | m |
| $`N,i,j,k`$ | Jumlah stasiun, indeks pasangan, dan indeks stasiun luar dalam rumus spiral. | Tanpa satuan |
| $`\theta_{ij},\theta_k`$ | Azimut dari stasiun $`i`$ ke $`j`$, dan sudut posisi stasiun luar ke-$`k`$. $`\mathrm{atan2}`$ memilih kuadran sudut berdasarkan kedua selisih koordinat. Offset sudut spiral 0,35 dinyatakan dalam radian. | rad, dengan 25° dikonversi ke rad |
| $`k,\Delta x,c,f`$ | Di rumus Nyquist, $`k`$ adalah bilangan gelombang horizontal, $`\Delta x`$ spasi receiver, $`c`$ kecepatan fase, dan $`f`$ frekuensi. | rad/m, m, m/s, Hz |

## Compliance P-SV

Untuk konvensi $`e^{-i\omega t}`$, $`\widetilde V_\alpha=V_\alpha(1-i/(2Q_\alpha))`$, $`\widetilde\mu=\rho\widetilde V_s^2`$, dan

```math
\widetilde\omega=\omega+i\eta,\qquad \nu_\alpha=\sqrt{k^2-\widetilde\omega^2/\widetilde V_\alpha^2},\qquad \gamma=2k^2-\widetilde\omega^2/\widetilde V_s^2.
```

Cabang akar memiliki bagian riil tidak negatif. Pada partisi $`U_{d/u}`$ dan $`T_{d/u}`$, tanda atas adalah gelombang turun dan tanda bawah gelombang naik.

```math
U_{d/u}=\begin{pmatrix}k&\pm\nu_s\\\pm\nu_p&k\end{pmatrix},\qquad T_{d/u}=-\widetilde\mu\begin{pmatrix}\pm2k\nu_p&\gamma\\\gamma&\pm2k\nu_s\end{pmatrix}.
```

Untuk half-space, $`K_N=T_dU_d^{-1}`$. Pada lapisan berhingga, $`E_j=\mathrm{diag}(e^{-\nu_ph_j},e^{-\nu_sh_j})`$, $`R_j=-(T_u-K_{j+1}U_u)^{-1}(T_d-K_{j+1}U_d)`$, dan $`X_j=E_jR_jE_j`$. Rekursi dan compliance permukaan adalah

```math
K_j=(T_d+T_uX_j)(U_d+U_uX_j)^{-1},\qquad G_{zz}(\omega,k)=-[K_1^{-1}]_{22}.
```

Invers 2×2 menggunakan adjugat dibagi determinan $`ad-bc`$. Jika magnitudo determinan di bawah $`10^{-30}`$, kode menambahkan $`10^{-30}(1+i)`$. Pergeseran numerik v6 adalah $`\eta=0.005\omega`$ pada sintesis standar. Model atenuasi sederhana dan nilai pergeseran adalah konfigurasi kode. Rujukan menjelaskan keluarga metode gelombang berlapis. [Kennett (1980), *Seismic waves in a stratified half space*](https://doi.org/10.1111/j.1365-246X.1980.tb04299.x) · [Rekursi Numba](../forward_mfj.py#L329) · [Rekursi PyTorch](../autodiff_likelihood.py#L20).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`t,\omega,\widetilde\omega,\eta`$ | Waktu, frekuensi sudut riil, frekuensi sudut kompleks, dan pergeseran imajiner yang dipakai dalam evaluasi numerik. | s untuk $`t`$, rad/s untuk ketiga frekuensi |
| $`\alpha,V_\alpha,\widetilde V_\alpha,Q_\alpha`$ | Jenis gelombang P atau S, kecepatan riil, kecepatan kompleks, dan faktor kualitas untuk jenis tersebut. | m/s untuk kecepatan, faktor kualitas tanpa satuan |
| $`\rho,\widetilde\mu,k`$ | Densitas, modulus geser kompleks, dan bilangan gelombang horizontal. | kg/m³, Pa, rad/m |
| $`\nu_\alpha,\gamma`$ | Akar vertikal yang mengatur peluruhan atau propagasi vertikal, dan faktor kuadratik yang muncul pada traksi di partisi P–SV. Simbol $`\nu_\alpha`$ pada bagian ini berbeda dari rasio Poisson. | m⁻¹ dan m⁻² |
| $`U_{d/u},T_{d/u},d,u`$ | Matriks partisi perpindahan dan traksi. Subskrip $`d`$ berarti gelombang turun dan $`u`$ gelombang naik. | Mengikuti basis gelombang pada rekursi |
| $`j,N,h_j`$ | Indeks lapisan, indeks half-space terakhir, dan tebal lapisan hingga ke-$`j`$. | Indeks tanpa satuan, tebal dalam m |
| $`E_j,R_j,X_j`$ | Matriks diagonal peluruhan satu lintasan, matriks refleksi di dasar lapisan, dan refleksi yang dipropagasikan kembali ke atas lapisan. | Tanpa satuan |
| $`K_j,K_{j+1},K_N,K_1`$ | Matriks impedansi yang menghubungkan perpindahan ke traksi di bagian atas lapisan $`j`$, lapisan di bawahnya, half-space, dan permukaan. | Pa/m |
| $`G_{zz}(\omega,k),[K_1^{-1}]_{22}`$ | Compliance vertikal terhadap pembebanan vertikal pada domain bilangan gelombang, dan elemen diagonal vertikal invers impedansi yang dipakai untuk menghitungnya. | m/Pa |
| $`a,b,c,d`$ | Empat elemen matriks 2×2 pada penjelasan invers. Huruf $`d`$ di sini adalah elemen kanan bawah. | Mengikuti matriks yang diinvers |

## Kuadratur Sommerfeld dan suku statik

Pada offset $`r>0`$, $`A=(1-\widetilde\nu_1)/\widetilde\mu_1`$. Keluaran Green dengan normalisasi internal proyek adalah

```math
G(\omega,r)\approx\sum_i w_i[k_iG_{zz}(\omega,k_i)-A]J_0(k_i r)+\frac Ar.
```

Suku $`A/r`$ berasal dari $`\int_0^\infty J_0(kr)\,dk=1/r`$. Untuk pembandingan perpindahan gaya titik satuan dalam [uji COMPOSTI](../diagnostics_review_v6/benchmark_composti/RINGKASAN.md), keluaran v6 dibagi $`2\pi`$. Grid menggabungkan 32 titik awal linier dari $`10^{-10}\,\mathrm m^{-1}`$ dengan segmen log yang mulai pada $`\omega/(5000\,\mathrm{m/s})`$. Segmen log memiliki $`\max(64,\lfloor1024d\rfloor)`$ titik dan ujung $`k_{\max}=\max[3\omega/\max(c_{\min}/2,40\,\mathrm{m/s}),8/r_{\min}]`$, dengan $`d`$ densitas kuadratur. Ujung bersama dihitung sekali. Untuk titik interior, $`w_i=(k_{i+1}-k_{i-1})/2`$. Pada ujung, bobot adalah setengah lebar interval tetangga. Pembantu trapesium juga digunakan pada jarak dan frekuensi. Untuk satu node pembantu mengembalikan 1. [Kennett (1980)](https://doi.org/10.1111/j.1365-246X.1980.tb04299.x) · [NIST DLMF §10.22](https://dlmf.nist.gov/10.22) · [Solusi statik Boussinesq sebagai limit homogen](https://pmc.ncbi.nlm.nih.gov/articles/PMC10302802/) · [Basis kuadratur v6](../forward_mfj.py#L803).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`r,r_{\min},\omega`$ | Offset sumber–receiver, offset terkecil dalam gather, dan frekuensi sudut pada baris spektrum. | m, m, rad/s |
| $`i,k_i,k_{\max},w_i`$ | Indeks node, bilangan gelombang pada node, batas atas integral, dan bobot trapesium pada node itu. Untuk integral terhadap $`k`$, bobot adalah lebar efektif dalam bilangan gelombang. | Indeks tanpa satuan, sisanya m⁻¹ |
| $`c_{\min},d`$ | Kecepatan fase terkecil pada grid citra dan pengali densitas kuadratur yang mengatur jumlah node logaritmik. | m/s dan tanpa satuan |
| $`J_0,G_{zz}(\omega,k_i)`$ | Fungsi Bessel jenis pertama orde nol dan compliance vertikal pada node $`i`$. Argumen $`k_i r`$ tidak bersatuan. | $`J_0`$ tanpa satuan, $`G_{zz}`$ dalam m/Pa |
| $`\widetilde\nu_1,\widetilde\mu_1,A`$ | Rasio Poisson kompleks lapisan teratas, modulus geser kompleks lapisan teratas, dan koefisien asimtot statik. Di sini $`\widetilde\nu_1`$ adalah rasio Poisson. | Tanpa satuan, Pa, Pa⁻¹ |
| $`G(\omega,r),A/r`$ | Fungsi Green vertikal pada offset dengan normalisasi internal proyek dan suku statik yang ditambahkan kembali setelah integrasi. | m/N, dengan skala internal $`2\pi`$ relatif terhadap respons gaya titik satuan |

## Sumber, gangguan dan FFT

Gather bersih adalah $`U_j(f,r)=S_j(f)G(f,r)`$. Sumber aktif memakai bentuk spektrum Ricker $`S_a=x e^{-x}`$ dengan $`x=(f/f_0)^2`$ dan $`f_0`$ default 25 Hz. Kanal pasif memakai $`S_p=[1+(f/(6\,\mathrm{Hz}))^{1.5}]^{-1}`$. Dalam observasi aktif, gain per offset adalah $`\exp[\mathcal N(0,0.35^2)]`$ dan pergeseran fase $`\mathcal N(0,0.08^2)`$. Pada pasif, fase tambahan adalah $`0.045a_r\omega r/(260\,\mathrm{m/s})`$ dengan $`a_r`$ faktor arah dari pasangan stasiun. Noise aditif per sel adalah $`b_j\,\mathrm{RMS}_r(U_j)(X+iY)/\sqrt2`$ setelah perturbasi multiplikatif, dengan $`X,Y`$ normal standar independen, $`b_a=0.08`$ dan $`b_p=0.045/\sqrt{\max(N_r,1)}`$. Ini model gangguan sintetik proyek. Kanal pasif merupakan proksi fungsi Green empiris. [Kode sintesis](../forward_mfj.py#L854) · [Rujukan Ricker pada sintesis aktif, Mreyen dkk.](https://doi.org/10.1093/gji/ggaf498) · [Wapenaar, hubungan korelasi dan fungsi Green](https://doi.org/10.1103/PhysRevLett.93.254301).

Untuk ilustrasi jejak, Green dihitung pada grid FFT berjarak $`\Delta f=1/(N_t\Delta t)`$ dengan $`N_t\ge2048`$. Bin dalam pita diisi dengan $`\overline U(f,r)e^{-i\omega t_0}`$, dan bin di luar pita diisi nol. Hasil inverse real FFT memberi $`u(t,r)`$, dengan $`t_0=0.045`$ s aktif dan $`0.12`$ s pasif. Opsi noise pada jejak ilustratif memberi skala RMS global 0,04, sehingga realisasi ini berbeda dari observasi inversi. [Konvensi NumPy FFT](https://numpy.org/doc/stable/reference/routines.fft.html#implementation-details) · [Kode jejak](../forward_mfj.py#L898).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`j,a,p`$ | Indeks kanal, dengan $`a`$ untuk aktif dan $`p`$ untuk pasif. | Tanpa satuan |
| $`U_j(f,r),S_j(f),G(f,r)`$ | Gather kompleks per kanal, bobot spektrum sumber, dan respons Green. Penulisan $`G(f,r)`$ mengevaluasi $`G(\omega,r)`$ pada $`\omega=2\pi f`$. | $`S_j`$ tanpa satuan. $`U_j,G`$ pada skala sintetik internal |
| $`S_a,S_p,x,f_0`$ | Bentuk sumber aktif, bentuk sumber pasif, argumen Ricker kuadrat, dan frekuensi pusat aktif. Angka 6 Hz adalah skala bentuk sumber pasif. | $`f_0`$ dalam Hz, sisanya tanpa satuan |
| $`\mathcal N(0,\sigma^2),X,Y`$ | Distribusi normal dengan mean nol dan ragam $`\sigma^2`$, serta dua pengambilan normal standar independen untuk noise kompleks. | Tanpa satuan untuk gain dan draw, fase diukur dalam rad |
| $`a_r,N_r,b_j,b_a,b_p`$ | Faktor azimut pada lag $`r`$, jumlah pasangan dalam lag, dan pengali noise aditif untuk kanal $`j`$, aktif, serta pasif. | Tanpa satuan |
| $`\mathrm{RMS}_r(U_j)`$ | Akar rata-rata kuadrat magnitudo gather sepanjang offset pada frekuensi tetap. | Sama dengan $`U_j`$ |
| $`N_t,\Delta t,\Delta f`$ | Jumlah sampel waktu, langkah waktu, dan jarak antarbin frekuensi FFT. | Hitungan, s, Hz |
| $`\overline U(f,r),t,t_0,u(t,r)`$ | Konjugat kompleks gather, waktu evaluasi, tunda sumber, dan jejak riil hasil inverse FFT. | $`t,t_0`$ dalam s. Amplitudo jejak mengikuti skala sintetik |

## Pencitraan MFJ

Argumen fungsi Bessel adalah $`x=\max(2\pi fr/c,10^{-4})`$. Identitas $`H_0^{(1,2)}=J_0\pm iY_0`$ memberi kernel unit-modulus $`H_0^{(1,2)}/(|H_0^{(1,2)}|+10^{-30})`$. Proyeksi menggunakan cabang $`H_0^{(2)}`$ untuk konvensi gelombang maju proyek. Pada pembantu Bessel, $`x`$ dibatasi minimal $`10^{-8}`$. Untuk $`x\le12`$, kode memakai suku deret $`t_0=1`$, $`t_j=-t_{j-1}x^2/(4j^2)`$ dan bilangan harmonik $`H_j=\sum_{m=1}^j1/m`$.

```math
J_0(x)\approx\sum_{j=0}^n t_j,\qquad
Y_0(x)\approx\frac{2}{\pi}\left[(\log(x/2)+\gamma_E)J_0(x)-\sum_{j=1}^n t_jH_j\right].
```

Penjumlahan berhenti paling lambat pada $`n=33`$ atau setelah suku di bawah $`10^{-16}`$. Untuk $`x>12`$, $`\chi=x-\pi/4`$ dan pasangan aproksimasi adalah $`\sqrt{2/(\pi x)}(P\cos\chi-Q\sin\chi,\;P\sin\chi+Q\cos\chi)`$ dengan polinomial terpotong $`P,Q`$ dalam $`1/x`$ yang tercantum di kode. [NIST DLMF §10.4](https://dlmf.nist.gov/10.4.E3) · [Deret §10.8](https://dlmf.nist.gov/10.8) · [Asimtotik §10.17](https://dlmf.nist.gov/10.17) · [Kernel Bessel v6](../forward_mfj.py#L33).

```math
\widehat U_m=\frac{U_m}{|U_m|+0.02\max_l |U_l|+10^{-30}},\qquad
I(f,c)=\left|\sum_m \widehat U_m\,T(r_m)\sqrt{r_m}\,\Delta r_m\,\frac{H_0^{(2)}(x_m)}{|H_0^{(2)}(x_m)|+10^{-30}}\right|.
```

Hasil dibagi maksimum per frekuensi jika puncaknya cukup besar. Dengan koordinat radial ternormalisasi $`u`$ dan fraksi tepi $`e`$, taper kiri adalah $`(1-\cos(\pi u/e))/2`$. Taper kanan menggunakan $`1-u`$, dan bagian tengah bernilai 1. Pembanding mengganti kernel dengan $`e^{-ix}`$ atau $`J_0(x)`$ dengan pra-pembobotan yang sama. Normalisasi kernel, taper, whitening, dan $`\sqrt r`$ adalah varian proyek. [Xi dkk. (2021), MFJ](https://doi.org/10.1093/gji/ggab008) · [Implementasi proyeksi v6](../forward_mfj.py#L616).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`f,r,c,x,x_m`$ | Frekuensi, offset, kecepatan fase yang dipindai, argumen Bessel, dan nilai argumen itu di offset ke-$`m`$. | Hz, m, m/s. Argumen tanpa satuan |
| $`J_0,Y_0,H_0^{(1)},H_0^{(2)}`$ | Fungsi Bessel jenis pertama dan kedua orde nol, serta fungsi Hankel jenis pertama dan kedua orde nol. | Tanpa satuan |
| $`t_j,H_j,\gamma_E,n`$ | Suku deret Bessel ke-$`j`$, bilangan harmonik, konstanta Euler–Mascheroni, dan indeks suku terakhir dalam aproksimasi. $`H_j`$ di deret berbeda dari fungsi Hankel. | Tanpa satuan |
| $`\chi,P,Q`$ | Fase asimtotik dan polinomial pendekatan yang mengalikan kosinus serta sinus. | rad untuk $`\chi`$, polinomial tanpa satuan |
| $`U_m,U_l,\widehat U_m`$ | Amplitudo kompleks pada jejak $`m`$ dan $`l`$ dalam baris frekuensi yang sama, serta amplitudo $`m`$ setelah whitening. Indeks $`l`$ dipakai untuk mencari puncak. | $`U_m,U_l`$ mengikuti gather, $`\widehat U_m`$ tanpa satuan |
| $`T(r_m),\Delta r_m,u,e`$ | Bobot taper, bobot jarak radial trapesium, koordinat radial ternormalisasi, dan fraksi tepi taper. | $`\Delta r_m`$ dalam m, sisanya tanpa satuan |
| $`I(f,c)`$ | Besar proyeksi per frekuensi dan kecepatan fase. Citra yang dikembalikan kode dinormalisasi terhadap puncak tiap frekuensi. | Citra keluaran tanpa satuan |

## Fusi frekuensi native

Normalisasi kontras memakai $`D'=[(D-D_{\min})/(D_{\max}-D_{\min})]^p`$ per irisan yang bervariasi. Irisan konstan menjadi nol. Benchmark A/B memakai $`p=1.35`$, dan komparasi rancangan memakai default $`p=1.4`$. Bagian dalam overlap memakai irisan native yang tersedia di kedua kanal. Dengan $`s=(f-f_{\rm lo})/(f_{\rm hi}-f_{\rm lo})`$ pada overlap, $`\alpha=\sin^2(\pi s/2)`$, sementara $`\alpha=0`$ di kiri dan 1 di kanan. Bentuk spektrum diukur dengan

```math
v_j=\max\left[\left(\frac{\mathrm{mean}_c D'_j}{\max_c D'_j+10^{-12}}\right)^2,10^{-4}\right],\quad
P_a=\frac{a\,\alpha}{v_a},\quad P_p=\frac{p_{\rm avail}(1-\alpha)}{v_p}.
```

$`a,p_{\rm avail}`$ adalah indikator ketersediaan. Saat hanya satu kanal tersedia, $`\alpha`$ disesuaikan untuk kanal itu. Bobot $`w_a=P_a/(P_a+P_p+10^{-30})`$ dan $`w_p=P_p/(P_a+P_p+10^{-30})`$. Irisan hasil $`w_aD'_a+w_pD'_p`$ kembali dinormalisasi min–max. Proksi $`v_j`$ menggunakan **rata-rata**, sesuai kode, dan tidak dikalibrasi sebagai ragam noise. [Aturan fusi proyek](../fusion.py#L91) · [Fusi dan grid native](../fusion.py#L116).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`D,D_{\min},D_{\max},D',D'_j`$ | Spektrum yang akan dinormalisasi, nilai minimum dan maksimum irisan, serta spektrum ternormalisasi umum dan untuk kanal $`j`$. | Skala citra ternormalisasi |
| $`p,j,a,p_{\rm avail}`$ | Eksponen kontras, indeks kanal, indikator ketersediaan aktif, dan indikator ketersediaan pasif. Indikator bernilai 0 atau 1. Subskrip $`p`$ pada $`v_p,P_p,w_p`$ berarti kanal pasif. | Tanpa satuan |
| $`f,f_{\rm lo},f_{\rm hi},s,\alpha`$ | Frekuensi, batas bawah dan atas overlap (default 8 dan 12 Hz), posisi relatif dalam overlap, serta bobot transisi awal. | Frekuensi dalam Hz. $`s,\alpha`$ tanpa satuan |
| $`c,\mathrm{mean}_c,\max_c`$ | Kecepatan fase dan operasi rata-rata serta maksimum sepanjang grid kecepatan pada satu frekuensi. | $`c`$ dalam m/s |
| $`v_j,v_a,v_p`$ | Proksi bentuk spektral kanal $`j`$, aktif, atau pasif, yang dihitung dari rasio rata-rata ke puncak. | Tanpa satuan |
| $`P_a,P_p,w_a,w_p`$ | Bobot mentah aktif dan pasif, serta bobot fusi setelah normalisasi. | Tanpa satuan |

## Transport 1D dan target inversi

Untuk irisan $`o,s`$ pada grid kecepatan $`c`$, $`\mathrm{NCC}=o^\mathsf{T}s/\max(\lVert o\rVert_2\lVert s\rVert_2,10^{-30})`$. Massa diskrit dibentuk dari $`\max(D,0)^\beta+10^{-12}`$, dengan $`\beta=2`$ pada likelihood. Setelah normalisasi massa, transport monoton memberikan

```math
W_2^2=\int_0^1|F_o^{-1}(u)-F_s^{-1}(u)|^2\,du=\sum_{ij}\gamma_{ij}(c_i-c_j)^2.
```

Numba mengalirkan sisa massa antar-support dengan dua penunjuk. PyTorch membangun gabungan breakpoint CDF dan mengintegrasikan biaya pada interval kuantil. Keduanya menggunakan posisi kecepatan aktual. [Peyré dan Cuturi (2019), §2.6](https://arxiv.org/abs/1803.00567) · [Numba](../smvi_inversion.py#L30) · [PyTorch](../autodiff_likelihood.py#L139).

```math
J=(1-\lambda)\overline{1-\mathrm{NCC}}+\lambda\,\overline{W_2^2}/s_W^2,\qquad
\log\widetilde\pi(z)=-J/(2\sigma_D^2)+\log|\det Dm(z)|.
```

Garis atas adalah rata-rata berbobot trapesium pada tiap pita kontigu yang dipertahankan. Mask $`c`$ untuk $`f\le12`$ Hz memakai batas default 680 m/s, dengan minimal dua bin bila tersedia. $`s_W`$ tetap diambil dari konfigurasi atau 0,15 kali span grid parent. Keluaran metrik transport adalah $`\sqrt{\overline{W_2^2}}`$ dalam m/s. Bentuk $`J`$, mask, dan skala adalah pilihan proyek. Target laten mengikuti pembaruan berbasis loss dan perubahan variabel prior. Gradien dihitung melalui autodiff atau beda hingga sentral $`(\ell(z+\delta e_i)-\ell(z-\delta e_i))/(2\delta)`$. [Bissiri, Holmes, dan Walker (2016)](https://doi.org/10.1111/rssb.12158) · [Target v6](../smvi_inversion.py#L262) · [Dokumentasi PyTorch autograd](https://docs.pytorch.org/docs/stable/generated/torch.autograd.grad.html).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`o,s,c_i,c_j`$ | Vektor intensitas observasi dan sintetik pada satu frekuensi, serta kecepatan fase pada bin $`i`$ dan $`j`$. | Intensitas ternormalisasi, kecepatan dalam m/s |
| $`\mathrm{NCC},o^\mathsf{T}s,\lVert\cdot\rVert_2`$ | Kesamaan kosinus irisan, hasil kali titik, dan norma Euclidean. Superskrip $`\mathsf T`$ adalah transpose. | NCC tanpa satuan. Norma mengikuti masukan |
| $`D,\beta`$ | Intensitas citra per bin dan pangkat penajaman untuk konversi ke massa sebelum dinormalisasi. | Tanpa satuan |
| $`F_o^{-1},F_s^{-1},u,\gamma_{ij}`$ | Fungsi kuantil distribusi massa observasi dan sintetik, level kuantil antara 0 dan 1, serta massa yang ditranspor dari bin $`i`$ ke bin $`j`$. | Kuantil dalam m/s. $`u,\gamma_{ij}`$ tanpa satuan |
| $`W_2,W_2^2,s_W`$ | Jarak Wasserstein orde dua, biaya kuadratik, dan skala kecepatan tetap untuk menormalisasi biaya. | m/s, (m/s)², m/s |
| $`J,\lambda,\overline{(\cdot)},\sigma_D`$ | Misfit spektral gabungan, bobot komponen Wasserstein, rerata berbobot terhadap frekuensi, dan skala yang mengatur pengaruh misfit pada log target. | $`J,\lambda,\sigma_D`$ tanpa satuan. Rerata mengikuti satuan argumennya |
| $`\widetilde\pi(z),Dm(z),\ell(z)`$ | Target laten tak ternormalisasi, Jacobian transformasi laten ke fisik, dan log target yang dibedakan. | Log target dihitung pada koordinat laten tanpa satuan |
| $`\delta,e_i`$ | Langkah beda hingga dan vektor satuan koordinat laten ke-$`i`$. | Tanpa satuan |

## Campuran Gaussian dan optimasi Stein

```math
q(z)=\sum_{k=1}^K w_k\mathcal N(z;\mu_k,L_k L_k^\mathsf{T}),\qquad z_k=\mu_k+L_k\epsilon,\quad\epsilon\sim\mathcal N(0,I),\quad w_k=1/K.
```

$`L_k`$ dibangun dari entri bawah dan diagonal positif $`\exp(\ell_k)`$. Untuk $`y_k=L_k^{-1}(z-\mu_k)`$, log kontribusi komponen adalah $`\ell_k^\ast=\log w_k-\frac d2\log(2\pi)-\sum_i\log L_{k,ii}-y_k^\mathsf{T}y_k/2`$. Tanggung jawab $`r_k=\exp(\ell_k^\ast-\log\sum_j\exp\ell_j^\ast)`$ memberi skor $`\nabla_z\log q=\sum_k r_k(-L_k^{-\mathsf T}y_k)`$.

ELBO adalah $`\mathbb E_q[\log\widetilde\pi(z)-\log q(z)]`$ hingga konstanta. Dengan $`v=\nabla_z(\log\widetilde\pi-\log q)`$, kontribusi gradien pathwise adalah $`v`$ untuk mean dan $`v\epsilon^\mathsf{T}`$ untuk Cholesky. Diagonal log memakai aturan rantai $`d(\exp\ell)/d\ell=\exp\ell`$. Vektor Stein $`\theta_k`$ berisi mean, log-diagonal, dan segitiga bawah ketat. Dengan $`k(\theta_j,\theta_i)=\exp(-\lVert \theta_j-\theta_i\rVert^2/h)`$, $`h=\max[\mathrm{median}_{j<i}\lVert \theta_j-\theta_i\rVert^2/\log(K+1),10^{-8}]`$, pembaruannya

```math
\Phi_i=K^{-1}\sum_j [k(\theta_j,\theta_i)\,g_j+\alpha\nabla_{\theta_j}k(\theta_j,\theta_i)].
```

$`g_j`$ adalah gradien per komponen dan $`\alpha`$ mengatur repulsi. Untuk $`K=1`$ transport langsung memakai gradien. Pembaruan Adam untuk ascent memakai $`m_t=0.9m_{t-1}+0.1g_t`$, $`v_t=0.999v_{t-1}+0.001g_t^2`$, dan $`\theta\gets\theta+\gamma\widehat m_t/(\sqrt{\widehat v_t}+10^{-8})`$, dengan koreksi bias standar. Implementasi membatasi gradien dan log-diagonal. Inisialisasi menempatkan sampel Latin hypercube $`u_{ij}=(\pi_j(i)+\epsilon_{ij})/n`$ dalam strata, memakai logit, menyeleksi pusat berdasarkan target, lalu dapat merapikannya dengan L-BFGS-B. [Zhao dan Curtis (2026), pracetak §2.4](https://arxiv.org/html/2610.23603v1) · [Kingma dan Ba, Adam](https://arxiv.org/abs/1412.6980) · [McKay dkk., Latin hypercube](https://doi.org/10.1080/00401706.1979.10489755) · [Byrd dkk., L-BFGS-B](https://doi.org/10.1137/0916069) · [Implementasi SMVI](../smvi_inversion.py#L338).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`z,z_k,K,k,d`$ | Vektor laten, draw dari komponen $`k`$, jumlah komponen, indeks komponen, dan dimensi laten. | Tanpa satuan |
| $`q(z),w_k,\mu_k,L_k,L_kL_k^\mathsf{T}`$ | Densitas campuran, bobot komponen, mean laten, faktor Cholesky, dan matriks kovariansi laten. | Mengikuti koordinat laten tanpa satuan |
| $`\mathcal N,\epsilon,I`$ | Distribusi Gaussian yang diparameterkan dengan mean dan kovariansi, draw normal standar multivariat, dan matriks identitas. | Tanpa satuan |
| $`\ell_k,L_{k,ii},y_k,\ell_k^\ast,r_k`$ | Vektor log-diagonal, entri diagonal ke-$`i`$ dari Cholesky, residual laten yang dipetakan dengan invers Cholesky, log kontribusi berbobot, dan tanggung jawab komponen pada $`z`$. | Tanpa satuan |
| $`\mathbb E_q,\widetilde\pi(z),v`$ | Ekspektasi terhadap campuran, target laten, dan selisih skor target terhadap skor campuran. | Mengikuti koordinat laten tanpa satuan |
| $`\theta_k,k(\theta_j,\theta_i),h`$ | Vektor parameter bebas komponen, nilai kernel RBF antara dua komponen, dan bandwidth kernel. Huruf $`k`$ pada nama fungsi adalah kernel. | Tanpa satuan |
| $`g_j,\Phi_i,\alpha`$ | Gradien estimasi ELBO komponen $`j`$, arah pembaruan Stein untuk komponen $`i`$, dan pengali repulsi. | Mengikuti parameter bebas dalam ruang laten |
| $`t,g_t,m_t,v_t,\widehat m_t,\widehat v_t,\gamma`$ | Iterasi Adam, gradien masuk ke Adam, dua momen bergerak, kedua momen setelah koreksi bias, dan laju pembaruan. $`t`$ di sini adalah nomor iterasi. | Tanpa satuan dalam koordinat laten |
| $`u_{ij},\pi_j(i),\epsilon_{ij},n`$ | Sampel ternormalisasi untuk kandidat $`i`$ pada koordinat $`j`$, indeks strata hasil permutasi, jitter dalam strata, dan jumlah kandidat awal. Jitter Latin hypercube berbeda dari $`\epsilon`$ Gaussian. | Tanpa satuan |

## Ensemble MCMC

Stretch move untuk walker $`x`$ dengan pasangan $`z_j`$ memakai $`y=z_j+\zeta(x-z_j)`$ dan $`\zeta=[1+(a-1)U]^2/a`$, $`U`$ seragam pada $`[0,1]`$. Log rasio penerimaan adalah $`(d-1)\log\zeta+\log\widetilde\pi(y)-\log\widetilde\pi(x)`$. Nilai tidak hingga ditolak. Setelah burn-in, perkiraan $`\tau=1+2\sum_{\ell>0}\rho_\ell`$ memakai ACF FFT, dengan penghentian pada autokorelasi negatif atau jendela $`5\tau`$. Kode melaporkan $`\mathrm{ESS}=n_{\rm draw}/\tau`$ dan memeriksa panjang per walker terhadap $`50\max\tau`$. Pemeriksaan panjang itu adalah diagnostik awal. [Goodman dan Weare (2010)](https://msp.org/camcos/2010/5-1/p04.xhtml) · [Implementasi](../smvi_inversion.py#L627).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`x,z_j,y`$ | Posisi walker yang sedang diperbarui, posisi walker pasangan berindeks $`j`$, dan posisi proposal baru di ruang laten. | Tanpa satuan |
| $`\zeta,a,U,d`$ | Faktor peregangan terambil, skala peregangan dengan $`a>1`$, draw seragam, dan dimensi vektor laten. | Tanpa satuan |
| $`\widetilde\pi(x),\widetilde\pi(y)`$ | Nilai target tak ternormalisasi pada posisi saat ini dan pada proposal. | Mengikuti koordinat laten tanpa satuan |
| $`\ell,\rho_\ell,\tau`$ | Jeda linimasa chain, koefisien autokorelasi pada jeda itu, dan estimasi waktu autokorelasi terintegrasi. | $`\ell,\tau`$ dalam langkah MCMC, $`\rho_\ell`$ tanpa satuan |
| $`n_{\rm draw},\mathrm{ESS}`$ | Jumlah draw gabungan setelah burn-in dan ukuran sampel efektif yang diperkirakan dari autokorelasi. | Hitungan sampel |

## Prediksi dan presisi Monte Carlo

Untuk $`n`$ draw seimbang pada setiap komponen aktif, $`\overline D_k=n^{-1}\sum_i D_{ki}`$. Estimator prediksi, ragam Monte Carlo atas mean, dan ragam respons campuran adalah

```math
\widehat D=\sum_k w_k\overline D_k,\qquad
\widehat{\mathrm{Var}}_{\rm MC}(\widehat D)=\sum_k w_k^2\widehat{\mathrm{Var}}(D_k)/n,\qquad
V_{\rm respons}=\left[\sum_k w_k\,\overline{D_k^2}-\widehat D^2\right]_+.
```

Bootstrap mengambil ulang draw dalam masing-masing komponen. Simpangan baku $`J(\widehat D)`$ antar-bootstrap memberi galat Monte Carlo misfit. Kode mengecek RMS galat mean citra, galat misfit, dan selisih misfit dengan setengah draw. Jumlah draw bertambah hingga syarat presisi terpenuhi atau anggaran habis. Besaran ini kondisional pada campuran yang sudah difit. [Efron (1979), bootstrap](https://doi.org/10.1214/aos/1176344552) · [Implementasi](../posterior_prediction.py#L11).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`k,i,n,w_k`$ | Indeks komponen, indeks draw, jumlah draw per komponen berbobot positif, dan bobot komponen yang dinormalisasi. | Tanpa satuan |
| $`D_{ki},\overline D_k,\widehat D`$ | Citra untuk draw $`i`$ dari komponen $`k`$, mean citra komponen, dan prediksi mean campuran berbobot. | Intensitas citra ternormalisasi |
| $`\widehat{\mathrm{Var}}(D_k),\widehat{\mathrm{Var}}_{\rm MC}(\widehat D)`$ | Ragam sampel respons dalam komponen dan ragam estimator mean prediksi akibat jumlah draw terbatas. Akar ragam MC memberi galat baku per sel. | Kuadrat satuan citra |
| $`\overline{D_k^2},V_{\rm respons},[x]_+`$ | Mean kuadrat respons dalam komponen, ragam respons campuran, dan operasi mengambil nilai maksimum antara $`x`$ dan nol. | Kuadrat satuan citra pada persamaan |
| $`J(\widehat D),\mathrm{RMS}`$ | Misfit yang dihitung pada citra prediksi mean, dan agregasi akar rata-rata kuadrat untuk galat baku citra. | Misfit tanpa satuan. RMS mengikuti masukan |

## Metrik diagnostik

- **Profil kedalaman.** Dengan $`e(z)=(\overline V_s(z)-V_s^\star(z))/V_s^\star(z)`$, $`\mathrm{RMSE}_{\rm rel}(\%)=100\sqrt{\mathrm{mean}_z e(z)^2}`$ dan $`\mathrm{MAPE}(\%)=100\mathrm{mean}_z|e(z)|`$. [Kode metrik profil](../ablation_experiments.py#L23).
- **Citra dan jejak.** Galat relatif L2 adalah $`\lVert x-x_{\rm ref}\rVert_2/\lVert x_{\rm ref}\rVert_2`$. Untuk jarak citra, MAE memakai $`\mathrm{mean}|D-D_{\rm ref}|`$. NCC pada uji geometri memakai cosine similarity atas seluruh jendela overlap terhadap citra bebas noise dari kombinasi yang sama. [Komparasi geometri](../acquisition_design.py#L225) · [Komparasi COMPOSTI](../diagnostics_review_v6/benchmark_composti/compare.py#L89).
- **Noise terkendali.** $`N'=N\lVert S\rVert_2/(\lVert N\rVert_2\,10^{\mathrm{SNR}_{\rm dB}/20})`$. Stress-test sumber lemah mengalikan gather aktif dengan $`\sigma[(f-15)/2]\exp(-r/26)`$ untuk $`f`$ dalam Hz dan $`r`$ dalam m. Gain kopling lognormal mempunyai log-mean $`-\sigma_A^2/2`$ agar rata-ratanya satu. Skala latar noise dibandingkan dengan sumber awal pada skenario lemah. [Kode stress-test](../ablation_experiments.py#L266).
- **Jejak dan visualisasi.** Seleksi jejak memakai node terdekat dari target radial seragam dan lebar setengah wiggle $`0{,}35\min\Delta r_{\rm terpilih}`$. Cadangan KDE memakai $`\widehat f(x)=\sum_i\exp(-\frac12\sum_j[(x_j-x_{ij})/h_j]^2)/(n(2\pi)^{d/2}\prod_j h_j)`$ dengan bandwidth per koordinat $`h_j=(\mathrm{sd}(x_{\cdot j})+10^{-12})\,b`$. SciPy digunakan ketika tersedia. [Kode visualisasi](../visualize.py#L44) · [Dokumentasi SciPy KDE](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html).
- **Pembanding eksternal.** PyProp8 menggunakan propagator berlapis. Uji compliance memasangkan moduli kompleks dan titik bilangan gelombang. Uji jejak pita aktif memakai kuadratur tersendiri. COMPOSTI memakai solver Levin eksternal. Konjugasi, normalisasi gaya titik, dan kontrol $`\eta`$ serta taper dicatat pada laporan komparasi. [Valentine dan Sambridge, PyProp8](https://doi.org/10.21105/joss.04217) · [Niskanen dan Lähivaara, COMPOSTI](https://doi.org/10.1016/j.softx.2022.101298) · [Laporan COMPOSTI](../diagnostics_review_v6/benchmark_composti/RINGKASAN.md). Kurva Rayleigh untuk overlay dihitung oleh [disba](https://github.com/keurfonluu/disba), dengan $`T=1/f`$ dan konversi satuan; hasil Guided-P kosong pada v6.
- **Identitas analitik untuk tes.** Pada half-space homogen, fungsi Rayleigh adalah $`R=\gamma^2-4k^2\nu_p\nu_s`$. Tes memakai $`G_{zz}=-\nu_p(\widetilde\omega^2/\widetilde V_s^2)/(\widetilde\mu R)`$. Untuk kasus homogen dengan $`V_p=2V_s`$, tes tiba Rayleigh menentukan $`x=(c_R/V_s)^2`$ dari akar fisik $`x^3-8x^2+20x-12=0`$ dan memeriksa $`t_0+r/c_R`$. [Dasar gelombang Rayleigh, Aki dan Richards](https://www.ldeo.columbia.edu/~richards/Aki_Richards.html) · [Tes compliance](../tests/test_scientific_regressions.py#L82) · [Tes waktu tiba](../tests/test_reassessment.py#L190).

**Dekomposisi residual dan ablasi picking.** Untuk dua irisan nonnol, identitas yang digunakan adalah

```math
1-\mathrm{NCC}=\frac{1}{2}\left\lVert \frac{o}{\lVert o\rVert_2}-\frac{s}{\lVert s\rVert_2}\right\rVert_2^2.
```

Skor residual B menjumlahkan suku kuadrat ini pada rentang kecepatan terpisah. Pembandingan noise baru memakai $`\Delta J=J_{\rm pred}-J_{\rm prior}`$ pada observasi yang sama. Ablasi picking minimisasi residual $`(c_{\rm Rayleigh}-c_{\rm pick})/(20\,\mathrm{m/s})`$, dengan substitusi penalti jika kurva tidak dapat dihitung. [Kode dekomposisi](../diagnostics_review_v6/benchmark_b_residual/analyze.py#L63) · [Evaluasi noise baru](../diagnostics_review_v6/benchmark_b_restarts/check_new_noise.py) · [Residual picking](../ablation_experiments.py#L158).

**Keterangan simbol**

| Simbol | Arti dalam bagian ini | Satuan atau skala |
| --- | --- | --- |
| $`z,e(z),\overline V_s(z),V_s^\star(z)`$ | Kedalaman, galat relatif pada kedalaman tersebut, profil Vs rerata sampel, dan profil Vs model acuan. | m untuk $`z`$, m/s untuk kedua profil. $`e(z)`$ tanpa satuan |
| $`\mathrm{RMSE}_{\rm rel},\mathrm{MAPE},\mathrm{MAE}`$ | Akar rerata kuadrat galat relatif, rerata absolut persen galat, dan rerata absolut selisih citra. | Dua metrik profil dalam persen. MAE mengikuti skala citra |
| $`x,x_{\rm ref},D,D_{\rm ref}`$ | Array keluaran dan referensinya pada galat L2, serta citra dan citra referensi pada MAE. | Sama dalam tiap pasangan. Galat L2 relatif tanpa satuan |
| $`S,N,N',\mathrm{SNR}_{\rm dB}`$ | Gather bersih, noise kompleks sebelum penskalaan, noise setelah penskalaan, dan target rasio sinyal terhadap noise. Huruf $`N`$ di sini adalah array noise. | $`N'`$ mengikuti skala amplitudo $`S`$. SNR dalam dB |
| $`\sigma,\sigma_A,f,r,\Delta r_{\rm terpilih}`$ | Fungsi sigmoid, simpangan baku log-gain kopling, frekuensi, offset, serta jarak antara dua jejak terpilih bertetangga. Konstanta pelemahan menggunakan 15 dan 2 Hz serta 26 m. | $`f`$ dalam Hz. $`r,\Delta r_{\rm terpilih}`$ dalam m. $`\sigma,\sigma_A`$ tanpa satuan |
| $`\widehat f(x),x_j,x_{ij},x_{\cdot j},i,j`$ | Densitas KDE pada titik $`x`$, koordinat evaluasi ke-$`j`$, koordinat ke-$`j`$ dari sampel $`i`$, dan seluruh sampel pada koordinat $`j`$. Huruf $`f`$ pada KDE berarti densitas. | Koordinat mengikuti parameter. Densitas per volume ruang parameter |
| $`n,d,h_j,b,\mathrm{sd}`$ | Jumlah sampel KDE, dimensi, bandwidth koordinat ke-$`j`$, pengali bandwidth, dan simpangan baku. | $`h_j`$ dan simpangan baku mengikuti satuan koordinat |
| $`T,c_R,R,x`$ | Periode gelombang, kecepatan fase Rayleigh, fungsi Rayleigh homogen, dan rasio kuadrat $`c_R/V_s`$ yang menjadi variabel polinomial. | s, m/s, m⁻⁴, tanpa satuan |
| $`o,s,\lVert o\rVert_2,\lVert s\rVert_2`$ | Irisan observasi dan sintetik nonnol, beserta norma Euclidean masing-masing yang dipakai untuk menormalisasi arah vektornya. | Skala citra ternormalisasi |
| $`\Delta J,J_{\rm pred},J_{\rm prior}`$ | Selisih misfit, misfit prediksi yang disimpan, dan misfit model tengah prior pada realisasi observasi yang sama. | Tanpa satuan |
| $`c_{\rm Rayleigh},c_{\rm pick}`$ | Kecepatan fase kurva Rayleigh pembanding dan kecepatan fase yang diperoleh dengan picking. Residualnya dibagi 20 m/s. | m/s, residual tanpa satuan |
