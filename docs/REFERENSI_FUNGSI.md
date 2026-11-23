# Referensi fungsi kode MFJ–SMVI v6

Gunakan nama fungsi untuk membuka definisi kode. Kolom persamaan mengarah ke rumus dan rujukan yang terkait dengan fungsi tersebut. Tanda pisah berarti fungsi menjalankan orkestrasi, validasi, format keluaran, atau tugas bantu tanpa persamaan ilmiah tersendiri.

Panduan mencakup fungsi, method, properti, fungsi lokal, dan skrip tingkat modul pada kode aktif serta tes. Snapshot di `results_v6/source_snapshot/` menyimpan salinan kode eksperimen yang ditelusuri lewat manifest. [Kumpulan persamaan dan rujukan](PERSAMAAN_KODE.md)

- [Alur utama dan model](#alur-utama-dan-model)
- [Sintesis, geometri, dan fusi](#sintesis-geometri-dan-fusi)
- [Target, inferensi, dan prediksi](#target-inferensi-dan-prediksi)
- [Eksperimen lanjutan dan pelaporan](#eksperimen-lanjutan-dan-pelaporan)
- [Pemeriksaan numerik dan diagnostik](#pemeriksaan-numerik-dan-diagnostik)
- [Tes](#tes)


## Alur utama dan model

### [run_benchmark.py](../run_benchmark.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`code_hashes`](../run_benchmark.py#L48) | Menginventarisasi kode Python, tes dan konfigurasi dependensi untuk menghitung SHA-256 tiap berkas. | — |
| [`model_summary`](../run_benchmark.py#L61) | Merangkum parameter acuan, rata-rata, kuantil, cakupan per parameter, metrik profil dan biaya dari hasil inversi. | [Profil](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`run_single_benchmark`](../run_benchmark.py#L98) | Dari model, prior, grid, konfigurasi dan seed, membuat observasi, menjalankan SMVI, membandingkan prior serta acuan, dan menyimpan NPZ, ringkasan serta gambar. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`quadrature_check`](../run_benchmark.py#L221) | Menjalankan ulang gather bersih pada densitas dasar dan ganda untuk beberapa frekuensi dan offset, lalu mengembalikan selisih L2 dan relatif per titik. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`benchmark_frequency_grids`](../run_benchmark.py#L261) | Membentuk grid frekuensi aktif dan pasif dengan node overlap yang sama, dengan kerapatan yang bergantung pada mode cepat. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`_benchmark_worker`](../run_benchmark.py#L274) | Menjalankan satu kasus dalam proses terpisah dan mengembalikan ringkasan serta pemeriksaan kuadratur. | — |
| [`main`](../run_benchmark.py#L280) | Membaca CLI, memilih anggaran serta cakupan eksperimen, memeriksa folder keluaran, dan mengatur eksekusi, manifest, snapshot serta laporan. | — |

### [models.py](../models.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`LayeredEarthModel.n_layers`](../models.py#L56) | Mengembalikan jumlah entri Vs, termasuk half-space. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`LayeredEarthModel.poisson_ratio`](../models.py#L60) | Menghitung rasio Poisson elastik per lapisan dari Vp dan Vs. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`LayeredEarthModel.vp_vs_ratio`](../models.py#L66) | Mengembalikan rasio Vp/Vs per lapisan. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`LayeredEarthModel.interface_depths`](../models.py#L70) | Menghitung kedalaman kumulatif antarmuka lapisan berhingga. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`LayeredEarthModel.to_flat_arrays`](../models.py#L74) | Mengembalikan tuple (h, Vs, Vp, rho, Qs, Qp) sebagai array float64 kontigu. | — |
| [`LayeredEarthModel.summary_table`](../models.py#L88) | Merender parameter model dan kedalaman antarmuka sebagai tabel teks. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.__post_init__`](../models.py#L135) | Menyusun nama parameter dengan urutan Vs, Vp, dan ketebalan berhingga. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.dim`](../models.py#L144) | Mengembalikan jumlah parameter inversi, 11 untuk A/B. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.lower_bounds`](../models.py#L148) | Menggabungkan batas bawah kotak Vs, Vp, dan ketebalan. | — |
| [`InversionParameterSpace.upper_bounds`](../models.py#L152) | Menggabungkan batas atas kotak Vs, Vp, dan ketebalan. | — |
| [`InversionParameterSpace.true_vector`](../models.py#L156) | Mengubah model acuan menjadi vektor numerik untuk evaluasi. | — |
| [`InversionParameterSpace.model_to_vector`](../models.py#L159) | Mengepak Vs, Vp dan ketebalan berhingga dari model ke vektor float64. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.vector_to_model`](../models.py#L169) | Memvalidasi vektor fisik dan membentuk model dengan properti tetap dari acuan. Menambahkan ketebalan nol sebagai penanda half-space. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace._conditional_lower`](../models.py#L197) | Membangun batas bawah yang bergantung pada Vs untuk koordinat Vp. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.physical_to_latent`](../models.py#L203) | Memetakan parameter fisik interior ke koordinat logit. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.latent_to_physical`](../models.py#L211) | Memetakan satu vektor atau batch koordinat laten ke model fisik dengan sigmoid kondisional. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.latent_jacobian_diag`](../models.py#L225) | Menghitung diagonal Jacobian transformasi laten ke fisik. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`InversionParameterSpace.latent_log_abs_det_jacobian`](../models.py#L232) | Menjumlahkan log diagonal Jacobian secara stabil pada sumbu parameter. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`get_benchmark_a_lvl`](../models.py#L241) | Mengembalikan model acuan A dengan low-velocity layer dan ruang prior A. | — |
| [`get_benchmark_b_water_table`](../models.py#L279) | Mengembalikan model acuan B dengan kontras Vp dangkal dan ruang prior B. | — |
| [`evaluate_depth_profiles`](../models.py#L318) | Menyampel model ensemble pada grid kedalaman, lalu mengembalikan profil Vs, Vp dan rasio Poisson serta ringkasan kuantilnya. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |


## Sintesis, geometri, dan fusi

### [forward_mfj.py](../forward_mfj.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`njit`](../forward_mfj.py#L23) | Menyediakan dekorator cadangan ketika Numba tidak bisa diimpor. | — |
| [`njit.decorator`](../forward_mfj.py#L24) | Mengembalikan fungsi Python tanpa kompilasi. | — |
| [`_bessel_j0_y0_tensor_numba`](../forward_mfj.py#L33) | Menghitung tensor J0 dan Y0 untuk argumen yang sama dengan aproksimasi rasional dan asimtotik. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`ActiveSurveyGeometry.offsets_primary`](../forward_mfj.py#L146) | Memberi jarak jejak untuk tembakan primer dengan offset dekat 4 m. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`ActiveSurveyGeometry.offsets_secondary`](../forward_mfj.py#L151) | Memberi jarak jejak untuk tembakan sekunder dengan offset dekat 10 m. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`ActiveSurveyGeometry.combined_unique_offsets`](../forward_mfj.py#L156) | Menghitung jarak sumber–receiver dari konfigurasi tembakan dua arah dan mengembalikan jarak unik terurut. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`PassiveSurveyGeometry.station_xy`](../forward_mfj.py#L198) | Mengembalikan koordinat yang diberikan atau membentuk geometri linear maupun inti nested dengan stasiun spiral. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`PassiveSurveyGeometry._compute_binned_2d_pairs`](../forward_mfj.py#L227) | Menghitung jarak, jumlah pasangan per lag atau bin, serta faktor azimut. Dapat mempertahankan lag aktual hingga pembulatan numerik. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`PassiveSurveyGeometry.interstation_offsets`](../forward_mfj.py#L292) | Mengakses vektor jarak radial hasil pengelompokan pasangan pasif. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`PassiveSurveyGeometry.pair_redundancy_weights`](../forward_mfj.py#L298) | Mengakses jumlah pasangan yang menyumbang ke setiap lag radial. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`PassiveSurveyGeometry.azimuthal_anisotropy_factors`](../forward_mfj.py#L304) | Mengakses faktor arah kosinus orde dua per lag untuk gangguan sintetik. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) · [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) |
| [`_invert_2x2_numba`](../forward_mfj.py#L317) | Mengembalikan empat entri invers matriks kompleks 2×2 dengan pelindung determinan. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`_psv_surface_compliance_single_wk`](../forward_mfj.py#L329) | Mengevaluasi satu compliance vertikal kompleks dari frekuensi sudut, bilangan gelombang dan properti lapisan melalui rekursi impedansi. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`_psv_surface_compliance_grid_numba`](../forward_mfj.py#L463) | Mengulang compliance skalar pada grid frekuensi–bilangan gelombang dengan pergeseran imajiner yang diberikan. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`PyProp8ForwardAdapter.is_available`](../forward_mfj.py#L510) | Melaporkan apakah pustaka PyProp8 berhasil diimpor. | — |
| [`PyProp8ForwardAdapter.compute_seismograms`](../forward_mfj.py#L514) | Menghasilkan seismogram elastik komponen vertikal melalui PyProp8 untuk sumber hampir di permukaan, offset, jumlah sampel dan dt. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`PyProp8ForwardAdapter.compute_seismograms.ricker`](../forward_mfj.py#L546) | Memberi bobot sumber aktif pada frekuensi sudut yang diminta PyProp8. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) |
| [`CPSHspec96Adapter.is_available`](../forward_mfj.py#L575) | Memeriksa keberadaan executable hspec96 pada PATH. | — |
| [`CPSHspec96Adapter.write_model96`](../forward_mfj.py#L581) | Menulis properti lapisan ke format MODEL96 dengan konversi satuan untuk CPS. | — |
| [`_mfj_hankel_transform_numba`](../forward_mfj.py#L616) | Memutihkan gather, menerapkan taper dan bobot radial, melakukan proyeksi kernel, lalu menormalkan citra menurut puncak per frekuensi. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`compute_trapezoidal_weights`](../forward_mfj.py#L689) | Menghitung vektor bobot trapesium pada koordinat terurut, termasuk aturan satu node. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) |
| [`build_modified_hanning_taper`](../forward_mfj.py#L703) | Membentuk taper kosinus tepi pada offset, dengan bobot satu di bagian tengah. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`precompute_hankel_kernels`](../forward_mfj.py#L722) | Menyusun dua kernel Hankel ternormalisasi dari frekuensi, kecepatan dan offset. Kernel orde nol cabang kedua dipakai dalam MFJ. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`FullWavefieldMFJOperator.__init__`](../forward_mfj.py#L772) | Menyimpan geometri dan grid, lalu menyiapkan offset, bobot jarak, taper, dan kernel setiap kanal. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) · [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`FullWavefieldMFJOperator.quadrature_grid`](../forward_mfj.py#L803) | Membentuk grid k dan basis Bessel berbobot per frekuensi serta offset. Dapat menyimpan hasil dalam cache. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) |
| [`FullWavefieldMFJOperator._green_gather`](../forward_mfj.py#L831) | Mengintegrasikan compliance, mengembalikan suku statik dan menghasilkan respons kompleks frekuensi–offset tanpa spektrum sumber. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) · [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) |
| [`FullWavefieldMFJOperator.synthesize_frequency_gather`](../forward_mfj.py#L854) | Mengalikan Green dengan spektrum sumber kanal dan menambah gangguan bila diminta. Mengembalikan gather kompleks dan placeholder None. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) |
| [`FullWavefieldMFJOperator.synthesize_time_gather`](../forward_mfj.py#L898) | Menghasilkan waktu, offset dan jejak riil dari evaluasi Green pada grid FFT dengan penundaan sumber. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) |
| [`FullWavefieldMFJOperator.compute_mfj_spectrum`](../forward_mfj.py#L927) | Memilih kernel serta bobot kanal, lalu mengubah gather kompleks menjadi citra MFJ. Argumen kompatibilitas CC-FJ tidak mengubah backend. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`FullWavefieldMFJOperator.forward_active_passive_spectra`](../forward_mfj.py#L953) | Mengembalikan dua citra MFJ dan dua gather kompleks pada grid native masing-masing kanal. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`compute_theoretical_dispersion_curves`](../forward_mfj.py#L990) | Menghitung pasangan frekuensi–kecepatan Rayleigh dari disba dan menyaringnya ke grid kecepatan yang diminta. Daftar guided P kosong. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [acquisition_design.py](../acquisition_design.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`FieldArrayDesign.n_stations`](../acquisition_design.py#L27) | Mengembalikan jumlah baris koordinat stasiun pasif. | — |
| [`FieldArrayDesign.pairwise_distances_and_azimuths`](../acquisition_design.py#L30) | Mengembalikan jarak radial dan sudut untuk setiap pasangan stasiun unik. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`build_four_candidate_field_designs`](../acquisition_design.py#L46) | Membentuk konfigurasi geometri linear, L-array, cincin dan nested 2D beserta bentangan aktif. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`evaluate_acquisition_metrics`](../acquisition_design.py#L125) | Merangkum jumlah pasangan, rentang jarak, persentase bin terisi, momen azimut, dan proksi frekuensi batas. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`make_operator`](../acquisition_design.py#L145) | Menerjemahkan satu rancangan akuisisi dan grid menjadi operator forward dan MFJ. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`plot_acquisition_design_comparison`](../acquisition_design.py#L164) | Menyimpan gambar tata letak, distribusi jarak serta arah, dan mengembalikan metrik rancangan. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`transform_spectrum`](../acquisition_design.py#L207) | Memproyeksikan gather dengan kernel phase-only Hankel, gelombang datar atau Bessel berdiri dengan pra-pembobotan yang konsisten. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`plot_dispersion_image_acquisition_comparison`](../acquisition_design.py#L225) | Membentuk citra bersih dan bising pada rancangan yang dibandingkan, menyimpan panel dan melaporkan NCC overlap serta MAE terhadap referensi masing-masing. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [fusion.py](../fusion.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`njit`](../fusion.py#L20) | Menyediakan dekorator cadangan ketika Numba tidak tersedia. | — |
| [`njit.decorator`](../fusion.py#L21) | Mengembalikan fungsi asli yang dibungkus. | — |
| [`normalize_frequency_slices_numba`](../fusion.py#L59) | Melakukan min–max per irisan, dengan eksponen kontras opsional. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`estimate_slice_variance_numba`](../fusion.py#L91) | Menghitung proksi bentuk spektral per irisan dari kuadrat rasio rata-rata ke puncak, dengan batas bawah. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`ActivePassiveSpectralFuser.__init__`](../fusion.py#L116) | Menyiapkan frekuensi fusi dari grid native dan preferensi apertur, memakai node bersama di bagian dalam overlap. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`ActivePassiveSpectralFuser._setup_broadband_interpolation`](../fusion.py#L152) | Menghitung bobot apertur transisi sinus kuadrat sebagai fungsi frekuensi. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`ActivePassiveSpectralFuser._interp_spectrum_to_broadband`](../fusion.py#L170) | Menyalin irisan yang punya node native cocok ke grid fusi, dengan toleransi numerik 1e-9 Hz. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`ActivePassiveSpectralFuser.fuse_spectra`](../fusion.py#L183) | Mengembalikan citra fusi, irisan sejajar tiap kanal, bobot dan proksi bentuknya dari dua citra native. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`ActivePassiveSpectralFuser.generate_fused_spectrum_for_model`](../fusion.py#L281) | Menjalankan sintesis dan MFJ dua kanal, lalu fusi. Menyertakan citra mentah serta gather dalam kamus keluaran. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |


## Target, inferensi, dan prediksi

### [smvi_inversion.py](../smvi_inversion.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`njit`](../smvi_inversion.py#L20) | Membungkus fungsi tanpa JIT ketika Numba tidak ada. | — |
| [`discrete_wasserstein2`](../smvi_inversion.py#L30) | Menghitung biaya transport kuadratik diskrit dengan memasangkan sisa massa secara monoton pada grid kecepatan aktual. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`compute_spectral_misfit_numba`](../smvi_inversion.py#L76) | Memvalidasi dimensi dan menghitung J, rata-rata 1−NCC serta akar rata-rata W2 kuadrat pada irisan yang dipertahankan. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`PickingFreeSpectralLikelihood.__init__`](../smvi_inversion.py#L160) | Memeriksa data dan mengikat fuser, prior, kanal, mask, skala transport dan konfigurasi target. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`PickingFreeSpectralLikelihood.set_frequency_weights`](../smvi_inversion.py#L233) | Mengintegrasikan bobot pada masing-masing pita kontigu yang didukung kanal dan mask. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`PickingFreeSpectralLikelihood.spectral_metrics`](../smvi_inversion.py#L262) | Menerapkan mask kecepatan pada dua citra dan mengembalikan tiga metrik kecocokan. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`PickingFreeSpectralLikelihood.evaluate_physical_model`](../smvi_inversion.py#L281) | Menghasilkan citra kandidat tanpa noise dari vektor fisik serta log likelihood dan metriknya. Menaikkan penghitung evaluasi. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`PickingFreeSpectralLikelihood.log_posterior_latent`](../smvi_inversion.py#L290) | Memetakan kandidat laten, menjumlahkan log likelihood dan Jacobian, dan menolak kandidat invalid. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`PickingFreeSpectralLikelihood.score_and_post_latent`](../smvi_inversion.py#L306) | Memilih autodiff atau beda hingga sentral untuk menghitung gradien target, bersama nilai target dan metrik. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`SteinMixtureVariationalInference.__init__`](../smvi_inversion.py#L339) | Memvalidasi konfigurasi dan menyiapkan mean laten, faktor Cholesky, bobot sama dan RNG. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference._cholesky`](../smvi_inversion.py#L354) | Membentuk faktor triangular bawah dengan diagonal eksponensial positif. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference._initialize_mixture_components`](../smvi_inversion.py#L357) | Menginisialisasi pusat dari Latin hypercube, menyaring berdasarkan target dan jarak, lalu dapat menjalankan refinemen L-BFGS-B. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference._initialize_mixture_components.objective`](../smvi_inversion.py#L390) | Mengembalikan negatif target laten dan negatif gradien untuk minimizer. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference._mixture_details`](../smvi_inversion.py#L415) | Mengembalikan log densitas campuran, skor terhadap z, dan tanggung jawab komponen. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference._compute_mixture_log_prob_and_score`](../smvi_inversion.py#L431) | Meneruskan log densitas dan skor, dengan membuang tanggung jawab yang tidak dibutuhkan caller. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference._stein_kernel_update`](../smvi_inversion.py#L435) | Menggabungkan gaya per komponen melalui kernel RBF pada mean dan parameter Cholesky dengan repulsi Stein. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference.run`](../smvi_inversion.py#L453) | Menjalankan optimasi dan sampling campuran, menilai pusat terbaik dan rerata parameter, lalu membangun prediksi spektrum beserta riwayat dan diagnostik. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`SteinMixtureVariationalInference.run.update`](../smvi_inversion.py#L462) | Melakukan satu langkah Adam ascent dengan koreksi bias dan pembatasan gradien. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`SteinMixtureVariationalInference.sample_posterior`](../smvi_inversion.py#L595) | Memilih label menurut bobot, menarik Gaussian laten, dan memetakan draw ke parameter fisik. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`AffineInvariantEnsembleMCMC.__init__`](../smvi_inversion.py#L607) | Mengikat likelihood dan jumlah walker minimum, serta memvalidasi burn-in dan skala stretch. | — |
| [`AffineInvariantEnsembleMCMC.run`](../smvi_inversion.py#L627) | Memperbarui walker dengan stretch move dan melaporkan chain, sampel, tingkat penerimaan, autokorelasi, ESS serta status panjang chain. | [MCMC](PERSAMAAN_KODE.md#ensemble-mcmc) |

### [autodiff_likelihood.py](../autodiff_likelihood.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`_inverse`](../autodiff_likelihood.py#L14) | Menghitung adjugat dibagi determinan pada tensor 2×2 dengan pelindung determinan. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`compliance`](../autodiff_likelihood.py#L20) | Menjalankan rekursi compliance P–SV dalam PyTorch untuk grid kandidat, dengan gradien yang dapat ditelusuri. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`normalize`](../autodiff_likelihood.py#L132) | Menormalkan irisan tensor dengan min–max dan pangkat kontras. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`wasserstein2`](../autodiff_likelihood.py#L139) | Mengintegrasikan biaya jarak kuadrat pada interval gabungan breakpoint CDF setelah normalisasi massa. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`AutodiffSpectralLikelihood.__init__`](../autodiff_likelihood.py#L166) | Membangun tensor tetap untuk prior, observasi, k, basis Bessel, kernel, sumber dan indeks native kedua kanal. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) · [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`AutodiffSpectralLikelihood._physical`](../autodiff_likelihood.py#L224) | Memetakan tensor z ke Vs, Vp, ketebalan dengan half-space dan log-determinan Jacobian. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`AutodiffSpectralLikelihood._spectrum`](../autodiff_likelihood.py#L244) | Menyintesis citra kanal bersih, menormalkannya dan mengembalikan irisan pada grid fusi serta mask ketersediaan. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`AutodiffSpectralLikelihood.spectrum_and_target`](../autodiff_likelihood.py#L268) | Menggabungkan kanal atau memilih kanal yang diminta, menerapkan mask, lalu mengembalikan spektrum, target dan tiga metrik sebagai tensor. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`AutodiffSpectralLikelihood.score_and_post_latent`](../autodiff_likelihood.py#L343) | Menghitung gradien PyTorch terhadap z, memeriksa nilai hingga, lalu mengembalikan gradien dan metrik ke caller NumPy. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |

### [posterior_prediction.py](../posterior_prediction.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`evaluate_mixture_prediction`](../posterior_prediction.py#L11) | Menarik sampel seimbang per komponen, mencitrakan setiap draw, dan mengembalikan rata-rata, spread, galat MC serta laporan presisi adaptif. | [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |


## Eksperimen lanjutan dan pelaporan

### [ablation_experiments.py](../ablation_experiments.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`profile_metrics`](../ablation_experiments.py#L23) | Membandingkan mean profil Vs ensemble dengan acuan pada grid kedalaman tetap dan menghasilkan RMSE relatif serta MAPE. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`fit_case`](../ablation_experiments.py#L34) | Menyiapkan likelihood dengan kanal, bobot pita dan mask tertentu, lalu mengembalikan hasil SMVI beserta konfigurasi efektif. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`save_inference`](../ablation_experiments.py#L56) | Mengekspor entri array hasil ke NPZ dan konfigurasi ke JSON bila ada. | — |
| [`summary_case`](../ablation_experiments.py#L65) | Merangkum statistik parameter, galat profil, metrik model mean, biaya evaluasi, dan konfigurasi. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`run_phd_ablation_suite`](../ablation_experiments.py#L78) | Menjalankan pembanding kanal, kurva picking jika disba tersedia, jendela Rayleigh, Gaussian tunggal, dan ensemble MCMC. Menyimpan metrik serta gambar. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) · [MCMC](PERSAMAAN_KODE.md#ensemble-mcmc) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`run_phd_ablation_suite.residual`](../ablation_experiments.py#L158) | Menghitung residual kurva picking terhadap referensi moda Rayleigh untuk least squares, dengan perlakuan jika dukungan kurva kurang. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`add_complex_noise`](../ablation_experiments.py#L266) | Menghasilkan noise kompleks terukur dan mengembalikan gather bising beserta SNR aktual dalam dB. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`stress_observation`](../ablation_experiments.py#L274) | Melemahkan sumber, menambah gangguan kopling serta noise pada gather, lalu mencitrakan dan menghitung ulang bobot fusi. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`run_noise_and_weak_source_stress_test`](../ablation_experiments.py#L308) | Membandingkan skenario sumber/noise/geometri melalui fit ulang dan menyimpan gambar serta metrik kondisi aktual. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [visualize.py](../visualize.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`gaussian_kde.__init__`](../visualize.py#L44) | Menyiapkan sampel dan bandwidth diagonal untuk implementasi cadangan KDE. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`gaussian_kde.__call__`](../visualize.py#L53) | Menghitung rerata kernel Gaussian ternormalisasi pada titik evaluasi. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`_overlay_dispersion_branches`](../visualize.py#L91) | Menggambar kurva Rayleigh referensi pada sumbu spektrum jika tersedia. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`select_radial_traces`](../visualize.py#L112) | Memilih jejak mendekati target radial seragam dan menghitung lebar wiggle dari jarak minimum yang dipilih. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`plot_waveforms_and_fusion`](../visualize.py#L123) | Membuat panel jejak waktu, citra tiap kanal, citra fusi, dan bobot kanal dari observasi serta geometri. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) · [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`plot_fc_spectrum_comparison`](../visualize.py#L310) | Membandingkan observasi dan prediksi spektrum, menampilkan residual dan riwayat optimasi, lalu menyimpan gambar. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`plot_posterior_profiles_and_marginals`](../visualize.py#L447) | Menyusun profil, kuantil kedalaman dan distribusi parameter dari sampel yang tersimpan. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [generate_pdf_manuscript.py](../generate_pdf_manuscript.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`build_manuscript`](../generate_pdf_manuscript.py#L13) | Merender teks naskah dari report tersimpan, angka benchmark, konfigurasi, dan jalur hasil. | — |
| [`build_diagnostic_report`](../generate_pdf_manuscript.py#L251) | Merender laporan diagnostik dengan metrik, peringatan terstruktur dan catatan editorial yang disimpan terpisah. | — |
| [`main`](../generate_pdf_manuscript.py#L312) | Memvalidasi manifest dan status sebelum menulis Markdown atau membangun PDF. Memiliki jalur khusus laporan diagnostik. | — |

### [replot_saved_results.py](../replot_saved_results.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`main`](../replot_saved_results.py#L28) | Memeriksa kompatibilitas sumber dengan hasil tersimpan, memuat array dan membuat ulang gambar yang diminta. | — |

### [research_quality.py](../research_quality.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`assess_research_readiness`](../research_quality.py#L4) | Mengevaluasi status run, verifikasi fisika, kalibrasi, konvergensi, presisi prediktif dan uji resolusi. Mengembalikan penilaian yang dapat diaudit. | — |
| [`require_research_readiness`](../research_quality.py#L45) | Menghentikan penerbitan ketika pemeriksaan kesiapan belum terpenuhi. | — |


## Pemeriksaan numerik dan diagnostik

### [validate_reassessment.py](../validate_reassessment.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`independent_compliance_check`](../validate_reassessment.py#L15) | Memasangkan moduli kompleks dan grid P–SV dengan PyProp8, lalu mengukur galat compliance pada titik uji A/B. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`main`](../validate_reassessment.py#L56) | Menjalankan pemeriksaan solver dan menyimpan ringkasan serta hash kode. | — |

### [diagnostics_review_v6/benchmark_b_residual/analyze.py](../diagnostics_review_v6/benchmark_b_residual/analyze.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`spectra`](../diagnostics_review_v6/benchmark_b_residual/analyze.py#L40) | Menghasilkan citra gabungan, aktif dan pasif untuk vektor fisik Benchmark B dalam urutan kanal yang dipakai audit. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`row_terms`](../diagnostics_review_v6/benchmark_b_residual/analyze.py#L63) | Memisahkan kontribusi NCC dan transport per frekuensi, lalu membagi kontribusi NCC ke tiga pita kecepatan. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [diagnostics_review_v6/benchmark_b_residual/plot.py](../diagnostics_review_v6/benchmark_b_residual/plot.py)

Membaca hasil analisis residual B dan merender kontribusi per pita dan kanal ke PNG. [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik).

### [diagnostics_review_v6/benchmark_b_restarts/check_new_noise.py](../diagnostics_review_v6/benchmark_b_restarts/check_new_noise.py)

Membuat atau memuat delapan realisasi noise baru B. Menilai prediksi tersimpan bersama prior pada observasi yang sama dan menulis hasil berpasangan. [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik).

### [diagnostics_review_v6/benchmark_b_restarts/check_resolution.py](../diagnostics_review_v6/benchmark_b_restarts/check_resolution.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`clean`](../diagnostics_review_v6/benchmark_b_restarts/check_resolution.py#L48) | Membangun citra MFJ aktif tanpa gangguan untuk model acuan atau kandidat pada grid diagnostik rapat. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |

### [diagnostics_review_v6/benchmark_b_restarts/check_waveform.py](../diagnostics_review_v6/benchmark_b_restarts/check_waveform.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`independent`](../diagnostics_review_v6/benchmark_b_restarts/check_waveform.py#L51) | Membentuk gather aktif B dengan propagator PyProp8, kuadratur linier dan pemulihan suku statik, dengan langkah dan batas k yang diterima sebagai masukan. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) · [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`to_time`](../diagnostics_review_v6/benchmark_b_restarts/check_waveform.py#L76) | Memberi bobot sumber dan menempatkan spektrum terkonjugasi pada bin FFT untuk mendapatkan jejak band-limited. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) |
| [`relative`](../diagnostics_review_v6/benchmark_b_restarts/check_waveform.py#L82) | Mengembalikan selisih L2 relatif terhadap array pembanding. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [diagnostics_review_v6/benchmark_b_restarts/run.py](../diagnostics_review_v6/benchmark_b_restarts/run.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`main`](../diagnostics_review_v6/benchmark_b_restarts/run.py#L30) | Memeriksa hash dan observasi asli B, mengambil seed optimasi, menjalankan fit ulang dengan pengaturan sepadan dan menyimpan hasil serta ringkasan. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |

### [diagnostics_review_v6/benchmark_b_restarts/summarize.py](../diagnostics_review_v6/benchmark_b_restarts/summarize.py)

Merangkum hasil restart B terhadap prior, memeriksa batas presisi, dan menulis laporan dari prediksi tersimpan. [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo).

### [diagnostics_review_v6/benchmark_composti/compare.py](../diagnostics_review_v6/benchmark_composti/compare.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`run_composti`](../diagnostics_review_v6/benchmark_composti/compare.py#L19) | Meneruskan lapisan, frekuensi, sumber, dan offset ke ekstensi Levin untuk mengembalikan perpindahan vertikal kompleks lewat berkas pertukaran. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`native_without_shift`](../diagnostics_review_v6/benchmark_composti/compare.py#L36) | Mengevaluasi compliance pada frekuensi sudut riil dengan grid v6 atau rentang slowness dan taper COMPOSTI. Menormalkan respons per gaya titik dengan pembagi 2π. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) · [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) |
| [`time_gather`](../diagnostics_review_v6/benchmark_composti/compare.py#L80) | Menempatkan spektrum pada bin FFT yang sesuai dan menghasilkan sumbu waktu serta jejak riil. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) |
| [`relative_l2`](../diagnostics_review_v6/benchmark_composti/compare.py#L89) | Mengukur selisih relatif terhadap referensi pada sumbu yang dipilih atau pada seluruh array. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`compare_case`](../diagnostics_review_v6/benchmark_composti/compare.py#L93) | Menyamakan model, offset dan wavelet A/B, menjalankan solver, mengukur selisih spektrum, jejak dan MFJ, lalu menyimpan hasil per benchmark. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) · [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`make_figure`](../diagnostics_review_v6/benchmark_composti/compare.py#L180) | Membuat panel jejak pada empat offset, tiga citra, dan residual v6 terhadap COMPOSTI. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`main`](../diagnostics_review_v6/benchmark_composti/compare.py#L211) | Memilih mode pekerja atau audit A/B, lalu mencatat API, commit upstream, normalisasi dan ringkasan pengujian. | — |

### [diagnostics_review_v6/check_followup_controls.py](../diagnostics_review_v6/check_followup_controls.py)

Menjalankan pemeriksaan massa transport, skala parent, bobot pita dan presisi prediksi. Menulis ringkasan kontrol. [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo).

### [diagnostics_review_v6/finalize_verification.py](../diagnostics_review_v6/finalize_verification.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`main`](../diagnostics_review_v6/finalize_verification.py#L8) | Memeriksa bahwa laporan dapat diregenerasi secara konsisten, catatan editorial dipertahankan, gerbang naskah bekerja, dan hash sumber serta artefak terlacak. | — |

### [diagnostics_review_v6/verify_saved_results.py](../diagnostics_review_v6/verify_saved_results.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`main`](../diagnostics_review_v6/verify_saved_results.py#L18) | Memuat observasi dan hasil tersimpan, menjalankan ulang forward serta diagnostik Monte Carlo, dan menguji apakah metrik serta prediksi bisa direproduksi. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) · [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |

### [output/pdf/figures/render_b_labels.py](../output/pdf/figures/render_b_labels.py)

Memuat manifest dan artefak B yang cocok, lalu membuat ulang tiga panel dengan label yang digunakan naskah. [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik).


## Tes

### [tests/test_followup_fixes.py](../tests/test_followup_fixes.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`fake_fuser`](../tests/test_followup_fixes.py#L18) | Membuat fuser tiruan untuk pengujian target dan bobot frekuensi. | — |
| [`test_transport_rejects_nonfinite_masses_without_hanging`](../tests/test_followup_fixes.py#L28) | Memastikan transport menolak massa tidak hingga tanpa terhenti pada loop. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_likelihood_rejects_invalid_observations_before_forward`](../tests/test_followup_fixes.py#L53) | Memastikan observasi invalid ditolak sebelum sintesis. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_fixed_wasserstein_scale_is_invariant_to_empty_velocity_extension`](../tests/test_followup_fixes.py#L61) | Menguji invariansi biaya ketika support kecepatan bertambah tanpa mengubah transport dan skala tetap. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_channel_frequency_integral_and_boundary_weights`](../tests/test_followup_fixes.py#L82) | Memastikan tiap kanal mengintegrasikan pita dan bobot ujung yang semestinya. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`prediction_target`](../tests/test_followup_fixes.py#L99) | Menyediakan konfigurasi serta target tiruan untuk prediksi Monte Carlo. | [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`prediction_target.Target.evaluate_physical_model`](../tests/test_followup_fixes.py#L103) | Menghasilkan citra deterministik tiruan dari parameter untuk menguji prediksi. | [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`prediction_target.Target.spectral_metrics`](../tests/test_followup_fixes.py#L112) | Mengukur metrik buatan dari mean citra dalam tes prediksi. | [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`test_balanced_prediction_uses_mixture_weights_and_separates_mc_error`](../tests/test_followup_fixes.py#L118) | Memeriksa bobot campuran, alokasi draw seimbang, dan pemisahan spread respons dari galat mean. | [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`test_predictive_precision_failure_is_reported_at_the_budget_limit`](../tests/test_followup_fixes.py#L146) | Memastikan flag presisi gagal ketika anggaran habis tanpa memenuhi toleransi. | [Prediksi dan MC](PERSAMAAN_KODE.md#prediksi-dan-presisi-monte-carlo) |
| [`test_report_regeneration_preserves_structured_warnings_and_separate_notes`](../tests/test_followup_fixes.py#L167) | Memeriksa determinisme laporan dan pelestarian peringatan serta catatan di berkas terpisah. | — |
| [`test_likelihood_rejects_complex_spectra_instead_of_discarding_imaginary_data`](../tests/test_followup_fixes.py#L223) | Memastikan input citra kompleks tidak diterima sebagai spektrum riil. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_likelihood_rejects_dc_frequency_before_green_evaluation`](../tests/test_followup_fixes.py#L230) | Memastikan frekuensi nol ditolak sebelum perhitungan Green. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |

### [tests/test_reassessment.py](../tests/test_reassessment.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`test_layered_compliance_against_independent_pyprop8`](../tests/test_reassessment.py#L22) | Menguji compliance berlapis A/B terhadap propagator PyProp8. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`test_fusion_uses_native_slices_without_inventing_frequencies`](../tests/test_reassessment.py#L30) | Memastikan proyeksi irisan dan fusi mengikuti support native yang sah. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |
| [`test_autodiff_matches_numba_spectra_target_and_finite_difference`](../tests/test_reassessment.py#L61) | Membandingkan spektrum, nilai target, dan gradien antara backend serta beda hingga. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_wiggles_include_aperture_ends_and_do_not_overlap`](../tests/test_reassessment.py#L97) | Memastikan tepi apertur disertakan dan lebar jejak sesuai jaraknya. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`test_stein_kernel_acts_on_covariances_and_has_correct_normalization`](../tests/test_reassessment.py#L106) | Memeriksa faktor normalisasi dan pengaruh parameter Cholesky pada kernel. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_multiple_components_recover_both_known_modes`](../tests/test_reassessment.py#L126) | Mengukur apakah campuran menangkap dua moda target buatan beserta statistik sampelnya. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_multiple_components_recover_both_known_modes.Target.log_posterior_latent`](../tests/test_reassessment.py#L141) | Menghitung log densitas untuk target dua moda analitik. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_multiple_components_recover_both_known_modes.Target.score_and_post_latent`](../tests/test_reassessment.py#L149) | Mengembalikan skor analitik target dua moda dan metrik dummy tes. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_multiple_components_recover_both_known_modes.Target.evaluate_physical_model`](../tests/test_reassessment.py#L157) | Mengembalikan citra dan metrik tiruan yang diperlukan oleh antarmuka SMVI. | — |
| [`test_completed_quick_run_cannot_be_published`](../tests/test_reassessment.py#L169) | Memastikan run cepat berstatus selesai belum otomatis melewati syarat penerbitan. | — |
| [`test_independent_waveform_adapter_runs_without_hidden_failure`](../tests/test_reassessment.py#L178) | Memastikan adapter PyProp8 menghasilkan jejak berdimensi dan bernilai sah. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`test_fft_convention_gives_expected_homogeneous_rayleigh_arrivals`](../tests/test_reassessment.py#L190) | Memeriksa waktu kedatangan Rayleigh homogen dan arah moveout setelah transformasi FFT. | [Sumber, noise, FFT](PERSAMAAN_KODE.md#sumber-gangguan-dan-fft) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`test_benchmark_overlap_is_synthesized_at_matching_native_frequencies`](../tests/test_reassessment.py#L222) | Memastikan grid aktif dan pasif memakai node native identik di overlap. | [Fusi](PERSAMAAN_KODE.md#fusi-frekuensi-native) |

### [tests/test_scientific_regressions.py](../tests/test_scientific_regressions.py)

| Definisi | Peran dan hasil | Persamaan dan rujukan |
| --- | --- | --- |
| [`test_conditional_transform_and_jacobian`](../tests/test_scientific_regressions.py#L15) | Menguji round-trip transformasi prior, syarat Vp/Vs dan Jacobian terhadap beda hingga. | [Model dan Jacobian](PERSAMAAN_KODE.md#model-berlapis-dan-jacobian) |
| [`test_discrete_wasserstein_translation_nonuniform_grid`](../tests/test_scientific_regressions.py#L25) | Menguji biaya kuadratik antara dua massa yang bergeser di grid tidak seragam. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_frequency_mask_uses_hz`](../tests/test_scientific_regressions.py#L32) | Memastikan mask pasif dikendalikan oleh nilai Hz. | [Transport dan target](PERSAMAAN_KODE.md#transport-1d-dan-target-inversi) |
| [`test_full_covariance_density_score_and_sampling_agree`](../tests/test_scientific_regressions.py#L40) | Memastikan faktor Cholesky konsisten dalam densitas, skor dan pengambilan sampel. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_initialization_does_not_read_truth`](../tests/test_scientific_regressions.py#L52) | Memastikan penyiapan komponen tidak mengakses parameter acuan. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_initialization_does_not_read_truth.NoTruth.true_vector`](../tests/test_scientific_regressions.py#L56) | Menggagalkan tes bila inisialisasi mencoba membaca model acuan. | — |
| [`test_mcmc_normal_has_no_clipping_atoms`](../tests/test_scientific_regressions.py#L62) | Memeriksa sebaran sampel MCMC pada target normal dan absennya massa akibat clipping. | [MCMC](PERSAMAAN_KODE.md#ensemble-mcmc) |
| [`test_mcmc_normal_has_no_clipping_atoms.Normal.log_posterior_latent`](../tests/test_scientific_regressions.py#L65) | Menyediakan target normal analitik dengan log densitas proporsional ke negatif setengah norma kuadrat. | [MCMC](PERSAMAAN_KODE.md#ensemble-mcmc) |
| [`test_geometry_preserves_actual_distances`](../tests/test_scientific_regressions.py#L73) | Membandingkan lag pasif dan jumlah pasangan terhadap jarak yang dihitung dari koordinat. | [Geometri](PERSAMAAN_KODE.md#geometri-dan-sampling-spasial) |
| [`test_homogeneous_compliance_matches_analytic_rayleigh_function`](../tests/test_scientific_regressions.py#L82) | Membandingkan compliance setengah ruang homogen terhadap fungsi Rayleigh analitik. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`small_operator`](../tests/test_scientific_regressions.py#L94) | Membangun operator geometri kecil untuk pengujian numerik cepat. | — |
| [`test_imaged_spectrum_depends_on_compliance`](../tests/test_scientific_regressions.py#L99) | Memastikan perubahan compliance memengaruhi citra permukaan yang dihasilkan. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) · [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`test_imaged_spectrum_depends_on_compliance.fake`](../tests/test_scientific_regressions.py#L104) | Mengganti kernel compliance selama tes dengan respons buatan agar pengaruhnya dapat diisolasi. | — |
| [`test_noise_snr_is_measured_before_imaging`](../tests/test_scientific_regressions.py#L110) | Memastikan SNR gather yang dihasilkan sesuai target sebelum transformasi ke citra. | [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`test_synthetic_conventional_transforms_are_actual_kernels`](../tests/test_scientific_regressions.py#L117) | Memeriksa operator pembanding gelombang datar dan Bessel pada masukan sintetik. | [MFJ](PERSAMAAN_KODE.md#pencitraan-mfj) |
| [`test_single_gaussian_optimizer_learns_correlated_target`](../tests/test_scientific_regressions.py#L125) | Menguji apakah satu komponen mempelajari target Gaussian dengan kovariansi berkorelasi. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_single_gaussian_optimizer_learns_correlated_target.GaussianTarget.log_posterior_latent`](../tests/test_scientific_regressions.py#L132) | Mengevaluasi log densitas target Gaussian dari kuadrat jarak Mahalanobis. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_single_gaussian_optimizer_learns_correlated_target.GaussianTarget.score_and_post_latent`](../tests/test_scientific_regressions.py#L134) | Mengembalikan gradien analitik log densitas Gaussian dan diagnostik tiruan. | [Campuran dan Stein](PERSAMAAN_KODE.md#campuran-gaussian-dan-optimasi-stein) |
| [`test_single_gaussian_optimizer_learns_correlated_target.GaussianTarget.evaluate_physical_model`](../tests/test_scientific_regressions.py#L138) | Menyediakan citra dan metrik tiruan agar optimizer dapat dievaluasi secara terisolasi. | — |
| [`test_layer_partition_does_not_change_homogeneous_compliance`](../tests/test_scientific_regressions.py#L145) | Memastikan pemisahan medium homogen menjadi beberapa lapisan identik tidak mengubah compliance. | [Compliance P–SV](PERSAMAAN_KODE.md#compliance-p-sv) |
| [`test_selected_quadrature_converges_on_active_and_passive_offsets`](../tests/test_scientific_regressions.py#L153) | Memeriksa bahwa perubahan densitas kuadratur memenuhi toleransi pada subset uji. | [Kuadratur](PERSAMAAN_KODE.md#kuadratur-sommerfeld-dan-suku-statik) · [Metrik dan verifikasi](PERSAMAAN_KODE.md#metrik-diagnostik) |
| [`test_render_refresh_rejects_changed_numerical_sources`](../tests/test_scientific_regressions.py#L159) | Memastikan pembuatan ulang gambar menolak sumber numerik berubah. Tes membutuhkan hasil historis v4. | — |
| [`test_manuscript_numbers_come_from_completed_report`](../tests/test_scientific_regressions.py#L172) | Memastikan teks laporan memuat angka, status MCMC, dan batas interpretasi yang berasal dari laporan historis v4. | — |
