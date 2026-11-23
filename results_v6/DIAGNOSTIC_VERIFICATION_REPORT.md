# Laporan diagnostik — bukan hasil penelitian tervalidasi

Versi: v6_followup_fixes. Mode: quick_verification.

Status penerbitan: ditahan

Pemeriksaan yang belum lolos: diagnostic_or_quick_run, benchmark_a_optimization_not_validated, benchmark_b_optimization_not_validated, independent_posterior_reference_missing, frequency_velocity_resolution_missing_or_failed, repeated_data_calibration_missing_or_failed, independent_waveform_validation_missing_or_failed

## Benchmark A: Low-Velocity Layer (Mode-Jumping & Osculation Test)

Metrik pada model mean: `{"misfit": 0.6374174160515531, "ncc_similarity": 0.780329846176737, "wasserstein_w2_mps": 201.14412881394347}`.
Metrik prediksi rata-rata sampel: `{"misfit": 0.2976356401120396, "ncc_similarity": 0.8680766579992989, "wasserstein_w2_mps": 131.64880531398052, "draws": 256, "component_draw_counts": [64, 64, 64, 64], "misfit_mc_se": 0.006747039117444793, "misfit_mc_resampling_interval95": [0.2872993389042141, 0.30910767699962677], "mean_spectrum_mc_se_rms": 0.008358127972966492, "precision_target": 0.01, "precision_target_met": true, "bootstrap_repeats": 32, "sampling": "equal draws per component, component means combined with mixture weights", "scope": "Monte Carlo precision conditional on fitted mixture; not posterior calibration", "precision_history": [{"draws": 128, "misfit": 0.3054627628425494, "misfit_mc_se": 0.008007545423982824, "mean_spectrum_mc_se_rms": 0.01178350197319035, "nested_half_misfit_difference": 0.012246461302031109, "precision_target_met": false}, {"draws": 256, "misfit": 0.2976356401120396, "misfit_mc_se": 0.006747039117444793, "mean_spectrum_mc_se_rms": 0.008358127972966492, "nested_half_misfit_difference": 0.007827122730509839, "precision_target_met": true}]}`.

Konvergensi posterior belum dinyatakan. Angka ini adalah diagnosis numerik.

### Perbandingan terhadap prior dan peringatan

Model tengah prior: `{"misfit": 0.3694319105504493, "ncc_similarity": 0.7913616833744423, "wasserstein_w2_mps": 138.29433264588351, "vs_profile_relative_rmse_pct": 15.431769026211137, "vs_profile_mape_pct": 6.113906485671191}`.
Evaluasi model sejati (validasi saja): `{"misfit": 0.17391184626967587, "ncc_similarity": 0.9436665362407017, "wasserstein_w2_mps": 105.97266282993117}`.

**Peringatan hasil:** mean_model_fit_worse_than_prior_midpoint

Galat Monte Carlo misfit: 0.006747039117444793. Jumlah draw: 256. Target presisi terpenuhi: True.
Galat Monte Carlo adalah ketelitian estimasi dari campuran yang dipasang, bukan kalibrasi posterior.

## Benchmark B: Shallow Groundwater Table (Guided-P Sensitivity Test)

Metrik pada model mean: `{"misfit": 0.378673352004771, "ncc_similarity": 0.8048020172474734, "wasserstein_w2_mps": 295.1674628510084}`.
Metrik prediksi rata-rata sampel: `{"misfit": 0.46370034630557255, "ncc_similarity": 0.8106725308655718, "wasserstein_w2_mps": 343.3139452648211, "draws": 512, "component_draw_counts": [128, 128, 128, 128], "misfit_mc_se": 0.007017714935311411, "misfit_mc_resampling_interval95": [0.4549038654133423, 0.4766546937850171], "mean_spectrum_mc_se_rms": 0.0046406844407953455, "precision_target": 0.01, "precision_target_met": true, "bootstrap_repeats": 32, "sampling": "equal draws per component, component means combined with mixture weights", "scope": "Monte Carlo precision conditional on fitted mixture; not posterior calibration", "precision_history": [{"draws": 128, "misfit": 0.44823033200724277, "misfit_mc_se": 0.011909651896918377, "mean_spectrum_mc_se_rms": 0.009206598001005026, "nested_half_misfit_difference": 0.01405231428965692, "precision_target_met": false}, {"draws": 256, "misfit": 0.4617529742759348, "misfit_mc_se": 0.010437450078430633, "mean_spectrum_mc_se_rms": 0.006618619885757127, "nested_half_misfit_difference": 0.013522642268692042, "precision_target_met": false}, {"draws": 512, "misfit": 0.46370034630557255, "misfit_mc_se": 0.007017714935311411, "mean_spectrum_mc_se_rms": 0.0046406844407953455, "nested_half_misfit_difference": 0.0019473720296377395, "precision_target_met": true}]}`.

Konvergensi posterior belum dinyatakan. Angka ini adalah diagnosis numerik.

### Perbandingan terhadap prior dan peringatan

Model tengah prior: `{"misfit": 0.3838779640174564, "ncc_similarity": 0.9019283241687677, "wasserstein_w2_mps": 332.8199427616656, "vs_profile_relative_rmse_pct": 9.819462898345733, "vs_profile_mape_pct": 3.9693555200659416}`.
Evaluasi model sejati (validasi saja): `{"misfit": 0.33542442087561697, "ncc_similarity": 0.9134234878133312, "wasserstein_w2_mps": 310.7876804416959}`.

**Peringatan hasil:** predictive_point_estimate_worse_than_prior_midpoint, posterior_profile_RMSE_worse_than_informative_prior_midpoint

Galat Monte Carlo misfit: 0.007017714935311411. Jumlah draw: 512. Target presisi terpenuhi: True.
Galat Monte Carlo adalah ketelitian estimasi dari campuran yang dipasang, bukan kalibrasi posterior.

## Catatan editorial terpisah

Galat Monte Carlo menyatakan ketelitian estimasi dari campuran yang dipasang. Galat tersebut berbeda dari spread respons model, noise pengukuran, dan kalibrasi posterior.

