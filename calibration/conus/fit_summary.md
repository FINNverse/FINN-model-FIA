# calibration/conus: fit summary

Fitted 2026-11-10 05:24 by `calibration/conus/02_fit.R` in 235 min (35.3 s per epoch, 1 thread) on cpu (FINN 0.2.0, torch 0.15.1, R 4.5.0).

Data: `conus-annual-v1/train` (1500 sites, 20 species, annual step, per-site intervals, 47892 tree intervals). Holdout for validation: `holdout`.

Processes: FINN::growth, FINN::mortality, FINN::regeneration_adult (inits K 50 stems/ha/yr, adult threshold 20 cm, half-saturation 15 adults/ha), FINN::competition; climate network hidden = c(20, 20), dropout 0.1, on `~tavg_c + prec_mm + cwd_z`; reg_floor 0.001; recruits_dbh 12.9; 4 patches of 0.01681 ha.

Fit: 400 epochs, batch 100, lr 0.01, weights "auto" with weight_scale "auto" (baseline-scaled losses), loss families dbh = mse, ba = mse, trees = poisson, growth = mse, mortality = binomial, regeneration = nbinom, growth and mortality aggregated per tree, env_autoscale TRUE, seed 1.

Final losses (baseline-scaled, 1 = no better than the mean): dbh 0.505, ba 0.087, trees 0.199, growth 0.629, mort 0.873, reg 0.283.

Fitted adult threshold x: min 3.3, median 7.8, max 62.4 cm. Half-saturation: min 5.6, median 84.6, max 325.3 adults/ha. Caps K: median 31.6 stems/ha/yr (adult_parameters.csv).

Environment scaling stored with the model: tavg_c (center 6.709, scale 3.634); prec_mm (center 478.305, scale 318.727); cwd_z (center 0.380, scale 0.936).

parameters.pt: 26 parameter tensors (nn_mortality.0.weight, nn_mortality.0.bias, nn_mortality.3.weight, nn_mortality.3.bias, nn_mortality.6.weight, nn_mortality.6.bias, nn_growth.0.weight, nn_growth.0.bias, nn_growth.3.weight, nn_growth.3.bias, nn_growth.6.weight, nn_growth.6.bias, nn_regeneration.0.weight, nn_regeneration.0.bias, nn_regeneration.3.weight, nn_regeneration.3.bias, nn_regeneration.6.weight, nn_regeneration.6.bias, par_mortality_unconstrained, par_growth_unconstrained, reg_logK, reg_adult_logx, reg_adult_logA, par_regeneration_unconstrained, par_theta_recruits_raw, par_competition_unconstrained); sha256 d40ce38f94bc307fd17fc172cf12d8c29f5b95415eb08e711d1f1a2838f2a9a9. build_model() reproduces the fitted simulation (max |diff| 0).

Validation report: `Rscript calibration/conus/03_validation.R`.
