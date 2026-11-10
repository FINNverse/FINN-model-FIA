# calibration/oregon: fit summary

Fitted 2026-11-10 01:47 by `calibration/oregon/02_fit.R` in 18 min (2.2 s per epoch, 1 thread) on cpu (FINN 0.2.0, torch 0.15.1, R 4.5.0).

Data: `oregon-10yr-v1/train_blocked` (590 sites, 28 species incl. other, 2 census steps of 10 years). Holdout for validation: `holdout_blocked`.

Processes: FINN::growth, FINN::mortality, FINN::regeneration_saturation (K init 500), FINN::competition; climate `~temp + prec` linear per species (hidden = NULL); reg_floor 0.01; recruits_dbh 12.9; 4 patches of 0.01681 ha.

Fit: 500 epochs, batch 50, lr 0.01, weights "auto" with weight_scale "auto" (baseline-scaled losses), loss families dbh = mse, ba = mse, trees = poisson, growth = mse, mortality = binomial, regeneration = nbinom, env_autoscale TRUE, seed 42.

Final losses (baseline-scaled, 1 = no better than the mean): dbh 0.972, ba 0.165, trees 0.376, growth 0.603, mort 1.023, reg 0.462.

Fitted recruitment caps K (stems per ha per decade): min 0.0, median 37.9, max 131271.7.

Environment scaling stored with the model: temp (center 6.872, scale 2.445); prec (center 1078.080, scale 673.151).

parameters.pt: 9 parameter tensors (nn_mortality.0.weight, nn_growth.0.weight, nn_regeneration.0.weight, par_mortality_unconstrained, par_growth_unconstrained, reg_logK, par_regeneration_unconstrained, par_theta_recruits_raw, par_competition_unconstrained); sha256 0d29d459a7ae10e692550c1ce64d9a94d42d8f03b7a971108ae12c57af32d815. build_model() reproduces the fitted simulation (max |diff| 0).

Validation report: `Rscript calibration/oregon/03_validation.R`.
