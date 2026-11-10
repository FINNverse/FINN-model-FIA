# v0.2 pilots: holdout comparison

Generated 2026-11-10 10:45 by `calibration/v02/compare_pilots.R` from 7 pilots. Oregon pilots: full training set (585-590 plots), holdout 118 plots, 300 epochs. CONUS pilots: 300 training plots / 200 holdout plots (fixed seeds), 250 epochs. One thread. `_persist` = the first inventory carried forward; `delta` = change against it; ratios = predicted / observed totals; bare ground = basal area (m2/ha) at the end of a 300-400 year run at three climate points and the largest drift over its last fifth.

## States against persistence

| pilot | ba_r2 | ba_r2_persist | trees_r2 | trees_r2_persist | ba_tot_r2 | ba_tot_r2_persist | trees_tot_r2 | trees_tot_r2_persist | dba_tot_r2 | dtrees_tot_r2 |
|---|---|---|---|---|---|---|---|---|---|---|
| conus_C0_v01 | 0.865 | 0.869 |  0.822 | 0.883 | 0.683 | 0.749 |  0.185 | 0.704 | -0.2680 |  -1.8400 |
| conus_C1_v2_cap_theta_cosine | 0.882 | 0.869 |  0.900 | 0.883 | 0.778 | 0.749 |  0.741 | 0.704 |  0.1150 |   0.1010 |
| conus_C2_v2_staged | 0.875 | 0.869 |  0.903 | 0.883 | 0.773 | 0.749 |  0.759 | 0.704 |  0.0926 |   0.1630 |
| oregon_O0_v01 | 0.855 | 0.864 |  0.473 | 0.765 | 0.587 | 0.677 | -2.590 | 0.552 | -0.4800 |  -7.1300 |
| oregon_O1_cap_cosine | 0.847 | 0.864 |  0.569 | 0.765 | 0.640 | 0.677 | -1.480 | 0.552 | -0.2910 |  -4.6100 |
| oregon_O2_v2_cap_theta | 0.800 | 0.864 | -0.683 | 0.765 | 0.453 | 0.677 | -5.860 | 0.552 | -0.9600 | -14.5000 |
| oregon_O3_v2_staged | 0.835 | 0.864 |  0.783 | 0.765 | 0.652 | 0.677 |  0.563 | 0.552 | -0.2480 |   0.0106 |

## Rates

| pilot | growth_r2 | growth_r2_tree | mort_auc | mort_auc_tree | reg_spearman | growth_ratio | deaths_ratio | reg_ratio | ba_ratio | trees_ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| conus_C0_v01 |  0.1010 | -0.203 | 0.733 | 0.732 | 0.2320 | 1.010 | 0.745 | 1.030 | 1.160 | 1.450 |
| conus_C1_v2_cap_theta_cosine | -0.0647 | -0.318 | 0.747 | 0.748 | 0.2600 | 0.984 | 0.750 | 0.679 | 1.030 | 1.020 |
| conus_C2_v2_staged | -0.0842 | -0.412 | 0.741 | 0.746 | 0.2680 | 0.943 | 0.665 | 0.240 | 1.040 | 0.988 |
| oregon_O0_v01 |  0.2090 | NA | 0.667 | NA | 0.3180 | 0.959 | 1.020 | 4.900 | 1.190 | 2.080 |
| oregon_O1_cap_cosine |  0.2700 | NA | 0.694 | NA | 0.3290 | 0.952 | 0.997 | 4.020 | 1.160 | 1.900 |
| oregon_O2_v2_cap_theta |  0.1420 | NA | 0.687 | NA | 0.3080 | 0.931 | 0.920 | 6.520 | 1.220 | 2.330 |
| oregon_O3_v2_staged |  0.1450 | NA | 0.658 | NA | 0.0524 | 0.981 | 0.974 | 0.269 | 0.959 | 0.909 |

## Change, bare ground, cost

| pilot | delta_ba_obs | delta_ba_pred | delta_trees_obs | delta_trees_pred | bare_ba_min | bare_ba_max | bare_drift_max_pct | minutes | epochs | K_max | theta_median |
|---|---|---|---|---|---|---|---|---|---|---|---|
| conus_C0_v01 | -0.00525 | 0.04810 | -0.365 |  2.410 |  62.5 | 163.0 |  5.61 | 53.2 | 250 |   175.00 | 5.1000 |
| conus_C1_v2_cap_theta_cosine | -0.00525 | 0.00583 | -0.365 | -0.217 |  55.6 |  79.0 | 25.90 | 52.9 | 250 |     5.95 | 0.0609 |
| conus_C2_v2_staged | -0.00525 | 0.00649 | -0.365 | -0.440 |  36.2 |  49.7 | 18.30 | 50.5 | 250 |     5.95 | 0.0609 |
| oregon_O0_v01 |  0.06750 | 0.17500 |  0.328 |  7.190 |  99.2 | 114.0 |  5.65 |  8.5 | 300 | 45900.00 | 1.6300 |
| oregon_O1_cap_cosine |  0.06750 | 0.15400 |  0.328 |  6.040 | 110.0 | 141.0 |  5.00 | 11.0 | 300 |   238.00 | 3.0000 |
| oregon_O2_v2_cap_theta |  0.06750 | 0.19100 |  0.328 |  8.810 | 136.0 | 213.0 |  4.32 | 10.4 | 300 |   238.00 | 0.0853 |
| oregon_O3_v2_staged |  0.06750 | 0.04500 |  0.328 | -0.253 |  39.3 |  63.3 | 20.50 | 11.7 | 300 |   238.00 | 0.0853 |

## Final losses (baseline-scaled)

| pilot | loss_last |
|---|---|
| conus_C0_v01 | dbh 0.531 ba 0.083 trees 0.324 growth 0.749 mort 0.875 reg 0.757 |
| conus_C1_v2_cap_theta_cosine | dbh 0.495 ba 0.116 trees 0.248 growth 0.645 mort 0.714 reg 0.893 |
| conus_C2_v2_staged | dbh 0.000 ba 0.000 trees 0.000 growth 0.000 mort 0.000 reg 0.854 |
| oregon_O0_v01 | dbh 0.861 ba 0.285 trees 0.609 growth 0.370 mort 1.076 reg 1.466 |
| oregon_O1_cap_cosine | dbh 0.707 ba 0.107 trees 0.515 growth 0.316 mort 0.963 reg 1.608 |
| oregon_O2_v2_cap_theta | dbh 0.567 ba 0.151 trees 0.705 growth 0.299 mort 0.978 reg 1.048 |
| oregon_O3_v2_staged | dbh 0.000 ba 0.000 trees 0.000 growth 0.000 mort 0.000 reg 0.860 |
