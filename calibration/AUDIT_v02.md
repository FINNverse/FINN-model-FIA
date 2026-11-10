# v0.2 audit: FIA data and the v0.1 fits (both variants)

Written 2026-11-10 from `calibration/v02/audit_v02.R` (tables in `calibration/v02/audit/`),
the v0.1 reports (`calibration/<variant>/report.md`) and the FIA design. Numbers are
holdout unless stated. The traps are those of the v0.2 brief (section 1).

## 1. Summary of the findings

| # | finding | evidence | consequence for v0.2 |
|---|---|---|---|
| A | Env rows and observation windows are aligned in both variants. | per site: env starts at year 1, one row per simulated year, first census year == first `period_length`, env covers the last census (590/590, 118/118, 1500/1500, 1000/1000). | nothing to change (trap 1 clear). |
| B | The v0.1 recruitment target is 92 % of the true ingrowth and contains re-recorded trees. | On the 2500 conus plots 8003 trees count as recruits (no `PREV_TRE_CN`); 720 more crossed 12.7 cm with a microplot record before and were called "alive" (never recruits); of the 8003, 261 carry RECONCILECD 3 ("missed live tree") and 761 no code. Oregon: the export has the same status rule. | data v2: a living tree is "new" at the first retained record of its chronology (`v02/prepare_conus_v2_annual.R`): train recruits +3.6 %, holdout -4.4 %, growth records +1115. Small; not the cause of the weak recruitment skill. |
| C | No harvest in the cohort. | 0 cut records (STATUSCD 3) and no land-use deaths on the model plots: the M0/M1 tiers select unmanaged plots, so the mortality target is natural + disturbance death. Selection, not censoring. | keep; state in the card that the models represent unmanaged-plot demography. |
| D | Macroplot thinning touches 1 % of the conus records (1336 of 142 778 alive records >= 12.7 cm, states OR/WA/CA) and 18 % of the Oregon trees. | audit D; Oregon: 5376 trees >= 61 cm, 1 in 6 kept. | keep (sample sizes count tallied trees; densities are right in expectation). The Oregon variant's big-tree demography rests on 899 trees. |
| E1 | **Mortality is under-predicted, most for big trees** (conus). | deaths predicted / observed 0.64; per-tree death probability ratio by size 0.72 (12.7-20 cm), 0.64, 0.51, 0.50, 0.53 (> 60 cm); Picea engelmannii 0.37, Pinus edulis 0.07, Acer rubrum 0.22, Juniperus 0.43. | the mortality loss competes with the state losses (BWI trap 9): stage the fit, mortality on its own loss. |
| E2 | **Growth of big trees is under-predicted** (conus). | per-tree growth predicted / observed by size 1.06, 1.06, 0.98, 0.78, 0.47 (> 60 cm); R2 0.22-0.25 below 60 cm, 0.06 above. | FINN's growth size decline is too steep at the top; staging and more epochs may help; a FINN-side question (growth kernel shape). |
| E3 | **States end 12 % above the observations** (conus) because too few trees die; the stand drift has the wrong sign. | basal area predicted / observed 1.12, stems 1.12; observed mean change per plot x census -0.030 m2 per patch (stands lose basal area on average over 4-12 years), predicted +0.010; change R2 ~ 0, sign agreement 0.60 (BA), 0.44 (stems); the loss of BA is largest on dense plots (> 20 m2/ha: observed -1.2 to -1.5 stems per patch, predicted -0.05 to -0.4). | same as E1: the model does not reproduce stand-level losses, which are the disturbance mortality FINN has no process for; a per-plot change metric must be in the v0.2 report. |
| E4 | **Oregon recruitment composition is wrong** although the total is only 1.8x. | recruits predicted / observed: `other` 27 % of the predicted recruits (1 % observed), Abies amabilis 5.5x, Abies concolor 4.6x, Abies lasiocarpa 2.4x; Douglas-fir 0.20x, western hemlock 0.21x, ponderosa 0.25x. Mortality ratio by species 0.25-1.95. | the joint fit traded the regeneration composition for the state fit (BWI: "drifted to fir+maple"): drop `other` (data v2), bound K, stage the regeneration fit on its own loss. |
| E5 | **Oregon stems grow instead of shrinking.** | observed mean change +0.33 stems per patch, predicted +2.30; stand delta-stems R2 -0.71; stems predicted / observed 1.31 after two decades. | recruitment cap and staging (E4). |
| F1 | **Neither fit had settled.** | Oregon: losses jumped in the last fifth (dbh 0.29 -> 0.97, ba 0.09 -> 0.17); mortality loss 1.02 = never below its intercept baseline. conus: the regeneration network's change per snapshot was still growing at the end (last/first 2.5-7.2), mortality nets 1.1-1.5, K 0.85. | cosine lr decay (brief trap 10); fewer trainable regeneration parameters (fixed K, frozen theta). |
| F2 | **Recruitment caps did not bind (Oregon).** | fitted K per ha per decade: `other` 131 272, Abies grandis 54 804, Abies concolor 26 022, Lithocarpus 7727, Quercus garryana 3639, Juniperus 2838, Abies procera 2746, Douglas-fir 2058; observed recruits per ha per decade: q99 60, max 521. conus: K max 450 per ha per year (Pinus contorta), median 32. | FINN's custom parameters can be fixed but not bounded (`createProcess(custom_parameters = list(reg_logK = list(init, trainable)))`): fix K per species at 2 x the 95th percentile of the observed recruits of that species (floor 20 per decade / 2 per year). |
| F3 | Recruitment dispersion is low but not collapsed. | theta per species: Oregon median 0.084 (min 0.056, max 2.8), conus median 0.061 (min 0.048, max 12). | freeze theta at the v0.1 values clipped to [0.05, 2] (brief trap 7) so a regeneration-only stage cannot flatten the NB loss. |
| G | **The `other` class takes the recruitment (Oregon).** | bare ground after 1000 years: 118-415 m2/ha, Tsuga heterophylla 61-73 % at the wet points, Thuja plicata 22-53 % at the dry ones, `other` 0-10 %; in the holdout `other` receives 27 % of the predicted recruits. | drop `other` (data v2, 27 species; 1 % of the observed basal area). |
| H | Interval spread (conus): 4-12 years, 79 % ten-year intervals (West), 12 % five-year (East). | 1836 train intervals: 4:11, 5:218, 6:70, 7:84, 8:1, 10:1449, 11:1, 12:2. Change R2 is negative in every class; the 10-year class carries the negative observed drift (-0.044 m2 per patch), the 5-7-year classes a positive one (+0.025 to +0.053). | nothing to change in the data; the per-site path handles it. Report change metrics per interval class. |
| I | Recruitment is rare and bursty in the data. | observed recruitment cells: 95-96 % zero; q99 60 (Oregon) / 45 (conus) per ha per period, max 521 / 803; v0.1 predicts < 0.5 in 61 % (Oregon) / 89 % (conus) of the cells. | the NB with low theta is the right family; judge recruitment by the ratio and by Spearman, not by R2. |

Not found: an env-offset error (A), expanded sample sizes (`n_at_risk`, `growth_n` count tallied trees), a censored-harvest bias (C), an all-NA `period_length` or variable patch counts.

## 2. Design facts checked

- **Tally thresholds.** Subplot (7.3 m radius, 0.01681 ha): trees >= 5 in (12.7 cm). Microplot (2.1 m, 0.00135 ha, 8 % of the subplot area): saplings 1-5 in. Macroplot (17.95 m, 0.1012 ha): trees >= 24 in (61 cm) in the Pacific-Northwest design (CA/OR/WA; `MACRO_BREAKPOINT_DIA`). The models tally >= 12.7 cm; FINN's `recruits_dbh = 12.9` is ingrowth into that threshold. The preparation drops every sapling record, so the microplot never enters; the only consequence is the status rule of finding B.
- **Expansion.** `TPA_UNADJ` identifies the geometry (6.018 subplot, 74.97 microplot, 0.999 macroplot per acre); densities per patch come from the subplot area, never from TPA; `n_at_risk` and `growth_n` are tallied trees. Macroplot records are thinned 1 in 6 (deterministic per tree) rather than weighted, so counts stay integers.
- **Fates.** STATUSCD 1 alive, 2 dead, 3 cut (absent in the cohort); AGENTCD 80 (land use) drops the plot; `RECONCILECD` 1 = ingrowth, 3 = missed live tree, 6-8 = other reasons; v0.1 and v0.2 both keep code-3 trees as recruits (261 of 8003) because the alternative, dropping them, would remove real trees from the stand.
- **Mortality target.** Closed cohort of the trees alive at the interval start (`n_at_risk`), deaths of any natural cause including fire, insects, weather; disturbed intervals are in the data (FINN-fia dropped them in its later paper fits; the tool keeps them, so the mortality it learns includes background disturbance).
- **Interval and env.** Annual conus: env rows 1..23 per plot, census at the simulated year of the plot's own interval; Oregon: one env row per 10-year step. Verified per site (A).

## 3. Candidate settings and the hypotheses they test

| pilot | data | change vs the previous pilot | hypothesis (audit finding) |
|---|---|---|---|
| O0 / C0 | v1 | v0.1 model at the pilot budget (300 / 250 epochs, pilot plots) | the comparison baseline |
| O1 | v1 | K fixed at the data-derived cap; cosine lr | F2 + F1: bounded recruitment and a settled fit stop the stems drift (E5) |
| O2 | v2 (no `other`) | + theta frozen at the v0.1 values [0.05, 2] | G + F3: the lumped class and a free dispersion are what the regeneration loss exploited |
| O3 | v2 | + staged: demography first with recruitment frozen, then recruitment alone | E4 (loss competition): regeneration fitted on its own likelihood gets the composition right |
| C1 | v2 (complete ingrowth) | K fixed cap; theta frozen; cosine | B + F1 + F2 + F3 |
| C2 | v2 | + staged as O3 | E1 + E3: mortality learned on its own loss, not traded against the states |

GLM warm start: not piloted. The per-species GLM warm start of FINN-bwi maps linear climate coefficients onto FINN's lazily created linear layers (`initEnv`); the conus variant has 20-20 networks and the Oregon variant two covariates and two censuses per plot, so the staged fit with the regeneration process trained alone (O3/C2) is the cheaper test of the same hypothesis. If O3/C2 still get the Oregon composition wrong, a GLM warm start for regeneration (intercept + climate per species from the observed recruits, frozen) is the next candidate.

## 4. Pilot results

Filled in by the pilot queue (`calibration/v02/pilots/<pilot>/summary.csv`, `metrics.csv`); see the
comparison table in `calibration/v02/PILOTS_v02.md` and the notes.
