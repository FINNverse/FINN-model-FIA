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
| J | **FINN's simulated annual mortality is floored at 1 %** (model side, not data). | `R/finn-class.R:1346` at 0848b60: `m_draw = torch_clamp(m, 0.01, 0.99)` before the death draw. The per-tree mortality loss sees the unclamped rate; every state loss and every simulation sees >= 1 %/yr. Where it acts here: Oregon's observed mortality is 0.066 per decade = 0.68 %/yr overall and lower for Juniperus occidentalis, Calocedrus, Quercus garryana, Alnus (no holdout deaths at all), so the simulated decade cannot lose fewer than 9.6 % of any cohort; conus: Juniperus osteosperma, Pinus edulis and Quercus alba have observed annual rates of 0.3-0.6 %, below the floor. In the v0.1 Oregon fit the mortality loss never left its baseline (1.02): the state losses were met with a floored death rate while the per-tree term asked for less. In bare-ground runs the floor sets a lower bound on turnover, i.e. an upper bound on the equilibrium stem number, which the Oregon run exceeded anyway through recruitment (G). | recorded, not compensated (coordinator's instruction): no weight or data change targets it. Max's side decides whether the clamp becomes an argument (branch `zoo-v02-mortfloor`). The v0.2 reports state the floor next to the mortality ratios. |
| I | Recruitment is rare and bursty in the data. | observed recruitment cells: 95-96 % zero; q99 60 (Oregon) / 45 (conus) per ha per period, max 521 / 803; v0.1 predicts < 0.5 in 61 % (Oregon) / 89 % (conus) of the cells. | the NB with low theta is the right family; judge recruitment by the ratio and by Spearman, not by R2. |

Not found: an env-offset error (A), expanded sample sizes (`n_at_risk`, `growth_n` count tallied trees), a censored-harvest bias (C), an all-NA `period_length` or variable patch counts.

## 2. Design facts checked

- **Tally thresholds.** Subplot (7.3 m radius, 0.01681 ha): trees >= 5 in (12.7 cm). Microplot (2.1 m, 0.00135 ha, 8 % of the subplot area): saplings 1-5 in. Macroplot (17.95 m, 0.1012 ha): trees >= 24 in (61 cm) in the Pacific-Northwest design (CA/OR/WA; `MACRO_BREAKPOINT_DIA`). The models tally >= 12.7 cm; FINN's `recruits_dbh = 12.9` is ingrowth into that threshold. The preparation drops every sapling record, so the microplot never enters; the only consequence is the status rule of finding B.
- **Expansion.** `TPA_UNADJ` identifies the geometry (6.018 subplot, 74.97 microplot, 0.999 macroplot per acre); densities per patch come from the subplot area, never from TPA; `n_at_risk` and `growth_n` are tallied trees. Macroplot records are thinned 1 in 6 (deterministic per tree) rather than weighted, so counts stay integers.
- **Fates.** STATUSCD 1 alive, 2 dead, 3 cut (absent in the cohort); AGENTCD 80 (land use) drops the plot; `RECONCILECD` 1 = ingrowth, 3 = missed live tree, 6-8 = other reasons; v0.1 and v0.2 both keep code-3 trees as recruits (261 of 8003) because the alternative, dropping them, would remove real trees from the stand.
- **Mortality target.** Closed cohort of the trees alive at the interval start (`n_at_risk`), deaths of any natural cause including fire, insects, weather; disturbed intervals are in the data (FINN-fia dropped them in its later paper fits; the tool keeps them, so the mortality it learns includes background disturbance).
- **Interval and env.** Annual conus: env rows 1..23 per plot, census at the simulated year of the plot's own interval; Oregon: one env row per 10-year step. Verified per site (A).

## 2b. Species set and site rule (agreed rare-species rule, 2026-11-10)

No lumped `other` class. A plot is kept when the species outside the modelled set S
hold < x % of its basal area at the START census of every interval (never the end
census, which would select on the outcome); the end-census share is recorded per plot
and the validation is stratified by it (bins < 1 %, 1-5 %, > 5 %). Trees outside S on
kept plots are dropped (option (a)). Option (b), one fixed shade-only out-of-S class,
is NOT available through the FINN interface: `optimizeSpecies` / `optimizeEnv` act on a
whole process, not on one species row, there is no per-species gradient mask, and a
class with NA observations still recruits through the per-species floor (no
`regen_mask` in 0848b60); it would need a FINN change (proposal in the notes), so x is
kept at 2 %, where the competition error of the dropped trees is small.

Trade-off (`calibration/v02/species_tradeoff.R`, start census only; blocks and climate
cells as in the v1 preparations; pool = plots with four consistently sampled subplots):

| variant | x | candidate S | plots kept (share of pool) | plots in held-out blocks | held-out blocks populated | climate cells covered | tavg / prec 5-95 % span | weakest species: plots / trees / held-out blocks / span | species below the coverage thresholds |
|---|---|---|---|---|---|---|---|---|---|
| conus | 1 % | top 20 (v0.1 set) | 6 666 (22 %) | 1 610 | 33 / 54 | 0.39 | 0.63 / 0.66 | 23 / 74 / 2 / 0.25 | 17 of 20 |
| conus | 2 % | top 20 | 6 982 (23 %) | 1 655 | 35 / 54 | 0.41 | 0.63 / 0.69 | 35 / 107 / 2 / 0.27 | 16 |
| conus | 5 % | top 20 | 7 748 (26 %) | 1 813 | 36 / 54 | 0.46 | 0.63 / 0.73 | 60 / 173 / 4 / 0.29 | 14 |
| conus | 1 % | top 40 | 12 596 (42 %) | 2 740 | 42 / 54 | 0.57 | 0.95 / 0.76 | 39 / 170 / 1 / 0.22 | 33 of 40 (strict thresholds) |
| conus | 2 % | top 40 | 13 193 (44 %) | 2 805 | 43 / 54 | 0.59 | 0.95 / 0.77 | 50 / 219 / 1 / 0.22 | 32 |
| conus | 2 % | top 65 | 18 048 (61 %) | 3 755 | 46 / 54 | 0.71 | 0.92 / 0.81 | 135 / 248 / 0 / 0.22 | 49 of 65 |
| **conus** | **1 %** | **S = 36 (chosen)** | 12 378 (41 %) | 2 694 | 37 / 54 | 0.54 | 0.95 / 0.75 | 116 / 287 / 1 / 0.22 | 2 |
| **conus** | **2 %** | **S = 36 (chosen)** | **12 950 (43 %)** | **2 756** | **39 / 54** | **0.56** | **0.93 / 0.76** | **146 / 392 / 1 / 0.22** | **1** |
| conus | 3 % | S = 36 | 13 436 (45 %) | 2 818 | 39 / 54 | 0.58 | 0.91 / 0.77 | 177 / 478 / 1 / 0.22 | 1 |
| conus | 5 % | S = 36 | 14 262 (48 %) | 2 945 | 39 / 54 | 0.61 | 0.90 / 0.79 | 240 / 661 / 1 / 0.22 | 0 |
| oregon | 2 % | all 27 export species | 694 (90 %) | 124 | 6 / 6 | 0.82 | 0.98 / 1.00 | 5 / 24 / 0 / 0.11 | 9 of 27 |
| oregon | 2 % | top 20 | 642 (83 %) | 122 | 6 / 6 | 0.77 | 0.94 / 1.00 | 15 / 47 / 0 / 0.27 | 3 |
| oregon | 2 % | top 16 | 590 (76 %) | 111 | 6 / 6 | 0.75 | 0.94 / 0.99 | 18 / 40 / 1 / 0.31 | 1 |
| **oregon** | **1 %** | **S = 18 (chosen)** | 601 (78 %) | 120 | 6 / 6 | 0.75 | 0.92 / 1.00 | 15 / 46 / 1 / 0.30 | 1 |
| **oregon** | **2 %** | **S = 18 (chosen)** | **620 (80 %)** | **122** | **6 / 6** | **0.75** | **0.93 / 1.00** | **15 / 47 / 1 / 0.30** | **1** |
| oregon | 5 % | S = 18 | 647 (83 %) | 124 | 6 / 6 | 0.75 | 0.94 / 1.00 | 16 / 47 / 1 / 0.31 | 1 |

Full table incl. x = 3 % and |S| = 10, 15, 25, 30, 50: `calibration/v02/species_tradeoff.csv`;
per-species counts: `species_identifiability_<variant>.csv`.

**Choice (a coverage choice, not a measured identifiability limit).** S = the species of
the candidate set that, on the plots kept at x = 2 %, have at least k plots, m trees at
the start census, presence in n held-out blocks and a temperature span of at least a
fraction of the pool's range: conus k = 100, m = 500, n = 1, span 0.2 applied to the
top-40 set -> 36 species (dropped: Ulmus americana 234 trees, Liquidambar styraciflua 52
plots, Quercus stellata 91 plots, Fraxinus pennsylvanica 50 plots); Oregon k = 10, m = 50,
n = 1, span 0.1 applied to the 27 export species -> 18 (dropped: nine species on 5-24
plots with no presence in a held-out block, among them Lithocarpus, Picea sitchensis,
Quercus garryana). Re-selecting the plots with exactly S leaves one weak species per
variant (conus: 392 start-census trees on 146 plots; Oregon: 15 plots / 47 trees), kept
because they sit in one held-out block and will get a holdout error of their own.

Why S = 36 and not the v0.1 20 for conus: the larger set doubles the plots that meet the
rule (12 950 vs 6 982 at 2 %), raises the populated held-out blocks from 35 to 39 of 54
and the climate cells covered from 0.41 to 0.56 (temperature span 0.63 -> 0.93: the
southern pines and oaks, the northern hardwoods and the Great Basin conifers enter), and
every one of the v0.1 species gets more plots (Quercus alba 655 -> 1 618 at |S| = 65).
Why x = 2 %: it is the smallest threshold at which the held-out block coverage is at its
plateau (37 -> 39 -> 39 -> 39 blocks for x = 1, 2, 3, 5 %) with the climate span
unchanged, and without option (b) the dropped out-of-S trees must stay a negligible
competitor (<= 2 % of basal area). Oregon: x makes little difference (601-647 plots);
2 % for the same reason. Population claim for the cards: conus 43 % of the four-subplot
cohort plots, Oregon 80 % of the export's plots meet the rule at the interval-start
censuses.

Data versions built under the rule: `oregon-10yr-v2` (561 plots after makeObsData and
the empty-start-plot filter: 450 train / 111 holdout; end-census out-of-S share of the
kept plots: 601 below 1 %, 21 at 1-5 %, 3 above 5 %) and `conus-annual-v2` via
`annual-v7-S36`: the rule at every interval-start census keeps 12 550 of the 29 847
four-subplot plots (end-census share: 11 511 below 1 %, 664 at 1-5 %, 375 above 5 %),
9 268 plots after makeObsData, the land-use and the span/climate filters (v0.1: 4 954),
173 two-degree blocks with 35 held out (seed 21; a new draw, the block set changed with
the plots), train 1 500 / holdout 1 000 plots (57 772 / 33 836 tree intervals). Both
`siteID_dt.csv` carry the shares (`share_out_first`, `share_out_start_max`,
`share_out_end`).

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
