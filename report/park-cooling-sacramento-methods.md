# Park spillover cooling in Sacramento: matched-pixel analysis (pilot)

Methods notes for the analysis in this repository. R 4.6.1, seed **20261006**. Scripts in `code/` run in numbered order, and each can be rerun on its own.

## Main result

Built neighbourhoods (≥ 2% building cover within 150 m), 790 urban parks. Each distance band is compared with matched pixels in the 400–500 m band around parks, on the same Landsat date. Effect on summer land surface temperature (LST), with 95% CIs from a 1 km block bootstrap (500 replicates):

| Distance from park edge | Matched difference | Doubly robust | Negative control |
|---|---|---|---|
| 0–50 m | −1.27 °C [−1.46, −1.04] | −1.26 | flagged |
| 50–100 m | −0.63 °C [−0.80, −0.46] | −0.67 | passed |
| 100–150 m | −0.39 °C [−0.55, −0.22] | −0.42 | passed |
| 150–200 m | −0.20 °C [−0.34, −0.06] | −0.23 | passed |
| 200–300 m | −0.03 °C [−0.14, 0.09] | −0.13 | passed |
| 300–400 m | −0.05 °C [−0.13, 0.03] | −0.05 | passed |
| Park interior (on-site) | −2.79 °C [−3.27, −2.30] | — | flagged |

All covariates are balanced (|SMD| < 0.1, variance ratio 0.5–2) in the 0–200 m bands. The main finding is spillover cooling of about 0.6 °C at 50–100 m, falling to about 0.2 °C at 150–200 m and to zero by about 200 m. These are lower bounds if the 400–500 m reference band is itself slightly cooled.

This design was added after the original specification showed a flat −0.65 °C offset out to 500 m (see below). It is better suited to this city, and the next step is to confirm it in a second city with the design fixed in advance.

## Original specification (reported alongside)

Controls: pixels > 800 m from any CPAD or OSM green space and from permanent water; same built-neighbourhood restriction. Matched difference [SE]:

| Distance | Exogenous covariates | + building density | Unrestricted, exogenous |
|---|---|---|---|
| 0–50 m | −1.63 [0.19] | −1.81 [0.19] | −1.91 [0.20] |
| 50–100 m | −1.04 [0.18] | −1.34 [0.16] | −1.32 [0.18] |
| 100–150 m | −0.74 [0.18] | −1.13 [0.16] | −1.05 [0.19] |
| 150–200 m | −0.65 [0.18] | −1.09 [0.17] | −0.99 [0.17] |
| 200–300 m | −0.73 [0.17] | −0.87 [0.17] | −0.90 [0.18] |
| 300–400 m | −0.68 [0.17] | −0.78 [0.16] | −0.81 [0.18] |
| 400–500 m | −0.63 [0.18] | −0.71 [0.17] | −0.78 [0.19] |
| Park interior | −3.21 [0.34] | −2.88 [0.32] | −3.49 [0.35] |

- **The effect levels off beyond about 150 m.** Combined with the main design, this indicates a level difference between park-adjacent and park-poor neighbourhoods, not cooling that reaches 500 m.
- **Distance to downtown and to the Delta was less closely balanced in the 50–500 m bands** (|SMD| up to 0.23). In the 200–500 m bands, 70–86% of treated observations were matched, because each date has fewer controls than treated pixels and matching is 1:1.
- **Adding present-day building density made the estimates larger.** Park-adjacent pixels in Sacramento are denser than the controls, so building density acts mainly as a confounder here (the controls are a different kind of land), not as a mediator.

## Sensitivity (decision rules written before running)

- **Negative control** (GHSL 1975 built-up as the outcome; pairs whose park dates from after 1975):
  - Passed at 50–500 m (|standardized difference| ≤ 0.08, all CIs include 0).
  - Flagged at 0–50 m (−3.45 pp, standardized difference −0.29) and inside parks (−6.15 pp, −0.51). GHSL's 100 m cells straddling the park edge are a likely cause; a finer built-up layer would resolve this.
- **Placebo parks:** null in every band, so the decision rule says continue.
  - 74 placebo parks placed inside control territory, smaller than the real ones so that they fit.
  - The median 95% CI half-width was 0.78 °C, so the check detects artefacts larger than about that size.
- **Rosenbaum bounds** (Wilcoxon, 1:1 pairs aggregated to treated locations):
  - Significance holds up to Γ = 2.9 at 0–50 m, 2.05 at 50–100 m, and about 1.5–1.6 at 100–500 m (4.3 for park interiors).
  - These are upper values, because neighbouring pairs are not independent.

## Heterogeneity (main design)

- **No park-size threshold was found.** Pooled cooling at 50–150 m is detectable in every size class, down to the 0.2–0.5 ha minimum: −0.74 °C [−1.18, −0.29] for 0.2–0.5 ha, and −0.43 to −0.66 °C in the larger classes. Pixels near small parks are often close to other green space as well (school fields, HOA land, street trees), which may add to the cooling measured around them.
- **Meta-regression** (652 parks, 50–150 m ATT):
  - Park tree canopy share is the clear predictor: −0.93 °C per unit share (SE 0.21, p < 0.001).
  - Log area (+0.05, p = 0.16), perimeter-to-area ratio (+4.75, p = 0.06) and internal water (−0.23, p = 0.27) showed no significant effect.
  - Adjusted R² = 0.05. Log area and perimeter-to-area ratio are correlated (r = −0.66).
- **Per-park cooling distance** (123 parks with ≥ 5 locations in every band): 50 m for 33 parks, 150 m for 19, 200 m for 10; 61 show no detectable effect at single-park scale.

## Data and design decisions

1. **Study area:** the Sacramento 2020 Census urban area (1,224 km²) as a pilot for a design first written for Chattogram. Grid: Landsat native 30 m, EPSG:32610.
2. **LST:** 20 Landsat 8/9 C2 L2 scenes, two per year 2016–2025 (25 Jun–20 Aug), path 44 / row 33, each ≥ 99.7% clear.
   - Scenes were also screened for smoke (≤ 10% high-aerosol pixels), which removed wildfire-smoke dates that the cloud mask had passed as clear.
3. **Parks:** CPAD 2026 holdings via the CNRA map service, dissolved to units, kept if ≥ 0.2 ha with summer median NDVI > 0.3.
   - Main: 790 urban parks (>50% owned by a city, county, park district, community services district or HOA; excluding golf courses, cemeteries and planned parks).
   - Secondary: all 850 green units.
   - Controls are kept away from all CPAD land and from 154 OSM-only green parks.
4. **Matching covariates:** elevation, slope, northness, eastness, distance to permanent water (JRC > 50%), distance to the CBD (State Capitol), and distance to the coast.
   - Natural Earth's coastline reaches into the Delta, so "coast" here is effectively distance to the estuary (32–97 km).
   - GHSL 1990/2000 were not used for matching because most parks predate them.
5. **Built-neighbourhood restriction:** ≥ 2% Microsoft building footprint cover within 150 m, fixed before any estimate; the unrestricted runs are reported too.
6. **Pools:** thinned to one pixel per 120 m.
   - Treated: subsampled to 40,000 locations.
   - Controls: 6,178 locations (8,851 unrestricted), since relatively little of Sacramento is 800 m from both green space and water.
7. **Matching:** coarsened exact matching retained only 3–5% of treated units, so every band used Mahalanobis nearest-neighbour matching: 1:1, no replacement, 0.2 SD caliper on the propensity logit, exact on scene date.
   - Common-support trimming removed 120 treated locations in the original specification and 14 in the main design.
8. **Bootstrap:** 1 km blocks resampled with the matching held fixed (cluster bootstrap of the matched sample).
9. **Negative control and canopy sources:** the negative-control outcome is GHSL 1975 built-up for post-1975 parks. Canopy share and water come from ESA WorldCover 2021 (10 m).

## Scope

- LST is ground-surface temperature at about 10:40 local time, not air temperature.
- Results describe built Sacramento neighbourhoods in summer; Chattogram's monsoon climate and dense cores are the next test.

## Files

- **Tables (`data/`):**
  - `results_table.csv` (all designs and bands, with negative-control results and Rosenbaum Γ)
  - `balance_table_headline.csv`, `balance_table.csv`
  - `sensitivity_*.csv`, `diagnostic_vs_400_500m_ring.csv`
  - `per_park_att.csv`, `cooling_beyond_50m_by_park_area.csv`, `heterogeneity_meta_regression.csv`
  - `park_screening_summary.csv`, `landsat_scene_screening.csv`, `landsat_scenes_selected.csv`
- **Figures (`images/`):**
  - `love-plot-headline.png`, `love-plot-original-spec.png`
  - `map-trimmed-common-support.png`
  - `distance-decay-curve.png`
  - `forest-plot-per-park.png`
  - `park-area-vs-cooling-distance.png`
