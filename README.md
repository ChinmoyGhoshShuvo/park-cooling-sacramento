# Do Parks Cool the Streets Around Them? Matched-Pixel Evidence from Sacramento, 2016–2025

**Research project (individual work)** · Chinmoy Ghosh Shuvo · October 2026

> **Status:** unpublished. This page shows the research question, study area, data and method overview only.
> **Full results and code will be released after publication.**

## Research question

How far does a park's cooling reach into the neighbourhood around it, and which parks cool their surroundings most?

Comparing temperatures near and far from parks mixes the park's effect with everything else that differs between those places, such as height above the river, nearness to water and distance from downtown. This project uses **matching**: each pixel near a park is paired with a comparable pixel farther away, with the same terrain and the same distances to water, downtown and the Delta, imaged on the same day. The study asks:

1. How much cooler are neighbourhoods at different distances from a park edge than comparable neighbourhoods farther away?
2. At what distance does the cooling fade?
3. Does a park's size, shape, tree canopy or internal water best explain how much it cools its surroundings?

The design was first written for Chattogram, Bangladesh. Sacramento serves as the pilot: its summer skies are reliably clear, the terrain is flat, one Landsat tile covers the whole city, and California keeps a complete public inventory of parks.

## Study area

**Sacramento, California:** the 2020 Census urban area "Sacramento, CA", 1,224 km², analysed on the native 30 m Landsat grid (WGS 84 / UTM 10N).

## Data

| Dataset | Used for | Resolution / period |
|---|---|---|
| Landsat 8/9 Collection 2 Level-2 via the Microsoft Planetary Computer STAC API | Land surface temperature (LST), NDVI, cloud and water masks, smoke screening | 30 m; summer scenes, 2016–2025 |
| California Protected Areas Database (CPAD) 2026a, CNRA map service | Parks and other protected open space | Vector |
| OpenStreetMap (`leisure=park`, `leisure=garden`, `landuse=recreation_ground`) | Cross-check for parks missing from CPAD | Vector |
| Copernicus DEM GLO-30 | Elevation, slope, aspect | 30 m |
| JRC Global Surface Water | Permanent water and distance to it | 30 m |
| Microsoft US Building Footprints | Building cover around each pixel | Vector |
| GHSL GHS-BUILT-S 1975 (R2023A) | Negative-control outcome | 100 m |
| ESA WorldCover 2021 | Tree canopy and water inside each park | 10 m |
| Natural Earth 1:10m coastline | Distance to the Sacramento–San Joaquin Delta | Vector |
| US Census TIGER 2023 urban areas | Study area boundary | Vector |

## Method overview

- A fully scripted, reproducible pipeline in **R**.
- **Scene choice:** clear summer Landsat 8/9 scenes from the single tile that covers the city, with an aerosol screen that removes wildfire-smoke days the cloud mask misses.
- **Parks:** CPAD holdings merged into park units and screened by size and greenness; urban parks are those run by cities, counties, park districts and homeowners' associations, excluding golf courses and cemeteries.
- **Pixel pools:** built-up neighbourhood pixels in distance bands up to 500 m from each park edge, plus park interiors, thinned to limit spatial autocorrelation.
- **Matching:** Mahalanobis nearest-neighbour matching on terrain and distances to water, downtown and the Delta, exact on the imaging date, with propensity-score trimming. Present-day building density is deliberately left out of the main specification, because parks shape how densely the land around them is built.
- **Estimation:** the average treatment effect on the treated (ATT) for each distance band, with a doubly robust check and block-bootstrap confidence intervals.
- **Robustness checks**, with decision rules set before running: a negative control (1975 built-up area), placebo parks and Rosenbaum bounds, plus a comparison of two control groups.
- **Park-level analysis:** per-park estimates, park-size classes, and a regression on area, shape, canopy share and internal water.

## Tools

R 4.6 (terra, sf, MatchIt, cobalt, fixest, rbounds, rstac, exactextractr, osmdata, data.table, ggplot2) · Microsoft Planetary Computer STAC API · QGIS 4

## Contact

Chinmoy Ghosh Shuvo · Open to collaboration and knowledge sharing. Feel free to reach out on [LinkedIn](https://www.linkedin.com/in/chinmoyghosh034).

## References

USGS Landsat 8-9 Collection 2 Level-2 Science Product Guide · Cook et al. (2014) *Remote Sensing* 6:11244–11266 · Microsoft Planetary Computer (Source et al. 2022) · GreenInfo Network, California Protected Areas Database 2026a · Pekel et al. (2016) *Nature* 540:418–422 · Pesaresi & Politis (2023) GHS-BUILT-S R2023A, JRC · Zanaga et al. (2022) ESA WorldCover 10 m 2021 v200 · Ho, Imai, King & Stuart (2011) *Journal of Statistical Software* 42(8) · Rosenbaum (2002) *Observational Studies*, 2nd ed. · Austin & Small (2014) *Statistics in Medicine* 33(24):4306–4319 · Bowler et al. (2010) *Landscape and Urban Planning* 97(3):147–155.

*Full results and code will be released after publication.*
