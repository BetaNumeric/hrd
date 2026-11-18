# Catalog stars and empirical Hubble profiles

The page constructs its spatial and sorted views from individual star sprites. It does not load, reconstruct, fade, or cut up a scene photograph. All datasets use the same rendering path and retain each star's identity and depth while sorting.

## Gaia: nearby and deeper distant samples

**Settings → Star data** offers three Gaia DR3 datasets. **5,000 pc (Gaia + BJ)** is selected by default; it uses published geometric distance posteriors. The two original inverse-parallax views remain available as more precise local references. The **100–2,000 pc** sample includes the faint main sequence and white-dwarf region as well as giants and very luminous candidates.

| Menu option | Asset | Displayed stars | Distance selection | Extra apparent-magnitude cutoff |
| --- | --- | ---: | --- | --- |
| Within 100 pc (Gaia) | `gaia-nearby-catalog.json` | 13,023 | Parallax > 10 mas; distance < 100 pc | None |
| 100–2,000 pc (Gaia) | `gaia-distant-catalog.json` | 24,000 | 0.5 ≤ parallax ≤ 10 mas; distance 100–2,000 pc | **None** |
| 5,000 pc (Gaia + BJ) | `gaia-bailer-jones-catalog.json` | 24,000 | Geometric posterior median ≤ 5,000 pc; fractional interval half-width ≤ 20% | **None** |

The nearby sample is preserved exactly. The distant sample is rebuilt from deeper Gaia data; the earlier G < 9 limit has been removed. Of its displayed stars, **19,716 have absolute G > 4**, below the old artificial cutoff. The displayed distant sample reaches apparent G = **20.602** and absolute G = **15.455**. These are sample extrema, not new selection thresholds. The distance ranges are disjoint; a star exactly at 100 pc belongs to the distant view.

The nearby sample has **9,413 main-sequence-region, 3,500 white-dwarf-region, 43 giant-region, and 67 other stars**. The deeper distant sample has **16,063 main-sequence-region, 4,578 white-dwarf-region, 1,982 giant-region, 221 very-bright-region, and 1,156 other stars**. Region names describe approximate positions in the color–magnitude diagram, not confirmed evolutionary or spectral types. In particular, absolute G < −4 does not establish a supergiant classification. Empty highlight regions are disabled.

- Sources: **ESA/Gaia/DPAC Gaia DR3**, using only `gaiadr3.gaia_source`. [Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3), [source column definitions](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html), and [photometric quality prescription: Riello et al. (2021)](https://doi.org/10.1051/0004-6361/202039587). FLAME luminosities and GSP-Phot temperatures are not required.
- **Quality criteria for the two inverse-parallax samples:** parallax/error > 10, RUWE < 1.4, visibility periods > 8, G flux SNR > 50, both BP/RP flux SNR > 10, and −1 ≤ BP−RP ≤ 6. The builder also requires absolute G between −10 and 18 and |C*| ≤ 3σ using Riello et al. Table 2 and Eq. 18 (σ evaluated at G ≥ 4). Stars are not selected on stellar-model estimates or luminosity uncertainty intervals. The query template is [gaia-expanded.adql](../scripts/gaia-expanded.adql); exact queries are embedded in each asset. The 5 kpc sample has its own posterior-distance criterion, described below.
- **Deeper distant field:** query 100–2,000 pc stars with `random_index < 5000000`, without a G-magnitude cutoff. Gaia's published `random_index` is a random permutation of source indices, giving an all-sky random cohort without sorting the full table. Of 206,173 downloaded sources, 199,635 pass the additional color-excess and magnitude-range checks. These sources form the field candidate pool for the smooth sampling procedure below.
- **Distant enrichment pools:** query the white-dwarf region and very-bright region across the entire **100–2,000 pc** range with the same quality criteria and no apparent-magnitude cut. The white-dwarf query returns 120,968 sources, of which 114,519 pass the additional checks; the very-bright query returns 299, of which 228 pass. Also use all distant sources from the cached G < 9 query to supply the rare upper main sequence continuously, not just a hard very-bright bin. **G < 9 is an enrichment-pool condition, never a filter on the combined sample.** All source pools have the same astrometric and photometric quality criteria. Queries must return fewer than the 500,000-row ceiling; a build fails on possible truncation.
- **Smooth distant sampling:** deduplicate the field and enrichment pools and choose 24,000 sources using deterministic weighted sampling without replacement. Hard per-region quotas would create artificial density steps at their boundaries, so the importance varies smoothly with measured absolute G and BP−RP. Its exact expression and the source-ID hash procedure are stored in the asset and builder. The input-cohort inclusion is 1 for the fully queried G < 9 and white-dwarf-region pools and `5000000 / 1811709771` otherwise; this factor reduces enrichment-pool overrepresentation before applying the visual importance weights. Keep the smallest `-log(U) × inclusion / importance`, with U a midpoint from the first 52 bits of SHA-256(source_id). The Gaia source count is documented in the [Gaia DR3 overview](https://gea.esac.esa.int/archive/documentation/GDR3/). This is a teaching selection, not a reconstruction of population frequencies or a complete density sample. No source is moved in the diagram or assigned a synthetic measurement.
- **Nearby reproducibility:** retain the established nearby sample from its legacy parent union of the local quality-filtered query and a G < 9 query. The parent selection takes 16,000 field stars in SHA-256 order and adds the original region quotas before partitioning by distance. This preserves all 13,023 existing nearby source records and their appearance. Its queries and parent-selection metadata remain in the nearby asset. The deeper distant sample uses the separate smooth method above.
- These are teaching samples with extra sources from sparse diagram regions, **not complete censuses or uniform density maps**. Gaia detection and quality limits still exclude some faint or uncertain stars. Faint white dwarfs naturally become harder to observe at large distances; this is not replaced with synthetic stars. The nearby spherical boundary is its distance cut, and the distant view deliberately excludes the inner 100 pc.
- Distances: `1000 / parallax_mas`, with displayed formal uncertainty `1000 × parallax_error / parallax²`. At the selected SNR the formal relative error is below 10%. No parallax zero-point correction is applied; systematic errors remain. These are distance estimates.
- Positions: published ICRS RA/Dec at **J2016.0** and inverse-parallax distances are transformed with Astropy into solar-barycentric Galactic Cartesian coordinates. X points toward the Galactic center, Y toward Galactic longitude 90°, and Z toward Galactic north. Coordinates in parsecs are stored to six decimal places; source IDs remain strings to preserve their 64-bit precision. Display X = Galactic X, vertical = Galactic Z, and depth = Galactic Y. **One common scale on all three spatial axes** preserves 3D shape and relative distances. The nearby view uses a 100 pc scale and the distant view uses 2,000 pc; axes show the corresponding parsec ticks. Both spatial views mark the Sun at the origin.
- Sorted view: both samples use the **same BP−RP color and absolute G magnitude ranges**, so diagram positions, color mapping, and sprite exposure remain comparable. Smaller magnitudes appear higher. `M_G = G + 5 log10(parallax_mas) − 10`. Absolute G is distance-corrected brightness in Gaia's G passband, not bolometric luminosity in solar units. Depth retains physical Galactic Y on each sample's spatial scale. All six camera views use the corresponding axes. Switching datasets resets zoom and pan while preserving camera direction and the spatial/sorted layout.
- **No extinction or reddening correction** is applied. Dust can make distant stars redder and fainter, broadening the branches. Temperature and bolometric luminosity are not inferred or fabricated. Star details show the exact Gaia ID, color, absolute G, position, and estimated distance.
- Appearance: empirical Hubble profiles model point sources; they are not photographs of these Gaia stars. Tint follows measured color and exposure follows absolute G with a display stretch. Every star retains its profile, color, exposure, and physical depth throughout sorting. Sprite size does not depict stellar radius, and brightness is not an Earth-view inverse-square image.

Data credit: ESA/Gaia/DPAC, Gaia Collaboration et al. (2023), Gaia DR3. DPAC receives support from national institutions participating in the Gaia Multilateral Agreement.

All three Gaia field views include the Sun as a selectable reference star at the origin. It sorts with the catalog stars using the published solar estimates BP−RP = 0.816 and absolute G = 4.67. The ring follows it in both layouts. [Solar twin study](https://doi.org/10.1093/mnras/stag1864).

## Extended Gaia / Bailer-Jones view: within 5,000 pc

The new cumulative view uses Gaia DR3 positions and photometry, joined by exact source ID to the **geometric** distance catalog of [Bailer-Jones et al. (2021), AJ 161, 147](https://doi.org/10.3847/1538-3881/abd806). It includes the inner 2 kpc; it is not a shell. The boundary is a cut on the posterior median, so some stars' upper distance bounds extend past 5,000 pc. [Catalog documentation and license](https://dc.g-vo.org/tableinfo/gedr3dist.main).

Of the 24,000 displayed stars, **6,515 have median distances beyond 2,000 pc**. The diagram regions contain 14,883 main-sequence-region, 3,793 white-dwarf-region, 2,804 giant-region, 598 very-bright-region, and 1,922 other sources. These labels describe diagram locations rather than confirmed classifications. The complete downloaded pools contain 94,299 field, 392,341 bright, and 226,445 white-dwarf-region records before the additional corrected-excess and magnitude-range checks. As an independent completeness check, the bright pool's 392,341 source IDs matched the GAVO Gaia DR3 mirror exactly; only ESA archive records are used in the asset.

- **Distance method:** use the published `r_med_geo`, `r_lo_geo`, and `r_hi_geo`, in parsecs. These are the median and 16th/84th percentiles of the geometric posterior. The method uses EDR3 parallaxes with their published correction and a direction-dependent spatial prior. It does not use the photogeometric color–magnitude prior. No extra parallax-SNR or positive-parallax cut is imposed by this builder. These remain inferred distances, with increasing uncertainty at depth.
- **Quality selection:** positive ordered percentiles, median ≤ 5,000 pc, and `(r_hi_geo − r_lo_geo) / (2 × r_med_geo) ≤ 0.2`. This is a posterior interval criterion, not a claim that all parallax measurements have 20% errors. Gaia RUWE, visibility-period, flux-SNR, color, corrected-excess, and diagram-range criteria otherwise match the original views; the old parallax/error > 10 condition does not apply here. The precise ADQL is in [gaia-bailer-jones.adql](../scripts/gaia-bailer-jones.adql) and the generated asset.
- **Source pools:** an all-sky field cohort with `random_index < 1000000`, a full G < 10 bright enrichment pool, and a full white-dwarf-region enrichment pool, each with the same distance-posterior and photometry quality criteria. The combined selection has no additional apparent-magnitude cutoff. Smooth importance sampling without replacement uses the same functional weights as the 2 kpc view; candidate inclusion is 1 for the complete bright/white-dwarf pools and `1000000 / 1811709771` otherwise. There are no hard display-group quotas. Exact pool sizes, retained region counts, and sampling expressions are recorded in the JSON asset.
- **Complete query results:** the builder uses asynchronous Gaia jobs over disjoint source-ID sky tiles covering `[0, 12 × 2^59)`. Only jobs reaching `COMPLETED` are read. Tiles that exceed the archive execution limit are bisected and retried; every leaf result must be below its row ceiling. This matters because the synchronous CSV endpoint can return a plausible partial table when a long join times out. Partial responses are never used as finished query pools. The `queries` metadata gives the untiled equivalent; `queryExecution` records the tiling procedure.
- **Coordinates and diagram:** RA/Dec at J2016.0 plus the geometric median distance are transformed into Galactic XYZ with Astropy. One common scale maps 5,000 pc to the spatial plot extent. Absolute G is `G − 5 log10(distance_pc) + 5`. Observed BP−RP and absolute G use the same plot ranges, color mapping, profiles, and exposure as the other Gaia samples. Dust/extinction is not corrected. Sorting changes only color and brightness coordinates; depth keeps estimated physical Galactic Y.
- **Selected-star uncertainty:** details show the exact Gaia source ID, median distance, and asymmetric percentile bounds. The orange line follows the star's line of sight between the two published distances. Sorting maps each bound through the same absolute-G formula while retaining its associated depth. In Front/Back it becomes the distance-only vertical magnitude interval. These endpoints exclude photometric and extinction uncertainties; they are not a full stellar-property error budget. Other datasets retain their existing formal inverse-parallax distance display.
- **Interpretation:** the sample emphasizes sparse diagram regions and is not a census, population-frequency estimate, or uniform spatial-density map. Faint nearby stars, luminous distant stars, dust, Gaia detection limits, quality cuts, and the geometric prior all affect its appearance. Very-bright-region candidates are not automatically classified as supergiants.

Distance-catalog credit: **C. A. L. Bailer-Jones, J. Rybizki, M. Fouesneau, M. Demleitner, and R. Andrae (2021)**, licensed **[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)**. Changes: quality/distance selection, sampling, DR3 photometry crossmatch, Cartesian transformation, and uncertainty visualization as described above. Astrometry/photometry credit: ESA/Gaia/DPAC. The published geometric distance intervals are preserved unchanged.

## Hyades: real cluster membership and real 3D structure

**Settings → Star data → Hyades cluster** selects `hyades-catalog.json`. Membership comes from the [GCNS open-cluster membership table](https://dc.zah.uni-heidelberg.de/tableinfo/gcns.hyacob), associated with [Gaia Collaboration (2021), A&A 649, A6](https://doi.org/10.1051/0004-6361/202039498). This is a published list of **probable** members, not a sky-cone selection or a claim of certain membership.

The table contains 920 Hyades candidates. We select the 630 whose **published** `dist_cluster` is at most 20 pc, then match their Gaia source IDs to DR3. Positive parallax with signal-to-noise > 50, finite G magnitude and BP−RP color, and −0.5 ≤ BP−RP ≤ 5.5 leave **607 displayed candidates**. Membership statistics and published center distances are preserved in the asset. We do not reject binaries by RUWE, select on luminosity uncertainty bounds, or invent missing stellar properties.

Only **189** displayed stars have published positive FLAME luminosity estimates. To include the other members, the whole Hyades diagram uses **absolute Gaia G magnitude**:

```text
M_G = G + 5 log10(parallax_mas) − 10
```

This is distance-corrected G-band brightness, not bolometric solar luminosity. No extinction correction is applied. Star details show published luminosity and temperature when available, and “Not published” otherwise. Inconsistent luminosity intervals remain flagged independently of membership and plotting. Tint follows BP−RP; sprite brightness uses absolute G with a visual stretch.

3D positions use the same inverse-parallax distance and Astropy Galactic transform as the Gaia distance samples. The display subtracts one common center (the mean XYZ of the displayed sample), with one common parsec scale on all axes. This is a rigid translation and uniform scaling: **all relative positions, pairwise distances, and orientations are preserved**. It is not a dynamical reconstruction. The published 20 pc membership cut uses the GCNS center; the display's sample center is not claimed to be the cluster's exact barycenter. “Cluster X/Y/Z” axes are offsets from this display center, aligned with the Galactic axes. The raw solar-barycentric coordinates remain in star details. Sorted depth retains the same cluster-relative Galactic Y coordinate.

Membership-table credit: Gaia Collaboration / GCNS, Richard Smart (2023), licensed **[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)**. Changes: Hyades-only inner-region and astrometric selection, DR3 crossmatch, coordinate transformation, and display sampling described above. Astrometry/photometry credit: ESA/Gaia/DPAC. The source table also lists Coma Berenices candidates; they are not included here.

## Legacy reference data (not in the dataset menu)

The visible selector contains the three Gaia distance samples and the Hyades. The earlier illustrative and Omega Centauri models remain in the source as internal references.

### Omega Centauri: observed photometry, simulated positions

`omega-centauri-catalog.json` contains 10,000 measured stars selected from the central-field proper-motion catalog of **Anderson & van der Marel (2010), ApJ 710, 1032**, distributed by CDS/VizieR as `J/ApJ/710/1032/pm_stars`.

- [Catalog and documentation](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/710/1032)
- [Persistent dataset identifier](https://doi.org/10.26093/cds/vizier.17101032)
- Records retain published reference-frame pixel positions and calibrated VEGAMAG F435W/F625W magnitudes at source precision. “Catalog star 0001” is a local sample index; the reference coordinates identify the source record.

The downloaded central field has 108,507 records. Valid display candidates have 16 ≤ F435W ≤ 25.5, 15.5 ≤ F625W ≤ 24, and −0.5 ≤ F435W−F625W ≤ 2.9. These cuts discard unmeasured-band sentinels and extreme color outliers. We retain all candidates brighter than F625W = 17 or bluer than F435W−F625W = 0.5, then select remaining records by a deterministic hash of their reference positions, to a total of 10,000. This favors visibility of sparse blue/bright sources; the sample is not a complete or statistically representative stellar census. The underlying proper-motion catalog also has selection limits, especially for saturated bright stars.

**Cluster coordinates are simulated**, using a fixed random seed and a truncated, centrally concentrated distribution. Published sky positions are retained only for identification in this first version. The sorted front view plots measured color index horizontally and F625W apparent magnitude vertically, with brighter (smaller) magnitudes higher. It does not infer temperature, bolometric luminosity, evolutionary type, membership, or line-of-sight distance. Side and top views preserve simulated depth.

## Star images

`hubble-acs-psf.json` is derived from the **STScI ACS/WFC SM3 empirical STDPSF library**, using its F435W and F625W files. These are measured point-spread-function models: representations of how unresolved stellar light appears through Hubble's optics, rather than a unique photograph of each catalog star.

- [STScI ACS point-spread-function library and documentation](https://www.stsci.edu/hst/instrumentation/acs/data-analysis/psf)
- Four representative detector-position profiles (zero-based FITS cube planes 20, 33, 57, 69) retain all **101×101 native samples** of the 4×-supersampled source profiles. They are encoded logarithmically as interleaved 8-bit band values. This is a visual asset, not a calibration product for scientific photometry.
- Each star keeps the same profile, color, and brightness throughout the transition. The palette maps color index to a blue-to-red scale; it is not a calibrated true-color composite. Stars are unresolved, so visible size comes from the profile and exposure, not a depicted stellar radius.
- WebGL2 decodes the profiles to a half-float texture and interpolates linear light intensities. It accumulates linear stellar light in a floating-point image, then applies one display stretch. Ordinary displays use 2× supersampling, and high-DPI displays use up to their native 3× resolution, within a six-million-pixel budget (at least native CSS resolution). Star footprints grow slowly with zoom to keep neighboring stars distinct. A Canvas2D sprite fallback uses the same native profiles with high-quality image resampling when floating-point WebGL2 is unavailable; it approximates compositing and quantizes display levels.

Zoom has **no fixed minimum or maximum**; only nonpositive or nonfinite numerical states are rejected. The wheel zooms around the pointer, Shift-drag pans, and two-finger gestures zoom and pan around their midpoint. Buttons and +/− keys also zoom. Reset view returns to the centered Top view at 1×. Projection combines camera offsets before magnifying to improve precision at deep zoom. Sprite footprints stop growing at the native profile width, independently of camera zoom, to avoid giant blurred point sources. The selected profiles already use all of their native 101×101 spatial samples; their optical diffraction and soft wings remain. Upscaling them does not resolve stellar surfaces or create new telescope detail.

The original 4,200-star illustrative dataset remains only as an internal reference in the source. Its temperatures, luminosities, population labels, and positions are synthetic; it uses the same empirical profiles as the observed dataset.

## Installing the app

Open the diagram on its HTTPS site and use your browser’s **Install app** or **Add to Home Screen** option. The app opens in its own window and works offline after the first successful load. The star data, textures, and font are embedded in `index.html`.

For a local installation test, serve the project folder with `python -m http.server --bind 127.0.0.1 8000` and open `http://127.0.0.1:8000/`. You can also open `index.html` directly for the portable browser version.

## Rebuilding

The JSON assets are also embedded in `index.html`, so it opens directly from disk. The dataset menu offers the three Gaia distance samples and the Hyades cluster.

Install NumPy and Astropy, then run:

```powershell
python scripts/build_star_assets.py --download
python scripts/build_expanded_gaia_assets.py --download
python scripts/build_gaia_bj_assets.py --download
python scripts/build_hyades_assets.py --download
```

Sources are cached in `.star-data-cache`. Subsequent builds can omit `--download`, or supply `--cache-dir <directory>`. The builders regenerate their respective JSON assets, synchronize embedded HTML copies, and update `source-manifest.json` with public source URLs/query references and SHA-256 hashes. `build_gaia_assets.py` remains a compatibility entry point for the Gaia builder, which writes both distance samples. The Gaia builder requires Astropy; the Hyades and Hubble builders also use NumPy. Hyades membership and Gaia queries are preserved in its JSON asset.

`build_gaia_bj_assets.py` builds only the 5 kpc asset and its embedded copy. It imports the shared color-excess/region definitions from the existing Gaia builder and uses [gaia_archive.py](../scripts/gaia_archive.py) for completed, tiled downloads. `--download --refresh` rebuilds the combined CSV caches from completed per-tile caches; a normal offline rebuild reuses the combined CSV files. The existing nearby, 2 kpc, Hyades, and Hubble-profile assets are not regenerated by this builder.

To rebuild only the star profiles, use `python scripts/build_star_assets.py --psf-only --download`; this does not require or change the legacy Omega Centauri catalog.

The interaction was inspired by [NASA/Hubble's H–R animation](https://science.nasa.gov/asset/hubble/constructing-the-hertzsprung-russell-diagram-for-globular-star-cluster/). This implementation does not claim to reproduce that animation's undocumented production method.
