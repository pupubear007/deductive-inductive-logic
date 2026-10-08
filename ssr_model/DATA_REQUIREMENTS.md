# Input data for `ssrforecast`

The model reads one **data directory** with up to four CSV tables plus the image and map files
they point to. Only `units.csv` and `observations.csv` are mandatory; every other input is
optional and the model trains on whatever subset a dataset has. Missing modalities are masked,
not imputed.

```
my_data/
  units.csv          one row per experimental unit
  observations.csv   disease assessments and assays (long format)
  timeseries.csv     environment readings (long format)
  images.csv         image index: growth chamber, greenhouse, field, drone, satellite
  images/ ...        the image files (PNG/JPG/TIF, or .npy for multispectral)
  maps/ ...          optional per-pixel disease maps (.npy or grayscale PNG, values 0-1)
```

Run the checker before training. It writes `my_data/data_report.md` with a prioritised
**DATA REQUEST** of what is missing:

```sh
python -m ssrforecast.validate --data my_data --config configs/default.yaml
```

`python -m ssrforecast.synthetic --out data/synthetic` writes a complete example in this format.

---

## 1. `units.csv` (required)

A *unit* is whatever one disease outcome belongs to: a pot or plant (growth chamber, greenhouse),
a plot, or a field-season.

| column | required | meaning |
|---|---|---|
| `unit_id` | yes | unique id |
| `experiment_type` | yes | `growth_chamber`, `greenhouse` or `field` |
| `season` | yes | year, or chamber run / greenhouse trial id. Test splits hold out whole seasons. |
| `site` | recommended | farm / location / chamber id |
| `cultivar` | recommended | soybean line or variety (e.g. Williams 82 and the check lines) |
| `isolate` | if inoculated | *S. sclerotiorum* isolate id |
| `inoculation_method` | recommended | `natural` (field), `cut_petiole`, `mycelial_plug`, `colonized_petal`, `wound`, `ascospore_spray` |
| `planting_date` | recommended | ISO date `YYYY-MM-DD` |
| `inoculation_date` | if inoculated | ISO date. Needed when observations are given in `dpi`. |
| `lat`, `lon` | field | decimal degrees |
| `station_id` | optional | links the unit to a weather station in `timeseries.csv` |
| any other column | optional | e.g. `row_spacing_cm`, `seeding_rate`, `fungicide`, `irrigation`, `prior_ssr_history`, canopy-architecture traits, RNA-seq marker scores. List them under `tabular.numeric` or `tabular.categorical` in the config to use them. |

`cut_petiole`, `mycelial_plug`, `colonized_petal` and `wound` are treated as an intervention
**do(court = 1)**: the infection court is created by the experimenter, so these units teach the
model the infection and lesion-growth steps, not the spore and petal steps.

## 2. `observations.csv` (required): disease assessment and assays

Long format, one row per reading.

| column | meaning |
|---|---|
| `unit_id` | which unit |
| `date` **or** `dpi` | calendar date, or days post inoculation |
| `assay` | assay name, which must appear under `assays:` in the config |
| `value` | the reading |

Assay names understood by the default config (rename or add your own in the YAML):

| `assay` | kind | mechanism level | value |
|---|---|---|---|
| `apothecia_scouting` | binary | apothecia | 1 if apothecia found on that visit, else 0 |
| `spore_trap_qpcr` | binary | spore | 1 if *S. sclerotiorum* DNA detected, else 0 (or give Cq and a threshold) |
| `petal_qpcr` | binary | court (petal colonization) | 1 / 0 |
| `stem_qpcr` | binary | infection | 1 / 0 |
| `incidence` | fraction | infection | proportion of plants with SSR lesions, 0-1 (use `scale: 100` if in percent) |
| `dsi` | continuous | (0 without infection) | disease severity index, 0-100 |
| `lesion_length_mm` | continuous | (0 without infection) | stem lesion length (growth chamber / greenhouse) |
| `yield_kg_ha` | continuous | none | yield (any unit; set `scale` accordingly) |

Each binary assay carries a **sensitivity `se` and specificity `sp`** for its level. The default
`1.0 / 1.0` is the paper's idealised, sound and complete assay. Put measured values in the config,
or set `learn: true` to estimate them. Scouting is the typical sound-but-incomplete assay
(`sp = 1`, `se < 1`): a negative visit is not proof of absence.

Repeated readings over time are expected and useful, for example lesion length at 3, 5, 7, 9, 11
and 14 dpi, or weekly scouting through R1-R3.

## 3. `timeseries.csv` (optional, strongly recommended): environment

| column | meaning |
|---|---|
| `unit_id` **or** `station_id` | whose reading |
| `date` **or** `dpi` | timestamp; sub-daily readings are averaged per day |
| `variable` | e.g. `air_temp`, `rh`, `rain`, `leaf_wetness`, `soil_moisture`, `canopy_temp`, `canopy_rh`, `wind` |
| `value` | the reading |

Only variables listed under `series_variables` in the config are used. Each sample sees the
`lookback_days` up to and including its issue day and **nothing after it**.

Canopy-level moisture (leaf wetness, canopy RH, soil moisture) identifies the microclimate gate.
With macro weather alone, that gate is only a proxy.

## 4. `images.csv` (optional): every kind of image

| column | meaning |
|---|---|
| `unit_id` | which unit |
| `date` **or** `dpi` | acquisition time |
| `modality` | `gc_rgb`, `greenhouse_rgb`, `field_rgb` (ground), `drone_rgb`, `drone_ms` (multispectral), `satellite` |
| `path` | path relative to the data directory |
| `target_map_path` | optional: co-registered disease map (for the map diffusion head) |

* RGB images: PNG / JPG / TIF of any size; they are resized to `size` in the config.
* Multispectral and satellite: `.npy` arrays shaped `[C, H, W]` or `[H, W, C]`. Set `channels` per
  modality in the config, with the band order fixed and documented.
* Disease maps: one channel, values 0-1 (e.g. incidence per grid cell from expert annotation or a
  segmentation model). The map head needs 50 or more of them to be worth training.
* Several images of one unit over time are fine. A sample uses the latest image of each modality
  taken on or before its issue day.

---

## What each component needs

| component | needs |
|---|---|
| gated chain head (levels and gates) | at least one binary or fraction assay with both positives and negatives (the checker requires ≥ 20 of each) |
| identifying *which* step failed | readings at more than one level. With only incidence, only the product of the steps is identified. |
| outcome regression | any continuous assay |
| vector diffusion (generative prior over outcomes) | continuous assays |
| map diffusion (spatial / dispersal structure) | ≥ 50 images with `target_map_path` |
| testing on an unseen season | ≥ 3 seasons (or chamber runs) per experiment type |

## Suggested request to collaborators (field data), in priority order

1. Field-season list with planting date, cultivar, location, management (row spacing, fungicide
   timing and product, irrigation) and SSR history.
2. End-of-season **incidence** and **DSI** (R6-R7), plus yield.
3. Daily weather for each field, or the station ID and distance.
4. Canopy sensors (leaf wetness, canopy RH, soil moisture) through R1-R3, where available.
5. **Apothecia scouting** and/or **spore-trap** results through flowering, where available.
6. Drone or satellite imagery around R3-R7, with any disease annotations.

From your own growth-chamber and greenhouse work: per-plant lesion length time-courses with
cultivar, isolate, inoculation method and date, chamber logs, plant images, and (as numeric
columns in `units.csv`) RNA-seq marker scores if you want to test them as predictors.
