# Inputs for `ssaggr`

Everything is set in one YAML config (see `configs/real_local.yaml`). Relative paths are
resolved from the config's folder. Data never go into git: `data/` and `runs/` are ignored.

## Required

| config key | file | format |
|---|---|---|
| `counts` | one pathogen-only count matrix per host | TSV, genes × libraries, first column gene id. Ids like `Ss_jgi\|Sclsc1\|8988\|SS1G_02403T0` are reduced to `SS1G_02403`. Library names must be `<hpi>_<rep>_<host>_<isolate>` (non-inoculated: `0_<rep>_<host>_NC`). |
| `phenotype_table` | mean sAUDPC per isolate and crop | long CSV `crop, isolate, mean, …` (CH4 `outputs_fig2/TableS1_long.csv`, the final Table S1), or a wide isolate × crop table |
| `host_to_crop` | mapping | host code in the library names → crop column, e.g. `{Gm: Soybean, Ha: Sunflower}` |

## Strongly recommended

| config key | file | why |
|---|---|---|
| `library_totals` | CSV `host,sample,host_total_counts` | gives pathogen transcript share = pathogen / (pathogen + host) reads, the colonization covariate. Column sums of the `*_HOST_only` matrices are enough. |
| `annotation` | TSV with `SS1G`, `signal_peptide`, `is_effector`, `cazyme`, `cazyme_family`, … | labels the candidate determinants |

## Optional

| config key | file | use |
|---|---|---|
| `per_plant_saudpc` | per-plant sAUDPC `crop, isolate, sAUDPC` (CH4 `outputs_fig2/sAUDPC_per_plant.csv`) | isolate × crop share of variance |
| `lesion_files` | raw lesion time courses, one per crop | fallback only: recomputes sAUDPC from every rated day, which can differ from the CH4 pipeline |
| `in_vitro` | any isolate-level table | each trait is tested as an assay (oxalic acid, appressoria, radial growth, later lesion pH, oxalate output and tissue buffering) |
| `isolate_aliases` | mapping | e.g. `Ss1980: "1980"`, `SsPotter: SSPotter` |
| `classes` | mapping | a-priori classes; `high` / `low` are used in the resolution table |

## Settings worth deciding on purpose

- `timepoint_sets`: each time point alone, then pooled (default). Report all of them; do not
  pick whichever gives the best result.
- `target_scale`: `z` (default, the proposal), `raw`, `rank` or `log`. `stability_scale`: `raw`
  (default) or `log`.
- `hypothesis_classes` / `stability_cutoffs`: your hypothesis, stated before the analysis.
- `n_permutations`: 720 enumerates every ordering of 6 isolates. Use a random subset for more
  isolates.
- `thin_to_equal_depth`: binomial thinning of pathogen reads before normalisation.

## For the proposal's design (12 isolates × 5 hosts)

The same files, with more host codes in `counts` and `host_to_crop`, work unchanged. The
stability regression needs at least 3 hosts for its deviation term. Botrytis panels with many
sequenced isolates fit the same format; that is where the leave-one-isolate-out tests gain power.
