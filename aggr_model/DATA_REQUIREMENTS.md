# Inputs for `ssaggr`

Everything is set in one YAML config (see `configs/real_local.yaml`). Relative paths are
resolved from the config's folder. Data never go into git: `data/` and `runs/` are ignored.

## Required

| config key | file | format |
|---|---|---|
| `counts` | one pathogen-only count matrix per host | TSV, genes × libraries, first column gene id. Ids like `Ss_jgi\|Sclsc1\|8988\|SS1G_02403T0` are reduced to `SS1G_02403`. Library names must be `<hpi>_<rep>_<host>_<isolate>` (non-inoculated: `0_<rep>_<host>_NC`). |
| `phenotype_table` | isolate × crop mean sAUDPC | CSV, first column isolate, one column per crop (rank columns are ignored). CH4 Table 2 works as is. |
| `host_to_crop` | mapping | host code in the library names → crop column, e.g. `{Gm: Soybean, Ha: Sunflower}` |

## Strongly recommended

| config key | file | why |
|---|---|---|
| `library_totals` | CSV `host,sample,host_total_counts` | gives pathogen transcript share = pathogen / (pathogen + host) reads, the colonization covariate. Column sums of the `*_HOST_only` matrices are enough. |
| `annotation` | TSV with `SS1G`, `signal_peptide`, `is_effector`, `cazyme`, `cazyme_family`, … | labels the candidate determinants |

## Optional

| config key | file | use |
|---|---|---|
| `lesion_files` | per-plant lesion time courses, one per crop (`Subject, Ss_isolate, 1dpi, 2dpi, …`) | isolate × crop share of variance; sAUDPC = AUDPC / [D(n−1)/n] as in CH4 |
| `in_vitro` | any isolate-level table | each trait is tested as an assay (oxalic acid, appressoria, radial growth, later lesion pH, oxalate output and tissue buffering) |
| `isolate_aliases` | mapping | e.g. `Ss1980: "1980"`, `SsPotter: SSPotter` |
| `classes` | mapping | a-priori classes; `high` / `low` are used in the resolution table |

## Settings worth deciding on purpose

- `timepoints`: `["48"]` follows the proposal. Run `["24"]` and `["96"]` as sensitivity checks;
  do not pick whichever gives the best result.
- `target_scale`: `log` (centred log sAUDPC) or `rank`.
- `n_permutations`: 720 enumerates every ordering of 6 isolates. Use a random subset for more
  isolates.
- `thin_to_equal_depth`: binomial thinning of pathogen reads before normalisation.

## For the proposal's design (12 isolates × 5 hosts)

The same files, with more host codes in `counts` and `host_to_crop`, work unchanged. The
stability regression needs at least 3 hosts for its deviation term. Botrytis panels with many
sequenced isolates fit the same format; that is where the leave-one-isolate-out tests gain power.
