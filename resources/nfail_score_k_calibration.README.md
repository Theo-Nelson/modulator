# NFail-SCORE k-ratio calibration

`nfail_score_k_calibration.tsv` gives calibrated **k** values for the site filter
`aggregation.zn.nfail_score_k` (config), applied in `aggregate_by_gene.row_pass_filter`.

## What the filter does

modkit computes a site's stoichiometry only from reads it confidently assigns to a canonical or a
modified base; reads where the modification basecaller cannot confidently decide either are counted
as **NFail**. Error-prone (false-positive) sites carry a large NFail; true modification sites carry a
small one. **NFail-SCORE** (Nelson et al., *NFail-SCORE: Spurious-Call Omission via Read-failure
Evidence*) scores each site by

```
k = Nmod / (NFail + 1)
```

and keeps sites with `k >= nfail_score_k`. Because it is computed from the reads modkit already
discards, it needs **no matched modification-free control** once k is calibrated for a given
modification basecaller + version.

## Columns

- `Model` — basecaller model (HAC/SUP) and modification (e.g. `SUP | m6A DRACH`).
- `Dorado_version` — Dorado release the calibration is for.
- `A [P(A)>]`, `B [P(mod)>]` — the modkit canonical / modified probability thresholds used.
- **`k (ratio)`** — the calibrated NFail-SCORE k to set as `aggregation.zn.nfail_score_k`.
- `Max_WT_HEK293T_sites_recovered`, `N_FP0_parameter_conditions`, `N_conditions_at_max`,
  `FP_stoichiometry` — supporting benchmark numbers.

## Usage

Look up your `(model, Dorado_version, modification)` row and set the k in the config, e.g. for
Dorado **SUP m6A-DRACH v2.0.0**:

```yaml
aggregation:
  zn:
    nfail_score_k: 0.4
```

`k = 0.0` disables the filter (keeps every site); the shipped **default is `k = 1.0`** (config.yaml
`filters.nfail_score_k`), which on a real bedMethyl subset keeps ~650 of ~223,130 rows — do not assume
the default table is unfiltered. Ranges (e.g. `1.8–2.2`) mean several conditions tied at the
maximum; pick within the range for your P(A)/P(mod) thresholds. The score is easily recalibrated for
each new basecaller as it is released.
