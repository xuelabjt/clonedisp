# CloneDisp

The method is named **Effective Clone Dispersion (ECD)**. **CloneDisp** is the
reference implementation. It calculates two linked metrics from cell-level
sample, state, and clone labels:

- **effective clone number**, `K_eff`: the inverse-Simpson, or Hill number
  with `q=2`, representation of the clones used by a cell state;
- **Clonal Dispersion Index**, `CDI`: the effective clone number divided by
  the number of valid clones available in that patient.

The method is deliberately **clone-agnostic**: it uses clone identity but not
the expression or CNV similarity between clones. It answers whether cells of
a state are concentrated in a few clones or spread across the valid clones.

## Definition

For patient `p` and state `s`, let:

- `c` index clones containing at least `min_clone_cells` cells;
- `K_valid` be the number of such valid clones in patient `p`;
- `n_c,s` be the number of cells in state `s` assigned to clone `c`;
- `q_c = n_c,s / sum_c(n_c,s)`.

Then:

```text
K_eff = 1 / sum(q_c^2)
CDI   = K_eff / K_valid
```

The state is evaluated only when it contains at least `min_state_cells` cells
within the valid clones.

`CDI=1` means the state is evenly dispersed across every valid clone. The
theoretical lower bound for a patient is `1/K_valid`, reached when all eligible
state cells occupy one clone. It is not generally 0; a two-clone patient has a
lower bound of 0.5.

## Install

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e ".[dev]"
pytest
```

## Command line

Input is a table with one row per cell:

```csv
sample,state,clone
P1,e0,C1
P1,e0,C2
P1,e1,C1
```

Run:

```bash
clonedisp \
  --input examples/example_labels.csv \
  --sample sample \
  --state state \
  --clone clone \
  --min-clone-cells 5 \
  --min-state-cells 20 \
  --output-prefix cdi_results
```

The command writes:

- `cdi_results_metrics.csv`: patient-state metrics;
- `cdi_results_clone_validity.csv`: clone filtering audit;
- `cdi_results_state_summary.csv`: patient-level state summary with bootstrap
  confidence intervals;
- `cdi_results_skipped.csv`: states omitted because a threshold was not met;
- `cdi_results_run.json`: formulas, parameters, version, and input checksum.

Use `clonedisp --help` for all options. By default, a CSV is comma-delimited
and a `.tsv` or `.txt` file is tab-delimited.

## Python API

```python
import pandas as pd
from clonedisp import analyze_clone_dispersion

labels = pd.read_csv("labels.csv")
result = analyze_clone_dispersion(
    labels,
    sample_col="sample",
    state_col="state",
    clone_col="clone",
    min_clone_cells=5,
    min_state_cells=20,
)

print(result.metrics)
print(result.clone_validity)
print(result.state_summary)
```

If a unique cell identifier is available, pass `cell_id_col="cell_id"`.
Duplicate identifiers are then rejected; otherwise each input row is treated
as one cell.

## Statistical notes

- The unit of analysis remains the patient-state pair. Do not treat cells as
  independent replicates.
- `K_valid` and `K_eff` depend on clone definitions and resolution. Report
  sensitivity across reasonable clone definitions where possible.
- The state-level bootstrap interval resamples patients, not cells. Wide and
  overlapping intervals are possible with a modest number of patients.
- CloneDisp measures distribution across clone identities. It does not measure
  phylogenetic distance between clones and is not a test of transcriptional
  plasticity, causality, or CNV neutrality.

## Citation

See [CITATION.cff](CITATION.cff). The software is released under the MIT
License.
