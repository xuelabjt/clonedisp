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

CloneDisp also supports **multi-definition analysis**. Each clone definition is
analyzed independently, then patient-level CDI values are compared across
definitions. Cell labels from different resolutions are never pooled.

The input table is the downstream product of a CNV workflow, not the raw
inferCNV matrix. See
[Generate CloneDisp labels from inferCNV](docs/input_from_infercnv.md) for the
R inferCNV pipeline, a reusable label-export script, and common input errors.

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

For plotting without the development tools:

```bash
python -m pip install ".[plot]"
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

## Multi-definition analysis

Create a manifest with one row per clone definition:

```csv
definition,path
res 0.5,labels/res0p5.csv.gz
res 0.6,labels/res0p6.csv.gz
res 0.7,labels/res0p7.csv.gz
```

Then run:

```bash
clonedisp-multi \
  --manifest definitions.csv \
  --output-prefix cdi_multi \
  --plot \
  --focus e1 \
  --focus e2
```

The command computes CDI separately for every definition and writes:

- `*_definition_patient_metrics.csv`: patient-state values within each
  definition;
- `*_definition_summary.csv`: patient-level mean CDI per definition and state;
- `*_state_summary.csv`: median, range, and rank across definitions;
- `*_skipped.csv` and `*_run.json`;
- `*.pdf` and `*.png` when `--plot` is enabled.

For tables without a sample column, `--cell-id cell --sample-regex
"^(T\\d+)_"` can derive sample identifiers from cell names.

## Python API

```python
import pandas as pd
from clonedisp import analyze_clone_dispersion, analyze_multi_definition

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

by_definition = {
    "res 0.5": labels_res05,
    "res 0.7": labels_res07,
    "res 1.0": labels_res10,
}
multi = analyze_multi_definition(by_definition)
print(multi.state_summary)
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
