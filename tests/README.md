# Testing StrainQC-nf

Three levels, following the testing levels used for DTU's yAMP pipeline. The first two run in CI on every push
and pull request (`.github/workflows/ci.yml`) with no HPC, containers or tools.

| Level | What it checks | Where | Runtime |
|---|---|---|---|
| **Unit + contract** (pytest) | `bin/` parsers, the `results.json` contract the frontend consumes, the CNV summariser; golden tests on a real breseq report and real CNVkit output | `tests/python/` | < 1 s |
| **Pipeline** (nf-test) | `-preview` dry runs: parameter switches and fail-fast input validation. `-stub` runs: channel wiring, per-process task counts, output layout for 2–3 samples | `tests/main.nf.test` | ~2.5 min |
| **Full run vs published truth** | real tools on real data; recall against Ottilie et al. (2022) supplementary tables, CRISPR-validated drivers, chr I aneuploidy | `tests/ottilie_benchmark/` | ~30 min on HPC |

## Run locally

```bash
# unit + contract (any Python ≥ 3.10)
pip install pytest pytest-cov
pytest --cov

# pipeline (Java 17 + Nextflow ≥ 23.10 + nf-test on PATH)
nf-test test tests/main.nf.test
```

The benchmark has its own instructions in `tests/ottilie_benchmark/README.md`.

## Adding a test

- A parser or a `results.json` field changes → add or adjust a case in `tests/python/`. If a
  field is added, also add it to `TOP_LEVEL_ALLOWED` in `test_contract.py` and to
  `docs/API_CONTRACT.md`. The contract test fails on any undocumented key.
- A new parameter, route or process → add a `-preview` case for its start-up logic and extend
  a `-stub` case with the task count you expect.
- Check that a new test can fail: break the code it guards once and confirm it goes red.
  This is how the CNVkit reference broadcast test was confirmed.
