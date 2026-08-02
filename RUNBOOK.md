# Runbook

End-to-end reproduction in ten commands or fewer.

> **Before running:** the package does not currently import - see the Status
> section in [README.md](README.md#status--what-runs-today). The sequence below
> is the intended flow once those imports are fixed.

## Complete run

```bash
# 1. Setup
conda env create -f environment.yml && conda activate syndrumnet
pip install -e .

# 2. Build all data (first time only, ~2-4 hours)
python scripts/build_all_data.py --config configs/default.yaml

# 3. Run the full pipeline (~1-2 hours)
python scripts/run_pipeline.py --config configs/default.yaml

# 4. Evaluate against known synergies
python scripts/evaluate.py --config configs/default.yaml

# 5. Generate all figures
python scripts/make_figures.py --config configs/default.yaml

# 6. View results
cat reports/tables/evaluation_summary.csv
ls reports/figures/
```

Or simply `make all`, which chains steps 2-5.

Optional:

```bash
pytest tests/ -v        # run the test suite
make appendix           # build Appendix_Code.pdf (requires pandoc + xelatex)
```

## Expected runtime

| Stage | Time |
|---|---|
| Data download | 1-2 hours (one-time) |
| Data preprocessing | 30-60 minutes (one-time) |
| Pipeline execution | 1-2 hours per run |
| Figure generation | 5-10 minutes |
| **First run, total** | **~3-5 hours** |
| **Subsequent runs** | **~1-2 hours** (data cached) |

## Expected outputs

| Path | Contents |
|---|---|
| `reports/tables/predictions_<disease>.csv` | One file per disease, all drug pair scores |
| `reports/tables/evaluation_summary.csv` | AUC-ROC / AUC-PR per disease |
| `reports/figures/*.png` | ~15-20 figures at 300 dpi |
| `logs/*.log` | Detailed execution logs |

## Troubleshooting

**Out of memory.** Reduce `n_cores` in the config, or process diseases one at a
time with repeated `--diseases` flags.

**Download failures.** Some sources need registration (PhosphoSitePlus). Retry
partial downloads with `python scripts/build_all_data.py --retry-failed`.

**Reproducibility.** Keep `random_seed` fixed in the config; data source versions
are logged to `data/raw/VERSIONS.txt` during the build.
