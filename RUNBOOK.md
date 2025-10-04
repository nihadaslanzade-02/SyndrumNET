\## RUNBOOK



\### Complete End-to-End Reproduction (≤10 commands)



```bash

\# 1. Setup

conda env create -f environment.yml \&\& conda activate syndrumnet

pip install -e .



\# 2. Build all data (first time only, ~2-4 hours)

python scripts/build\_all\_data.py --config configs/default.yaml



\# 3. Run full pipeline (~1-2 hours)

python scripts/run\_pipeline.py --config configs/default.yaml



\# 4. Evaluate against known synergies

python scripts/evaluate.py --config configs/default.yaml



\# 5. Generate all figures

python scripts/make\_figures.py --config configs/default.yaml



\# 6. View results

cat reports/tables/evaluation\_summary.csv

ls reports/figures/



\# Optional: Run tests

pytest tests/ -v



\# Optional: Generate code appendix

make appendix  # Creates Appendix\_Code.pdf

```



\### Expected Runtime

\- Data download: 1-2 hours (one-time)

\- Data preprocessing: 30-60 minutes (one-time)  

\- Pipeline execution: 1-2 hours per run

\- Figure generation: 5-10 minutes

\- \*\*Total first run: ~3-5 hours\*\*

\- \*\*Subsequent runs: ~1-2 hours\*\* (data cached)



\### Expected Outputs

\- `reports/tables/predictions\_\*.csv`: One file per disease with all drug pair scores

\- `reports/tables/evaluation\_summary.csv`: AUC/PR metrics

\- `reports/figures/\*.png`: ~15-20 publication-quality figures

\- `logs/\*.log`: Detailed execution logs

"""

