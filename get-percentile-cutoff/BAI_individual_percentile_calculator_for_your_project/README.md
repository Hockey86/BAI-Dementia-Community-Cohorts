# BAI individual percentile calculator

Converts each person's brain age index (BAI, in years) into a percentile among
similar healthy reference people from five community cohorts (ARIC, FHS, MESA,
MrOS, SOF).

BAI can be computed from sleep EEG (EDF format) using [Luna](https://zzz-luna.org/luna/ref/predict/).


## Usage

Requires Python 3 with `numpy` and `pandas`.

```bash
python compute_BAI_percentile.py example_input.csv output.csv
```

The input needs columns `BAI, age, sexM, race_Asian, race_Black, race_White, BMI, AHI`.
People are matched on sex, race, age, BMI and AHI. See `example_input.csv` and the top of
`compute_BAI_percentile.py` for definitions. The output adds `BAIPercentile` (0–100)
and `n_similar`. The percentile is left empty when fewer than 10 similar
reference people exist.

## Files

- `compute_BAI_percentile.py` — the calculator
- `dataset_healthy_reference.csv` — reference people (n = 4,875)
- `example_input.csv` — example input
