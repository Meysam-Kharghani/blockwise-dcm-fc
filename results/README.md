# Derived result tables

`tables_csv/` contains detailed numerical summaries for DCM, PEB/BMA, raw FC, residual FC, model fit, posterior prevalence, and block timing. The publication-oriented figures are stored separately under `figures/`.

Conventions:

- DCM and PEB matrices use target rows and source columns.
- Paired frequentist block contrasts are computed on participant-level run-averaged values unless otherwise indicated.
- Frequentist contrast tables report FDR-adjusted q values.
- PEB inference is summarized with posterior means and posterior probabilities.
- DCM posterior prevalence uses `Pp > .95` as the high-confidence threshold.
