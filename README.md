# Detection of Structural Breaks in Variance — ICSS Algorithm

## Overview

This repository contains a full SAS/IML implementation of the **Iterative Cumulative Sum of Squares (ICSS)** algorithm introduced by Inclán & Tiao (1994), along with an empirical application to the **USD/RUB monthly exchange rate (1995–2026)**.

The ICSS algorithm detects multiple structural breaks in the unconditional variance of a time series without requiring the number of breaks to be specified in advance.

---

## Repository Structure
```
.
├── Code_1_2_ICSS.sas          # ICSS module definition + validation on simulated data
├── ROUBLE_RUSSE_ICSS.sas      # Empirical application to the USD/RUB exchange rate
└── ICSS_paper.pdf             # Full report (methodology, results, interpretation)
```

---

## How It Works

The algorithm is built around three steps:

1. **Phase A — Forward & Backward search**: the series is scanned left-to-right and right-to-left to collect candidate breakpoints using the `Dk` statistic.
2. **Deduplication**: candidate points from both passes are merged and sorted.
3. **Phase B — Refinement loop**: each segment between detected breaks is re-tested until convergence (stable breakpoints, max shift ≤ 2).

The core test statistic is:
```
Dk = Ck / CT − k / T
```

where `Ck` is the cumulative sum of squares up to observation `k`. Under constant variance, `Dk` fluctuates around zero. A break is flagged when:
```
M = sqrt(T/2) × max|Dk|  >  1.358   (5% critical value)
```

---

## Files

### `Code_1_2_ICSS.sas`
- Defines and **stores** the reusable `icss_test` module (`store module=icss_test`)
- Validates the implementation on **simulated data**: T = 600 observations with true breaks at t = 200 and t = 400 (σ² = 1 → 3 → 1)
- Produces three diagnostic plots: observed series, cumulative sum of squares Ck, and Dk statistic — all with detected breakpoints marked in red

### `ROUBLE_RUSSE_ICSS.sas`
- Imports and cleans monthly USD/RUB data from a CSV file
- Computes log-returns: `r_t = ln(P_t / P_{t-1})`
- Runs an ADF stationarity test and a Ljung-Box autocorrelation test
- Fits an **ARMA(1,1) filter** on returns to satisfy the iid assumption required by ICSS
- Applies the ICSS algorithm on the ARMA residuals
- Outputs a table of breakpoint dates with corresponding exchange rates and variance estimates

---

## Results Summary

Four structural breaks were detected in the USD/RUB series:

| Regime | Break obs. | Date         | Variance |
|--------|-----------|--------------|----------|
| 1      | 8         | Jan 1995     | 0.0055   |
| 2      | 51        | Jan 1999     | 0.0101   |
| 3      | 257       | Jan 2016     | 0.0015   |
| 4      | 333       | Jan 2022     | 0.0035   |
| 5      | —         | —            | 0.0027   |

Each break aligns with a documented macroeconomic or geopolitical shock: post-Soviet transition (1995), Russian sovereign default (1998), Crimea annexation + oil price collapse (2014–2016), and invasion of Ukraine (2022).

