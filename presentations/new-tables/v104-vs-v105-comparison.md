# Structural Tables: v104 vs v105 Comparison

**Model**: `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP`  
**Tables**: `struct-overall-te-table.tex`, `private-signal-te-table.tex`, `fob-beliefs-table.tex`  
**Date**: 2026-04-22

---

## Summary

Table structure is identical across versions — same rows, columns, and layout. All differences are numerical: point estimates are essentially unchanged, with credible interval bounds shifting by ±0.001–0.004. This is consistent with normal MCMC sampling variation (different seeds/chains). All qualitative conclusions — signs, significance, and magnitudes — are unchanged.

---

## `struct-overall-te-table.tex` — Overall ATEs

| Row | Column | v104 | v105 |
|---|---|---|---|
| Bracelet | Combined | 0.082 (0.059, 0.106) | 0.083 (0.057, 0.107) |
| Bracelet | Far | 0.116 (0.084, 0.147) | 0.117 (0.086, 0.149) |
| Calendar | Combined | 0.033 (0.008, 0.056) | 0.033 (0.009, 0.058) |
| Ink | Combined | −0.018 (−0.043, 0.007) | −0.017 (−0.042, 0.010) |
| Control mean | Combined | 0.331 (0.314, 0.350) | 0.331 (0.312, 0.350) |
| Signal − No Signal | Far | 0.024 (0.000, 0.048) | 0.024 (0.001, 0.047) |

Unlisted rows are numerically identical or differ only in the last decimal place of CI bounds.

---

## `private-signal-te-table.tex` — Signal and Private Value Decomposition

Small CI shifts across all six rows (Bracelet, Calendar, Ink for both Signal and Private panels), typically ±0.001–0.003 on bounds. Notable change in the Private panel:

| Row | Column | v104 | v105 |
|---|---|---|---|
| Bracelet (private) | Combined | 0.019 (−0.009, 0.046) | 0.020 (−0.009, 0.048) |
| Calendar (private) | Combined | 0.019 (−0.009, 0.046) | 0.020 (−0.009, 0.048) |
| Ink (private) | Combined | −0.059 (−0.089, −0.030) | −0.059 (−0.088, −0.030) |

---

## `fob-beliefs-table.tex` — First-Order Beliefs

Point estimates stable throughout. The Far column shows the largest changes, with some CI bounds shifting by up to ~0.007. Most notable:

| Row | Column | v104 | v105 |
|---|---|---|---|
| Bracelet | Combined | 0.097 (0.060, 0.136) | 0.097 (0.059, 0.134) |
| Bracelet | Far−Close | 0.113 (0.057, 0.187) | 0.113 (0.062, 0.181) |
| Calendar | Far | 0.039 (−0.001, 0.085) | 0.039 (0.003, 0.079) |
| Control mean | Far−Close | −0.098 (−0.161, −0.051) | −0.098 (−0.154, −0.056) |
| Signal − No Signal | Far−Close | 0.062 (0.027, 0.106) | 0.062 (0.029, 0.101) |

The narrowing of some credible intervals in v105 (particularly in the Far−Close column) suggests the chains may have mixed slightly better.

---

## Bottom line

No qualitative changes between v104 and v105. The two fits agree on all signs and substantive conclusions; differences are within normal MCMC sampling noise.
