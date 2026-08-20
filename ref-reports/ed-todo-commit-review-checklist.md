# Ed ToDo Commit Review Checklist

Use this to track whether each atomic change is ready to keep, needs revision, or should be dropped.

Code PR: https://github.com/EdJeeOnGitHub/takeup/pull/1  
Analysis branch: `ed-refine-todos` into `w-beliefs`  
Paper branch: local Overleaf branch `../overleaf/overleaf-takeup:ed-refine-todos`

| ToDo | Scope | Commit(s) | Review status | Notes / requested changes |
|---|---|---|---|---|
| 1 Prior deworming | Analysis | `54fa50d`, `7705024`, `7cfa8f3` | Pending |  |
| 1 Prior deworming | Paper | `b78969a`, `9baa834` | Pending |  |
| 2 Balance control means | Analysis | `f9d733e` | Accept | Cherry-picked to `w-beliefs` as `690360f` and pushed to `origin/w-beliefs`. |
| 4 `\psi` USD units | Analysis | `219dc73` | Accept | Cherry-picked to `w-beliefs` as `572018c` and pushed to `origin/w-beliefs`. |
| 4 `\psi` USD units | Paper | `15fc692` | Accept | Applied to Overleaf `main` as `7fe05b8` and pushed to `overleaf/master`; omitted the “estimate implies...” clause per review. |
| 5/8 Potential-distance text | Paper | `9f5df01` | Accept | Applied to Overleaf `main` as `905c279` and pushed to `overleaf/master`; includes user wording tweak removing the pooling/weighting sentence in Appendix E. |
| 6 `\sigma_u` identification | Analysis | `fe9c8d8`, `7705024` | Accept | Applied to `w-beliefs` as `7b58c6b` and pushed to `origin/w-beliefs`; kept the sigma-u prior/posterior plotting script with variable-filtered Stan reads and Canva palette. |
| 6 `\sigma_u` identification | Paper | `10cfa9c`, `9baa834` | Accept | Applied to Overleaf `main` as `68712c0` and pushed to `overleaf/master`; includes final user wording on the identification paragraph and priors paragraph, plus the prior/posterior figure. |
| 10 Figure 3(b) ordering | Paper | None | Accept | No change needed; reviewer concern is not borne out by the figure/data, so do not cherry-pick `c013bdb`. |
| 11 HMC diagnostics | Analysis | `e7e2ac1` | Accept | Applied to `w-beliefs` as `7e7a545` and pushed to `origin/w-beliefs`; kept only the diagnostics helper script and removed treedepth from the script output. |
| 11 HMC diagnostics | Paper | `2a3ca83`, `9831b83` | Accept | Applied to Overleaf `main` as `9126723` and `4d1722e`, then pushed to `overleaf/master`; notes report \(\hat R\), ESS, and divergences, with no treedepth text. |
| 12 E3 feasible-distance RI | Analysis | `41873f6` | Accept | Cherry-picked to `w-beliefs` as `a3e78dc` and pushed to `origin/w-beliefs`. |
| 12 E3 feasible-distance RI | Paper | `90a16d5` | Accept | Cherry-picked/rebased to Overleaf `main` as `da36a88` and pushed to `overleaf/master`. |
| 17 Figure M1 parameters | Paper | `7df4781` | Accept | Applied to Overleaf `main` as `b05d100` and pushed to `overleaf/master`; parameters verified in `scratch/quick-vstar.jl`, with “type index” wording used without the hyphen. |
| 18 Closest PoT / ITT | Analysis | `fe26587`, `a22d7eb`, `88d3540` | Accept | Accepted analysis audit. Added guarded cluster-ID parsing, source-to-analysis merge checks, cluster-level audit output, documented switch-count assertion, and comparison to `docs/randomization-mismatches.csv`. |
| 18 Closest PoT / ITT | Paper | `9195391`; Overleaf `e9a6578`, `497dd78` | Accept | Added original-distance ITT take-up and FOB robustness tables to Overleaf; SOB table intentionally excluded. Remaining paper wording edits incorporated in `497dd78`. |
| 24 Normal/convolution wording | Paper | `d08cf63` | Accept | User edit replaces “normal-sum specification” with “Since \(w_i\) is the sum of two independent normal variables...”. |
| 25 SMS SE / table style | Analysis | `b59c86b` | Accept | Cherry-picked to `w-beliefs` as `51d3c4c` and pushed to `origin/w-beliefs`. |
| Tracker / commit map | Analysis | `b1a6458`, `faff25f` | Pending |  |

Review status suggestions: `Accept`, `Revise`, `Drop`, `Question`.
