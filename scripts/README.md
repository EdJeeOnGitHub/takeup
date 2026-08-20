# Executable workflows

Scripts are grouped by analysis domain. Run them from the repository root so
project-relative data and output paths resolve consistently. Shared functions
live under `R/`; side-effectful setup files remain here until separately
refactored into explicit functions.

- `workflow/`: staging, artifact registries, and replication packaging.
- `checks/`: validation and audits.
- `reduced-form/`, `balance/`, `structural/`, `policy/`: main analyses.
- `appendix/`: supported robustness and referee-response analyses.
- `design/`: data preparation and design workflows.
