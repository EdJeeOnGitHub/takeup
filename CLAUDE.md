# Agent Container Environment

You are running inside the agent container, not directly on the host.
Prefer these container-installed language binaries:

- Julia: `/usr/local/bin/julia` (`julia` is on `PATH`)
- R: `/usr/bin/R` (`R` is on `PATH`)
- Rscript: `/usr/bin/Rscript` (`Rscript` is on `PATH`)

Container CPU allocation: `8` cores. For expensive Julia or R jobs, use multiple threads where it is appropriate for the workload, for example `JULIA_NUM_THREADS=8 julia ...` or R/OpenMP/BLAS thread settings capped at `8`.

Useful container paths:

- Primary workdir: `/workspace/takeup`
- Julia depot: `/julia_depot`
- Scratch: `/scratch`
- renv cache: `/renv-cache`

R package installs:

- The R base image is configured for Posit Package Manager Linux binaries for Ubuntu noble via `https://packagemanager.posit.co/cran/__linux__/noble/latest`.
- Do not conclude binaries are disabled just because an activated `renv` project reports `getOption("repos")` as the lockfile repositories or `getOption("pkgType")` as `source`. Linux RSPM binaries are selected by the `__linux__/noble` repository plus the configured `HTTPUserAgent`, and `renv.config.repos.override` points restores/installs at that repository.
- Prefer `renv::restore(prompt = FALSE)` for locked projects. For ad hoc package installs, pass the override explicitly if needed: `repos <- getOption("renv.config.repos.override", getOption("repos")); install.packages(pkgs, repos = repos)`.
- If a compiled package fails to load with a missing shared library such as `libicui18n.so.*`, treat it as a stale compiled package/cache entry; reinstall or rebuild that package for the current image rather than changing the container user or installing as root.

Mounted repos:

- takeup: `/workspace/takeup` (rw)
- overleaf-takeup-slides: `/home/agent/projects/overleaf/overleaf-takeup-slides` (rw)

Other mounted data paths:

- none

If project-local instructions mention host-only paths such as `/home/ed/.juliaup/bin`, translate them to these container paths before running commands.
