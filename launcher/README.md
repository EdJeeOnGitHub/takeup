# Quest launcher setup

Backend: `quest-general` (account `p53238`). These presets are configured, not executed.

## Available experiments

`environment-check` loads the required packages without fitting a model.

- `models`: 4 CPUs, 24G, 48:00:00; commands are recorded in `launcher/jobs.json`.

## Inputs and results

Bind the input named `data` to a registered immutable artifact. Its root must contain: `analysis.RData`. The launcher stages it at `/inputs/data`; it is never fetched from GitHub by Quest.

The driver copies committed source to per-task writable scratch for legacy relative
paths and Stan compilation. Fits, logs, and execution receipts are retained below
`/outputs`; source and input mounts stay read-only. It does not submit nested Slurm
jobs, copy manuscript files, push results, or run a devcontainer.

## Environment

The launcher builds the digest-pinned Containerfile on limiting-factor and transfers
the resulting verified SIF. It restores the committed R package lock during image
build, never on Quest compute nodes. The base is the existing R 4.5.3/Julia runtime
on limiting-factor. Older R lockfiles may require compatibility fixes on the first
build; this setup is not a claim of bit-for-bit reproduction on a newer R version.
CmdStan 2.36.0 is installed when required. All image builds and scientific execution
remain untested until explicitly requested. Request an environment check first.

## Planning and launch

Commit the launcher files before planning. A local plan performs no image build or
submission:

```sh
research-run plan environment-check --repo . --backend quest-general
```

For server-side launching, use this project's full published commit SHA and one of
the experiment names above. The host registration authorizes the project/backend;
it does not automatically run any workflow. Quest's host-owned profile chooses
`short` for requests up to 4 hours, `normal` up to 48 hours and `long` up to 168 hours.

The launcher branch externalizes the former `data` and `multilvlr` Git submodules. Private data is an artifact; `multilvlr` is a registered code dependency pinned to the original submodule commit. Standalone checkouts must supply those paths explicitly. The main branch is unchanged.
