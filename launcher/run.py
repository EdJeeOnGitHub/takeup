#!/usr/bin/env python3
"""Generate launcher tasks or execute a committed preset in isolated writable scratch."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "launcher/jobs.json").read_text())

def generate(name):
    if name == "environment-check":
        resources = dict(cpus=1, memory="4G", time_limit="00:10:00")
    else:
        job = SETTINGS["jobs"][name]
        resources = {key: job[key] for key in ("cpus", "memory", "time_limit")}
    return {"tasks": [dict(id=name, name=name,
             command="python3 launcher/run.py execute " + name, **resources)]}

def execute(name):
    work = Path(os.environ["RESEARCH_RUN_WORK"])
    output = Path(os.environ["RESEARCH_RUN_OUTPUT"]) / name
    output.mkdir(parents=True, exist_ok=False)
    project = work / "project"
    shutil.copytree(ROOT, project, ignore=lambda directory, names: set(names) & ({".git", "__pycache__", "logs"} | ({"output"} if Path(directory) == ROOT / "data" else set())),
                    copy_function=shutil.copyfile)
    # Source is mounted read-only; compilation and legacy relative writes happen here.
    env = dict(os.environ, R_PROFILE_USER="/dev/null", R_ENVIRON_USER="/dev/null",
               RENV_CONFIG_AUTOLOADER_ENABLED="FALSE", R_LIBS_USER="/opt/research-library",
               R_LIBS_SITE="/opt/research-library:/usr/local/lib/R/site-library:/usr/lib/R/site-library", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", MAKEFLAGS="-j2")
    if SETTINGS.get("cmdstan"):
        env["CMDSTAN"] = "/opt/cmdstan/cmdstan-2.36.0"
    values = dict(output=str(output), input="/inputs/data")
    packages = SETTINGS["packages"]
    check = "p <- " + "c(" + ",".join(json.dumps(x) for x in packages) + "); " + \
            "stopifnot(all(vapply(p, requireNamespace, logical(1), quietly=TRUE))); sessionInfo()"
    checks = [["Rscript", "--vanilla", "-e", check]]
    if SETTINGS.get("cmdstan"):
        checks.append(["Rscript", "--vanilla", "-e", "stopifnot(dir.exists(cmdstanr::cmdstan_path())); print(cmdstanr::cmdstan_version())"])
    if SETTINGS.get("julia"):
        (work / "PCNMCMC").symlink_to("/workspace/PCNMCMC", target_is_directory=True)
        checks.append(["julia", "--startup-file=no", "--project=julia/full-covariance-gp", "-e", "using PCNMCMC, Distributions; println(VERSION)"])
    if SETTINGS["project"] == "takeup":
        (project / "multilvlr").symlink_to("/workspace/multilvlr", target_is_directory=True)
    commands = checks
    if name != "environment-check":
        job = SETTINGS["jobs"][name]
        for entry in job.get("requires", []):
            if not (Path(values["input"]) / entry).exists():
                raise SystemExit("Missing input artifact file/directory: " + entry)
        env.update({key: value.format(**values) for key, value in job.get("env", {}).items()})
        if SETTINGS["project"] == "takeup":
            data = project / "data"
            data.mkdir(exist_ok=True)
            for file in Path(values["input"]).iterdir():
                (data / file.name).symlink_to(file, target_is_directory=file.is_dir())
        # Preserve outputs of entry points that use project-relative paths.
        for relative in ("data/output", "logs"):
            destination = project / relative
            retained = output / relative.replace("/", "-")
            if destination.exists():
                shutil.move(str(destination), str(retained))
            else:
                retained.mkdir()
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(retained, target_is_directory=True)
        commands += [[arg.format(**values) for arg in command] for command in job["commands"]]
    record = dict(project=SETTINGS["project"], experiment=name,
                  commit=env.get("RESEARCH_RUN_COMMIT"), commands=commands, status="running")
    receipt = output / "execution.json"
    receipt.write_text(json.dumps(record, indent=2) + "\n")
    try:
        for index, command in enumerate(commands):
            print("Running:", command, flush=True)
            with (output / (str(index) + ".log")).open("w") as log:
                subprocess.run(command, cwd=project, env=env, stdout=log,
                               stderr=subprocess.STDOUT, check=True)
        record["status"] = "completed"
    except BaseException:
        record["status"] = "failed"
        raise
    finally:
        receipt.write_text(json.dumps(record, indent=2) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["generate", "execute"])
    parser.add_argument("experiment", choices=["environment-check"] + list(SETTINGS["jobs"]))
    args = parser.parse_args()
    if args.action == "generate":
        print(json.dumps(generate(args.experiment)))
    else:
        execute(args.experiment)
