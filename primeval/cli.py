"""primeval command-line wrapper around the Snakemake workflow."""
import argparse
import subprocess
import sys
from pathlib import Path

from .rundir import resolve_run_dir

# <repo>/primeval/cli.py -> <repo>/workflow/Snakefile
WORKFLOW = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"


def build_parser():
    p = argparse.ArgumentParser(
        prog="primeval",
        description="Run the primeval in silico PCR pipeline on a set of assemblies.",
    )
    p.add_argument("--run-name", default="results",
                   help="Name for this run's results directory (default: results).")
    p.add_argument("--directory", default=".",
                   help="Analysis directory holding assemblies/, config.yaml, and outputs "
                        "(default: current directory).")
    p.add_argument("--configfile", default=None,
                   help="Pipeline config (default: <directory>/config.yaml).")
    p.add_argument("--force", action="store_true",
                   help="Reuse the existing dated run directory and resume unfinished work.")
    p.add_argument("--cores", type=int, default=8,
                   help="CPU cores for Snakemake (default: 8).")
    p.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                   help="Arguments after -- are passed through to Snakemake.")
    return p


def build_snakemake_cmd(*, snakefile, workdir, configfile, results_dir_rel,
                        cores, quiet=True, extra=None):
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--directory", str(workdir),
        "--configfile", str(configfile),
        "--config", f"results_dir={results_dir_rel}",
        "--cores", str(cores),
    ]
    if quiet:
        cmd += ["--quiet", "rules"]
    cmd += list(extra or [])
    return cmd


def _count_assemblies(configfile, workdir):
    try:
        import yaml
        cfg = yaml.safe_load(Path(configfile).read_text()) or {}
        adir = Path(workdir) / cfg.get("assembly_dir", "assemblies")
        return len(list(adir.glob("*.fna")))
    except Exception:
        return None


DRY_RUN_FLAGS = {"-n", "--dry-run", "--dryrun"}


def _is_dry_run(extra):
    return any(a in DRY_RUN_FLAGS for a in extra)


def _tee_run(cmd, logf):
    lf = open(logf, "a") if logf else None
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            sys.stdout.write(line)
            if lf:
                lf.write(line)
        proc.wait()
    finally:
        if lf:
            lf.close()
    return proc.returncode


def _print_header(run_name, rel_results, n_assemblies):
    n = "unknown" if n_assemblies is None else str(n_assemblies)
    print(f"primeval  |  run: {Path(rel_results).name}")
    print(f"  workflow : {WORKFLOW}")
    print(f"  input    : assemblies/  ({n} assemblies)")
    print(f"  results  : {rel_results}/")
    print()


def _print_footer(rel_results):
    print()
    print(f"Done. Results in {rel_results}/")
    print(f"  {rel_results}/reports/species_summary.csv")
    print(f"  {rel_results}/reports/assay_summary_long.csv")
    print(f"  {rel_results}/reports/detection_by_assembly.csv")
    print(f"  {rel_results}/reports/figures/species_detection_heatmap.pdf")


def main(argv=None, runner=None):
    args = build_parser().parse_args(argv)
    workdir = Path(args.directory).resolve()
    configfile = Path(args.configfile).resolve() if args.configfile else workdir / "config.yaml"
    if not configfile.exists():
        print(f"error: config file not found: {configfile}", file=sys.stderr)
        print(f"  copy the template:  cp <primeval-repo>/config/config.yaml "
              f"{workdir}/config.yaml", file=sys.stderr)
        return 2

    extra = list(args.snakemake_args)
    if extra and extra[0] == "--":
        extra = extra[1:]
    dry = _is_dry_run(extra)

    run_dir, mode = resolve_run_dir(workdir / "results", args.run_name, args.force)
    rel_results = run_dir.relative_to(workdir).as_posix()
    # A dry run must not create the run directory or a run.log (it would burn the
    # run name and litter empty dated dirs). Only real runs materialize the dir.
    if dry:
        logf = None
    else:
        run_dir.mkdir(parents=True, exist_ok=True)
        logf = run_dir / "run.log"

    _print_header(args.run_name, rel_results, _count_assemblies(configfile, workdir))

    cmd = build_snakemake_cmd(
        snakefile=WORKFLOW, workdir=workdir, configfile=configfile,
        results_dir_rel=rel_results, cores=args.cores, extra=extra,
    )
    runner = runner or _tee_run
    rc = runner(cmd, logf)
    if rc != 0:
        print(f"error: run failed (exit {rc}); see {logf}", file=sys.stderr)
    elif dry:
        print(f"\nDry run complete (no files written). Real run would write to {rel_results}/")
    else:
        _print_footer(rel_results)
    return rc


if __name__ == "__main__":
    sys.exit(main())
