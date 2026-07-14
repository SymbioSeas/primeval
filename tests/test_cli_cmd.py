from primeval.cli import build_parser, build_snakemake_cmd, filter_terminal_line


def test_parser_defaults():
    args = build_parser().parse_args([])
    assert args.run_name == "results"
    assert args.directory == "."
    assert args.cores == 8
    assert args.force is False
    assert args.configfile is None


def test_parser_passthrough_after_ddash():
    args = build_parser().parse_args(["--run-name", "Vpop", "--", "-n", "--rerun-incomplete"])
    assert args.run_name == "Vpop"
    # REMAINDER keeps the leading "--"; main() strips it
    assert args.snakemake_args[-2:] == ["-n", "--rerun-incomplete"]


def test_build_cmd_core_flags():
    cmd = build_snakemake_cmd(
        snakefile="/repo/workflow/Snakefile", workdir="/wd",
        configfile="/wd/config.yaml", results_dir_rel="results/Vpop_2026-07-07",
        cores=8,
    )
    assert cmd[0] == "snakemake"
    assert "--snakefile" in cmd and "/repo/workflow/Snakefile" in cmd
    assert "--directory" in cmd and "/wd" in cmd
    assert "--configfile" in cmd and "/wd/config.yaml" in cmd
    assert "--config" in cmd and "results_dir=results/Vpop_2026-07-07" in cmd
    assert cmd[cmd.index("--cores") + 1] == "8"
    assert cmd[cmd.index("--quiet") + 1] == "rules"


def test_build_cmd_extra_appended():
    cmd = build_snakemake_cmd(
        snakefile="s", workdir="w", configfile="c",
        results_dir_rel="results/x", cores=4, extra=["-n"],
    )
    assert cmd[-1] == "-n"


def test_filter_progress_milestones():
    state = {}
    assert filter_terminal_line("5 of 100 steps (5%) done", state)
    assert filter_terminal_line("8 of 100 steps (8%) done", state) == []
    assert filter_terminal_line("12 of 100 steps (12%) done", state)


def test_filter_alerts_pass_through():
    state = {}
    assert filter_terminal_line("Error in rule run_blast:", state)
    assert filter_terminal_line("WARNING: something", state)
    assert filter_terminal_line("Traceback (most recent call last):", state)


def test_filter_suppresses_noise():
    state = {}
    assert filter_terminal_line("Building DAG of jobs...", state) == []
    assert filter_terminal_line("Select jobs to execute...", state) == []
    assert filter_terminal_line("Finished job 42.", state) == []
