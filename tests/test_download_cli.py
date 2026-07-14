"""Tests for the download-assemblies console wrapper."""
import shutil
import subprocess

import primeval.download_cli as dc


def test_wrapper_resolves_to_download_script():
    assert dc.SCRIPT.name == "download_assemblies.sh"
    assert dc.SCRIPT.exists()


def test_download_command_on_path_shows_usage():
    exe = shutil.which("download-assemblies")
    assert exe, "download-assemblies not on PATH (is the package installed?)"
    r = subprocess.run([exe, "-h"], capture_output=True, text=True)
    assert "TAXON" in (r.stdout + r.stderr)


def test_help_lists_assembly_source_flag():
    exe = shutil.which("download-assemblies")
    assert exe
    out = subprocess.run([exe, "-h"], capture_output=True, text=True)
    combined = out.stdout + out.stderr
    assert "-s" in combined and "assembly-source" in combined.lower()
    assert "refseq" in combined.lower()


def test_invalid_assembly_source_rejected(tmp_path):
    exe = shutil.which("download-assemblies")
    assert exe
    r = subprocess.run(
        [exe, "-s", "bogus", "-t", "Foo", "-o", str(tmp_path / "asm")],
        capture_output=True, text=True,
    )
    assert r.returncode != 0
    assert "refseq" in (r.stdout + r.stderr).lower()
