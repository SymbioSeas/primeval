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
