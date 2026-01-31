import pytest
import subprocess
import os
import shutil
import pandas as pd
import numpy as np
import bz2
import sys

@pytest.fixture(scope="module")
def demo_env():
    """
    Sets up the test environment
    """
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    demo_dir = os.path.join(project_root, "demo_data")
    out_dir = os.path.join(project_root, "tests/test_out_pytest")

    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    files_to_decompress = [
        "SRR12996245.1pct.kraken",
        "SRR12996245.1pct_1.fastq",
        "SRR12996245.1pct_2.fastq"
    ]

    print(f"\n[Setup] Preparing test data in {out_dir}...")

    for filename in files_to_decompress:
        source_bz2 = os.path.join(demo_dir, filename + ".bz2")
        dest_file = os.path.join(out_dir, filename)
        
        if os.path.exists(source_bz2):
            with bz2.open(source_bz2, "rb") as f_in, open(dest_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            pytest.fail(f"Required source file missing: {source_bz2}")

    yield {'out': out_dir, 'demo': demo_dir}

def test_guacamole_metrics(demo_env):
    """
    Runs guacamole and asserts efficiency and abundance ranges.
    """
    out_dir = demo_env['out']
    demo_dir = demo_env['demo']
    output_file_name = "SRR12996245.1pct.guaca"

    cmd = [
        sys.executable, "-m", "guacamole.guacamole",
        "--output", output_file_name,
        "--kraken_report", os.path.join(demo_dir, "SRR12996245.1pct_report.txt"),
        "--kraken_file", os.path.join(out_dir, "SRR12996245.1pct.kraken"),
        "--read_files", os.path.join(out_dir, "SRR12996245.1pct_1.fastq"), os.path.join(out_dir, "SRR12996245.1pct_2.fastq"),
        "--kraken_db", os.path.join(demo_dir, "demo_db"),
        "--read_len", "150",
        "--fragment_len", "400",
        "--length_correction", "True",
        "--threshold", "5",
        "--plot", "False"
    ]

    result = subprocess.run(
        cmd, cwd=out_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    
    # Assert successful run
    assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"

    # Check Efficiencies
    eff_path = os.path.join(out_dir, "efficiencies.txt")
    assert os.path.exists(eff_path), "efficiencies.txt was not created"
    
    efficiencies = np.loadtxt(eff_path)
    max_eff = np.argmax(efficiencies)
    
    assert 20 <= max_eff <= 30, f"Max efficiency {max_eff:.2f} is out of range [20, 30]"

    # Check Abundances
    guaca_path = os.path.join(out_dir, output_file_name)
    assert os.path.exists(guaca_path), f"{output_file_name} was not created"

    df = pd.read_csv(guaca_path, sep="\t")
    abundances = df['GuaCAMOLE_estimate']

    assert (np.sum(abundances) == pytest.approx(1.0, abs=0.01)), f"Sum of abundances is not 1: {np.sum(abundances)}"
    assert (abundances >= 0.03).all(), f"Some abundances are < 4%:\n{df[abundances < 0.04]}"
    assert (abundances <= 0.12).all(), f"Some abundances are > 12%:\n{df[abundances > 0.12]}"