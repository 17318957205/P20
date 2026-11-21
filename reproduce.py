"""Reproduce the regular or blind G20 run directly from a supplied data folder."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter

import matplotlib
import numpy as np

from analysis import BOOTSTRAP_COUNT, BOOTSTRAP_SEED, NOISE_COUNT, NOISE_SEED, analyse_l1, analyse_l2, analyse_l3
from figures import l1_figures, l2_figures, l3_figures, l3_spectrum_figure, write_l2_table, write_notes

ROOT = Path(__file__).resolve().parent
THRESHOLD = 0.030  # L3 specification in the supplied G20 blind package's README_FIRST.txt.
INPUTS = (
    "L1_data.csv", "L1_params.json", "L2_data.csv", "L2_meta.json",
    "L3_A.npy", "L3_d.csv", "L3_meta.json",
)
CODE_FILES = ("analysis.py", "figures.py", "reproduce.py", "requirements.txt")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def git_command(*args):
    return subprocess.run(
        ["git", "-C", str(ROOT), *args], check=True, text=True, stdout=subprocess.PIPE,
    ).stdout.strip()


def pipeline_commit():
    commit = git_command("rev-parse", "HEAD")
    if git_command("status", "--porcelain", "--", *CODE_FILES):
        raise ValueError("Commit changes to the pipeline files before generating marked results.")
    return commit


def verify_input(data):
    manifest = read_json(data / "manifest.json")
    group, dataset = manifest["group"], manifest["dataset_set"]
    if dataset not in ("regular", "blind"):
        raise ValueError(f"Unsupported dataset set {dataset!r}.")
    template_name = "results_template.json" if dataset == "regular" else "blind_results.json"
    entries = {entry["name"]: entry for entry in manifest["files"]}
    for name in (*INPUTS, template_name, "README_FIRST.txt"):
        entry = entries[name]
        contents = (data / name).read_bytes()
        if entry["bytes"] != len(contents) or entry["sha256"] != hashlib.sha256(contents).hexdigest():
            raise ValueError(f"Input differs from the supplied manifest: {name}.")
    documents = [read_json(data / name) for name in ("L1_params.json", "L2_meta.json", "L3_meta.json", template_name)]
    if any(doc["group"] != group or doc["dataset_set"] != dataset for doc in documents):
        raise ValueError("The input manifest, metadata and template do not identify the same group and set.")
    return documents


def validate_results(template, results, expected_trajectory_length):
    if template.keys() != results.keys():
        raise ValueError("Result root keys do not match the supplied template.")
    for stage in ("L1", "L2", "L3"):
        if template[stage].keys() != results[stage].keys():
            raise ValueError(f"{stage} result fields do not match the supplied template.")
    if len(results["L2"]["EF_trajectory_eV"]) != expected_trajectory_length:
        raise ValueError("L2 trajectory length does not match the data.")
    l3 = results["L3"]
    if len(l3["map"]) != 100 or abs(math.fsum(l3["map"]) / 100) > 1e-12:
        raise ValueError("L3 must have 100 row-major values and zero mean.")
    if type(l3["pass"]) is not bool or not l3["confidence"]:
        raise ValueError("L3 must include a Boolean pass field and a confidence statement.")
    # allow_nan=False also validates all numerical values nested in the output.
    return json.dumps(results, indent=2, allow_nan=False) + "\n"


def write_coursework_metrics(output, l1, l2, l3):
    """Write report diagnostics without adding fields to the marked JSON schema."""
    amplitude = l3["amplitude"]
    limit = l3["results"]["detection_limit_fraction"]
    coordinates = [
        (feature["x_probe_pitch"], feature["y_probe_pitch"])
        for feature in l3["results"]["features"]
    ]
    # Current-related Hall sensitivity 1/(e*n_si) equals mobility_si times sheet resistance.
    mobility_to_sensitivity = l1["rho_ohm_sq"] * 1e-4
    trajectory = l2["results"]["EF_trajectory_eV"]
    metrics = {
        "L1": {
            **{key: l1[key] for key in ("ef_op", "temperature", "dos_scale_c", "rho_ohm_sq", "n_band_op", "mobility_band",
                                       "closed_form_relative_difference", "checks", "crossings", "band_max_rel_diff")},
            "sweep_ef_eV": l1["ef"].tolist(),
            "sweep_band_density_cm2": l1["n_band_sweep_cm2"].tolist(),
            "sweep_fermi_density_cm2": l1["n_fermi_sweep_cm2"].tolist(),
            "sweep_sheet_resistance_ohm_sq": l1["rho_sweep_ohm_sq"].tolist(),
            "sweep_band_mobility_cm2_Vs": l1["mobility_band_sweep_cm2_Vs"].tolist(),
            "sweep_fermi_mobility_cm2_Vs": l1["mobility_fermi_sweep_cm2_Vs"].tolist(),
            "sensitivity_band_V_per_A_per_T": l1["mobility_band"] * mobility_to_sensitivity,
            "sensitivity_fermi_V_per_A_per_T": l1["results"]["mobility_corrected_cm2_Vs"] * mobility_to_sensitivity,
        },
        "L2": {
            **{key: l2[key] for key in ("temperature", "cone_rms", "gap_rms", "best_gap_kT",
                                       "bootstrap_median_kT", "gap_fraction")},
            "relative_residual_reduction": 1 - l2["gap_rms"] / l2["cone_rms"],
            "measurement_count": len(l2["indices"]),
            "mean_hall_brackets": [[int(l2["indices"][left]), int(l2["indices"][right])]
                                   for left, right in l2["mean_brackets"]],
            "ef_endpoints_eV": [trajectory[0], trajectory[-1]],
            "bootstrap_count": BOOTSTRAP_COUNT,
            "bootstrap_seed": BOOTSTRAP_SEED,
        },
        "L3": {
            **{key: l3[key] for key in ("eta", "threshold", "amplitude", "row_sum", "noise_sd_rms", "numerical_rank")},
            "singular_values": l3["singular_values"].tolist(),
            "map_mean": math.fsum(l3["results"]["map"]) / 100,
            "eta_factors": l3["factors"],
            "roughness_residual_pairs": l3["sensitivity"].tolist(),
            "point_spread_fwhm_probe_pitch": [width for _, _, width in l3["spread"]],
            "extrema_separation_probe_pitch": math.dist(*coordinates),
            "noise_adjusted_bounds_fraction": [amplitude - limit, amplitude + limit],
            "noise_simulations": NOISE_COUNT,
            "noise_seed": NOISE_SEED,
        },
    }
    (output / "coursework_metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )


def run(data, output, threshold):
    start = perf_counter()
    commit = pipeline_commit()
    l1_params, l2_meta, _, template = verify_input(data)
    group, dataset = template["group"], template["dataset_set"]
    print(f"{group} {dataset} | pipeline commit {commit}", flush=True)
    print(f"Python {platform.python_version()} | NumPy {np.__version__} | Matplotlib {matplotlib.__version__}")
    print(f"Data folder {data} | input manifest hashes verified")

    l1 = analyse_l1(data / "L1_data.csv", l1_params)
    l2 = analyse_l2(data / "L2_data.csv", float(l2_meta["T_K"]))
    l3 = analyse_l3(data / "L3_A.npy", data / "L3_d.csv", threshold)
    results = dict(template, pipeline_commit=commit, L1=l1["results"], L2=l2["results"], L3=l3["results"])
    payload = validate_results(template, results, len(l2["indices"]))
    print(f"L1 300 K reference factors {[round(value, 5) for value in l1['checks'].values()]}")
    print(f"L1 EF={l1['ef_op']:.9f} eV | C={results['L1']['correction_factor']:.6f} "
          f"| mobility={results['L1']['mobility_corrected_cm2_Vs']:.3f} cm^2/(V s)")
    print(f"L2 {len(l2['indices'])} indices, 3 repeats per index | "
          f"cone RMS={l2['cone_rms']:.6f}, best hard-gap RMS={l2['gap_rms']:.6f} "
          f"({l2['best_gap_kT']:g} kBT)")
    print(f"L2 selected {results['L2']['model']}, Eg={results['L2']['Eg_kT']:g} kBT | "
          f"hard-gap fraction={l2['gap_fraction']:.3f} ({BOOTSTRAP_COUNT} draws, seed {BOOTSTRAP_SEED})")
    print(f"L3 max row sum={l3['row_sum']:.3e} | "
          f"100 finite map values, mean={math.fsum(results['L3']['map']) / 100:.3e} | eta={l3['eta']:.8g}")
    print(f"L3 resolution={results['L3']['resolution_probe_pitch']:.6f} pitches | "
          f"delta95={results['L3']['detection_limit_fraction']:.6%} ({NOISE_COUNT} draws, seed {NOISE_SEED})")
    print(f"L3 peak={l3['amplitude']:.6%}, threshold={threshold:.3%} | {results['L3']['confidence']}")

    results_name = "results.json" if dataset == "regular" else "blind_results.json"
    hash_name = "commit_hash.txt" if dataset == "regular" else "pipeline_commit.txt"
    # Generate into a sibling temp directory so partial calculations don't appear as a complete run.
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".g20-", dir=output.parent) as temp:
        staging = Path(temp) / "run"
        staging.mkdir()
        (staging / results_name).write_text(payload, encoding="utf-8")
        (staging / hash_name).write_text(commit + "\n", encoding="utf-8")
        if dataset == "regular":
            l1_figures(group, l1, staging)
            l2_figures(group, l2, staging)
            l3_figures(group, l3, staging)
            l3_spectrum_figure(group, l3, staging)
            write_notes(group, l1, l2, l3, staging)
            write_l2_table(l2, staging)
            write_coursework_metrics(staging, l1, l2, l3)
        staging.rename(output)
    print(f"Wrote {results_name} and {hash_name} to {output} in {perf_counter() - start:.2f} s")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="Extracted regular or blind group data folder")
    parser.add_argument("--output", type=Path, required=True, help="New output folder")
    parser.add_argument("--threshold", type=float, default=THRESHOLD, help="L3 limit as a fraction, default 0.030")
    args = parser.parse_args()
    data, output = args.data.resolve(), args.output.resolve()
    if not data.is_dir():
        parser.error(f"data directory not found: {data}")
    if data == output or data in output.parents:
        parser.error("--output must be outside the raw data folder")
    if output.exists():
        parser.error("--output must be a new folder")
    if not math.isfinite(args.threshold) or args.threshold <= 0:
        parser.error("--threshold must be positive and finite")
    run(data, output, args.threshold)


if __name__ == "__main__":
    main()
