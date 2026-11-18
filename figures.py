"""Figures and short checkpoint notes generated from the computed results."""

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

from analysis import BOOTSTRAP_COUNT, BOOTSTRAP_SEED, KB, NOISE_COUNT, NOISE_SEED, normalise_pairs

BLUE = "#185a82"
ORANGE = "#ba510d"
GREEN = "#16806a"

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.19,
    "axes.axisbelow": True,
    "figure.dpi": 120,
    "legend.frameon": False,
    "savefig.dpi": 220,
})


def save(fig, output, filename):
    fig.savefig(output / filename, facecolor="white")
    plt.close(fig)


def l1_figure(group, result, output):
    fig, ax = plt.subplots(figsize=(7.8, 4.65), layout="constrained")
    ax.plot(result["ef"], result["ratios"], "o-", color=BLUE, ms=4, label="Supplied sweep")
    for level, energy in result["crossings"].items():
        ax.axhline(level, color="#757575", linewidth=0.8, linestyle="--", zorder=1)
        ax.scatter(energy, level, marker="s", s=28, color=GREEN, zorder=3)
        ax.annotate(f"C = {level} at {energy:.4f} eV", (energy, level), xytext=(7, -16),
                    textcoords="offset points", color="#205c50")
    ef, ratio = result["ef_op"], result["results"]["correction_factor"]
    ax.scatter(ef, ratio, color=ORANGE, marker="D", s=55, zorder=4,
               label=f"Exact operating point ({ef:.4f} eV, {ratio:.3f})")
    ax.set(xlabel=r"Fermi level $E_F$ (eV)", ylabel=r"Correction factor $C=n_{\mathrm{band}}/n_F$",
           title=f"{group} L1 | Carrier-count correction")
    ax.legend(loc="upper left")
    save(fig, output, "L1_correction.png")


def l2_figures(group, result, output):
    fig = plt.figure(figsize=(11.4, 5.55), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=(1.1, 1))
    pair_ax = fig.add_subplot(grid[:, 0])
    hall_ax = fig.add_subplot(grid[0, 1])
    sheet_ax = fig.add_subplot(grid[1, 1], sharex=hall_ax)
    curve = result["model_curve"]
    pairs = normalise_pairs(result["means"])
    gap = result["results"]["Eg_kT"]
    name = f"Hard-gap curve ({gap:g} kBT)" if gap else "Gapless cone"
    pair_ax.plot(curve[:, 0], curve[:, 1], color=ORANGE, lw=1.45, label=name)
    pair_ax.scatter(pairs[:, 0], pairs[:, 1], color=BLUE, s=18, alpha=0.76, label="Repeat means")
    pair_ax.axvline(0, color="#777777", lw=0.8, linestyle="--")
    pair_ax.set(title="Fit in normalised pair space", xlabel=r"$R_H/\max |R_H|$", ylabel=r"$R_s/\max R_s$")
    pair_ax.legend(loc="lower left")
    pair_ax.set_aspect("equal", adjustable="box")

    indices = result["indices"]
    for rep, (colour, marker) in enumerate(zip((BLUE, ORANGE, GREEN), ("o", "s", "^"))):
        for channel, ax in enumerate((hall_ax, sheet_ax)):
            ax.plot(indices, result["repeats"][:, rep, channel], color=colour, lw=0.8,
                    marker=marker, ms=2.1, alpha=0.79, label=f"Repeat {rep + 1}")
    for i, (a, b) in enumerate(result["mean_brackets"]):
        for ax in (hall_ax, sheet_ax):
            ax.axvspan(indices[a], indices[b], alpha=0.13, color="#413281",
                       label="Mean Hall sign change" if i == 0 else "_nolegend_")
    hall_ax.axhline(0, color="#666666", linestyle="--", lw=0.8)
    hall_ax.set(title="Raw repeat measurements", ylabel=r"Raw $R_H$")
    sheet_ax.set(xlabel="Measurement index", ylabel=r"Raw $R_{xx}$")
    hall_ax.tick_params(labelbottom=False)
    for ax in (hall_ax, sheet_ax):
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    hall_ax.legend(ncol=2, loc="lower left", fontsize=8)
    fig.suptitle(f"{group} L2 | Density-of-states fit and transport measurements", fontsize=13)
    save(fig, output, "L2_model.png")

    fig, ax = plt.subplots(figsize=(8.8, 4.7), layout="constrained")
    lower, upper = result["ef_interval"]
    ax.fill_between(indices, lower, upper, color=BLUE, alpha=0.23, label="16th–84th bootstrap percentiles")
    ax.plot(indices, result["results"]["EF_trajectory_eV"], "o-", color=BLUE, ms=2.8,
            lw=1.25, label="Fit to repeat means")
    for a, b in result["mean_brackets"]:
        ax.axvspan(indices[a], indices[b], color="#413281", alpha=0.10)
    ax.axhline(0, color="#777777", lw=0.8, linestyle="--")
    ax.set(title=f"{group} L2 | Fermi-level trajectory", xlabel="Measurement index", ylabel=r"$E_F$ (eV)")
    ax.legend()
    save(fig, output, "L2_trajectory.png")


def map_axes(ax, values, limit, title):
    image = ax.imshow(values, cmap="RdBu_r", vmin=-limit, vmax=limit, origin="lower",
                      extent=(-3, 3, -3, 3), interpolation="nearest")
    ax.grid(False)
    ax.set(title=title, xlabel="x (probe pitches)", ylabel="y (probe pitches)",
           xticks=(-3, -2, -1, 0, 1, 2, 3), yticks=(-3, -2, -1, 0, 1, 2, 3))
    return image


def l3_figures(group, result, output):
    values = np.asarray(result["results"]["map"]).reshape(10, 10)
    fig, ax = plt.subplots(figsize=(7.4, 5.8), layout="constrained")
    image = map_axes(ax, values, result["amplitude"], f"{group} L3 | Zero-mean non-uniformity")
    fig.colorbar(image, ax=ax, label="Deviation from reference conductivity", format=PercentFormatter(1.0))
    for feature, marker in zip(result["results"]["features"], ("^", "v")):
        x, y = feature["x_probe_pitch"], feature["y_probe_pitch"]
        value = feature["amplitude_fraction"]
        name = "Maximum" if feature["type"] == "positive" else "Minimum"
        ax.scatter(x, y, marker=marker, facecolor="white", edgecolor="#333333", linewidth=0.85, s=66,
                   label=f"{name} ({x:+.1f}, {y:+.1f})  {value:+.2%}")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=9)
    save(fig, output, "L3_map.png")

    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.45), layout="constrained",
                             gridspec_kw={"width_ratios": (1.23, 1, 1)})
    tradeoff = axes[0]
    roughness, residual = result["sensitivity"].T
    tradeoff.loglog(roughness, residual, "o-", color=BLUE, ms=4)
    tradeoff.scatter(roughness[2], residual[2], s=50, color=ORANGE, zorder=4)
    offsets = [(9, -12), (-11, -16), (-34, 10), (6, 7), (7, 4)]
    for x, y, factor, offset in zip(roughness, residual, result["factors"], offsets):
        tradeoff.annotate(f"{factor:g}×", (x, y), xytext=offset, textcoords="offset points", fontsize=8)
    tradeoff.margins(x=0.23, y=0.28)
    tradeoff.set(title=r"Regularisation ($\eta/\eta_0$)", xlabel=r"Roughness $\|Lm\|_2$",
                 ylabel=r"Residual $\|Am-\bar{d}\|_2$")
    limit = max(float(np.max(np.abs(values))) for _, values, _ in result["spread"])
    centres = result["centres"]
    for ax, (pixel, values, width), title in zip(axes[1:], result["spread"], ("Centre", "Off-centre")):
        image = map_axes(ax, values, limit, f"{title} | FWHM {width:.2f} pitches")
        row, col = pixel
        ax.scatter(centres[col], centres[row], marker="+", color="white", s=52, linewidth=1.25)
        ax.set_xticks((-2, 0, 2))
        ax.set_yticks((-2, 0, 2))
    fig.colorbar(image, ax=axes[1:], label="Recovered unit-impulse amplitude", shrink=0.8)
    fig.suptitle(f"{group} L3 | Regularisation and point-spread tests", fontsize=13)
    save(fig, output, "L3_method.png")


def write_notes(group, l1, l2, l3, output):
    a, b, c = (item["results"] for item in (l1, l2, l3))
    checks = ", ".join(f"{ef:g} eV → {value:.4f}" for ef, value in l1["checks"].items())
    (output / "L1_note.txt").write_text(
        f"{group} L1 - Carrier-count correction\n\n"
        f"At the operating point of {l1['ef_op']:.9f} eV, the band-edge count is "
        f"{l1['n_band_op']:.6e} cm⁻² and the Fermi-referenced count is {a['n_fermi_op_cm2']:.6e} cm⁻². "
        f"The correction factor is {a['correction_factor']:.4f}. Applying it raises the mobility "
        f"from {l1['mobility_band']:.1f} to {a['mobility_corrected_cm2_Vs']:.1f} cm² V⁻¹ s⁻¹. "
        "As doping increases, the band-edge count includes more occupied states below the Fermi level. "
        "The Fermi-referenced count mainly samples the thermal tail above it, so the ratio grows.\n\n"
        f"C = 2 and C = 4 occur at {l1['crossings'][2]:.5f} eV and {l1['crossings'][4]:.5f} eV. "
        f"The 300 K checks give {checks}. The maximum relative difference between "
        f"the supplied and recalculated band-edge counts is {l1['band_max_rel_diff']:.2e}.\n",
        encoding="utf-8",
    )
    lo, hi = b["Eg_uncertainty_kT"]
    kt = KB * l2["temperature"]
    decision = (
        f"The hard-gap model passes the residual and gap-size rules. The fitted gap is "
        f"{b['Eg_kT']:g} kBT ({b['Eg_kT'] * kt:.5f} eV). The gap reduces the carrier population near "
        "neutrality in the fitted model and raises its sheet resistance there."
        if b["model"] == "hard_gap" else
        "The gapless cone is selected by the residual and gap-size rules, so the reported gap is zero. "
        "At this temperature the fitted cone can explain the transport through neutrality."
    )
    brackets = ", ".join(f"{l2['indices'][left]} to {l2['indices'][right]}"
                         for left, right in l2["mean_brackets"])
    improvement = (1 - l2["gap_rms"] / l2["cone_rms"]) * 100
    (output / "L2_note.txt").write_text(
        f"{group} L2 - Density of states from Hall and sheet-resistance pairs\n\n"
        f"The cone RMS residual is {l2['cone_rms']:.5f}. The best hard-gap candidate "
        f"({l2['best_gap_kT']:g} kBT) gives {l2['gap_rms']:.5f}, "
        f"a reduction of {improvement:.1f}%. {decision}\n\n"
        f"The mean Hall sign change is bracketed by measurement indices {brackets}. "
        f"Across {BOOTSTRAP_COUNT} draws (seed {BOOTSTRAP_SEED}), "
        f"{l2['gap_fraction']:.1%} select a hard gap. The selected-gap median is "
        f"{l2['bootstrap_median_kT']:g} kBT and the 16th–84th percentile interval is "
        f"[{lo:g}, {hi:g}] kBT. The trajectory figure shows the percentile band at each of the "
        f"{len(l2['indices'])} measurement indices.\n",
        encoding="utf-8",
    )
    peaks = " and ".join(
        f"{feature['amplitude_fraction']:+.2%} at "
        f"({feature['x_probe_pitch']:+.1f}, {feature['y_probe_pitch']:+.1f})"
        for feature in c["features"]
    )
    amplitude, detection, threshold = l3["amplitude"], c["detection_limit_fraction"], l3["threshold"]
    (output / "L3_note.txt").write_text(
        f"{group} L3 - Conductivity non-uniformity\n\n"
        f"The three scans give an RMS channel standard deviation of {l3['noise_sd_rms']:.4e}. "
        f"The map uses the guide's reference η₀ = {l3['eta']:.6g}. The method figure shows "
        "how changing η changes the residual and roughness. The strongest reconstructed "
        f"positive and negative features are {peaks}, with coordinates in probe pitches.\n\n"
        f"The positive-lobe moment gives FWHM estimates of {l3['spread'][0][2]:.3f} and "
        f"{l3['spread'][1][2]:.3f} pitches for the central and off-centre impulses, respectively. "
        f"The conservative resolution is {c['resolution_probe_pitch']:.3f} pitches. The 95th "
        f"percentile of the largest absolute noise-map values is {detection:.3%} "
        f"({NOISE_COUNT} simulations, seed {NOISE_SEED}).\n\n"
        f"The largest absolute deviation is {amplitude:.3%} against the {threshold:.3%} specification. "
        f"The nominal decision is {'pass' if c['pass'] else 'fail'}. The noise-adjusted bounds are "
        f"{amplitude - detection:.3%} and {amplitude + detection:.3%}, giving a {c['confidence']} result.\n",
        encoding="utf-8",
    )


def write_l2_table(result, output):
    """Record exact per-index repeat standard deviations and bootstrap intervals."""
    with (output / "L2_trajectory_and_repeats.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["measurement_index", "RH_mean_raw", "RH_sd_raw", "Rxx_mean_raw", "Rxx_sd_raw",
                         "EF_eV", "EF_p16_eV", "EF_p84_eV"])
        for index, mean, sd, ef, low, high in zip(
            result["indices"], result["means"], result["sd"], result["results"]["EF_trajectory_eV"],
            *result["ef_interval"]
        ):
            writer.writerow([index, mean[0], sd[0], mean[1], sd[1], ef, low, high])
