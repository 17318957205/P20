"""Numerical methods from the QXU6028 Virtual Lab Practical Methods Guide."""

import csv
from pathlib import Path

import numpy as np

KB = 8.617333e-5  # eV/K
NA = 3.816e15  # atoms/cm^2
CHARGE = 1.602176634e-19  # C
EMAX = 2.5  # eV
BOOTSTRAP_SEED = 6028
BOOTSTRAP_COUNT = 200
NOISE_SEED = 6029
NOISE_COUNT = 500

# NumPy 1.26 uses the earlier name for the same trapezoidal rule.
trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def fermi(energy, ef, temperature):
    exponent = np.clip((energy - ef) / (KB * temperature), -60, 60)
    return 1.0 / (np.exp(exponent) + 1.0)


def count_carriers(ef, temperature, c, lower):
    energy = np.linspace(lower, EMAX, 20001)
    return float(NA * trapezoid(c * energy * fermi(energy, ef, temperature), energy))


def correction_at(ef, temperature, c):
    return count_carriers(ef, temperature, c, 0) / count_carriers(ef, temperature, c, ef)


def bracketed_crossing(target, ef, temperature, c):
    """Bisect the model curve inside the observed sweep."""
    lower, upper = float(ef[0]), float(ef[-1])
    if not correction_at(lower, temperature, c) <= target <= correction_at(upper, temperature, c):
        raise ValueError(f"The sweep does not bracket C={target}.")
    for _ in range(30):
        middle = (lower + upper) / 2
        if correction_at(middle, temperature, c) < target:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def analyse_l1(path: Path, params):
    data = np.genfromtxt(path, delimiter=",", names=True, ndmin=1)
    ef = data["EF_eV"]
    supplied = data["n_band_cm2"]
    c, temp, ef_op, rho = (float(params[key]) for key in ("c", "T", "EF_op_eV", "rho_ohm_sq"))
    if not (c > 0 and temp > 0 and 0 <= ef_op < EMAX and rho > 0):
        raise ValueError("Invalid L1 operating parameters.")
    if not (len(ef) >= 2 and np.all(np.diff(ef) > 0) and np.isfinite(supplied).all()):
        raise ValueError("Expected a finite, increasing L1 sweep.")

    references = {0.0: 1.00, 0.073: 2.00, 0.100: 2.60, 0.200: 5.10}
    checks = {ref: correction_at(ref, 300.0, c) for ref in references}
    for ref, expected in references.items():
        if abs(checks[ref] - expected) > 0.01:
            raise ValueError(f"L1 300 K reference failed at {ref:g} eV.")

    n_fermi = np.array([count_carriers(value, temp, c, value) for value in ef])
    n_band = np.array([count_carriers(value, temp, c, 0) for value in ef])
    if not np.allclose(supplied, n_band, rtol=5e-6, atol=0):
        raise ValueError("L1 supplied band counts differ from the band-edge quadrature.")
    ratios = supplied / n_fermi
    if not np.all(np.diff(ratios) > 0):
        raise ValueError("L1 correction is not increasing across the sweep.")
    n_f_op = count_carriers(ef_op, temp, c, ef_op)
    n_b_op = count_carriers(ef_op, temp, c, 0)
    kt = KB * temp
    closed = NA * c * kt * (ef_op * np.log(2) + (np.pi**2 / 12) * kt)
    if not np.isclose(n_f_op, closed, rtol=1e-5, atol=0):
        raise ValueError("L1 quadrature disagrees with the closed-form operating-point count.")
    return {
        "results": {
            "correction_factor": n_b_op / n_f_op,
            "n_fermi_op_cm2": n_f_op,
            "mobility_corrected_cm2_Vs": 1 / (CHARGE * n_f_op * rho),
        },
        "ef": ef,
        "ratios": ratios,
        "n_band_sweep_cm2": supplied,
        "n_fermi_sweep_cm2": n_fermi,
        "ef_op": ef_op,
        "temperature": temp,
        "dos_scale_c": c,
        "rho_ohm_sq": rho,
        "n_band_op": n_b_op,
        "closed_form_relative_difference": abs(n_f_op / closed - 1),
        "mobility_band": 1 / (CHARGE * n_b_op * rho),
        "checks": checks,
        "crossings": {level: bracketed_crossing(level, ef, temp, c) for level in (2, 4)},
        "band_max_rel_diff": float(np.max(np.abs(n_band / supplied - 1))),
    }


def load_pairs(path: Path):
    """Group pairs in first-appearance measurement order and sort repeat IDs within each group."""
    groups = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            index = int(row["measurement_index"])
            repeat = int(row["repeat"])
            pairs = groups.setdefault(index, {})
            if repeat in pairs:
                raise ValueError(f"Duplicate repeat {repeat} at measurement {index}.")
            pairs[repeat] = (float(row["R_H_raw"]), float(row["R_xx_raw"]))
    if not groups or any(set(pairs) != {1, 2, 3} for pairs in groups.values()):
        raise ValueError("Expected repeats 1, 2 and 3 at every measurement index.")
    repeats = np.asarray([[pairs[r] for r in (1, 2, 3)] for pairs in groups.values()])
    if not np.isfinite(repeats).all() or not np.all(repeats[:, :, 1] > 0):
        raise ValueError("Invalid Hall or sheet-resistance values.")
    return np.array(list(groups)), repeats


def normalise_pairs(pairs):
    scales = np.stack((np.max(np.abs(pairs[..., 0]), axis=-1), pairs[..., 1].max(axis=-1)), axis=-1)
    if np.any(scales <= 0):
        raise ValueError("Cannot normalise a channel with a zero or negative scale.")
    return pairs / scales[..., None, :]


def model_library(temperature):
    if temperature <= 0:
        raise ValueError("L2 temperature must be positive.")
    kt = KB * temperature
    energy = np.linspace(-EMAX, EMAX, 4000)
    ef = np.linspace(-8 * kt, 9 * kt, 400)
    gaps = np.concatenate(([0.0], np.arange(0.5, 8.01, 0.25)))
    # c, e and mobility are fixed positive scales that cancel in the normalised model curves.
    f = fermi(energy, ef[:, None], temperature)
    electrons = np.where(energy >= ef[:, None], f, 0.0)
    holes = np.where(energy <= ef[:, None], 1 - f, 0.0)
    curves = []
    for gap in gaps:
        dos = np.maximum(np.abs(energy) - gap * kt / 2, 0.0)
        n = NA * trapezoid(dos * electrons, energy, axis=1)
        p = NA * trapezoid(dos * holes, energy, axis=1)
        curves.append(np.column_stack(((p - n) / (p + n)**2, 1 / (p + n))))
    return ef, gaps, normalise_pairs(np.asarray(curves))


def nearest_model(pairs, curves, gaps):
    """Find the nearest model point for each measured pair before comparing model RMS values."""
    observed = normalise_pairs(pairs)
    hall_distance = observed[None, :, None, 0] - curves[:, None, :, 0]
    sheet_distance = observed[None, :, None, 1] - curves[:, None, :, 1]
    distance_sq = hall_distance**2 + sheet_distance**2
    nearest = distance_sq.argmin(axis=2)
    residual = np.sqrt(np.take_along_axis(distance_sq, nearest[..., None], axis=2).mean(axis=(1, 2)))
    best_gap = int(residual[1:].argmin() + 1)
    selected = best_gap if residual[best_gap] < 0.85 * residual[0] and gaps[best_gap] >= 1 else 0
    return selected, best_gap, residual, nearest[selected]


def hall_brackets(values):
    """Find adjacent nonzero readings of opposite sign, including any intervening zeros."""
    nz = np.flatnonzero(values)
    change = np.sign(values[nz[:-1]]) != np.sign(values[nz[1:]])
    return [(int(left), int(right)) for left, right in zip(nz[:-1][change], nz[1:][change])]


def analyse_l2(path: Path, temperature):
    indices, repeats = load_pairs(path)
    means = repeats.mean(axis=1)
    sd = repeats.std(axis=1, ddof=1)
    ef, gaps, curves = model_library(temperature)
    selected, best_gap, residual, nearest = nearest_model(means, curves, gaps)
    trajectory = ef[nearest]
    mean_brackets = hall_brackets(means[:, 0])
    if not mean_brackets or not hall_brackets(curves[selected, :, 0]):
        raise ValueError("The L2 means and selected model should span a Hall sign change.")

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    sample_gaps = np.empty(BOOTSTRAP_COUNT)
    sample_trajectories = np.empty((BOOTSTRAP_COUNT, len(indices)))
    rows = np.arange(len(indices))
    for sample in range(BOOTSTRAP_COUNT):
        draw = repeats[rows, rng.integers(0, 3, size=len(indices))]
        fit, _, _, assignment = nearest_model(draw, curves, gaps)
        sample_gaps[sample] = gaps[fit]
        sample_trajectories[sample] = ef[assignment]
    gap_quantiles = np.percentile(sample_gaps, [16, 50, 84])
    ef_interval = np.percentile(sample_trajectories, [16, 84], axis=0)
    return {
        "results": {
            "model": "hard_gap" if selected else "cone",
            "Eg_kT": float(gaps[selected]),
            "Eg_uncertainty_kT": gap_quantiles[[0, 2]].tolist(),
            "EF_trajectory_eV": trajectory.tolist(),
        },
        "indices": indices,
        "repeats": repeats,
        "means": means,
        "sd": sd,
        "mean_brackets": mean_brackets,
        "model_ef": ef,
        "model_curve": curves[selected],
        "cone_rms": float(residual[0]),
        "gap_rms": float(residual[best_gap]),
        "best_gap_kT": float(gaps[best_gap]),
        "bootstrap_median_kT": float(gap_quantiles[1]),
        "gap_fraction": float(np.mean(sample_gaps > 0)),
        "ef_interval": ef_interval,
        "temperature": temperature,
    }


def gradient_matrix(size):
    grid = np.arange(size**2).reshape(size, size)
    pairs = [
        (grid[:, :-1].ravel(), grid[:, 1:].ravel()),
        (grid[:-1, :].ravel(), grid[1:, :].ravel()),
    ]
    first = np.concatenate([start for start, _ in pairs])
    second = np.concatenate([end for _, end in pairs])
    rows = np.arange(len(first))
    gradient = np.zeros((len(first), size**2))
    gradient[rows, first] = -1
    gradient[rows, second] = 1
    return gradient


def reconstruction_operator(matrix, gradient, project, eta):
    """Solve the projected system together for all measurement-channel basis vectors."""
    gram = matrix.T @ matrix + eta * (gradient.T @ gradient)
    system = project @ gram @ project + 1e-9 * np.eye(matrix.shape[1])
    weights = np.linalg.lstsq(system, project @ matrix.T, rcond=None)[0]
    return weights - weights.mean(axis=0)


def point_spread(matrix, operator, pixel, centres):
    test = np.zeros(matrix.shape[1])
    test[np.ravel_multi_index(pixel, (len(centres), len(centres)))] = 1
    test -= test.mean()
    reconstructed = (operator @ (matrix @ test)).reshape(len(centres), len(centres))
    row, col = pixel
    x, y = np.meshgrid(centres, centres)
    distance_sq = (x - centres[col])**2 + (y - centres[row])**2
    weights = np.maximum(reconstructed, 0)
    sigma_sq = np.sum(weights * distance_sq) / (2 * np.sum(weights))
    return pixel, reconstructed, float(2.355 * np.sqrt(sigma_sq))


def analyse_l3(matrix_path: Path, scans_path: Path, threshold):
    size = 10
    matrix = np.load(matrix_path, allow_pickle=False)
    scans = np.loadtxt(scans_path, delimiter=",")
    if matrix.shape != (128, 100) or scans.shape != (3, 128):
        raise ValueError(f"Invalid L3 shapes: {matrix.shape}, {scans.shape}.")
    if not np.isfinite(matrix).all() or not np.isfinite(scans).all():
        raise ValueError("Non-finite L3 inputs.")
    row_sum = float(np.max(np.abs(matrix.sum(axis=1))))
    if row_sum >= 1e-11:
        raise ValueError("L3 forward matrix fails the constant-mode check.")

    mean_scan = scans.mean(axis=0)
    gradient = gradient_matrix(size)
    project = np.eye(size**2) - np.ones((size**2, size**2)) / size**2
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    eta = 1e-2 * singular_values[0]**2 / np.linalg.norm(gradient, 2)**2
    operator = reconstruction_operator(matrix, gradient, project, eta)
    flat_map = operator @ mean_scan
    flat_map -= flat_map.mean()
    factors = (0.01, 0.1, 1, 10, 100)
    sensitivity = []
    for factor in factors:
        trial = flat_map if factor == 1 else (
            reconstruction_operator(matrix, gradient, project, eta * factor) @ mean_scan
        )
        sensitivity.append((np.linalg.norm(gradient @ trial), np.linalg.norm(matrix @ trial - mean_scan)))

    centres = np.linspace(-2.7, 2.7, size)
    features = []
    # One peak per sign represents its neighbouring lobe. Coordinates are column x and row y.
    for label, peak in (("positive", flat_map.argmax()), ("negative", flat_map.argmin())):
        row, col = divmod(int(peak), size)
        features.append({
            "type": label,
            "x_probe_pitch": float(centres[col]),
            "y_probe_pitch": float(centres[row]),
            "amplitude_fraction": float(flat_map[peak]),
        })
    spread = [point_spread(matrix, operator, pixel, centres) for pixel in ((4, 4), (2, 7))]
    resolution = max(width for _, _, width in spread)
    scan_sd = scans.std(axis=0, ddof=1)
    rng = np.random.default_rng(NOISE_SEED)
    noise = rng.normal(0, scan_sd, size=(NOISE_COUNT, len(scan_sd)))
    simulated = noise @ operator.T
    detection = float(np.percentile(np.max(np.abs(simulated), axis=1), 95))

    amplitude = float(np.max(np.abs(flat_map)))
    nominal_pass = amplitude <= threshold
    confidence = "borderline"
    if amplitude + detection < threshold:
        confidence = "clear pass"
    elif amplitude - detection > threshold:
        confidence = "clear fail"
    return {
        "results": {
            "map": flat_map.tolist(),
            "features": features,
            "resolution_probe_pitch": resolution,
            "detection_limit_fraction": detection,
            "pass": nominal_pass,
            "confidence": confidence,
        },
        "centres": centres,
        "eta": float(eta),
        "threshold": threshold,
        "spread": spread,
        "factors": factors,
        "sensitivity": np.asarray(sensitivity),
        "singular_values": singular_values,
        "numerical_rank": int(np.count_nonzero(
            singular_values > max(matrix.shape) * np.finfo(matrix.dtype).eps * singular_values[0]
        )),
        "amplitude": amplitude,
        "row_sum": row_sum,
        "noise_sd_rms": float(np.sqrt(np.mean(scan_sd**2))),
    }
