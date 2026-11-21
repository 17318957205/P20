# QXU6028 virtual lab · Group G20

This repository contains the L1–L3 Python analysis, the original regular group
data, the generated checkpoints and the run records. The `g20-coursework-v3`
tag identifies the version that regenerates the final report figures and
numbers. The `g20-coursework` and `g20-coursework-v2` tags preserve earlier report runs. The
`g20-pipeline` tag preserves the frozen version recorded
in the phase checkpoints and the local blind run. The blind input ZIP is held
separately so it can be extracted into a clean directory.

## Reproduce a run

Install Git and Python 3.10 or later. Create a virtual environment if needed and
install the dependencies once with `python -m pip install -r requirements.txt`.
Run these commands from the repository root. Each output path must be new.

```text
git switch --detach g20-coursework-v3
python reproduce.py --data data/regular_G20 --output runs/coursework
```

The pipeline verifies the input hashes and writes `results.json`, the commit
hash, the eight report figures (including the separate L1 density and mobility sweeps, the two L2
channel overlays and the L3 singular-value spectrum),
three brief notes, the L2 per-index data table and `coursework_metrics.json`.
The latter records the report diagnostics separately from the marked schema. The
default L3 specification is 0.030, as stated in `README_FIRST.txt` in the
supplied blind archive.

For the blind run, extract only `blind_G20` into a clean input directory.
Then supply that directory in the same command and save its terminal output.
For example, when the extracted directory is next to this repository, run:

```text
git switch --detach g20-pipeline
python reproduce.py --data ../clean_input/blind_G20 --output runs/blind > run_log.txt 2>&1
```

The blind output directory contains `blind_results.json` and
`pipeline_commit.txt`. Upload them with `run_log.txt` as three separate files
in the course platform's live window. Check the command exit status first.

## Result conventions

- L2 model names are `cone` and `hard_gap`. The gap interval in
  `Eg_uncertainty_kT` is [16th percentile, 84th percentile] across selected
  bootstrap models. The bootstrap median and hard-gap fraction are in the L2
  note. `L2_trajectory_and_repeats.csv` gives the per-index interval and sample
  standard deviations.
- The L3 map is flat, row-major and zero-mean. Column index sets x and row index
  sets y, both increasing from −2.7 to +2.7 probe pitches. Each feature is
  represented by its signed extremum and coordinates. Map entries, feature
  amplitudes and the detection limit are fractions of reference conductivity.
- The L3 `pass` Boolean uses the nominal 3% specification. The confidence
  string adds the detection-limit margin from the guide, so a borderline
  result may have either nominal Boolean value.
- The same reference regularisation, point-spread inputs and random seeds are
  used in both runs. The terminal log records the run version and parameters.
