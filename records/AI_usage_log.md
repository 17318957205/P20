# G20 AI usage record

Tool used: OpenAI Codex. The requester provided the QXU6028 methods guide
and the two G20 archives and then asked Codex to start the scoped work
(“开始做”). Codex implemented the analysis, generated the checkpoint
material and checked the results.

The methods below were the implementation briefs supplied to Codex by the
assignment guide, especially its starter prompts on page 15.

| Stage | Guide instruction used as an implementation brief | Checks run |
|---|---|---|
| L1 | Integrate the linear DoS from 0 and from the exact Fermi level, then use the Fermi-referenced count for the mobility. | Tested the four 300 K reference points. Compared the supplied band-edge sweep to quadrature and checked the exact operating-point counts with adaptive integration. |
| L2 | Fit the signed Hall and sheet-resistance pairs to separately normalised model curves by nearest distance. Apply the gap rule and resample repeats. | Confirmed three paired repeats per index and checked the Hall sign change. Recomputed the residuals and nearest assignments with SciPy distances, then checked the 200 draws and the trajectory intervals. |
| L3 | Project out the constant mode, use the specified regularisation and reconstruct both test impulses and noise-only scans. | Checked input hashes, dimensions, row sums and map mean. Recomputed the map, two widths and 95th percentile noise limit using an independent zero-mean basis. |

The regular analysis was reviewed and committed before the blind numerical
run. Both runs used tag `g20-pipeline` and commit
`686d0010d58429684a9d37c9b245994d816bcdf4`. The two terminal logs record
the environment and results. See the error ledger for corrections.

## Final coursework preparation

The requester supplied the full assessment brief after the checkpoint
packages were completed. Codex compared Sections 5 and 7 with the existing
deliverables, added the L3 spectral diagnostic and supplementary report
metrics, reran the regular inputs from a clean Git checkout, and drafted
the short report, critical self-assessment and blank contribution form.
Codex also checked the rendered Word pages and prepared the local Git
bundle. The final analytical commit is the `g20-coursework` tag. The previous
`g20-pipeline` tag still identifies the checkpoint and blind-run evidence.

For the Transfer section, Codex checked the methods guide's Chapter 5
pointer and verified the P3HT transistor example against the source paper
by Sirringhaus et al. The group should confirm the relevant QXU6028
lecture slide and the device covered there. Group members will complete
their personal details, evidence and signatures themselves.

## Revision after receiving the course report guidance

Codex checked slides 6–10 and 14 against the report. It added the L1
band-edge and Fermi-referenced density sweeps and the separate L2 channel
overlays to the reproducible figures. It revised Transfer to use the
Chapter 3 mobility derivation identified on slide 10 and verified the
inkjet-printed P3HT device against Park and Wong (2022,
doi: 10.1021/acsaelm.1c00985). This resolved the lecture-reference note in
the earlier draft. Codex also checked the rendered report and independently
recomputed the marked results.

The revised report is reproduced from tag `g20-coursework-v2`. The previous
report and the phase and blind-run records remain traceable in Git.
The archived blind run was performed locally with the original frozen
code, and no QMplus receipt was supplied.
