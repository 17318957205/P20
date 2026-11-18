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
