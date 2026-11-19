# Nuclear forensics with gamma-ray spectroscopy: supplementary information

Jerry Chen, Webb Wu, Aika Saito and Yuheng Wu
Department of Physics and Astronomy, University College London

This repository holds the data, code, results and figures behind our interim report on the PHAS0051 experiment *Nuclear Forensics with Gamma-Ray Spectroscopy* (autumn 2026). Every number, table and figure in the report can be traced to a file here and regenerated from the raw spectra.

The report asks three questions of the spectra we had by mid-November: which fit gives honest peak areas and errors, how the errors fall with counting time, and why the third-floor detector performs worse than the basement one. It also compares the Cs-137/K-40 activity ratios of two Norwegian lake sediments.

## Contents

```
data/raw/          the six spectra as recorded (plus a second export of our own run and the Maestro text printout)
data/reference/    nuclear data and attenuation coefficients used by the analysis
data/interspec/    spectra exported for InterSpec and the peak tables read back from it
code/              the analysis module gamma_spec.py, its tests, and hpge_analysis.py (used only for a comparison in notebook 1)
notebooks/         the three analysis notebooks, executed, in the order they should be run
results/           every table the notebooks write, one folder per notebook (session1 holds two earlier inputs)
figures/report/    the figures of the report, named by their figure number
figures/session2/, figures/session3/, figures/interim/
                   the working figures each notebook draws as it runs
scripts/           make_figures.py, which draws the report figures from notebook 3
```

## Data

| File | System | Sample | Recorded | Live time | Recorded by | SHA-256 (first 12) |
|---|---|---|---|---|---|---|
| `BkgRun5_Long.Chn` | third floor, ORTEC DSPEC-LF | empty shield | 3 to 7 Sep 2026 | 3.88 d | the demonstrator | `40e52f417ef1` |
| `BkgRun5_Long_txt.Txt` | third floor | same run, Maestro text export | | | the demonstrator | `5dc3205d549c` |
| `backround_run_hpge_12_10.csv` | basement, Mirion Lynx | empty shield | 12 to 15 Oct 2026 | 2.94 d | another laboratory group | `e2c0784a6e61` |
| `thurs 15 oct - calibration sample 2 .csv` | basement, Lynx | IAEA-385 reference sediment | 15 to 19 Oct 2026 | 4.05 d | another laboratory group | `b0f17c64b25b` |
| `Norway1_22_10.csv` | basement, Lynx | Norway 1 sediment | 19 to 22 Oct 2026 | 2.75 d | another laboratory group | `8a3859ab06b9` |
| `norway sample 2 - 22 oct 2026.csv` | basement, Lynx | Norway 2 sediment | 22 to 26 Oct 2026 | 3.97 d | another laboratory group | `459790ea1bea` |
| `run_hpge_Tuesday_testrun.csv` | basement, Lynx | Norway 2 sediment, second count | 27 to 29 Oct 2026 | 1.96 d | us | `6f3fc583b166` |
| `run_hpge_Tuesday_testrun 2.csv` | | second download of the same run, identical | | | us | `6f3fc583b166` |

The file names are kept exactly as the files were received, typos included, so that they match the checksums and the laboratory records.

Our own third-floor detector was out of action from 20 October, so four of the basement spectra were shared with us by another laboratory group. Their names will be added here with their permission, which should also be confirmed before this repository is made public. The Lynx files carry no calibration and their header start times are not reliable (notebook 1, "The time interval of each run, checked against an independent clock"), so the dates above come from the file names and download times.

The reference tables in `data/reference/` are taken from public sources: gamma lines, half-lives and emission probabilities from the ENSDF evaluations through NuDat 3 (Brookhaven National Laboratory), photon attenuation coefficients for Al, Ge and Pb from NIST XCOM, crustal abundances from Rudnick and Gao (Treatise on Geochemistry 3, 2003) and median soil activities from the UNSCEAR 2000 Report, Annex B.

## Running the analysis

The notebooks were run with Python 3.12 and the package versions in `requirements.txt`. From the repository root:

```
pip install -r requirements.txt
python code/test_gamma_spec.py
jupyter nbconvert --to notebook --execute --inplace notebooks/1_lynx_spectra.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/2_counting_time.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/3_interim_analysis.ipynb
python scripts/make_figures.py
```

On a laptop the three notebooks take about 35 s, 100 s and 80 s. Each notebook finds the repository root by itself, so it can also be opened and run cell by cell in Jupyter. The simulations use fixed seeds, so a rerun reproduces the numbers in `results/` exactly. The one exception is the run timeline (`s2_timeline.csv`, `s3_timeline.csv`), which checks each file's header start time against its download time, read from the file's modification time. Git does not keep modification times, so after cloning that column reads "unknown". The shipped timeline files hold the values from the original downloads.

1. `1_lynx_spectra.ipynb` reads the four basement spectra, calibrates them without a stored calibration, measures the resolution and peak tail of the Lynx detector, and compares our fits with InterSpec. It writes `results/session2/`.
2. `2_counting_time.ipynb` identifies our own run as a second count of the Norway 2 pot and compares the 1.96-day and 3.97-day counts: errors, detection limits, the radon change and the activity ratios. It writes `results/session3/`.
3. `3_interim_analysis.ipynb` is the analysis behind the report: the K-40 peak fitted from first principles on both detectors, its statistical and systematic errors, simulated peaks, calibration and resolution, every background line, the comparison of the two backgrounds, the choice of fit statistic and peak shape, InterSpec, the counting-time summary, the plans, the projected significance, the lead estimate and the escape peaks. It writes `results/interim/`.

## Where each part of the report comes from

| Report | Notebook 3 section | Files |
|---|---|---|
| Fig. 1, decay series | | `figures/report/tikz/fig1_decay_chains.tex` |
| Fig. 2, background spectrum and photon interactions | 1, 24 | `fig2_background_spectrum`, `results/interim/ir_escape_peaks.csv` |
| Fig. 3, set-up and analysis chain | | `figures/report/tikz/fig3_setup_and_analysis_chain.tex` |
| Fig. 4, pull widths of simulated peaks | 5 | `fig4_pull_widths`, `ir_k40_toys.csv` |
| Fig. 5, choice of peak shape | 5, 12 | `fig5_k40_peak_shape`, `ir_k40_tilt_models.csv` |
| Fig. 6, plan and projected significance | 22 | `fig6_plan`, `ir_projection.csv`, `ir_projection_curve.csv` |
| Table 1, key results | 3 to 10, 19 | `ir_k40_gaussian_fits.csv`, `ir_k40_tilt_models.csv`, `ir_resolution.csv`, `ir_background_bands.csv`, `ir_cs137_background.csv`, `results/session3/s3_precision.csv`, `s3_agreement_lines.csv`, `s3_activity_ratios.csv` |
| Table 2, remaining weeks | 17, 22 | `ir_projection.csv` |
| Table A1, the two systems | 11 | `ir_file_formats.csv`, `ir_bin_normalisation.csv` |
| Table A2, the six spectra | | `results/session3/s3_timeline.csv` |
| Tables A3 and A4, K-40 fits and peak models | 3, 4, 12 | `ir_k40_gaussian_fits.csv`, `ir_k40_systematics.csv`, `ir_k40_tilt_models.csv` |
| Table A5, background lines | 7, 21 | `ir_background_lines.csv`, `ir_background_lines_simp.csv` |
| Fig. A1, basement K-40 peak | 3, 5 | `figA1_lynx_k40` |
| Fig. A2, fit statistic | 13 | `figA2_fit_statistic`, `ir_statistic_toys.csv`, `ir_statistic_on_data.csv` |
| Fig. A3, calibration and resolution | 6 | `figA3_calibration_resolution`, `ir_calibrations.json`, `ir_resolution_core.csv` |
| Fig. A4, the two backgrounds | 7, 8, 23 | `figA4_backgrounds`, `ir_background_lines.csv`, `ir_lead_thickness.csv` |
| Fig. A5, the two sediments | 10 | `figA5_sediments`, `results/session3/` |
| Fig. C1, counting check of the errors | 19 | `figC1_counting_check`, `ir_error_prediction.csv` |
| Fig. C2, fit window | 3, 4 | `figC2_window_scan`, `ir_k40_window_scan.csv`, `ir_k40_window_choice.csv` |

Figure files are in `figures/report/` as PDF and PNG, and result files in `results/interim/` unless another folder is given. Notebook 3 also holds analyses that the report does not use, among them a radon ingrowth table, the split of counting time between sample and background, self-absorption estimates and a limit on neutrinoless double beta decay, and its planning section predates the final schedule. The schedule in the report is the one to follow.

## Contact

Jerry Chen, Department of Physics and Astronomy, University College London.
