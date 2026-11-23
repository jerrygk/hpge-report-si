# Supplementary information for Nuclear Forensics with Gamma-Ray Spectroscopy

Jerry Chen, Aika Saito, Webb Wu and Yuheng Wu (in alphabetical order of surname)

Department of Physics and Astronomy, University College London

This repository holds the data, code, results and figures behind our interim report on the PHAS0051 experiment *Nuclear Forensics with Gamma-Ray Spectroscopy* (autumn 2026). The submitted report cites version [v1.0-interim](https://github.com/jerrygk/hpge-report-si/tree/v1.0-interim), and every number, table, figure and code listing in it can be traced to that version and regenerated from the raw spectra. The current version, v1.1, corrects errors that were found after the report was submitted, and they are listed under [Bug corrections](#bug-corrections-of-24-november-2026-version-11) below.

The report asks three questions of the spectra we had by mid-November. The first is which fit gives honest peak areas and errors, the second is how those errors fall as the counting time grows, and the third is why the third-floor detector performs worse than the basement one. It also compares the Cs-137/K-40 activity ratios of two Norwegian lake sediments.

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
| `BkgRun5_Long.Chn` | third floor, ORTEC DSPEC-LF | empty shield | 3 to 7 Sep 2026 | 3.88 d | the demonstrator, William Quinn | `40e52f417ef1` |
| `BkgRun5_Long_txt.Txt` | third floor | same run, Maestro text export | | | the demonstrator, William Quinn | `5dc3205d549c` |
| `backround_run_hpge_12_10.csv` | basement, Mirion Lynx | empty shield | 12 to 15 Oct 2026 | 2.94 d | the Thursday group | `e2c0784a6e61` |
| `thurs 15 oct - calibration sample 2 .csv` | basement, Lynx | IAEA-385 reference sediment | 15 to 19 Oct 2026 | 4.05 d | the Thursday group | `b0f17c64b25b` |
| `Norway1_22_10.csv` | basement, Lynx | Norway 1 sediment | 19 to 22 Oct 2026 | 2.75 d | the Thursday group | `8a3859ab06b9` |
| `norway sample 2 - 22 oct 2026.csv` | basement, Lynx | Norway 2 sediment | 22 to 26 Oct 2026 | 3.97 d | the Thursday group | `459790ea1bea` |
| `run_hpge_Tuesday_testrun.csv` | basement, Lynx | Norway 2 sediment, second count | 27 to 29 Oct 2026 | 1.96 d | us | `6f3fc583b166` |
| `run_hpge_Tuesday_testrun 2.csv` | | second download of the same run, identical | | | us | `6f3fc583b166` |

The file names are kept exactly as the files were received, typos included, in order that they match the checksums and the laboratory records.

Our own third-floor detector was out of action from 20 October, and four of the basement spectra were therefore shared with us by the Thursday group, the laboratory group that shares the basement detector with us. They have agreed to the publication of these spectra here. The Lynx files carry no calibration, and their header start times are not reliable, as notebook 1 shows in the section "The time interval of each run, checked against an independent clock". The dates above were hence taken from the file names and the download times.

The reference tables in `data/reference/` were taken from public sources. The gamma lines, half-lives and emission probabilities come from the ENSDF evaluations through NuDat 3 (Brookhaven National Laboratory), and the photon attenuation coefficients for aluminium, germanium and lead from NIST XCOM. The crustal abundances come from Rudnick and Gao (Treatise on Geochemistry 3, 2003) and the median soil activities from the UNSCEAR 2000 Report, Annex B.

## Running the analysis

The notebooks and the figure script were run with Python 3.11 and the package versions listed in `requirements.txt`. The commands below, run from the repository root, install those packages, run the tests of the analysis module, execute the three notebooks in order and redraw the report figures.

```
pip install -r requirements.txt
python code/test_gamma_spec.py
jupyter nbconvert --to notebook --execute --inplace notebooks/1_lynx_spectra.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/2_counting_time.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/3_interim_analysis.ipynb
python scripts/make_figures.py
```

On a laptop the three notebooks take about 35 s, 100 s and 80 s. Each notebook finds the repository root by itself, which means that it can also be opened and run cell by cell in Jupyter. Since the simulations use fixed seeds, a rerun reproduces the numbers in `results/` exactly. The one exception is the run timeline in `s2_timeline.csv` and `s3_timeline.csv`, which checks the start time in each file header against the time the file was downloaded, read from its modification time. Git does not keep modification times, and after cloning that column therefore reads "unknown". The timeline files shipped here hold the values from the original downloads.

The first notebook, `1_lynx_spectra.ipynb`, reads the four basement spectra and calibrates them without a stored calibration. It then measures the resolution and the peak tail of the Lynx detector and compares our fits with those of InterSpec, writing its tables to `results/session2/`.

The second notebook, `2_counting_time.ipynb`, identifies our own run as a second count of the Norway 2 pot. It compares the 1.96-day and 3.97-day counts through their errors, detection limits, radon change and activity ratios, and writes its tables to `results/session3/`.

The third notebook, `3_interim_analysis.ipynb`, holds the analysis behind the report and writes its tables to `results/interim/`. It fits the K-40 peak from first principles on both detectors and estimates the statistical and systematic errors of that fit, which are then tested on simulated peaks. It goes on to the calibration and resolution, every background line, the comparison of the two backgrounds and the choice of fit statistic and peak shape, and then compares the fits with InterSpec. The later sections summarise the counting-time results and cover the plans, the projected significance, the lead estimate and the escape peaks.

## Where each part of the report comes from

The figure, table, appendix and section numbers below are those of the submitted report, and they hold for both versions of this repository. The code listings of Appendix C.1 are copied from [v1.0-interim](https://github.com/jerrygk/hpge-report-si/tree/v1.0-interim). In v1.1 they differ only in the last line of `fit_peak`, which now also returns the number of channels where the model fell to zero or below, and the functions and notebook sections named in Table C1 are unchanged.

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
| Appendix C.1, peak model and fit | 3 | `code/gamma_spec.py`, `fit_peak` in Section 3 |
| Table C1, the remaining code | 2 to 23 | `code/gamma_spec.py` and the sections named in the table |

The figure files are in `figures/report/` as PDF and PNG, and the result files are in `results/interim/` unless another folder is given. In v1.1 the figures are drawn from the corrected results, and the figures exactly as submitted are those of v1.0-interim. Notebook 3 also holds analyses that the report does not use, among them a radon ingrowth table, the split of counting time between sample and background, self-absorption estimates and a limit on neutrinoless double beta decay. Its planning section was written before the final schedule was agreed, and the schedule in the report is the one to follow. Section 14 compares the K-40 rate from the single-peak fit with the rate from the line table, which appear in the report as Tables A3 and A5.

## Bug corrections of 24 November 2026 (version 1.1)

Version 1.1 corrects the following errors in version 1.0-interim. The corrections were made in our Tuesday laboratory session of 24 November 2026. Every number below comes from rerunning notebooks 2 and 3, and the figures in `figures/report` were redrawn from the corrected results with the layout and wording of the interim report.

**Projected significance** (notebook 3, Section 22). The error of the combined Norway 2 rate now follows the live time weights of the rate, $s^2 = w_1^2 s_1^2 + w_2^2 s_2^2$, from the two fitted errors, in place of the Poisson error of the combined counts. The difference after the planned runs changes from 6.6 to 6.3 standard deviations with the background drawn separately for each sediment, and from 7.6 to 7.1 with one background for both.

**Detection limits** (notebook 2, Section B4). The region of interest is now centred on the top of the tailed peak, through the new function `emg_window_fraction` in `gamma_spec.py` and its test. The limits fall by 1.2 to 1.7 %.

**Activity ratios** (notebook 2, Section A7). The error of the background run, which is subtracted from both the sample and IAEA-385, now enters once. Cs-137/K-40 is $0.0462 \pm 0.0040$ for Norway 1 and $0.0277 \pm 0.0024$ for Norway 2, in place of errors of 0.0044 and 0.0027.

**Counting check of the fitted errors** (notebook 3, Section 19). The peak is now removed from its window with the fraction of the tailed peak, about 98 %, in place of the Gaussian 98.8 %. The strong lines now exceed the upper bound by 0.4 to 0.7 %, in place of about 1.5 %.

**Tests.** Section 5 of notebook 3 now also refits simulated tailed peaks with the tailed model, which gives area pull widths of 1.01. A simulated goodness of fit value with no exceedance in 400 simulations is now reported as $p < 0.0075$ at 95 % confidence, with the number of exceedances in `n_exceed`. Every fit now reports the number of channels where its model is zero or below, in `n_nonpositive`, and none of the 32 line fits has any.

**Notebook text.** The notebooks now give the fall of the total rate in the recount of Norway 2, $1.45 \pm 0.19$ %, and the ratio of the two backgrounds over all recorded channels, 3.0, beside 3.8 between 40 and 2700 keV.

## Contact

Jerry Chen, Department of Physics and Astronomy, University College London.
