# Reference data

These tables hold nuclear and material constants, not measurements. The analysis module reads them from this folder.

`gamma_lines_ENSDF.csv` lists 152 gamma lines of the nuclides expected in natural, anthropogenic and cosmogenic backgrounds, with energy, emission probability, half-life and the Nuclear Data Sheets evaluation each comes from (ENSDF, through NuDat 3, Brookhaven National Laboratory). The uncertainty columns are in units of the last quoted digit. As the energies were stored as floating-point numbers, trailing zeros were dropped, so where an energy uncertainty matters the notebooks type the value in from the evaluation.

`primordial_nuclides.csv` gives half-lives and isotopic abundances (ENSDF, as distributed in InterSpec's `sandia.decay.xml`, Sandia National Laboratories), crustal concentrations (R. L. Rudnick and S. Gao, Treatise on Geochemistry 3, 1-64, 2003), median soil activities (UNSCEAR 2000 Report, Annex B) and the strongest useful gamma line of each parent, with the branching to the emitting daughter included.

`xcom_Z13.csv`, `xcom_Z32.csv` and `xcom_Z82.csv` are photon mass attenuation coefficients in cm2/g for aluminium, germanium and lead from NIST XCOM: coherent and incoherent scattering, photoelectric absorption, pair production in the nuclear and electron fields, and the totals. An energy that appears twice marks an absorption edge, the first row lying just below it.

`alpha_emitters_ENSDF.csv` lists the alpha emitters of the uranium and thorium series and Am-241 with half-life, alpha branching and Q-value (ENSDF).

`non_nuclear_features.csv` lists the non-nuclear features of a germanium spectrum, such as X-ray lines, escape peaks and the annihilation line, used to name peaks that are not gamma lines.
