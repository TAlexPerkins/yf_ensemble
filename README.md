# yf_ensemble

This repository contains code to reproduce results from the following preprint.

TA Perkins, JH Huber, QM Tran, RJ Oidtman, MK Walters, AS Siraj, SM Moore (2021) **Burden is in the eye of the beholder: Sensitivity of yellow fever disease burden estimates to modeling assumptions**. *medRxiv* doi:[10.1101/2021.01.06.21249311](https://www.medrxiv.org/content/10.1101/2021.01.06.21249311v1)

### License

This repository is being released under the MIT License.

### Organization

Contents of this repository are organized according to code, data, output, and figures. The data folder contains inputs to the analysis, and the output folder contains outputs of the analysis steps. This corresponds to the left and right columns, respectively, of Figure 1 in the preprint.

### Code

Within the code folder, file names indicate the step in the analysis that the code pertains to, consistent with the steps enumerated in the preprint. The only exception is prior_iceberg.R, which pertains to parameters for a prior distribution described in the Supplemental Appendix. This script was run before Step 1.

All code for the analysis was written in the R language (version 3.5.2) and was executed in part on a MacBook Pro (macOS 10.13.6) and in part on a Linux cluster. At the time the research was done, all R packages used in this analysis were available on CRAN and were straightforward to install and load. Some multi-panel figures were assembled using Microsoft Powerpoint.

### Data

Estimates of population and vaccination coverage, both stratified by age, year, and first administrative level, come from a study by Hamlet et al., the results of which can be accessed at https://shiny.dide.imperial.ac.uk/polici/

Serological data derive from published estimates we obtained from the literature. The file containing these data (data/yf_sero_data_with_coverage.csv) contains all data gathered by our search, with some data being filtered out prior to analysis. Filtering occurred if data fell outside our geographic scope (34 countries), outside our temporal scope (1980-2014), or did not result from a neutralization assay.

Data on reported cases and deaths were compiled from publicly available sources and shared with us by Tini Garske from Imperial College. While our analysis makes use of a version of those records shared with us in August 2018 (data/outbreaks_1969_2014.csv), an updated version of those data is now available at https://github.com/kjean/YF_outbreak_PMVC/blob/main/formatted_data/outbreaks_1980s-2018.csv
