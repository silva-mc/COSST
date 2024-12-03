1. Files
cosst_code: R code to run the climate-oriented seed sourcing tool (COSST). Please set up the species name on l. 27, the target restoration site coordinates on l. 30, the provenancing strategy on l. 33, the shared socioeconomic pathway on l. 36, and the timeframe on l. 39.
cosst_shiny: R Shiny App to visualize COSST
cosst_visuals: R Markdown HTML tutorial to reproduce the figures from the original publication.
seed.supp: Coodinates of the municipalities where the seed suppliers are active, please see Silva et al (2022) for more info: https://doi.org/10.3389/fevo.2022.1045591
sessionInfo: Information about the R session including the version of the packages used during the code development.

2. General notes
The R codes fetch occurrence and climatic data online but the user can use other data sources
Note that the code requires Java to be installed in the same computer operating system (32 or 64-bit): https://www.java.com/en/download/help/java_win64bit.html
We recommend the code to be run using R latest version to date (R version 4.4.1 or above)