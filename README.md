![](./www/logo.svg)

### A web service for Visually supervised protein Inference and protein Quantification.

### What is VIQoR?

VIQoR, is a user-friendly web service for Visually supervised protein Inference and protein Quantification implemented in R. The Shiny web interface integrates the post-identification processes involved in protein inference and relative protein abundance summarization, along with smart and novel interactive visualization modules to support the common researchers with a straight-forward tool for protein quantification, data browsing and data inspection.

#### Input

The input should be a peptide/PSM quantitative report in .csv format and a protein sequence database in .fasta format. The imported reports should be already filtered according to the user's preferences and based on the experiment type (FDR, number of missed cleavages etc.) but should not be preprocessed regarding protein inference. Modifications can be also included. Some of the software that can identify spectra, extract/assign intensities and export peptide/PSM reports are: [SearchGUI](https://compomics.github.io/projects/searchgui) or [MaxQuant](https://www.maxquant.org/).

#### Preprocessing and filtering

- PSM to peptide aggregation by intensity summation.
- Missing value filtering.
- Filtering of post-translational modifications (PTMs).
- Intensity Log2 transformation.
- Per sample zero center normalization (median, average, quantile).

#### Protein Inference

- Two algorithms infer protein groups by "strict" or "soft" parsimony.
- Protein inference by direct peptide to protein mapping is available too. (could support cases when input data is a result of a prior inference process).
- The 2-peptide rule per protein group is applied to reduce the protein identification FDR.

#### Protein Quantification

- Protein quantification by a factor analysis coupled to a weighted average summarization function.
- Factor analysis parameter is optimized by the Global Correlation Index (GCI).
- Signal-to-noise ratios to reduce the FQR.
- Additional summarization methods: weighted or total summation of peptide intensities corresponding to a protein group.

#### Interactive Data Inspection and Visualization features

- Interactive tables.
- Connected component graph of inferred protein groups.
- Lineplot of protein group quantitative profiles.
- Heatmaps of protein expressions and peptide abundances.
- VIQoR plot.

More information about the functionality of the tool can be found in the [Manual](https://github.com/vtsiamis88/viqor/blob/main/www/Manual.pdf).


### Installation

VIQoR works only for R version 4 or higher.

Download VIQoR's repository to your computer and install the following R packages:
```R
install.packages(c("shiny", "shinydashboard", "DT", "protr", "shinyjs", "igraph", "networkD3", 
                    "plotly", "shinyBS", "shinycssloaders", "ggplot2", "reshape2", "remotes",
                    "webshot2", "htmlwidgets", "V8", "seqinr", "dplyr", "seqinr", "heatmaply"))
                    
remotes::install_github("rstudio/webshot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "preprocessCore"))
                    
```

The dependencies installation is also included in app.R.

You also need to have the Chrome browser installed on your system or any other browser based on Chromium.

Load the app.R file in [Rstudio](http://rstudio.com) and run the shiny app.

### Docker container

The latest version of VIQoR is available as docker container veitveit/viqor:latest

Assuming that you have a working docker environment, run it via
```
docker run -it -p3838:3838 veitveit/viqor
```
and access the tool in your web browser via `localhost:3838` (the address might change when running docker in a virtual machine).


### Web-service

You can access the online service by any web browser (preferably Chrome or Mozilla) [here](https://computproteomics.bmb.sdu.dk/app_direct/VIQoR/).

### Contact
For software issues and general questions, please submit an issue.

### License
Apache-2.0