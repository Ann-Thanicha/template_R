

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Test R repro 

<!-- badges: start -->


[![CRAN_Status_Badge](https://img.shields.io/badge/GitBook-%23000000.svg?style=for-the-badge&logo=gitbook&logoColor=white)](https://github.com/Ann-Thanicha/template_R)

# Welcome to my unruly R project
# Be prepare if you want to use my code
<a href="https://amices.org/mice/"><img src="docs/tvcnd6e86cx51.jpg" align="center" height="300" /></a>

# Are you ready?
## First step: Installation required package

All the required packages can be installed as follows:

``` r
## First specify the packages of interest
packages = c("ggplot2", "sf", "sp","raster", "readxl",
             "dplyr", "RandomFields","RColorBrewer","ggridges", "cowplot", 
             "lintr", "styler", "docstring")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
```



## Project Structure

The project structure distinguishes three kinds of folders:
- read-only (RO): not edited by either code or researcher
- human-writeable (HW): edited by the researcher only.
- project-generated (PG): folders generated when running the code; these folders can be deleted or emptied and will be completely reconstituted as the project is run.


```
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── requirements.txt
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── R                  <- Source code for this project (HW)

```

## Add a citation file
Create a citation file for your repository using [cffinit](https://citation-file-format.github.io/cff-initializer-javascript/#/)

## License

This project is licensed under the terms of the [MIT License](/LICENSE).
