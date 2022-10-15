BCB
================
Jireh Huang
(<jirehhuang@ucla.edu>)

# Bayesian Causal Bandits with Backdoor Adjustment Prior

Jireh Huang and Qing Zhou

## Installation

On R version 3.6.0 or more recent, first install the
[phsl](https://github.com/jirehhuang/phsl) package according to the
instructions provided by the author. Unfortunately, some functionality
of the bcb package is not compatible with Windows due to the integration
of [BIDA](https://github.com/jopensar/BIDA) and
[modular-dag-sampling](https://github.com/ttalvitie/modular-dag-sampling),
copies of which are included in the `inst` folder with permission from
the respective authors.

Then, run the following code in R to install this package from GitHub
with its dependencies.

``` r
devtools::install_github("jirehhuang/bcb", dependencies = TRUE)
```

Alternatively, download the source package from [\[Google
Drive\]](https://drive.google.com/drive/folders/10F-s6wyE_RQ3TKUnDe2y3r-9nVHCB6Ml)
and run the following code in R to install the package from source.

``` r
install.packages("bcb_1.0.tar.gz", repos = NULL, dependencies = TRUE, type = "source")
```

## Reproducing Numerical Results

This section contains instructions for reproducing the experimental
results in the paper. The following instructions assume that the current
working directory is writable and contains the scripts in
`inst/scripts`.

``` r
setwd("inst/scripts")
```

### Main Experiments

The main experiments are described and presented in Section 6 and
Appendix C.

#### 1. Generating Networks and Data

``` r
source("generate_dkpar.R")
```

`generate_dkpar.R` creates two folders, randomly generating 100 discrete
Bayesian networks in `dkpar_3_6_3-d_0` and 100 Gaussian Bayesian
networks in `dkpar_3_6_3-g_0`, as well as 10 observational datasets for
each network. The script executed on a i5-9600k using a single CPU core
in about 40 minutes. The generated networks and data are available at
[\[Google
Drive\]](https://drive.google.com/drive/folders/1aAtAxcLeztgWNFPy2KceRJ5V9gUt79u8).

#### 2. Executing Algorithms

Remove `_0` from the generated directories from `generate_dkpar.R`,
resulting in folder names `dkpar_3_6_3-d` and `dkpar_3_6_3-g`.

``` r
source("run_dkpar.R")
```

`run_dkpar.R` executes 44 algorithms on the networks and datasets
`dkpar_3_6_3-{d, g}`, creating a folder for each algorithm and saving
the results for each execution. The script executed in about 7045 days
of single thread CPU time using the UCLA Hoffman2 computing cluster,
creating over 400 GB of output. Since these results are not practical to
store, only the compiled results are provided, discussed in the
following step.

#### 3. Compiling Results

``` r
source("compile_dkpar.R")
```

`compile_dkpar.R` creates the `concise` folder in `dkpar_3_6_3-{d, g}`
and, for each algorithm, compiles and saves the results of all
executions. Then, it averages the results of each algorithm, normalizing
where appropriate, and combines all results into `df_2.rds`. The script
executed in about 3 hours using 4 CPU cores on the UCLA Hoffman2
computing cluster. The compiled results are available at [\[Google
Drive\]](https://drive.google.com/drive/folders/1VP-WoJ5wDQM4LjOv_XZvXIcFvteezTII).

#### 4. Analyzing Results

Append `_done` to the compiled directories from `compile_dkpar.R`,
resulting in folder names `dkpar_3_6_3-d_done` and `dkpar_3_6_3-g_done`.
The compiled results for each algorithm are not necessary – only
`df_2.rds` and `method_grid.txt` in the `concise` folder.

``` r
source("analyze_dkpar.R")
```

`analyze_dkpar.R` creates the cumulative regret
(`cumulative.{eps, png}`), head start for competing algorithms in the
discrete setting (`hs-d.{eps, png}`), and edge support sum of absolute
errors (`essae.{eps, png}`) figures. The analysis results contain the
necessary files to execute `analyze_dkpar.R` and are available at
[\[Google
Drive\]](https://drive.google.com/drive/folders/1BNCybuaKkQZkNtUkVc67KTqkz5FbE1eH)
and are most manageable in terms of file size.

### Additional Experiments

The additional experiments evaluating the proposed backdoor adjustment
methodology are provided in Appendix D.

``` r
source("test_bda.R")
source("analyze_bda.R")
```

`test_bda.R` creates folders `test_bda-{d, g}`, investigating 24000
randomly generated scenarios for each distributional setting. In each
folder, an RDS file is compiled containing the simulation results. The
script executed in about 167 days of single thread CPU time using the
UCLA Hoffman2 computing cluster. `analyze_bda.R` creates the backdoor
adjustment coverage probabilities for the discrete
(`bda_coverage-d.{eps, png}`) and Gaussian (`bda_coverage-g.{eps, png}`)
simulations. The results are available at [\[Google
Drive\]](https://drive.google.com/drive/folders/1BvlbSjHGmVEo4wim4tHMzCj4fVBf2gOA).
