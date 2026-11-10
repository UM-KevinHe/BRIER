
# Description

`BRIER` is a comprehensive Bregman-divergence based integration framework for **improving genetic risk prediction** of complex traits in a **local target** cohort 
by leveraging information **from external sources** under diverse data-availability settings. 
It flexibly accommodates a broad range of outcome types, including continuous, binary, count, and time-to-event traits. 
The framework has three modules, `BRIER.FULL`, `BRIER.I`, and `BRIER.S`, fitted by the functions
`BRIERfull()`, `BRIERi()`, and `BRIERs()`, which address the major data-access constraints encountered in biomedical and statistical-genetics research.

# Using BRIER with an AI assistant?

We publish an AI-optimized reference at **<https://um-kevinhe.github.io/BRIER/llms.txt>** following the [llms.txt convention](https://llmstxt.org/). It gives the assistant a compact map of the package: decision tree, parameter semantics, worked examples drawn from the verified vignettes, and a list of common pitfalls, without forcing it to ingest the full pkgdown site. Point your assistant at the URL or paste the contents into the chat at the start of a session.

# Models

<table style="
    width:100%;
    text-align:center;
    border-collapse:collapse;
    border-top:2px solid black;        /* top rule */
    border-bottom:2px solid black;     /* bottom rule */
">

<tr style="border-bottom:2px solid black;">   
    <th style="text-align:center;">Model</th>
    <th style="text-align:center;">Target Local Data</th>
    <th style="text-align:center;">External Data</th>
    <th style="text-align:center;">Application Examples</th>
</tr>

<tr style="border-bottom:1px solid black;">
    <td><code>BRIER.FULL</code><br><code>BRIERfull()</code></td>    
    <td rowspan="2">
        Individual-level data <br>
        (e.g. outcome y and design matrix X)
    </td>
    <td>Individual-level data</td>
    <td rowspan="2">
        Modeling the polygenic risk of cholesterol level (continuous) [1] <br>
        Predicting drug response (binary) using gene expression profiles [2] <br>
        Constructing genetic prognostic models for survival data (time-to-event) [3]
    </td>
</tr>

<tr style="border-bottom:1px solid black;">  
    <td><code>BRIER.I</code><br><code>BRIERi()</code></td>
    <td>
        Pretrained models <br>
        (e.g., penalized regression coefficients)
    </td>
</tr>

<tr style="border-bottom:1px solid black;">   
    <td><code>BRIER.S</code><br><code>BRIERs()</code></td>
    <td>
        Summary-level data <br>
        (e.g., GWAS summary statistics, local LD matrix)
    </td>
    <td>Summary-level data or pretrained models</td>
    <td>
        Constructing PRS from GWAS summary statistics [4] <br>
        Building GReX models from eQTL summary statistics [5]
    </td>
</tr>

</table>


# Features

`BRIER` improves risk prediction in a target cohort by borrowing from external information, with an integration weight $\eta$ that sets how much to borrow. Predictors can be SNPs, gene expression, proteins or other features.

**Integrating external information**

- **Three modules for different data access and privacy constraints**: `BRIERfull` pools raw external cohorts; `BRIERi` integrates pretrained external models with individual-level target data; `BRIERs` needs only target summary statistics and an LD matrix.
- **Multiple external models**: one integration weight per source (`multi.method = "ind"`) or a weighted ensemble (`"stacking"`, `"stacking.c"`), with tools to screen, de-duplicate and reduce a large panel of external models (`screenExternals()`, `dedupExternals()`, `reduceExternalsPCA()`).

**Modeling and tuning**

- **Distribution families** for complex traits: Gaussian, binomial and Poisson.
- **Penalties**: LASSO, MCP and SCAD, the elastic net, and unpenalized covariates.
- **Hyperparameter tuning of $(\eta, \lambda)$**: on a validation set (MSPE, $R^2$, deviance, AUC), by cross-validation, by Bayesian optimization over $\eta$, or without held-out data by information criteria ($C_p$ for `BRIERi` and GIC for `BRIERs` came close to validation-set tuning in simulations).

**Practical genetic applications**

- **Computational efficiency**: a C++ coordinate-descent core with sequential strong-rule screening, and parallel fitting across $\eta$.
- **Genotype input**: genotype dosages, including imputed data; allele alignment between target and external data (`preprocessI()`, `preprocessS()`); GWAS p-values to correlations (`p2cor()`).
- **Public LD panels**: block LD from 1000 Genomes and UK Biobank for African, European and East Asian ancestry, on Zenodo and ready for `BRIERs`; `convertLDref()` also reads the PRS-CS panels directly, and `calLD()` with `getLDB()` computes LD from your own reference genotypes.
- **Downstream analysis**: `coef()`, `predict()` and `evalMetric()`; test performance along $\eta$ (`plot.eta()`); bootstrap comparisons with target-only and external-only models (`plot.box()`); and variable importance from bootstrap selection frequencies (`plot.importance()`).

# Installation

> **Note:** This package is currently under active development. Please report any issues or unexpected behavior.  
> Requires **R ≥ 4.0.0**.

Install from CRAN:

    install.packages("BRIER")

Or install the development version from GitHub:

    require("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/BRIER")

<br>


# Detailed Tutorial

Full package documentation and parameter explanations: [here](https://um-kevinhe.github.io/BRIER/)


# Real Data Examples

We provide two analyses on public data to demonstrate the usage of BRIER:

- [Modeling tumour mutational burden using gene expression from TCGA](https://um-kevinhe.github.io/BRIER/articles/Example_TCGA.html): `BRIERfull`, `BRIERi` and `BRIERs` predict tumour mutational burden in TCGA lung squamous cell carcinoma, borrowing from lung adenocarcinoma and bladder cancer.
- [Cross-ancestry height polygenic risk modeling from summary statistics](https://um-kevinhe.github.io/BRIER/articles/Example_panAFR.html): `BRIERs` builds a height PRS for African-ancestry individuals from the Pan-UKB AFR GWAS, the PGS Catalog score PGS004999 and a public UKBB LD panel, tuned by GIC without individual-level data.


# Data and Resources

We deposited the data used in this tutorial on Zenodo (DOI: to be added):

- public LD reference panels for `BRIERs` modeling of polygenic risk: built from **1000 Genomes** and **UK Biobank** (the [PRS-CSx](https://github.com/getian107/PRScsx) reference panels) on HapMap3 SNPs, for African, European and East Asian ancestry;
- the data, scripts and fitted models of the two real data examples and the one simulated data example used in the tutorial.


# Getting Help

If you encounter problems or bugs, please contact us:

- Ruiwen Zhou (zrayw@umich.edu)
- Kevin He (kevinhe@umich.edu)


# Reference

Li Q, Patrick MT, Zhang H, Khunsriraksakul C, Stuart PE, Gudjonsson JE, Nair R, Elder JT, Liu DJ, Kang J, Tsoi LC, He K. Bregman Divergence-Based Data Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment. arXiv preprint arXiv:2210.06025. 2022 Oct 12. (https://arxiv.org/abs/2210.06025)
