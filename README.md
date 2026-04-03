<style>
h1, h2, h3 {
    border-bottom: none !important;
    box-shadow: none !important;
}
</style>


# Description

`BRIER2` is a comprehensive Bregman-divergence–based integration framework for **improving genetic risk prediction** of complex traits in a **local target** cohort 
by leveraging information **from external sources** under diverse data-availability settings. 
It flexibly accommodates a broad range of outcome types, including continuous, binary, count, and time-to-event traits. 
The framework is implemented through three principled functions, `BRIER.FULL`, `BRIER.I`, and `BRIER.S`, 
which address the major data-access constraints encountered in biomedical and statistical-genetics research.


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
    <td><code>BRIER.FULL</code></td>    
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
    <td><code>BRIER.I</code></td>
    <td>
        Pretrained models <br>
        (e.g., penalized regression coefficients)
    </td>
</tr>

<tr style="border-bottom:1px solid black;">   
    <td><code>BRIER.S</code></td>
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

`BRIER2` extends the flexibility and scalability of risk-prediction modeling, with a particular emphasis on high-dimensional genetic applications:

- **Flexibility to integrate external information according to different data access and privacy constraints**
- **Support various distribution families for complex traits**  
- **A suite of hyperparameter tuning tools.**  
- **Ensemble weighting scheme integration of multiple external models.**  

Several features to improve practical concerns in genetic applications:

- **Improved computational efficiency and parallel computing.**  
- **Support input for imputed genotype data.**  
- **Precomputed ancestry-specific LD panels.**  
- **Various downstream analysis tools for interpretability**  


# Installation

> **Note:** This package is currently under active development. Please report any issues or unexpected behavior.  
> Requires **R ≥ 4.0.0**.

Install from CRAN:

    install.packages("BRIER2")

Or install the development version from GitHub:

    require("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/BRIER2")

<br>


# Detailed Tutorial

Full package documentation and parameter explanations: [here](https://um-kevinhe.github.io/BRIER2/)


# Getting Help

If you encounter problems or bugs, please contact us:

- zrayw@umich.edu
- kevinhe@umich.edu


# Reference

Li Q, Patrick MT, Zhang H, Khunsriraksakul C, Stuart PE, Gudjonsson JE, Nair R, Elder JT, Liu DJ, Kang J, Tsoi LC, He K. Bregman Divergence-Based Data Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment. arXiv preprint arXiv:2210.06025. 2022 Oct 12. (https://arxiv.org/abs/2210.06025)

Li Q, Patrick MT, Fritsche L, Wang D, Kang J, Tsoi LC, He K. Bregman Divergence-Based Data Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment. arXiv preprint arXiv:2210.06025. 2022 Oct 12. (https://arxiv.org/abs/2210.06025)