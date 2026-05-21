#' Simulated example data for BRIERfull
#'
#' A simulated Gaussian outcome dataset with one target cohort and three
#' external cohorts. The dataset was generated to illustrate the use of
#' `BRIERfull()` when individual-level data are available from both the target
#' and external cohorts.
#'
#' Each cohort contains uncorrelated predictors, a continuous Gaussian outcome,
#' and marginal summary statistics for the training split. The predictors were
#' simulated independently from a standard normal distribution with 200
#' predictors in each cohort. The true signal variables are located at columns
#' 20, 40, 60, ..., 200. The outcome was generated from a sparse linear model
#' with Gaussian noise. The noise variance was chosen so that the proportion of
#' phenotype variance explained by the linear predictor is approximately
#' `hsq = 0.2`.
#'
#' The target cohort and external cohorts share the same set of signal
#' variables but use different true coefficient values. This creates related
#' but non-identical genetic-effect structures across cohorts, mimicking a
#' setting where external cohorts provide useful but imperfectly transferable
#' information.
#'
#' @format A list named `Data_BRIERfull` with four components:
#' \describe{
#'   \item{target}{
#'     Target cohort. Contains training, validation, and testing data.
#'     The training sample size is 150, and the validation and testing sample
#'     sizes are both 150.
#'   }
#'   \item{external1}{
#'     First external cohort. Contains training and validation data.
#'     The training sample size is 300, and the validation sample size is 150.
#'   }
#'   \item{external2}{
#'     Second external cohort. Contains training and validation data.
#'     The training sample size is 300, and the validation sample size is 150.
#'   }
#'   \item{external3}{
#'     Third external cohort. Contains training and validation data.
#'     The training sample size is 300, and the validation sample size is 150.
#'   }
#' }
#'
#' Each cohort split contains some or all of the following elements:
#' \describe{
#'   \item{X}{Predictor matrix with 200 columns.}
#'   \item{y}{Continuous Gaussian outcome vector.}
#'   \item{XB}{True linear predictor used to generate the outcome.}
#'   \item{var_epsilon}{Noise variance used in the Gaussian outcome model.}
#'   \item{true_id}{Indices of the true signal variables: `20, 40, ..., 200`.}
#'   \item{true_beta}{True coefficient values used for the corresponding signal variables.}
#'   \item{hsq}{Target proportion of variance explained by the linear predictor; set to `0.2`.}
#'   \item{sumstats}{
#'     Marginal summary statistics for the training split, including correlation,
#'     test statistic, degrees of freedom, p-value, and sample size for each predictor.
#'   }
#' }
#'
#' The true coefficient values are:
#' \describe{
#'   \item{target}{
#'     `c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)`.
#'   }
#'   \item{external1}{
#'     `c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)`.
#'   }
#'   \item{external2}{
#'     `c(1.2, 1.7, 0.8, 1.3, 1.8, -1.2, -1.7, -0.8, -1.3, -1.8)`.
#'   }
#'   \item{external3}{
#'     `c(1.9, 2.4, 0.1, 0.6, 1.1, -1.9, -2.4, -0.1, -0.6, -1.1)`.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' data(Data_BRIERfull)
#' names(Data_BRIERfull)
#' }
#' @name Data_BRIERfull
#' @docType data
NULL


#' Simulated example data for BRIERi
#'
#' A simulated Gaussian outcome dataset with one target cohort and three
#' external cohorts. The dataset was generated to illustrate the use of
#' `BRIERi()` when individual-level data are available from both the target
#' and external cohorts.
#'
#' Each cohort contains uncorrelated predictors, a continuous Gaussian outcome,
#' and marginal summary statistics for the training split. The predictors were
#' simulated independently from a standard normal distribution with 200
#' predictors in each cohort. The true signal variables are located at columns
#' 20, 40, 60, ..., 200. The outcome was generated from a sparse linear model
#' with Gaussian noise. The noise variance was chosen so that the proportion of
#' phenotype variance explained by the linear predictor is approximately
#' `hsq = 0.2`.
#'
#' The target cohort and external cohorts share the same set of signal
#' variables but use different true coefficient values. This creates related
#' but non-identical genetic-effect structures across cohorts, mimicking a
#' setting where external cohorts provide useful but imperfectly transferable
#' information.
#'
#' @format A list named `Data_BRIERi` with two components:
#' \describe{
#'   \item{target}{
#'     Target cohort. Contains training, validation, and testing data.
#'     The training sample size is 150, and the validation and testing sample
#'     sizes are both 150.
#'   }
#'   \item{beta.external}{
#'     A matrix of trained coefficients for the three external cohorts.
#'   }
#' }
#'
#' Each cohort split contains some or all of the following elements:
#' \describe{
#'   \item{X}{Predictor matrix with 200 columns.}
#'   \item{y}{Continuous Gaussian outcome vector.}
#'   \item{XB}{True linear predictor used to generate the outcome.}
#'   \item{var_epsilon}{Noise variance used in the Gaussian outcome model.}
#'   \item{true_id}{Indices of the true signal variables: `20, 40, ..., 200`.}
#'   \item{true_beta}{True coefficient values used for the corresponding signal variables.}
#'   \item{hsq}{Target proportion of variance explained by the linear predictor; set to `0.2`.}
#'   \item{sumstats}{
#'     Marginal summary statistics for the training split, including correlation,
#'     test statistic, degrees of freedom, p-value, and sample size for each predictor.
#'   }
#' }
#'
#' The true coefficient values used to simulate the source cohorts are:
#' \describe{
#'   \item{target}{
#'     `c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)`.
#'   }
#'   \item{external1}{
#'     `c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)`.
#'   }
#'   \item{external2}{
#'     `c(1.2, 1.7, 0.8, 1.3, 1.8, -1.2, -1.7, -0.8, -1.3, -1.8)`.
#'   }
#'   \item{external3}{
#'     `c(1.9, 2.4, 0.1, 0.6, 1.1, -1.9, -2.4, -0.1, -0.6, -1.1)`.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' data(Data_BRIERi)
#' names(Data_BRIERi)
#' }
#' @name Data_BRIERi
#' @docType data
NULL


#' Simulated example data for BRIERs
#'
#' A simulated Gaussian outcome dataset with one target cohort and three
#' external cohorts. The dataset was generated to illustrate the use of
#' `BRIERs()` when individual-level data are available from both the target
#' and external cohorts.
#'
#' Each cohort contains uncorrelated predictors, a continuous Gaussian outcome,
#' and marginal summary statistics for the training split. The predictors were
#' simulated independently from a standard normal distribution with 200
#' predictors in each cohort. The true signal variables are located at columns
#' 20, 40, 60, ..., 200. The outcome was generated from a sparse linear model
#' with Gaussian noise. The noise variance was chosen so that the proportion of
#' phenotype variance explained by the linear predictor is approximately
#' `hsq = 0.2`.
#'
#' The target cohort and external cohorts share the same set of signal
#' variables but use different true coefficient values. This creates related
#' but non-identical genetic-effect structures across cohorts, mimicking a
#' setting where external cohorts provide useful but imperfectly transferable
#' information.
#'
#' @format A list named `Data_BRIERs` with two components:
#' \describe{
#'   \item{target}{
#'     Target cohort. Contains training, validation, and testing data.
#'     The training sample size is 150, and the validation and testing sample
#'     sizes are both 150.
#'   }
#'   \item{beta.external}{
#'     A matrix of trained coefficients for the three external cohorts.
#'   }
#' }
#'
#' Each cohort split contains some or all of the following elements:
#' \describe{
#'   \item{X}{Predictor matrix with 200 columns.}
#'   \item{y}{Continuous Gaussian outcome vector.}
#'   \item{XB}{True linear predictor used to generate the outcome.}
#'   \item{var_epsilon}{Noise variance used in the Gaussian outcome model.}
#'   \item{true_id}{Indices of the true signal variables: `20, 40, ..., 200`.}
#'   \item{true_beta}{True coefficient values used for the corresponding signal variables.}
#'   \item{hsq}{Target proportion of variance explained by the linear predictor; set to `0.2`.}
#'   \item{sumstats}{
#'     Marginal summary statistics for the training split, including correlation,
#'     test statistic, degrees of freedom, p-value, and sample size for each predictor.
#'   }
#' }
#'
#' The true coefficient values used to simulate the source cohorts are:
#' \describe{
#'   \item{target}{
#'     `c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)`.
#'   }
#'   \item{external1}{
#'     `c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)`.
#'   }
#'   \item{external2}{
#'     `c(1.2, 1.7, 0.8, 1.3, 1.8, -1.2, -1.7, -0.8, -1.3, -1.8)`.
#'   }
#'   \item{external3}{
#'     `c(1.9, 2.4, 0.1, 0.6, 1.1, -1.9, -2.4, -0.1, -0.6, -1.1)`.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' data(Data_BRIERs)
#' names(Data_BRIERs)
#' }
#' @name Data_BRIERs
#' @docType data
NULL


#' Example data for preprocessing individual-level genotype inputs
#'
#' A small example dataset used to demonstrate genotype preprocessing with
#' \code{\link{preprocessI}} before fitting \code{\link{BRIERi}}. The example
#' is based on a height PRS prediction analysis in the Michigan Genomics
#' Initiative (MGI), where the target cohort contains African-ancestry
#' individuals and the external models are derived from European-ancestry
#' resources.
#'
#' The dataset contains SNP annotation information for the target genotype
#' matrix and two external SNP-level coefficient files. Raw individual-level
#' genotype and phenotype data are not included due to privacy constraints.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{target.info}{A data.frame with 10,000 rows describing the SNPs in
#'     the target genotype matrix. Each row corresponds to one target genotype
#'     column. The data.frame includes columns \code{CHR} (chromosome),
#'     \code{BP} (base-pair position), \code{A2} (reference allele), and
#'     \code{A1} (alternative/effect allele), following PLINK conventions
#'     where \code{A1} is the effect allele.}
#'   \item{external.coef1}{A data.frame with 10,000 rows containing SNP-level
#'     coefficients from the first external model. The required columns are
#'     \code{CHR}, \code{BP}, \code{REF}, \code{ALT}, and \code{coef}, where
#'     \code{coef} is the external model coefficient for the corresponding
#'     SNP.}
#'   \item{external.coef2}{A data.frame with 163,140 rows containing SNP-level
#'     coefficients from the second external model. The required columns are
#'     \code{CHR}, \code{BP}, \code{REF}, \code{ALT}, and \code{coef}. This
#'     external model has substantially more SNPs than \code{external.coef1},
#'     allowing users to demonstrate merging and alignment of multiple
#'     external models.}
#' }
#'
#' @details
#' This object is intended for illustrating the preprocessing workflow:
#' first merge multiple external coefficient files with
#' \code{\link{mergeExternals}}, then align the merged external coefficients
#' to the target SNP set using \code{\link{preprocessI}}. The output of
#' \code{\link{preprocessI}} can be used to subset the target genotype matrix
#' and construct the \code{beta.external} matrix required by
#' \code{\link{BRIERi}}.
#'
#' Because \code{\link{BRIERi}} expects an intercept slot as the first row
#' of \code{beta.external}, users should prepend a zero row to the aligned
#' external coefficient matrix when no external intercept is used.
#'
#' @source
#' Michigan Genomics Initiative (MGI), African-ancestry subcohort
#' (\code{target.info}); MGI European-ancestry subcohort
#' (\code{external.coef1}); published European-ancestry height GWAS summary
#' statistics (\code{external.coef2}). Raw individual-level data are not
#' redistributed.
#'
#' @seealso \code{\link{preprocessI}}, \code{\link{mergeExternals}},
#'   \code{\link{BRIERi}}
#'
#' @examples
#' \dontrun{
#' data(Data_preprocessI)
#' names(Data_preprocessI)
#' }
#' @name Data_preprocessI
#' @docType data
#' @keywords datasets
NULL


#' Example data for preprocessing summary-level genotype inputs
#'
#' A small example dataset used to demonstrate summary-statistic preprocessing
#' with \code{\link{preprocessS}} before fitting \code{\link{BRIERs}}. The
#' example is based on a height PRS prediction analysis in the Michigan
#' Genomics Initiative (MGI), where the target data are provided as GWAS
#' summary statistics rather than individual-level genotypes.
#'
#' The dataset contains target GWAS summary statistics, target LD reference
#' SNP information, a precomputed target LD matrix, and two external
#' SNP-level coefficient files. The target LD reference panel defines the
#' canonical SNP order and allele orientation for downstream
#' \code{\link{BRIERs}} modeling.
#'
#' @format A list with five elements:
#' \describe{
#'   \item{target.ss}{A data.frame with 10,000 rows of target GWAS summary
#'     statistics, with columns \code{SNP}, \code{CHR} (chromosome),
#'     \code{BP} (base-pair position), \code{A2} (reference allele),
#'     \code{A1} (alternative/effect allele), \code{BETA} (effect size),
#'     \code{STAT} (z-statistic), \code{P} (p-value), and \code{NMISS}
#'     (sample size). PLINK conventions are followed where \code{A1} is the
#'     effect allele.}
#'   \item{target.ld}{A data.frame with 10,000 rows describing the SNPs in
#'     the target LD reference panel, with columns \code{SNP}, \code{CHR},
#'     \code{BP}, \code{A2}, and \code{A1}. Rows correspond to the rows and
#'     columns of \code{target.ld.mat}.}
#'   \item{target.ld.mat}{A 10,000 x 10,000 sparse LD matrix
#'     (\code{Matrix::dgCMatrix}) for the SNPs described in
#'     \code{target.ld}. After running \code{\link{preprocessS}}, this
#'     matrix should be subset using \code{processed$target.ld.keep} before
#'     being passed to \code{\link{BRIERs}}.}
#'   \item{external.coef1}{A data.frame with 10,000 rows containing SNP-level
#'     coefficients from the first external model, with required columns
#'     \code{CHR}, \code{BP}, \code{REF}, \code{ALT}, and \code{coef}.}
#'   \item{external.coef2}{A data.frame with 163,140 rows containing SNP-level
#'     coefficients from the second external model, with required columns
#'     \code{CHR}, \code{BP}, \code{REF}, \code{ALT}, and \code{coef}. This
#'     model has substantially more SNPs than \code{external.coef1},
#'     allowing users to demonstrate external-model merging, allele
#'     alignment, and zero-padding.}
#' }
#'
#' @details
#' This object is intended for illustrating the \code{\link{BRIERs}}
#' preprocessing workflow. Multiple external coefficient files can first be
#' combined using \code{\link{mergeExternals}}. Then,
#' \code{\link{preprocessS}} aligns the target summary statistics and
#' external coefficients to the target LD reference panel.
#'
#' When \code{target.ind = "gwas"}, \code{\link{preprocessS}} converts GWAS
#' summary statistics into marginal SNP-outcome correlations using p-values,
#' sample sizes, and effect directions. The processed \code{target.ss} object
#' contains the aligned correlation vector, and \code{target.ld.keep}
#' provides the indices needed to subset \code{target.ld.mat}.
#'
#' Unlike \code{\link{BRIERi}}, \code{\link{BRIERs}} does not require an
#' intercept row in \code{beta.external}.
#'
#' @source
#' Michigan Genomics Initiative (MGI), African-ancestry subcohort
#' (\code{target.ss}, \code{target.ld}, \code{target.ld.mat}); MGI
#' European-ancestry subcohort (\code{external.coef1}); published
#' European-ancestry height GWAS summary statistics
#' (\code{external.coef2}). Raw individual-level data are not redistributed.
#'
#' @seealso \code{\link{preprocessS}}, \code{\link{mergeExternals}},
#'   \code{\link{BRIERs}}, \code{\link{p2cor}}
#'
#' @examples
#' \dontrun{
#' data(Data_preprocessS)
#' names(Data_preprocessS)
#' }
#' @name Data_preprocessS
#' @docType data
#' @keywords datasets
NULL