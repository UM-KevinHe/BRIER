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