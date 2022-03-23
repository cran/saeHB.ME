#' saeHB.ME: Small Area Estimation with Measurement Error using Hierarchical Bayesian Method
#'
#' Implementation of small area estimation using Hierarchical Bayesian (HB) Method when auxiliary variable measured with error.
#' The 'rjags' package is employed to obtain parameter estimates.
#'
#' @section Authors:
#' Azka Ubaidillah, Muhammad Rifqi Mubarak
#'
#' @section Email:
#' Muhammad Rifqi Mubarak \email{rifqi.mubarak@bps.go.id}
#'
#' @section Functions:
#' \describe{
#'     \item{\code{\link{meHBNormal}}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under normal distribution.}
#'     \item{\code{\link{meHBt}}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under student-t distribution.}
#'     }
#'
#' @references Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New York: John Wiley and Sons, Inc <doi:10.1002/9781118735855>.
#' @references Ybarra, L.M. and Lohr, S. L. (2008). Small area estimation when auxiliary information is measured with error. Biometrika 95, 919-931 <doi:10.1093/biomet/asn048>.
#' @references Ntzoufras, I. (2009), Bayesian Modeling Using WinBUGS. 1st Edn., Wiley, New Jersey, ISBN-10: 1118210352.
#' @docType package
#' @name saeHB.ME
#' @importFrom stringr str_replace_all
#' @importFrom grDevices graphics.off
#' @importFrom graphics par
#' @import rjags
#' @import coda
#' @import stats
NULL

