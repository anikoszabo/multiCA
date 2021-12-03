
#'Multinomial Cochran-Armitage trend test
#'
#'The \code{multiCA.test} function performs a multinomial generalization of the 
#' Cochran-Armitage trend test.
#'
#'
#'@export
#'@param x a two-dimensional matrix of event counts with the outcomes as rows and ordered groups as columns.
#'@param \dots other arguments 
#'@return a list with two components
#' \item{overall}{an object of class "htest" with the results of the overall test}
#' \item{individual}{a vector with adjusted p-values for individual outcomes}
#'@author Aniko Szabo
#'@references Szabo, A. (2016) Test for trend with a multinomial outcome.  
#'@keywords nonparametric 
#'@examples
#'
#'data(stroke)
#'## using formula interface
#'multiCA.test(Type ~ Year, weights=Freq, data=stroke)
#'
#'##using Westfall's multiple testing adjustment
#'multiCA.test(Type ~ Year, weights=Freq, data=stroke, p.adjust.method="Westfall")
#'
#'## using matrix interface and testing only the first 3 outcomes
#'strk.mat <- xtabs(Freq ~ Type + Year, data=stroke)
#'multiCA.test(strk.mat, outcomes=1:3)
#'
#'@name multiCA.test

multiCA.test <- function(x,...) UseMethod("multiCA.test")
 
