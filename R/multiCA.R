
#' @keywords internal
#' @importFrom stats terms xtabs


.multiCA.test <- function(x, scores, outcomes){
  K <- nrow(x)
  full <- length(outcomes) == K  #full test
  
  nidot <- apply(x, 2, sum)
  n <- sum(nidot)
  
  cbar <- sum(nidot * scores)/n
  
  s2 <- sum(nidot * (scores - cbar)^2)
  pdot <- prop.table(rowSums(x))[outcomes]
  nonz <- (pdot > 0)
  
  if (!any(nonz)) return(1)
  
  X <- x[outcomes, ,drop=FALSE] %*% (scores - cbar)

  #individual tests
  Tt <-  X[nonz] / sqrt(pdot[nonz] * (1-pdot[nonz])* s2)
  CAT <- Tt^2
  CAT.p.value <- pchisq(CAT, df=1, lower.tail=FALSE)
  
  #overall test
  if (full || sum(pdot) >= 1){
    W <- ( sum(X[nonz]^2 / pdot[nonz])) / s2
  } else {
    W <- (sum(X)^2 / (1-sum(pdot)) + sum(X[nonz]^2 / pdot[nonz])) / s2
  }

  df <- length(outcomes) - full
  p.value <- pchisq(W, df=df, lower.tail=FALSE)
  
  
    # correlation of T
    sqrt.or.pdot <- sqrt(pdot[nonz]/(1-pdot[nonz]))
    Sigma0 <- -outer(sqrt.or.pdot, sqrt.or.pdot)
    diag(Sigma0) <- 1

    # contrast matrix
    if (full){
      coefs <- sqrt(pdot[nonz] * (1-pdot[nonz]))
      C <- rbind(coefs[-1], diag(K-1))
    } else {
      C <- diag(length(nonz))
    }
  

  res <- list(statistic = W, parameter = df, p.value = p.value, 
              indiv.statistics = Tt, indiv.p.value = CAT.p.value,
              sigma0 = Sigma0, contrast = C)
  return(res)
}

#'@rdname multiCA.test
#'@method multiCA.test default
#'@param scores non-decreasing numeric vector of the same length as the number of ordered groups. Defaults to linearly increasing values
#'@param outcomes integer or character vector defining the set of outcomes (by row index or row name) over which the trend should be tested. Defaults to all outcomes.
#'@param p.adjust.method character string defining the correction method for individual outcome p-values. Defaults to "closed.set" when \code{length(outcomes)<=3}, and "Holm-Shaffer" otherwise.
#'@export
#' @importFrom utils str
#' @importFrom multcomp glht adjusted parm

multiCA.test.default <- function(x, scores=1:ncol(x), outcomes=1:nrow(x),
  p.adjust.method=c("none","closed.set","Holm-Shaffer", "single-step", "Westfall"),...){
  if (!is.matrix(x)) {
    cat(str(x))
    stop("x should be a two-dimensional matrix")
}
  if (length(scores) != ncol(x)) stop("The length of the score vector should equal the number of columns of x")

  testres <- .multiCA.test(x=x, scores=scores, outcomes=outcomes)
 
  W <- c(W = testres$statistic)
  df <- c(df = testres$parameter)

  p.value <- testres$p.value
  null.value <- 0
  names(null.value) <- sprintf("slope for outcomes %s", deparse(substitute(outcomes)))

  res <- list(statistic = W, parameter = df, p.value = p.value, 
              method="Multinomial Cochran-Armitage trend test",
              alternative="two.sided",
              null.value=null.value,
              data.name = deparse(substitute(x)))
  class(res) <- "htest"

  
    if (missing(p.adjust.method)){
      if (length(outcomes)<=3) p.adjust.method <- "closed.set"
      else p.adjust.method <- "Holm-Shaffer"
    } else {
      p.adjust.method <- match.arg(p.adjust.method)
    }

    full.set <- (length(outcomes) == nrow(x)) 
    if (p.adjust.method=="none") {
      indiv.res <- testres$indiv.p.value
    } else if (p.adjust.method=="closed.set") {
      
        mytest <- function(hypotheses){
          .multiCA.test(x, scores, hypotheses)$p.value
        }
        indiv.res <- .p.adjust.closed(mytest, outcomes, remove=full.set)  
      
    } else if (p.adjust.method=="Holm-Shaffer") {
      
          s <- seq_along(testres$indiv.p.value)
          if (full.set) s[2] <- 3
          o <- order(testres$indiv.p.value)
          ro <- order(o)
          indiv.res <- pmin(1, cummax((length(outcomes) - s + 1L) * testres$indiv.p.value[o]))[ro]
      
    } else if (p.adjust.method %in% c("single-step", "Westfall")) {
      
          if (full.set) {
              testparm <- parm(testres$indiv.statistics[-1], testres$sigma0[-1,-1])
          } else {
              testparm <-  parm(testres$indiv.statistics, testres$sigma0)                            
              }
          g1 <- glht(model = testparm,  linfct = testres$contrast)
          indiv.res <- summary(g1, test=adjusted(type=p.adjust.method,...))$test$pvalues
      
    }
    attr(indiv.res, "method") <- p.adjust.method
  

  return(list(overall = res, individual = indiv.res))  
}

#'@rdname multiCA.test
#'@method multiCA.test formula
#'@param formula a formula of the form \code{outcome ~ group} where \code{outcome} is a factor representing the categorical outcome and \code{group} is the grouping variable over which the trend is tested.
#'@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@param subset an optional vector specifying a subset of observations to be used.
#'@param na.action      a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#'@param weights an integer-valued variable representing the number of times each \code{outcome} - \code{group} combination was observed.
#'@export
#' @importFrom stats terms xtabs

multiCA.test.formula <- function(formula, data, subset, na.action,  weights, ...){
    if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
        "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    responsevar <- attr(attr(mf, "terms"), "response")
    response <- mf[[responsevar]]
    weightvar <- which(names(mf)=="(weights)")
    w <- if(length(weightvar) > 0)  mf[[weightvar]] else rep(1L, nrow(mf))
    g <- factor(mf[,-c(responsevar, weightvar)])

    tab <- xtabs(w ~ response + g)
    multiCA.test(tab, ...)
}

#'Internal functions
#'
#' These internal functions perform the closed set p-value adjustment calculation
#' for the multivariate Cochran-Armitage trend test. The logical constraint
#' on the possible number of true null hypotheses is incorporated.
#'
#'
#' @name internal
#' @importFrom bitops bitAnd
#' @keywords internal
.bit2boolean <- function (x, N) 
{
  base <- 2^(1:N - 1)
  bitAnd(x, base) != 0
}

#' @param test function that performs the local test. The function should accept a subvector of the hypotheses argument as input, and return a p-value.
#' @param  hypotheses identifiers of the collection of elementary hypotheses. 
#' @param remove logical indicator of whether hypotheses of length N-1 should be removed
#' @param  ...  additional parameters to the 'test' function
#' @return  numeric vector of adjusted p-values for each hypothesis
#' @keywords internal
#' @name internal
.p.adjust.closed <- function (test, hypotheses, remove=FALSE, ...) 
{
  N <- length(hypotheses)
  Nmax <- log2(.Machine$integer.max + 1)
  if (N > Nmax) 
    stop("no more than ", Nmax, " hypotheses supported in full closed testing.\n Use a shortcut-based test.")
  closure <- 1:(2^N - 1)
  base <- 2^(1:N - 1)
  offspring <- function(x) {
    res <- bitAnd(x, closure)
    res[res != 0]
  }
  lengths <- rowSums(sapply(base, function(bs) bitAnd(closure, bs) != 0))
  
  idx <- sort.list(lengths, decreasing = TRUE)
  closure <- closure[idx]
  lengths <- lengths[idx]
  if (remove)  closure <- closure[lengths != (N-1)]
  
  adjusted <- numeric(2^N - 1)
  for (i in closure) {
    if (adjusted[i] < 1) {
      localtest <- test(hypotheses[.bit2boolean(i,N)], ...)
      if (localtest > adjusted[i]) {
        offs <- offspring(i)
        adjusted[offs] <- pmax(adjusted[offs], localtest)
      }
    }
  }
  
  out <- adjusted[base]
  names(out) <- hypotheses
  return(out)
}

#' Non-centrality parameter for chi-square distribution
#'
#' Calculates the non-centrality parameter for a chi-square distribution for a given 
#' quantile. This is often needed for sample size calculation for chi-square based tests.
#'
#'@details The function is modeled after the SAS function CNONCT. If \code{p} is larger 
#' than the cumulative probability of the central chi-square distribution at \code{x}, then
#' there is no solution and NA is returned.
#'
#'@param x a numeric value at which the distribution was evaluated
#'@param p a numeric value giving the cumulative probability at \code{x}
#'@param df an integer giving the degrees of freedom of the chi-square variable
#'@examples
#' (ncp <- cnonct(qchisq(0.95, df=10), 0.8, df=10))
#' ## check
#' pchisq(qchisq(0.95, df=10), df=10, ncp=ncp)  ## 0.8
#'@export
#'@importFrom stats pchisq uniroot

cnonct <- function(x, p, df){
 
  if (pchisq(x, df=df) < p) return(NA)

  f <- function(ncp){pchisq(x, df=df, ncp=pmax(0,ncp)) - p}

  res <- uniroot(f, interval=c(0, 100), extendInt="downX", tol=.Machine$double.eps^0.5)
  res$root
}

#' Power calculations for the multinomial Cochran-Armitage trend test
#'
#' Given the probabilities of outcomes, compute the power of the overall multinomial 
#' Cochran-Armitage trend test or determine the sample size to obtain a target power. 
#'
#'@details 
#' The sample size calculation depends only on \code{p.ave} - the weighted average probability of 
#' each outcome, and \code{slopes} - the weighted regression slope of each outcome.
#'
#' The values of these two key inputs can be specified in three ways: 
#'
#' 1. directly passing \code{p.ave} and  \code{slopes}, or 
#'
#' 2. specifying exactly two of the parameters \code{p.ave}, \code{slopes}, \code{p.start}, and \code{p.end}. 
#' In this case the full matrix of outcome probabilities will be inferred
#' assuming linearity within each outcome.
#'
#' 3. specifying the full matrix of outcome probabilities \code{pmatrix}.  
#' 
#' The calculation is based on approximating the distribution of the test statistic
#' under the alternative with a non-central chi-squared distribution instead of the correct
#' weighted mixture of multiple non-central chi-squares. This results in bias in the power
#' away from 50% - values above it are somewhat overestimated, while values under it are
#' underestimated.
#'
#' @param N integer, the total sample size of the study. If \code{NULL} then \code{power} needs to be specified.
#' @param power target power. If \code{NULL} then \code{N} needs to be specified.
#' @param pmatrix numeric matrix of hypothesized outcome probabilities in each group,  with 
#' the outcomes as rows and ordered groups as columns. The columns should add up to 1. 
#' @param p.ave numeric vector of average probability of each outcome over the groups  
#' weighted by \code{n.prop}.
#' @param p.start,p.end numeric vectors of the probability of each outcome for the  
#' first / last ordered group
#' @param slopes numeric vector of the hypothesized slope of each outcome when regressed  
#' against the column \code{scores} with weights \code{n.prop}. The values should add up to zero,
#' as the total probability is always 1 and has no trend.
#' @param scores non-decreasing numeric vector of the same length as the number of ordered groups  
#' giving the trend test scores. Defaults to linearly increasing values.
#' @param n.prop numeric vector describing relative sample sizes of the ordered groups.  
#' Will be normalized to sum to 1. Defaults to equal sample sizes.
#' @param G integer, number of ordered groups
#' @param sig.level significance level
#' @return object of class "power.htest"
#'
#' @examples
#' power.multiCA.test(power=0.8, p.start=c(0.1,0.2,0.3,0.4), p.end=c(0.4, 0.3, 0.2, 0.1), 
#'                      G=5, n.prop=c(3,2,1,2,3))
#'
#' ## Power of stroke study with 100 subjects per year and observed trends
#' data(stroke)
#' strk.mat <- xtabs(Freq ~ Type + Year, data=stroke)
#' power.multiCA.test(N=900, pmatrix=prop.table(strk.mat, margin=2))
#' @seealso \code{\link{power.CA.test}} for simpler (and more precise) power calculation 
#' with a binomial outcome 
#' @references Szabo, A. (2018). Test for Trend With a Multinomial Outcome. The American Statistician, 73(4), 313–320. 
#' @export
#' @importFrom stats pchisq qchisq weighted.mean


power.multiCA.test <- function(N=NULL, power=NULL, pmatrix=NULL, p.ave=NULL, p.start=NULL, 
                               p.end=NULL, slopes=NULL, scores=1:G, n.prop=rep(1, G),
                               G=length(p.ave), sig.level=0.05){
  if (sum(sapply(list(N, power), is.null)) != 1) 
        stop("exactly one of 'N',  and 'power' must be NULL")
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")


  if (!is.null(pmatrix)){
    K <- nrow(pmatrix)
    G <- ncol(pmatrix)
    if (!isTRUE(all.equal(colSums(pmatrix), rep(1, G), 
                          check.attributes=FALSE, use.names=FALSE))) 
      stop("pmatrix should have column sums of 1.")
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- as.vector(pmatrix %*% (n.prop * (scores-cbar))) / s2
    p.ave <- as.vector(pmatrix %*% n.prop)
  }
  else {
   if (sum(sapply(list(p.ave, slopes, p.start, p.end), is.null)) != 2) 
        stop("Either  pmatrix, or exactly two of 'p.ave', 'slopes', 'p.start', and 'p.end' must be specified (ie not NULL)")

  if (!is.null(p.ave) & !is.null(slopes)){
    if (length(p.ave) != length(slopes)) 
      stop("p.ave and slopes should have the same length")
    K <- length(p.ave)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
  }
  else if (!is.null(p.ave) & !is.null(p.start)){
    if (length(p.ave) != length(p.start)) 
      stop("p.ave and p.start should have the same length")
    K <- length(p.ave)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- (p.start - p.ave) / (scores[1] - cbar)
  }
  else if (!is.null(p.ave) & !is.null(p.end)){
    if (length(p.ave) != length(p.end)) 
      stop("p.ave and p.end should have the same length")
    K <- length(p.ave)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- (p.end - p.ave) / (scores[G] - cbar)
  }
  else if (!is.null(p.start) & !is.null(p.end)){
    if (length(p.start) != length(p.end)) 
      stop("p.start and p.end should have the same length")
    K <- length(p.start)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- (p.end - p.start) / (scores[G] - scores[1])
    p.ave <- p.start - slopes * (scores[1] - cbar)
  }
  else if (!is.null(p.start) & !is.null(slopes)){
    if (length(p.start) != length(slopes)) 
      stop("p.start and slopes should have the same length")
    K <- length(p.start)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    p.ave <- p.start - slopes * (scores[1] - cbar)
  }
  else if (!is.null(p.end) & !is.null(slopes)){
    if (length(p.end) != length(slopes)) 
      stop("p.end and slopes should have the same length")
    K <- length(p.end)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    p.ave <- p.end - slopes * (scores[G] - cbar)
  }
    
  
    if (!isTRUE(all.equal(sum(slopes), 0, check.attributes=FALSE, use.names=FALSE))) 
        stop("Implied or specified values of slopes should sum to 0.")
    if (!isTRUE(all.equal(sum(p.ave), 1, check.attributes=FALSE, use.names=FALSE))) 
        stop("Implied or specified values of p.ave should sum to 1.")
    check <- outer(1:K, 1:G, function(j,i)p.ave[j] + slopes[j]*(scores[i]-cbar))
    if (!all(check >= 0) || !(all(check <=1)))
      stop("The parameters do not define a valid probability matrix")  
  
}

  
  df <- K - 1
  crit <- qchisq(sig.level, df=df, lower.tail=FALSE)
  ncp0 <- sum(slopes^2 / p.ave) * s2 
  if (missing(power)){
    ncp <- ncp0 * N
    power <- pchisq(crit, df=df, ncp=ncp, lower.tail=FALSE)
   } 
   else {
     ncp <- cnonct(crit, p=1-power, df=df)
     N <-  ncp / ncp0
   }

   res <- structure(list(n = N, n.prop = n.prop, p.ave=p.ave, slopes = slopes, G = G, 
                         sig.level = sig.level, power = power,  
                         method = "Multinomial Cochran-Armitage trend test"), 
                     class = "power.htest")
   res
   }

#' Power calculations for the Cochran-Armitage trend test
#'
#' @param N integer, the total sample size of the study. If \code{NULL} then \code{power} needs to be specified.
#' @param power target power. If \code{NULL} then \code{N} needs to be specified.
#' @param pvec numeric vector of hypothesized outcome probabilities in each group. 
#' @param scores non-decreasing numeric vector of the same length as the number of ordered groups  
#' giving the trend test scores. Defaults to linearly increasing values.
#' @param n.prop numeric vector describing relative sample sizes of the ordered groups.  
#' Will be normalized to sum to 1. Defaults to equal sample sizes.
#' @param sig.level significance level
#' @param alternative character string specifying the alternative hypothesis
#' @return object of class "power.htest"
#'
#' @examples
#' # sample size required to detect with 80% power a decreasing trend over 4 groups 
#' # with 3:2:1:2 sample-size distribution at a 2.5% significance level
#' power.CA.test(power=0.8, pvec=c(0.4, 0.3, 0.2, 0.1), n.prop=c(3,2,1,2),
#'               alternative = "less", sig.level=0.025)
#'
#' # power of a 2-sided test to detect a logistic increase with slope 0.2 over 5 groups
#' # with groups of size 10 with unequal dose spacing
#' doses <- c(0,1,2,4,8)
#' p0 <- 0.05 # event probability at lowest dose
#' logit.props <- log(p0/(1-p0)) + doses * 0.2
#' p <- 1 / (1 + exp(-logit.props)) # hypothesized probabilities at each dose
#' power.CA.test(N = 10 * 5, pvec=p, scores = doses)
#'
#' @export
#' @references Nam, J. (1987). A Simple Approximation for Calculating Sample Sizes for Detecting Linear Trend in Proportions.
#' Biometrics, 43(3), 701-705. 
#' @importFrom stats pnorm qnorm 
#'

power.CA.test <- function(N=NULL, power=NULL, pvec=NULL, scores=seq_along(pvec), 
                          n.prop=rep(1, length(pvec)), sig.level=0.05,
                          alternative = c("two.sided", "less", "greater")){
  if (sum(sapply(list(N, power), is.null)) != 1) 
        stop("exactly one of 'N',  and 'power' must be NULL")
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")
  if (any(pvec < 0) | any(pvec > 1))
        stop("All probabilities in 'pvec' should be between 0 and 1")
  if (length(pvec) != length(scores) | length(pvec) != length(n.prop))
      stop("Vectors 'pvec', 'scores', and 'n.prop', if specified, should have the same lengths.")
  
        
  alternative <- match.arg(alternative)        
        
  n.prop <- n.prop / sum(n.prop)
  
  
  sbar <- sum(scores * n.prop)
  pbar <- sum(pvec * n.prop)
  v.nu <- sum(n.prop * (scores-sbar)^2)
  v0 <- pbar * (1-pbar) * v.nu
  v <- sum(n.prop * pvec * (1-pvec) * (scores-sbar)^2)
  mu <- sum(n.prop * pvec * (scores - sbar))
  
  if (alternative == "two.sided"){
    crit <- qchisq(sig.level, df=1, lower.tail=FALSE)
    ncp0 <- mu^2 / v
  } else {
    crit <- qnorm(sig.level, lower.tail=FALSE)
    if (alternative == "less") crit <- (-1) * crit
  }
  
  
  if (missing(power)){
    if (alternative == "two.sided"){
      ncp <- ncp0 * N
      power <- pchisq(v0/v * crit, df=1, ncp=ncp, lower.tail=FALSE)
    } else {
      term <- (sqrt(v0)*crit - mu * sqrt(N)) / sqrt(v)
      power <- pnorm(term, lower.tail = (alternative == "less"))
    }
   } 
   else {
    if (alternative == "two.sided"){
      ncp <- cnonct(v0/v *crit, p=1-power, df=1)
      N <-  ncp / ncp0
    } else {
      zb <- qnorm(power, lower.tail=TRUE)
      N <- (sqrt(v0) * crit + sqrt(v) * zb)^2 / mu^2
    }
   }

   res <- structure(list(n = N, n.prop = n.prop, p = pvec, 
                        alternative = alternative,  
                        sig.level = sig.level, power = power,  
                         method = "Cochran-Armitage trend test"), 
                     class = "power.htest")
   res
   }
        
