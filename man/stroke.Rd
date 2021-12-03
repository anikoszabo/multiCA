\name{stroke}
\alias{stroke}
\docType{data}
\title{
Stroke types over time
}
\description{
Nakajima et al. (2014) collected information on stroke patients over a 9-year period. For each patient, the type of stroke was classified into one of 5 categories by etiology.
}
\usage{data("stroke")}
\format{
  A data frame with 45 observations on the following 3 variables.
  \describe{
    \item{\code{Type}}{a factor with levels \code{Small vessel occlusion}, \code{Large artery atherosclerosis}, \code{Cardioembolism}, \code{Other determined aetiology}, and \code{Undetermined aetiology} giving the etiology of the stroke}
    \item{\code{Year}}{a numeric vector with the year of the observation}
    \item{\code{Freq}}{a numeric vector with the number of patients with a stroke of the given etiology that year}
  }
}

\source{
Nakajima, M., Y. Inatomi, T. Yonehara, Y. Hashimoto, T. Hirano, and Y. Ando (2014).
Temporal trends in oral intake ability 3 months after acute ischaemic stroke: analysis of
a single-centre database from 2003 to 2011. J Rehabil Med 46 (3), 200--205.
}

\examples{
data(stroke)
xtabs(Freq ~ Type + Year, data=stroke)

strk.props <- prop.table(xtabs(Freq ~ Year+Type, data=stroke), margin=1)
matplot(strk.props, type="l")
}
\keyword{datasets}
