\name{getCharSTM2u}

\alias{getCharSTM2u}

\title{Computing Some Overall Characteristics in \code{compContourM2u}}

\description{
The function computes some overall characteristics 
of directional regression quantiles in the output 
of \code{\link{compContourM2u}}, namely the list 
\code{COutST$CharST}.  It makes possible 
to obtain some useful information without saving 
any file on the disk, and it can be easily modified 
by the users according to their wishes. 
}

\usage{
getCharSTM2u(Tau, N, M, P, BriefDQMat, CharST, IsFirst)
}

\arguments{
\item{Tau}{the quantile level in (0, 0.5).}

\item{N}{the number of observations.}

\item{M}{the dimension of responses.}

\item{P}{the dimension of regressors including the intercept.}

\item{BriefDQMat}{the method-specific matrix containing 
                  the rows of a potential individual output 
                  file corresponding to 
                  \code{CTechST$BriefOutputI} = 1. See 
                  the details below.}

\item{CharST}{the output list, updated with each run 
              of the function.}

\item{IsFirst}{the indicator equal to one in the 
                         first run of \code{getCharSTM2u} 
                         (when \code{CharST} is initialized) 
                         and equal to zero otherwise.}
}

\details{
This function is called inside \code{\link{compContourM2u}}.  
First, it is called with 
\code{BriefDQMat} = \code{NULL}, 
\code{CharST} = \code{NULL} and 
\code{IsFirst} = 1 to initialize 
the output list \code{CharST}, and then 
it is called with \code{IsFirst} = 0 
successively for the content of each potential 
output file corresponding to 
\code{CTechST$BriefOutputI} = 1, i.e., even if 
the output file(s) are not stored on the disk owing to 
\code{CTechST$OutSaveI} = 0.  

It still remains to describe in detail the content of possible 
output files, describing the optimal conic segmentation of 
the directional space that lies behind the optimization 
problem involved.  

If \code{CTechST$BriefOutputI} = 1, then the rows of such 
files are vectors of length \code{1+1+M+P*M+M} of the form 
\code{c(ConeID, Nu, UVec, vec(ACOMat), MuBRow)} where 
\describe{
\item{ConeID}{is the number/order of the cone related to the 
              line.  If \code{M} > 2, then a cone can appear in 
              the output repeatedly (under different numbers).}

\item{Nu}{is the number of negative residuals corresponding 
          to the interior directions of the cone.}

\item{UVec}{is a normalized vector of the cone.  
            It is usually its vertex direction but it may 
            also be its interior vector pointing to a vertex 
            of the artificial intersection of the cone with 
            the bounding box \code{[-1,1]^M}.  The 
            max normalization is used if the breadth-first 
            search algorithm is employed and the L2 
            normalization is used in the other case (when 
            \code{M} = 2 and \code{CTechST$D2SpecI} = 1).}

\item{ACOMat}{is the matrix describing \code{AVec}, 
              \code{AVec} = \code{ACOMat\%*\%UVec}.}

\item{MuBRow}{is the constant vector of the Lagrange 
              multipliers corresponding to \code{BVec}. 
              Its inner product with \code{UVec} is equal 
              to the optimal value \code{Psi} of the objective 
              function for that direction.}
}

Recall that \code{c(BVec, AVec)} stands for the coefficients 
of the regression quantile hyperplane associated with 
\code{UVec} and always \code{BVec} = \code{UVec}.  

If \code{CTechST$BriefOutputI} = 0, then the rows of the 
potential output file(s) are longer 
(of length \code{1+1+P*M+M+P+P}) 
because they contain two more vectors appended at the end.  
The rows are of the form 
\code{c(ConeID, Nu, UVec, vec(ACOMat), MuBRow, MuR0Row, IZ)} 
where 
\describe{
\item{MuRORow}{is the constant vector of the Lagrange 
              multipliers corresponding to zero residuals 
              associated with the interior of the cone.  
              That is to say that all directions 
              from the interior of the cone result in the 
              regression \code{Tau}-quantile hyperplanes 
              containing the same \code{P} observations.} 

\item{IZ}{is the vector containing original indices of the 
          \code{P} observations with zero residuals
          for all directions from the interior of the cone.}
}

}
\value{
\code{getCharSTM2u} returns a list with the following 
components:

\item{NUESkip}{the number of (skipped) directions (and 
               corresponding hyperplanes) artificially 
               induced by intersecting the cones with 
               the \code{[-1,1]^M} bounding box}

\item{NAZSkip}{the number of (skipped) hyperplanes (and 
               corresponding directions) not counted in NUESkip 
               and with at least one coordinate of 
               \code{AVec} zero.}

\item{NBZSkip}{the number of (skipped) hyperplanes (and 
               corresponding directions) not counted in NUESkip 
               and with at least one coordinate of 
               \code{BVec} zero.}

\item{HypMat}{(for \code{M > 4}) the component is missing  \cr
              (for \code{M <= 4}) the matrix with \code{M + P} 
               columns containing (in rows) all the distinct 
               regression Tau-quantile hyperplane 
               coefficients \code{c(BVec, AVec)} 
               rounded to the eighth decimal digit and 
               sorted lexicographically.  
               This matrix can be used for the computation 
               of the regression \code{Tau}-quantile contour.}

\item{CharMaxMat}{the matrix with the (slightly rounded) maxima 
                  of certain directional regression Tau-quantile 
                  characteristics over all remaining vertex 
                  directions.  \cr
                  If \code{P} = 1, then \code{CharMaxMat} has 
                  only two rows: 
                  \describe{
                  \item{}{\code{c(UVec, max(Psi))}, and}
                  \item{}{\code{c(UVec, max(|MuBRow|))},}
                  }
                  respectively.  \cr
                  If \code{P} > 1, then the rows of 
                  \code{CharMaxMat} are as follows: 
                  \describe{
                  \item{}{\code{c(UVec, max(|Psi|))},}
                  \item{}{\code{c(UVec, max(MuBRow))},}
                  \item{}{\code{c(UVec, max(|c(a_2,...,a_P)|))},}
                  \item{}{\code{c(UVec, max(|a_2|))},}
                  \item{}{\code{...},}
                  \item{}{\code{c(UVec, max(|a_P|))},}
                  }
                  respectively.  If \code{P} = 2, 
                  then the last row is missing for not 
                  being included twice. 
                  }

\item{CharMinMat}{the matrix with the (slightly rounded) minima 
                  of certain directional regression Tau-quantile 
                  characteristics over all remaining vertex 
                  directions.  \cr
                  If \code{P} = 1, then \code{CharMinMat} has 
                  only two rows: 
                  \describe{
                  \item{}{\code{c(UVec, min(Psi))}, and}
                  \item{}{\code{c(UVec, min(|MuBRow|))},}
                  }
                  respectively.  \cr
                  If \code{P} > 1, then \code{CharMinMat} 
                  has three rows: 
                  \describe{
                  \item{}{\code{c(UVec, min(Psi))},}
                  \item{}{\code{c(UVec, min(|MuBRow|))}, and}
                  \item{}{\code{c(UVec, min(|c(a_2,...,a_P)|))},}
                  }
                  respectively.}

Note that \code{||} symbolizes the Euclidean norm, 
and that the vertices (\code{UVec}) in the rows of 
\code{CharMaxMat} and \code{CharMinMat} are generally 
different and denote (one of) the direction(s) where 
the row maximum or minimum is attained.  
}

\examples{
##Run print(getCharSTM2u) to examine the default setting.
}
