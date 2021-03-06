\name{compContourM1/2u}

\alias{compContourM1u}
\alias{compContourM2u}
\alias{compContourM1/2u}

\title{Directional Regression Quantile Computation}

\description{
The functions \code{compContourM1u} and \code{compContourM2u} 
may be used to obtain not only directional regression quantiles 
for \emph{all} directions, but also some related overall 
statistics.  Their output may also be used for the evaluation 
of the corresponding regression quantile regions by means of 
\code{\link{evalContour}}.  The functions use different 
methods and algorithms, namely \code{compContourM1u} is based 
on [01] and [06] and \code{compContourM2u} results from [03] 
and [07].  The corresponding regression quantile regions are 
nevertheless virtually the same.  See all the references below 
for further details and possible applications.
}

\usage{
compContourM1u(Tau = 0.2, YMat = NULL, XMat = NULL, CTechST = NULL)
compContourM2u(Tau = 0.2, YMat = NULL, XMat = NULL, CTechST = NULL)
}

\arguments{
\item{Tau}{the quantile level in (0, 0.5).}
\item{YMat}{the N x M response matrix with two to six 
            columns, \code{N > M+P-1}.  
            Each row corresponds to one observation.}
\item{XMat}{the N x P design matrix including the (first) 
            intercept column. The default NULL value 
            corresponds to the unit vector of the right 
            length.  
            Each row corresponds to one observation.}
\item{CTechST}{the (optional) list with some parameters 
               influencing the computation and its output.  
               Its default value can be generated by 
               method-dependent \code{\link{getCTechSTM1/2u}} 
               and then modified by the user before its use 
               in \code{compContourM1/2u}.}
}

\details{
Generally, the performance of the functions deteriorates with 
increasing Tau, N, M, and P as for their reliability and 
time requirements. Nevertheless, they should work fine 
at least for two-dimensional problems 
up to N = 10000 and P = 10, 
for three-dimensional problems
 up to N = 500 and P = 5,
and for four-dimensional problems
 up to N = 150 and P = 3. 

Furthemore, common problems related to the computation 
can fortunately be prevented or overcome easily.  

\bold{Bad data} - the computation may fail if the processed 
data points are in a bad configuration (i.e., if they are not 
in general position or if they would lead to a quantile 
hyperplane with at least one zero coefficient), which mostly 
happens when discrete-valued/rounded/repeated observations, 
dummy variables or bad random number generators are employed.  
Such problems can often be prevented if one perturbs the data 
with a random noise of a reasonably small magnitude before the 
computation, splits the model into separate or independent 
submodels, cleverly uses affine equivariance, or replaces 
a few identical observations with a copy of them weighted 
by the total number of their occurrences.  

\bold{Bad Tau} - the computation may fail for a finite number of 
problematic quantile levels, e.g., if Tau is an integer multiple 
of 1/N in the location case with unit weights (when the sample 
quantiles are not uniquely defined).  Such a situation may occur 
easily for Tau's with only a few decimal digits or in a 
fractional form, especially when the number of observations 
changes automatically during the computation.  The problem 
can be fixed easily by perturbing Tau with a sufficiently small 
number in the right direction, which should not affect the 
resulting regression quantile contours although it may slightly 
change the other output.  The strategy is also adopted 
by \code{compContourM1/2u}, but only in the location case and 
with a warning output message explaining it.  

\bold{Bad scale} - the computation may fail easily for badly 
scaled data.  That is to say that the functionality has been 
heavily tested only for the observations coming from 
a centered unit hypercube.  Nevertheless, you can always 
change the units of measurements or employ full affine 
equivariance to avoid all the troubles.  Similar problems 
may also arise when properly scaled data are used with highly 
non-uniform weights, which frequently happens in local(ly) 
polynomial regression.  Then the weights can be rescaled 
in a suitable way and the observations with virtually zero 
weights can be excluded from the computation.  

\bold{Bad expectations} - the computation and its output need 
not meet false expectations.  Every user should be aware of 
the facts that the computation may take a long time or fail 
even for moderately sized three-dimensional data sets, that the 
\code{HypMat} component is not always present in the list 
\code{COutST$CharST} by default, and that the sample regression
quantile contours can be not only empty, but also unbounded and 
crossing one another in the general regression case.  

\bold{Bad interpretation} - the output results may be easily 
interpreted misleadingly or erroneously.  That is to say that 
the quantile level Tau is not linked to the probability 
content of the sample (regression) Tau-quantile region in 
any straightforward way.  Furthermore, any meaningful 
parametric quantile regression model should include 
as regressors not only the variables influencing the 
trend, but also all those affecting the dispersion 
of the multivariate responses.  Even then the cuts 
of the resulting regression quantile contours parallel 
to the response space cannot be safely interpreted 
as conditional multivariate quantiles except for 
some very special cases.  Nevertheless, such a 
conclusion could somehow be warranted in case of 
nonparametric multiple-output quantile regression; see [09].  
}

\value{
Both compContourM1u and compContourM2u may display some 
auxiliary information regarding the computation on the screen 
(if \code{CTechST$ReportI} = 1) or store their in-depth
output (determined by \code{CTechST$BriefOutputI}) in the output
files (if \code{CTechST$OutSaveI} = 1) with
the filenames beginning with the string contained in 
\code{CTechST$OutFilePrefS}, followed by the file number
padded with zeros to form six digits 
and by the extension \file{.dqo}, respectively.  The first 
output file produced by \code{compContourM1u} would 
thus be named \file{DQOutputM1_000001.dqo}.

Both compContourM1u and compContourM2u always return 
a list with the same components. Their interpretation is 
also the same (except for CharST that itself contains some 
components that are method-specific):

\item{CharST}{the list with some default or 
              user-defined output.  
              The default one is provided 
              by function \code{\link{getCharSTM1u}} 
              for \code{compContourM1u} and by 
              function \code{\link{getCharSTM2u}} 
              for \code{compContourM2u}.  
              A user-defined function generating 
              its own output can be employed instead 
              by changing \code{CTechST$getCharST}.}

\item{CTechSTMsgS}{the (possibly empty) string that informs 
                   about the problems with input \code{CTechST}.}

\item{ProbSizeMsgS}{the (possibly empty) string that warns 
                    if the input problem is very large.}

\item{TauMsgS}{the (possibly empty) string that announces 
               an internal perturbation of \code{Tau}.}

\item{CompErrMsgS}{the (possibly empty) string that decribes 
                   the error interrupting the computation.}

\item{NDQFiles}{the counter of (possible) output files, 
                i.e., as if \code{CTechST$OutSaveI} = 1.}

\item{NumB}{the counter of (not necessarily distinct) optimal 
            bases considered.}

\item{PosVec}{the vector of length N that desribes 
              the position of individual (regression) 
              observations with respect to the 
              exact (regression) Tau-quantile 
              contour.  
              The identification is reliable only after a 
              successful computation.  
              \code{PosVec[i]} = 0/1/2 if the \code{i}-th 
              observation is in/on/out of the contour.  
              If \code{compContourM2u} is used with 
              \code{CTechST$SkipRedI} = 1, then \code{PosVec}
              correctly detects only all the outer observations.}

\item{MaxLWidth}{the maximum width of one layer of the 
                 internal algorithm.}

\item{NIniNone}{the number of trials when the initial 
                solution could not be found at all.}

\item{NIniBad}{the number of trials when the found 
               initial solution did not have 
               the right number of clearly nonzero 
               coordinates.}

\item{NSkipCone}{the number of skipped cones (where 
                 an interior point could not be found).}

If \code{CTechST.CubRegWiseI} = 1, then the last four 
components are calculated over all the individual 
orthants.
}

\references{
[01] Hallin, M., Paindaveine, D. and Šiman, M. (2010)
       Multivariate quantiles and multiple-output regression quantiles:
        from L1 optimization to halfspace depth.
       \emph{Annals of Statistics} \bold{38}, 635--669.

[02] Hallin, M., Paindaveine, D. and Šiman, M. (2010)
       Rejoinder (to [01]).
       \emph{Annals of Statistics} \bold{38}, 694--703.

[03] Paindaveine, D. and Šiman, M. (2011)
       On directional multiple-output quantile regression.
       \emph{Journal of Multivariate Analysis} \bold{102}, 193--212.

[04] Šiman, M. (2011)
       On exact computation of some statistics based on projection
        pursuit in a general regression context.
       \emph{Communications in Statistics - Simulation and Computation} \bold{40}, 948--956.

[05] McKeague, I. W., López-Pintado, S., Hallin, M. and Šiman, M. (2011)
       Analyzing growth trajectories.
       \emph{Journal of Developmental Origins of Health and Disease} \bold{2}, 322--329.

[06] Paindaveine, D. and Šiman, M. (2012)
       Computing multiple-output regression quantile regions.
       \emph{Computational Statistics & Data Analysis} \bold{56}, 840--853.

[07] Paindaveine, D. and Šiman, M. (2012)
       Computing multiple-output regression quantile regions
        from projection quantiles.
       \emph{Computational Statistics} \bold{27}, 29--49.

[08] Šiman, M. (2014)
       Precision index in the multivariate context.
       \emph{Communications in Statistics - Theory and Methods} \bold{43}, 377--387.

[09] Hallin, M., Lu, Z., Paindaveine, D. and Šiman, M. (2015)
       Local bilinear multiple-output quantile/depth regression.
       \emph{Bernoulli} \bold{21}, 1435--1466.
}

\examples{
##computing all directional 0.15-quantiles of 199 random points
##uniformly distributed in the unit square centered at zero
##- preparing the input
Tau  <- 0.15
XMat <- matrix(1, 199, 1)
YMat <- matrix(runif(2*199, -0.5, 0.5), 199, 2)
##- Method 1:
COutST <- compContourM1u(Tau, YMat, XMat)
##- Method 2:
COutST <- compContourM2u(Tau, YMat, XMat)
}
