\name{od.opt.param}
\alias{od.opt.param}
\title{Optimal Parameter Values In RaPKod}
\description{
Uses a heuristic formula to set optimal values for gamma and p.}
\usage{
od.opt.param(X, K1 = 6, K2 = 50, which.estim = "Gauss", RATIO = 0.1, 
            randomize = TRUE, sub.n = floor(nrow(X)))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame or an n x d matrix.}
  \item{K1}{universal constant used in the heuristic formula of the optimal parameter gamma.}
  \item{K2}{universal constant used in the heuristic formula of the optimal parameter p.}
  \item{which.estim}{specifies the estimation method of the parameters: either "Gauss"(default) or "general".}
  \item{RATIO}{optional parameter used in estimation method "Gauss"}
  \item{randomize}{optional parameter used in the estimation method "general".}
  \item{sub.n}{optional parameter used in the estimation method "general" if randomize=TRUE.}
}
\details{
  This function uses a heuristic formula to determine the optimal parameter values gamma and p, in the case when a Gaussian kernel is used. This formula is of the form  \eqn{gamma = K1 * |f|_2^{2/(d+2)} * n^{1/(d+2)}} and 
  \eqn{p = ceil(K2 * |f|_2^{2/(d+2)} * n^{2/(d+2)} )}, where \eqn{|f|_2} is the L2-norm of the density function of non-outliers \eqn{f} and \eqn{ceil(x)} denotes the smallest integer larger than \eqn{x}.

Two methods are proposed to estimate \eqn{|f|_2} and are specified by the argument which.estim: "Gauss" and "general". 

If which.estim="Gauss", the estimation is done as though \eqn{f} was a Gaussian density, which yields \eqn{|f|_2^{2/(d+2)} ) = (4*pi)^{-0.5}*exp(0.5*mean(log(1/ev)))}, where \eqn{ev} are the covariance eigenvalues of the non-outlier distribution. Note that the eigenvalues smaller than \eqn{ev[1]*RATIO}  (where \eqn{ev[1]} is the largest eigenvalue) are discarded to avoid numerical issues.

If which.estim="general", \eqn{|f|_2} is estimated without any assumption on \eqn{f}. However this method may fail in very high dimensions because of the dimensionality curse, since it relies on an estimation of the derivative of \eqn{F} at \eqn{0} where \eqn{F} is the cdf of the pairwise distance between two non-outliers. . Besides, to shorten the computation time, the optional argument 'randomize' can be set as TRUE, so that only a subset of size sub.n of the data is considered to estimate the cdf \eqn{F}.


}
\value{
\item{gamma.opt}{optimal value for gamma.}
\item{p.opt}{optimal value for p.}
\item{est.f2.pw}{estimation of \eqn{|f|_2^{2/(d+2)} }.} 
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{rapkod}}
}
\examples{
data(iris)

##Define data frame with non-outliers
inliers = iris[sample(which(iris$Species!="setosa"), 100, replace=FALSE),
                                              -which(names(iris)=="Species")]
                                              
param <- od.opt.param(inliers)

#display optimal gamma
param$gamma.opt
#display optimal p
param$p.opt

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
