\name{RaPKod-package}
\alias{RaPKod-package}
\alias{RaPKod}
\docType{package}
\title{
\packageTitle{RaPKod}
}
\description{
The RaPKod package implements a kernel method made for outlier detection. Namely, given a data set of reference typical observation (non-outliers or inliers), it tests each new observation in an online way to determine whether it is an outlier or not. This method uses random low-dimensional projections in a kernel space to build a test statistic whose asymptotic null-distribution (ie when the tested observation is not an outlier) is known. The RaPKod method has two parameters: gamma - the hyperparameter of the (Gaussian) kernel used - and p - the dimensionality of the random projection in the kernel space. %The optimal values of these parameters are provided by the function "od.opt.param" from this package. 
}
\details{
The package consists of two functions: the main function "rapkod" and the auxilary function "od.opt.param" which computes optimal parameters values in RaPKod.
%The DESCRIPTION file:
%\packageDESCRIPTION{RaPKod}
%\packageIndices{RaPKod}
%
}
\author{
\packageAuthor{RaPKod}

Maintainer: \packageMaintainer{RaPKod}
}
\references{
Kellner J., "Gaussian Models and Kernel Methods", PhD thesis, Universite des Sciences et Technologies de Lille (defended on December 1st, 2016)
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
%\keyword{ package }
\seealso{
\code{\link{rapkod}}, \code{\link{od.opt.param}} 
}
%\examples{
%~~ simple examples of the most important functions ~~
%}
