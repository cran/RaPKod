\name{rapkod}
\alias{rapkod}
\title{RaPKod: Random Projections Kernel Outlier Detection}
\description{RaPKod is a kernel method for detecting outliers in a given dataset on the basis of a reference set of non-outliers. To do so, it 'transforms' a tested observation into some kernel space (through a 'feature map') and then projects it onto a random low-dimensional subspace of this kernel space. Since the distribution of this projection is known in the case of a non-outlier, it allows RaPKod to control the probability of false alarm error (ie labelling a non-outlier as an outlier).}
\usage{
  rapkod(X, given.kern = FALSE, ref.n=NULL, gamma=NULL,  p=NULL, alpha = 0.05, 
	use.tested.inlier = FALSE, lowrank = "No", r.lowrk = ceiling(sqrt(nrow(X))), 
	K1 = 6, K2 = 50)
}
\arguments{
  \item{X}{either a data frame or an n x d matrix (if given.kern=FALSE), otherwise an n x n kernel matrix (if given.kern=TRUE). In the former case, a Gaussian kernel is used by default.}
  
  \item{given.kern}{If FALSE (default), each row of X is an observation. Otherwise X is a kernel matrix (in this case, gamma and p must be user-specified).}
  
  \item{ref.n}{the size of the reference non-outlier dataset. Must be smaller than n.}
  
  \item{gamma}{the hyperparameter of the Gaussian kernel \eqn{k(x, y) = exp( - gamma * || x - y ||^2)}. Set automatically by the program if not specified and given.kern=FALSE.}
  
  \item{p}{the number of dimensions of the projection made in the kernel space. Set automatically by the program if not specified and given.kern=FALSE.}
  
  \item{alpha}{the prescribed probability of false alarm error.}
  
  \item{use.tested.inlier}{If TRUE, each tested observation that is labelled as a non-outlier is appended to the reference dataset of non-outliers (the 'oldest' reference non-outlier is discarded). Set to FALSE by default.}
  
  \item{lowrank}{if lowrank="No" (default), the full kernel matrix is used. Otherwise, a low-rank approximation of the kernel matrix is computed: if "Nyst", it is approximated through Nystrom method; if "RKS", it is approximated by random Kitchen Sinks (in this case, X must be a dataset matrix, not a kernel matrix)}
  
  \item{r.lowrk}{if lowrank="Nyst" or "RKS", specifies the (low) rank of the approximated kernel matrix.}
  
  \item{K1}{universal constant used in the heuristic formula of the optimal parameter gamma.}
  
  \item{K2}{universal constant used in the heuristic formula of the optimal parameter p.}
  
 }
\details{
%%  ~~ If necessary, more details than the description above ~~
  If given.kern = FALSE, X is a dataset matrix whose first ref.n rows corresponds to the reference dataset of non-outliers. The (n - ref.n) other observations will be tested one by one by RaPKod to determine whether they are outliers or not. 
  
  If given.kern = TRUE, X must be a n x n Gram matrix. The kernel used to compute this Gram matrix should be of the form \eqn{k(x, y) = K(gamma * || x - y ||^2)} where \eqn{K} is a positive function. Also note that in this case, the parameters gamma and p must be specified by the user.%Also note that the default values of the universal constants K1 and K2 used to determine the optimal gamma and p have been calibrated for the use of a Gaussian kernel.
  }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{stats }{ a vector of length (n - ref.n) containing the test statistics for each tested observation.}
  
  \item{flag }{ a vector of length (n - ref.n) indicating which observations have been labelled as an outlier (TRUE in this case).}
  
  \item{pv }{ a vector of length (n - ref.n) containing p-values for each tested observation.}
  
  \item{gamma }{ the optimal value of gamma determined by the program (or the value provided by the user if it was user-specified).}
  
  \item{p }{ the optimal value of p determined by the program (or the value provided by the user if it was user-specified).}
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
\code{\link{od.opt.param}}
}
\examples{
data(iris)

##Define data frame with non-outliers
inliers = iris[sample(which(iris$Species!="setosa"), 100, replace=FALSE),
                                              -which(names(iris)=="Species")]
##Define data frame with outliers
outliers = iris[which(iris$Species=="setosa"),-which(names(iris)=="Species")]


X = rbind(inliers, outliers)

ref.n = 50
result <- rapkod(X, ref.n = ref.n, use.tested.inlier = FALSE, alpha = 0.05)


##False alarm error ratio obtained on tested non-outliers (should be close to 0.05)
mean(result$pv[1:(nrow(inliers)-ref.n)]<0.05, na.rm = TRUE)
##Missed detection error ratio obtained on tested outliers (should be close to 0)
mean(result$pv[-(1:(nrow(inliers)-ref.n))]>0.05, na.rm = TRUE)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
