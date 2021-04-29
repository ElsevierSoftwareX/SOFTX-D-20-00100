\name{meptransf}
\alias{meptransf}
\title{Multiscale Ensemble Patch Transforms of a Signal}
\description{
This function performs multiscale ensemble patch transforms of a signal for a sequence of size parameters.
}
\usage{
meptransf(tindex = NULL, signal, type = "rectangle", taus, 
    process = c("average", "average"), pquantile = c(0, 1), equantile = c(0, 1), 
    gamma = 1, boundary = "symmetric") 
}
\arguments{
\item{tindex}{time index at which a signal is observed. When it is \code{NULL}, the signal is supposed to be equally spaced.}
\item{signal}{a set of data or a signal observed at time \code{tindex}.}
\item{type}{patch type of \code{"rectangle"} or \code{"oval"}.}
\item{taus}{a sequence of size parameters for ensemble patch transform.}
\item{process}{specifies transform types for patch and ensemble processes:
           \code{process[1]} for patch process and process[2] for ensemble process.
           Each process has options of \code{"average"}, \code{"median"}, or \code{"envelope"}.
           Note that when \code{process[1]} is \code{"average"} or \code{"median"}, \code{process[2]} must be \code{"average"} or \code{"median"}. 
           When \code{process[1]} is \code{"envelope"}, lower and upper envelopes are obtained by \eqn{\code{pquantile[1]} \times 100 \%}-quantile 
           and \eqn{\code{pquantile[2]} \times 100 \%}-quantile of patches, respectively.
           When \code{process[2]} is \code{"envelope"}, ensemble lower and upper envelopes are obtained as 
           \eqn{\code{equantile[1]} \times 100 \%}-quantile and \eqn{\code{equantile[2]} \times 100 \%}-quantile  
           of lower and upper envelopes of shifted patches, respectively.}
\item{pquantile}{quantiles for lower and upper envelopes of patch transform. When it is \code{c(0, 1)},
   minimum and maximum of a patch are used for lower and upper envelopes, respectively.}
\item{equantile}{quantiles for lower and upper envelopes of ensemble patch transform.}            
\item{gamma}{controls the amount of envelope magnitude.}
\item{boundary}{specifies boundary condition from \code{"symmetric"}, \code{"periodic"} or \code{"none"}.}
}
\details{
This function performs multiscale ensemble patch transforms of a signal for a sequence of size parameters \code{taus}, and 
produces statistics and envelopes for ensemble patch transform.
When \code{process[1]} is \code{"average"} or \code{"median"}, outputs related to envelopes are defined as \code{NULL}.
When \code{process[2]} is \code{"envelope"}, outputs, \code{pstat} and \code{Epstat}, are defined as \code{NULL}.
}
\value{ 
\item{tindex}{time index at which a signal is observed.}
\item{signal}{a set of data or a signal observed at time \code{tindex}.}
\item{pstat}{matrix of centrality of patch transform for a sequence of size parameters \code{taus}.}
\item{Epstat}{matrix of centrality of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{psd}{matrix of standard deviation of patch transform for a sequence of size parameters \code{taus}.}
\item{Epsd}{matrix of standard deviation of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{pL}{matrix of lower envelope of patch transform for a sequence of size parameters \code{taus}.}
\item{pU}{matrix of upper envelope of patch transform for a sequence of size parameters \code{taus}.}
\item{pM}{matrix of mean envelope, \code{(pL + pU) / 2}, of patch transform for a sequence of size parameters \code{taus}.}
\item{pR}{matrix of distance between lower and upper envelopes, \code{(pU - pL)}, of patch transform for a sequence of size parameters \code{taus}.}
\item{EpL}{matrix of lower envelope of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{EpU}{matrix of upper envelope of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{EpM}{matrix of mean envelope, \code{(EpL + EpU) / 2}, of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{EpR}{matrix of distance between lower and upper envelopes, \code{(EpU - EpL)}, of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{rho}{vector of correlations between \code{(signal - ept)} component and \code{ept} component for a sequence of size parameters \code{taus}.    
       The \code{ept} component is component obtained by ensemble patch transform.}
\item{parameters}{a list of input parameters of \code{type}, \code{taus}, \code{process}, \code{pquantile}, \code{equantile}, \code{gamma}, and \code{boundary}.}
\item{nlevel}{the number of size parameters \code{taus}.}
}
\seealso{
\code{\link{eptransf}}, \code{\link{eptdecomp}}.
}
\examples{
#### example : composite of two components having different frequencies
ndata <- 1000 
tindex <- seq(0, 1, length=ndata)
comp1 <- cos(45*pi*tindex)
comp2 <- cos(6*pi*tindex)
f <- comp1 + comp2 

op <- par(mfcol=c(3,1), mar=c(2,2,2,1))
plot(tindex, f, main="a signal", xlab="", ylab="", type='l')
abline(h=0, lty=3)
plot(tindex, comp1, main="high-frequency component", xlab="", ylab="", type='l')
abline(h=0, lty=3)
plot(tindex, comp2, main="low-frequency component", xlab="", ylab="", type='l')
abline(h=0, lty=3)

\donttest{
#### Multiscale Ensemble Patch Transform according to tau's 
taus1 <- seq(20, 60, by=2)
outmulti <- meptransf(signal=f, taus=taus1, process=c("envelope", "average"),
    pquantile=c(0, 1)) 

#### To continue, click the plot in case of "locator(1)".
par(mfrow=c(2,2), mar=c(2,2,2,1))
for (i in 1:length(taus1)) {
  plot(f - outmulti$EpM[,i], type='l', main="", xlab="", ylab=""); abline(h=0, lty=3)
  title(paste0("Remaining component for tau=", taus1[i]))
  lines(comp1, col="red", lty=2, lwd=0.5)
  plot(outmulti$EpM[,i], type="l", main=, xlab="", ylab=""); abline(h=0, lty=3)
  title(paste0("Mean envelope of ensemble patch transform for tau=", taus1[i]))
  lines(comp2, col="red", lty=2, lwd=0.5); locator(1)  
}
}
par(op)
}
\keyword{nonparametric}
