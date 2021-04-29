\name{eptdecomp}
\alias{eptdecomp}
\title{Decomposition of a signal by Ensemble Patch Transform}
\description{
This function decomposes a signal into frequency component and residue of ensemble patch transform by sifting process.
}
\usage{
eptdecomp(tindex = NULL, signal, type = "rectangle", tau, 
    process = c("average", "average"), pquantile = c(0, 1), equantile = c(0, 1), 
    gamma = 1, boundary = "symmetric", stoprule = "type1", 
    tol = sd(signal, na.rm = TRUE) * 0.1^2, maxiter = 10, check = FALSE) 
}
\arguments{
\item{tindex}{time index at which a signal is observed. When it is \code{NULL}, the signal is supposed to be equally spaced.}
\item{signal}{a set of data or a signal observed at time \code{tindex}.}
\item{type}{patch type of \code{"rectangle"} or \code{"oval"}.}
\item{tau}{a size parameter for ensemble patch transform.}
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
\item{stoprule}{stopping rule of sifting. The \code{type1} stopping rule indicates that absolute values of ept (ensemble patch transform) component must be less than the user-specified tolerance level.
The \code{type2} is the stopping rule that the difference between ept components at the current and previous iterative steps must be less than the user-specified tolerance level.}
\item{tol}{tolerance for stopping rule of sifting.}
\item{maxiter}{the maximum number of sifting.}
\item{check}{specifies whether the sifting process is displayed. When \code{check = TRUE}, ept component and residue for each sifting step are displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
}
\details{
This function decomposes a signal into frequency component and residue of ensemble patch transform by sifting process for a size parameter.
}
\value{ 
\item{eptcomp}{matrix of ept (ensemble patch transform) component at each sifting step.}
\item{FC}{frequency component of ensemble patch transform by sifting process.}
\item{residue}{residue of ensemble patch transform by sifting process.}
\item{parameters}{a list of input parameters of \code{type}, \code{tau}, \code{process}, \code{pquantile}, \code{equantile}, \code{gamma}, \code{boundary}, 
   and output parameter \code{niter}, the number of sifting.}
}
\seealso{
\code{\link{eptransf}}, \code{\link{meptransf}}.
}
\examples{
#### example : composite of two components having different frequencies
ndata <- 1000 
tindex <- seq(0, 1, length=ndata)
comp1 <- cos(90*pi*tindex)
comp2 <- cos(10*pi*tindex)
f <- comp1 + comp2 

op <- par(mfrow=c(3,1), mar=c(2,2,2,1))
plot(tindex, f, main="a signal", xlab="", ylab="", type='l')
abline(h=0, lty=3)
plot(tindex, comp1, main="high-frequency component", xlab="", ylab="", type='l')
abline(h=0, lty=3)
plot(tindex, comp2, main="low-frequency component", xlab="", ylab="", type='l')
abline(h=0, lty=3)

#### Decomposition by Ensemble Patch Transform  
outdecom <- eptdecomp(signal=f, tau=21, process=c("envelope", "average"), pquantile=c(0, 1))

#### Decomposition result
plot(tindex, f, main="a signal", xlab="", ylab="", type='l'); abline(h=0, lty=3)
plot(outdecom$FC, type='l', main="", xlab="", ylab=""); abline(h=0, lty=3)
title(paste0("high-frequency component, tau=", 21))
lines(comp1, col="red", lty=2, lwd=0.5)
plot(outdecom$residue, type="l", main="", xlab="", ylab=""); abline(h=0, lty=3)
title(paste0("low-frequency component, tau=", 21))
lines(comp2, col="red", lty=2, lwd=0.5)
par(op)
}
\keyword{nonparametric}
