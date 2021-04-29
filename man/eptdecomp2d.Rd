\name{eptdecomp2d}
\alias{eptdecomp2d}
\title{Decomposition of an Image by Two-dimensional Ensemble Patch Transform}
\description{
This function decomposes an image into frequency component and residue of two-dimensional ensemble patch transform by sifting process.
}
\usage{
eptdecomp2d(x = NULL, y = NULL, z, type = "rectangle", tau, theta = 0, 
    process = c("average", "average"), pquantile = c(0, 1), equantile = c(0, 1), 
    gamma = 1, boundary = "reflexive", 
    stoprule = "type2", tol = 0.1^2, maxiter = 10, check = FALSE) 
}
\arguments{
\item{x, y}{locations of regular grid at which the values in image \code{z} are measured. When it is \code{NULL}, the image is supposed to be equally spaced.}
\item{z}{matrix of an image observed at location \code{(x, y)}.}
\item{type}{patch type of \code{"rectangle"} or \code{"oval"}.}
\item{tau}{a size parameter for two-dimensional ensemble patch transform: \code{tau[1]} for horizontal size and tau[2] for vertical size of a two-dimensional patch.
     When \code{length(tau)} is 1, the horizontal and vertical size are the same.}
\item{theta}{a degree of clockwise rotation of a patch.}
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
\item{boundary}{specifies boundary condition from \code{"reflexive"}, \code{"periodic"} or \code{"none"}.}
\item{stoprule}{stopping rule of sifting. The \code{type1} stopping rule indicates that absolute values of ept  (ensemble patch transform) component must be less than the user-specified tolerance level.
The \code{type2} is the stopping rule that the difference between ept components at the current and previous iterative steps must be less than the user-specified tolerance level.}
\item{tol}{tolerance for stopping rule of sifting.}
\item{maxiter}{the maximum number of sifting.}
\item{check}{specifies whether the sifting process is displayed. When \code{check = TRUE}, ept component and residue for each sifting step are displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
}
\details{
This function decomposes an image into frequency component and residue of two-dimensional ensemble patch transform by sifting process for a size parameter.
}
\value{ 
\item{eptcomp}{list of ept (ensemble patch transform) component at each sifting step when \code{check=TRUE}.}
\item{FC}{frequency component of ensemble patch transform by sifting process.}
\item{residue}{residue of ensemble patch transform by sifting process.}
\item{parameters}{a list of input parameters of \code{type}, \code{tau}, \code{theta}, \code{process}, \code{pquantile}, \code{equantile}, \code{gamma}, \code{boundary}, 
   and output parameter \code{niter}, the number of sifting.}
}
\seealso{
\code{\link{eptransf2d}}, \code{\link{meptransf2d}}.
}
\examples{
#### example : composite of two components having different frequencies
nr <- nc <- 128; x <- seq(0, 1, length=nr); y <- seq(0, 1, length=nc)

coscomp1 <- outer(cos(20 * pi * x), cos(20 * pi * y))
coscomp2 <- outer(cos(5* pi * x), cos(5 * pi * y))
cosmeanf <- coscomp1 + coscomp2

op <- par(mfcol=c(3,1), mar=c(0,0.5,2,0.5))
image(cosmeanf, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, main="a composite image")
image(coscomp1, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, main="high-frequency component")
image(coscomp2, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, main="low-frequency component")

\donttest{
#### Decomposition by Ensemble Patch Transform  
outcossift <- eptdecomp2d(z=cosmeanf, tau=8)

#### Decomposition Result
par(mfrow=c(2,2), mar=c(2,2,2,1))
image(outcossift$FC, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, main="Decomposed HF")
persp(outcossift$FC, theta = -30, phi = 45, col = "white", xlab="X", ylab="Y", main="Decomposed HF")
image(outcossift$residue, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, main="Residue")
persp(outcossift$residue, theta = -30, phi = 45, col = "white", xlab="X", ylab="Y", main="Residue")
}
par(op)
}
\keyword{nonparametric}
