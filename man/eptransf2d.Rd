\name{eptransf2d}
\alias{eptransf2d}
\title{Two-dimensional Ensemble Patch Transform of an Image} 
\description{
This function performs two-dimensional ensemble patch transform of an image for a size parameter.
}
\usage{ 
eptransf2d(x = NULL, y = NULL, z, type = "rectangle", tau, theta = 0, 
    process = c("average", "average"), pquantile = c(0, 1), equantile = c(0, 1), 
    gamma = 1, boundary = "reflexive") 
}
\arguments{
\item{x, y}{locations of regular grid at which the values in image \code{z} are measured. When those are \code{NULL}, the image is supposed to be equally spaced.}
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
}
\details{
This function performs two-dimensional ensemble patch transform of an image for a size parameter \code{tau}, and 
produces statistics and envelopes for two-dimensional ensemble patch transform.
When \code{process[1]} is \code{"average"} or \code{"median"}, outputs related to envelopes are defined as \code{NULL}.
When \code{process[2]} is \code{"envelope"}, outputs, \code{pstat} and \code{Epstat}, are defined as \code{NULL}.
}
\value{ 
\item{x, y}{locations of regular grid at which the values in image \code{z} are measured. When it is \code{NULL}, the image is supposed to be equally spaced.}
\item{z}{matrix of an image observed at location \code{(x, y)}.}
\item{pstat}{centrality of patch transform for size parameter \code{tau}.}
\item{Epstat}{centrality of ensemble patch transform for size parameter \code{tau}.}
\item{psd}{standard deviation of patch transform for size parameter \code{tau}.}
\item{Epsd}{standard deviation of ensemble patch transform for size parameter \code{tau}.}
\item{pL}{lower envelope of patch transform for size parameter \code{tau}.}
\item{pU}{upper envelope of patch transform for size parameter \code{tau}.}
\item{pM}{mean envelope, \code{(pL + pU) / 2}, of patch transform for size parameter \code{tau}.}
\item{pR}{distance between lower and upper envelopes, \code{(pU - pL)}, of patch transform for size parameter \code{tau}.}
\item{EpL}{lower envelope of ensemble patch transform for size parameter \code{tau}.}
\item{EpU}{upper envelope of ensemble patch transform for size parameter \code{tau}.}
\item{EpM}{mean envelope, \code{(EpL + EpU) / 2}, of ensemble patch transform for size parameter \code{tau}.}
\item{EpR}{distance between lower and upper envelopes, \code{(EpU - EpL)}, of ensemble patch transform for size parameter \code{tau}.}
\item{parameters}{a list of input parameters of \code{type}, \code{tau}, \code{theta}, \code{process}, \code{pquantile}, \code{equantile}, \code{gamma}, and \code{boundary}.}
\item{nlevel}{the number of size parameter \code{tau}. For \code{eptransf2d()}, \code{nlevel} is 1.}
}
\seealso{
\code{\link{meptransf2d}}, \code{\link{eptdecomp2d}}.
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
#### Ensemble average of Ensemble Patch Transform
outcos <- eptransf2d(z=cosmeanf, tau=12)
rangez <- range(cosmeanf) 

par(mfrow=c(3,1), mar=c(2,2,2,1))
image(outcos$Epstat, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, zlim=rangez, 
      main="ensemble average of patch mean")
persp(outcos$Epstat, theta = -30, phi = 45, col = "white", xlab="X", ylab="Y", 
      main="ensemble average of patch mean")
image(outcos$Epsd, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, 
      main="ensemble average of standard deviation")

#### Ensemble Envelope of Ensemble Patch Transform
outcos2 <- eptransf2d(z=cosmeanf, tau=12, process = c("envelope", "average"))
par(mfrow=c(2,2), mar=c(2,2,2,1))
image(outcos2$EpL, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, 
      main="ensemble average of lower envelope")
image(outcos2$EpU, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, 
      main="ensemble average of upper envelope")
image(outcos2$EpM, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, 
      main="ensemble average of mean envelope")
image(outcos2$Epsd, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, 
      main="ensemble average of standard deviation")
}
par(op)
}
\keyword{nonparametric}
