\name{meptransf2d}
\alias{meptransf2d}
\title{Multiscale Two-dimensional Ensemble Patch Transforms of an Image}
\description{
This function performs multiscale two-dimensional ensemble patch transforms of an image for a sequence of size parameters.
}
\usage{
meptransf2d(x = NULL, y = NULL, z, type = "rectangle", taus, theta = 0, 
    process = c("average", "average"), pquantile = c(0, 1), equantile = c(0, 1),
    gamma = 1, boundary = "reflexive") 
}
\arguments{
\item{x, y}{locations of regular grid at which the values in image \code{z} are measured. When those are \code{NULL}, the image is supposed to be equally spaced.}
\item{z}{matrix of an image observed at location \code{(x, y)}.}
\item{type}{patch type of \code{"rectangle"} or \code{"oval"}.}
\item{taus}{a matrix or vector of size parameters for two-dimensional ensemble patch transform. When it is a matrix, the first and second columns specify
     the horizontal and vertical sizes of a two-dimensional patch, respectively.
     When it is a vector, the horizontal and vertical size of a two-dimensional patch are the same.}
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
This function performs multiscale two-dimensional ensemble patch transforms of an image for a sequence of size parameters \code{taus}, and 
produces statistics and envelopes for two-dimensional ensemble patch transform.
When \code{process[1]} is \code{"average"} or \code{"median"}, outputs related to envelopes are defined as \code{NULL}.
When \code{process[2]} is \code{"envelope"}, outputs, \code{pstat} and \code{Epstat}, are defined as \code{NULL}.
}
\value{ 
\item{x, y}{locations of regular grid at which the values in \code{z} are measured. When those are \code{NULL}, image is supposed to be equally spaced.}
\item{z}{matrix of an image observed at \code{(x, y)}.}
\item{pstat}{list of centrality of patch transform for a sequence of size parameters \code{taus}.}
\item{Epstat}{list of centrality of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{psd}{list of standard deviation of patch transform for a sequence of size parameters \code{taus}.}
\item{Epsd}{list of standard deviation of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{pL}{list of lower envelope of patch transform for a sequence of size parameters \code{taus}.}
\item{pU}{list of upper envelope of patch transform for a sequence of size parameters \code{taus}.}
\item{pM}{list of mean envelope, \code{(pL + pU) / 2}, of patch transform for a sequence of size parameters \code{taus}.}
\item{pR}{list of distance between lower and upper envelopes, \code{(pU - pL)}, of patch transform for a sequence of size parameters \code{taus}.}
\item{EpL}{list of lower envelope of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{EpU}{list of upper envelope of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{EpM}{list of mean envelope, \code{(EpL + EpU) / 2}, of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{EpR}{list of distance between lower and upper envelopes, \code{(EpU - EpL)}, of ensemble patch transform for a sequence of size parameters \code{taus}.}
\item{rho}{vector of correlations between \code{(z - ept)} component and \code{ept} component for a sequence of size parameters \code{taus}.    
       The \code{ept} component is component obtained by ensemble patch transform.}
\item{parameters}{a list of input parameters of \code{type}, \code{taus}, \code{theta}, \code{process}, \code{pquantile}, \code{equantile}, \code{gamma}, and \code{boundary}.}
\item{nlevel}{the number of size parameters \code{taus}.}
}
\seealso{
\code{\link{eptransf2d}}, \code{\link{eptdecomp2d}}.
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
#### Multiscale Ensemble Patch Transform according to tau's 
taus1 <- seq(6, 12, by=2)
outcosmulti <- meptransf2d(z=cosmeanf, taus=taus1)

par(mfrow=c(length(taus1), 2), mar=c(2,2,2,1))
for (i in 1:length(taus1)) {
    estlowfreq <- outcosmulti$Epstat[[i]]
    image(estlowfreq, xlab="", ylab="", col=gray(0:100/100), axes=FALSE, 
        main=paste0("ensemble average of patch mean, tau=", taus1[i]))
    persp(estlowfreq, theta = -30, phi = 45, col = "white", xlab="X", ylab="Y", 
        main=paste0("ensemble average of patch mean, tau=", taus1[i]))
}
}
par(op)
}
\keyword{nonparametric}
