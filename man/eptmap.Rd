\name{eptmap}
\alias{eptmap}
\title{Multiscale Visualization of Ensemble Patch Transform of a Signal}
\description{
This function displays time-scale representation of ensemble patch transform of a signal for a sequence of size parameters.  
}
\usage{
eptmap(eptransf, taus = eptransf$parameters$tau, maptype = c("C", "D", "DC", "DD"), 
    stat = c("pstat", "Epstat", "pM", "EpM", "psd", "Epsd"),
    der = c("time", "tau"), ncolor = 100, ...) 
}
\arguments{
\item{eptransf}{R object of ensemble patch transform by \code{eptransf()} or \code{meptransf()}.}
\item{taus}{specifies size parameters for time-scale visualization.}
\item{maptype}{specifies \code{"C"} for centrality map, \code{"D"} for dispersion map, \code{"DC"} for derivative of centrality map and \code{"DD"} for derivative of dispersion map.}
\item{stat}{\code{"pstat"} for centrality of patch transform, \code{"Epstat"} for centrality of ensemble patch transform, 
\code{"pM"} for mean envelope of patch transform, \code{"EpM"} for mean envelope of ensemble patch transform, 
\code{"psd"} for standard deviation of patch transform and \code{"Epsd"} for standard deviation of ensemble patch transform.}
\item{der}{specifies derivative with respect to \code{"time"} or \code{"tau"}.}
\item{ncolor}{the number of colors (\eqn{\ge 1}) to be in the palette.}
\item{\dots}{graphical parameters for image.}
}
\details{
This function performs multiscale visualization of ensemble patch transform of a signal for a sequence of size parameters. 
This function creates images with \code{heat.colors(ncolor)} colors.
}
\value{
image
}
\seealso{
\code{\link{eptransf}}, \code{\link{meptransf}}, \code{\link{eptplot}}.
}
\examples{
# a doppler signal
n <- 1000
tindex <- seq(0, 1, length=n)
j <- 5
f <- 10 * sqrt(tindex*(1-tindex)) * sin((2*pi*(1+2^((9-4*j)/5))) / (tindex+2^((9-4*j)/5)))

set.seed(7)
fnoise <- f + 0.4 * rnorm(n)

op <- par(mar=c(2,2,2,1))
plot(f, type="l", , xlab="", ylab="", ylim=range(fnoise))
points(fnoise, cex=0.3)

\donttest{
taus <- seq(4, 64, by=4)

# try1 : Multiscale EPT by average patch transform and average ensemble transform 
try1 <- meptransf(tindex=tindex, signal=fnoise, taus=taus, process=c("average", "average"))

par(mfrow=c(2,2))
eptmap(try1, maptype="C", stat="pstat", main="centrality of patch transform")
eptmap(try1, maptype="D", stat="psd", main="standard deviation of patch transform")
eptmap(try1, maptype="C", stat="Epstat", main="centrality of ensemble patch transform")
eptmap(try1, maptype="D", stat="Epsd", main="standard deviation of ensemble patch transform")

eptmap(try1, maptype="DC", stat="Epstat", der="time", 
    main="derivative of centrality w.r.t time")
eptmap(try1, maptype="DC", stat="Epstat", der="tau", 
    main="derivative of centrality w.r.t tau")
eptmap(try1, maptype="DD", stat="Epsd", der="time", 
    main="derivative of standard deviation w.r.t time")
eptmap(try1, maptype="DD", stat="Epsd", der="tau", 
    main="derivative of standard deviation w.r.t tau", ncolor=70)

# try2 : Multiscale EPT by envelope patch transform and average ensemble transform 
try2 <- meptransf(tindex=tindex, signal=fnoise, taus=taus, process=c("envelope", "average"),
    pquantile=c(0, 1))

eptmap(try2, maptype="C", stat="pM", main="mean envelope of patch transform")
eptmap(try2, maptype="C", stat="EpM", main="mean envelope of ensemble patch transform")
eptmap(try2, maptype="DC", stat="EpM", der="time", 
    main="derivative of mean envelope w.r.t time")
eptmap(try2, maptype="DC", stat="EpM", der="tau", 
    main="derivative of mean envelope w.r.t time")
}
par(op)
}
\keyword{nonparametric}


