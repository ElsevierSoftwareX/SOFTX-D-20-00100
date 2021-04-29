\name{eptplot}
\alias{eptplot}
\title{Plot of Components by Ensemble Patch Transform of a Signal}
\description{
This function plots ensemble patch transform of a signal for a sequence of size parameters tau's.  
}
\usage{
eptplot(eptransf, taus = eptransf$parameters$tau)
}
\arguments{
\item{eptransf}{R object of ensemble patch transform by \code{eptransf()} or \code{meptransf()}.}
\item{taus}{specifies size parameters for which ensemble patch transform of a signal is displayed.}
}
\details{
This function plots ensemble patch transform of a signal for a sequence of size parameters \code{taus}.  
}
\value{
plot
}
\seealso{
\code{\link{eptransf}}, \code{\link{meptransf}}, \code{\link{eptmap}}.
}
\examples{
n <- 500                                                                      
set.seed(1)
x <- c(rnorm(n), arima.sim(list(order = c(1,0,0), ar = 0.9), n = n, sd=sqrt(1-0.9^2)))

taus <- seq(10, 100, by=10)

# eptr1 : Multiscale EPT by average patch transform and average ensemble transform 
eptr1 <- meptransf(tindex=1:(2*n), signal=x, taus=taus, process=c("average", "average"), 
    boundary="none")
names(eptr1)

op <- par(mfcol=c(4,1), mar=c(4,2,2,0.1))
plot(x, xlab="", type="l", main="signal")

eptplot(eptr1)
eptplot(eptr1, taus=20)
eptplot(eptr1, taus=c(20, 30))
lines(eptr1$Epstat[, 2], col="blue")
lines(eptr1$Epstat[, 3], col="red")

# eptr2 : Multiscale EPT by envelope patch transform and average ensemble transform 
eptr2 <- meptransf(tindex=1:(2*n), signal=x, type="oval", taus=taus, 
    process=c("envelope", "average"), pquantile=c(0,1), gamma=0.06, boundary="none")
names(eptr2)

plot(x, xlab="", type="l")
eptplot(eptr2)
eptplot(eptr2, taus=20)
eptplot(eptr2, taus=c(20, 30))
lines(eptr2$EpM[, 2], col="blue")
lines(eptr2$EpM[, 3], col="red")
par(op)
}
\keyword{nonparametric}
