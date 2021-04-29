\name{eptransf}
\alias{eptransf}
\title{Ensemble Patch Transform of a Signal}
\description{
This function performs ensemble patch transform of a signal for a size parameter.
}
\usage{
eptransf(tindex = NULL, signal, type = "rectangle", tau,
    process = c("average", "average"), pquantile = c(0, 1), equantile = c(0, 1), 
    gamma = 1, boundary = "symmetric") 
}
\arguments{
\item{tindex}{time index at which a signal is observed. When it is \code{NULL}, the signal is supposed to be equally spaced.}
\item{signal}{a set of data or a signal observed at time \code{tindex}.}
\item{type}{patch type of \code{"rectangle"} or \code{"oval"}.}
\item{tau}{size parameter for ensemble patch transform.}
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
This function performs ensemble patch transform of a signal for a size parameter \code{tau}, and 
produces statistics and envelopes for ensemble patch transform.
When \code{process[1]} is \code{"average"} or \code{"median"}, outputs related to envelopes are defined as \code{NULL}.
When \code{process[2]} is \code{"envelope"}, outputs, \code{pstat} and \code{Epstat}, are defined as \code{NULL}.
}
\value{ 
\item{tindex}{time index at which a signal is observed.}
\item{signal}{a set of data or a signal observed at time \code{tindex}.}
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
\item{parameters}{a list of input parameters of \code{type}, \code{tau}, \code{process}, \code{pquantile}, \code{equantile}, \code{gamma}, and \code{boundary}.}
\item{nlevel}{the number of size parameter \code{tau}. For \code{eptransf()} function, \code{nlevel} is 1.}
}
\seealso{
\code{\link{meptransf}}, \code{\link{eptdecomp}}.
}
\examples{
# a doppler signal
n <- 256
tindex <- seq(0, 1, length=n)
j <- 5
f <- 10 * sqrt(tindex*(1-tindex)) * sin((2*pi*(1+2^((9-4*j)/5))) / (tindex+2^((9-4*j)/5)))
fnoise <- f + 0.4 * rnorm(n)
  
#### Ensemble statistics
op <- par(mfrow=c(5,3), mar=c(2,2,2,1))
layout(matrix(c(1, 1, 1, 2:13), 5, 3, byrow = TRUE))

plot(f, main="a doppler signal", xlab="", ylab="", type='l', ylim=range(fnoise))
points(fnoise); abline(h=0, lty=3)

#### Ensemble Patch Transform  
taus <- c(5, 10, 20)

out <- list()
for (i in 1:length(taus)) 
    out[[i]] <- eptransf(signal=fnoise, tau=taus[i], process=c("average", "average"))
    
for (i in 1:length(taus)) {
    plot(out[[i]]$Epstat, type="l", xlab="", ylab="",
        main=paste0("ensemble average of patch mean, tau=", taus[i]))
    abline(h=0, lty=3)   
}

for (i in 1:length(taus)) 
    plot(out[[i]]$Epsd, type='l', xlab="", ylab="", 
         main=paste0("ensemble average of standard deviation, tau=", taus[i]))  

out2 <- list()
for (i in 1:length(taus)) 
    out2[[i]] <- eptransf(signal=fnoise, tau=taus[i], process=c("envelope", "average"))
    
for (i in 1:length(taus)) {
    plot(out2[[i]]$EpM, type="l", col="red", xlab="", ylab="", 
         ylim=range(c(out2[[i]]$EpU,out2[[i]]$EpL)),
         main=paste0("ensemble average of mean envelope, tau=", taus[i]))
    points(fnoise, cex=0.1)
    abline(h=0, lty=3); lines(out2[[i]]$EpU); lines(out2[[i]]$EpL)    
}

for (i in 1:length(taus))
    plot(out2[[i]]$EpR, type='l', xlab="", ylab="",
         main=paste0("ensemble average of envelope distance, tau=", taus[i]))
         
par(op)
}
\keyword{nonparametric}
