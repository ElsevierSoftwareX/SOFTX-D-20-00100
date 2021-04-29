\name{LOD}
\docType{data}
\alias{LOD}
\title{Length of Day Data}
\description{
The length-of-day dataset \code{comb2000_daily.eop} is available in the ftp site \cr
\code{ftp://euler.jpl.nasa.gov/keof/combinations/2000/}. The length-of-day (LOD) data was analyzed in Huang et al. (2003).
}
\usage{
data(LOD)
}
\format{
A list of LOD, YEAR, MONTH and DATE 
}
\references{
Huang, N.~E., Wu, M.~C., Long, S.~R., Shen, S., Qu, W., Gloerson, P. and Fan, K.~L. (2003)  
A confidence limit for the empirical mode decomposition and Hilbert spectral analysis. 
\emph{Proceedings of the Royal Society London A.}, \bold{459}, 2317--2345. 
}
\examples{
data(LOD)
names(LOD)

xt <- LOD$LOD[LOD$YEAR >= 1981 & LOD$YEAR <= 2000] # From 1981/1/1 to 2000/12/31
xt <- xt/10^4 # measured in millisecond

# EP transform for LOD
outLOD <- eptransf(signal=xt, tau=15, process=c("envelope", "average"), boundary="none")

# outLOD$EpM : candidate of remaining component 
eptplot(outLOD)

op <- par(mfcol=c(3,1), mar=c(2,2,2,1))
plot(xt, type='l', main="LOD", xlab="", ylab="", ylim=range(xt))
plot(xt - outLOD$EpM, type='l', main="candidate of frequency component 
    with half month period", xlab="", ylab=""); abline(h=0, lty=3)
plot(outLOD$EpM, type='l', main="candidate of remaining component", 
    xlab="", ylab="", ylim=range(xt))

\donttest{
# sifting
LODdecom1 <- eptdecomp(signal=xt, tau=15, process=c("envelope", "average"),  
    boundary="none", tol=sd(xt)*0.1^3, maxiter = 30)

# extraction of frequency component with half month period 
plot(xt, type='l', main="LOD", xlab="", ylab="", ylim=range(xt))
plot(LODdecom1$FC, type='l', main="frequency component 
    with half month period", xlab="", ylab=""); abline(h=0, lty=3)
plot(LODdecom1$residue, type='l', main="remaining component", 
    xlab="", ylab="", ylim=range(xt))

# EP transform for remaining signal from LODdecom1
outLOD2 <- eptransf(signal=LODdecom1$residue, tau=30, process=c("envelope", "average"),
    boundary="none")

# outLOD2$EpM : candidate of remaining component for residue signal from LODdecom1
plot(LODdecom1$residue, type='l', main="remaining component from LODdecom1", 
    xlab="", ylab="", ylim=range(xt))
plot(LODdecom1$residue - outLOD2$EpM, type='l', main="candidate of frequency component
    with one month period", xlab="", ylab=""); abline(h=0, lty=3)
plot(outLOD2$EpM, type='l', main="candidate of remaining component", 
    xlab="", ylab="", ylim=range(xt))

# sifting
LODdecom2 <- eptdecomp(signal=LODdecom1$residue, tau=30, process=c("envelope", "average"),
    boundary="none", tol=sd(xt)*0.1^3, maxiter = 50)

# extraction of frequency component with one month period 
plot(LODdecom1$residue, type='l', main="remaining component from LODdecom1", 
    xlab="", ylab="", ylim=range(xt))
plot(LODdecom2$FC, type='l', main="frequency component with one month period", 
    xlab="", ylab=""); abline(h=0, lty=3)
plot(LODdecom2$residue, type='l', main="remaining component", xlab="", ylab="",
    ylim=range(xt))

### Decomposition Result
ttt <- paste(LOD$YEAR, LOD$MONTH, LOD$DATE, sep="/")
ttt <- ttt[LOD$YEAR >= 1981 & LOD$YEAR <= 2000]
ttt <- as.Date(ttt)

att <- as.Date(c("1981/1/1", "1982/1/1", "1983/1/1", "1984/1/1", "1985/1/1", "1986/1/1",
                 "1987/1/1", "1988/1/1", "1989/1/1", "1990/1/1", "1991/1/1", "1992/1/1", 
                 "1993/1/1", "1994/1/1", "1995/1/1", "1996/1/1", "1997/1/1", "1998/1/1", 
                 "1999/1/1", "2000/1/1", "2001/1/1"))

plot(ttt, xt, type='l', main="LOD", xlab="", ylab="", xaxt = "n")
axis(1, at=att, labels=seq(1981, 2001)) 

plot(ttt, LODdecom1$FC, type="l", main="component with half month period by EPT", 
    xlab="", ylab="", xaxt = "n")
axis(1, at=att, labels=seq(1981, 2001)); abline(h=0, lty=3)

plot(ttt, LODdecom2$FC, type="l", main="component with one month period by EPT", 
    xlab="", ylab="", xaxt = "n")
axis(1, at=att, labels=seq(1981, 2001)); abline(h=0, lty=3)
}
par(op)
}
\keyword{datasets}


 
