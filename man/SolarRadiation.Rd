\name{SolarRadiation}
\docType{data}
\alias{SolarRadiation}
\title{Solar Radiation}
\description{
The solar radiations were hourly observed at Seoul, Daegu, and Busan in South Korea from September 1, 2003 to September 29, 2003. 
The data are available from Korea Meteorological Administration (\url{https://data.kma.go.kr}). 
Daegu and Busan, located in the southeast of the Korean Peninsula, are close to each other geographically, whereas Seoul is located in the middle of the Peninsula.
In addition, note that Daegu and Busan were severely damaged by a typhoon named ``MAEMI" at that time, while Seoul was hardly affected by the typhoon. 
It is expected that the climatic characteristics of Daegu and Busan are similar, and the pattern of Seoul seems to be different from the other two cities. 
}
\usage{
data(SolarRadiation)
}
\format{
A daraframe of Date, Seoul, Daegu and Busan.
}
\examples{
data(SolarRadiation)
names(SolarRadiation)

# ensemble patch transform
SolarRadiationEpU <- SolarRadiationEpL <- NULL

for(j in 1:3) {
    tmp <- eptransf(signal=SolarRadiation[,j+1], tau=24, process=c("envelope", "average"), 
        pquantile=c(0, 1), gamma=0)
    SolarRadiationEpU <- cbind(SolarRadiationEpU, tmp$EpU)
    SolarRadiationEpL <- cbind(SolarRadiationEpL, tmp$EpL)
}

# Correlation of the solar radiations at Seoul, Daegu, and Busan
cor(SolarRadiation[, 2:4])

# Correlation of ensemble average of upper envelope
cor(SolarRadiationEpU)

op <- par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(SolarRadiation[,2], type='l', main="(a) solar-radiation in Seoul and upper envelope",
    ylim=c(0, 3.3), xaxt="n"); axis(1, at=seq(1, 30*24, by=24), labels=seq(1, 30, by=1)) 
lines(SolarRadiationEpU[,1], lty=2); lines(SolarRadiationEpL[,1], lty=2)
plot(SolarRadiation[,3], type='l', main="(b) solar-radiation in Daegu and upper envelope",
    ylim=c(0, 3.3), xaxt="n"); axis(1, at=seq(1, 30*24, by=24), labels=seq(1, 30, by=1)) 
lines(SolarRadiationEpU[,2], lty=2); lines(SolarRadiationEpL[,2], lty=2)
plot(SolarRadiation[,4], type='l', main="(c) solar-radiation in Busan and upper envelope",
    ylim=c(0, 3.3), xaxt="n"); axis(1, at=seq(1, 30*24, by=24), labels=seq(1, 30, by=1)) 
lines(SolarRadiationEpU[,3], lty=2); lines(SolarRadiationEpL[,3], lty=2)
par(op)
}
\keyword{datasets}


 
