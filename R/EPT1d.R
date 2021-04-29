
#### The fourth version

eptdecomp <- function(tindex=NULL, signal, type="rectangle", tau, process=c("average", "average"), pquantile=c(0, 1), equantile=c(0, 1), gamma=1, boundary="symmetric",   
                           stoprule = "type1", tol=sd(signal, na.rm=TRUE)*0.1^2, maxiter=10, check=FALSE)
{       
    if (!(any(type == c("rectangle", "oval")))) stop("Patch type should be 'rectangle', or 'oval'." )
    if (floor(tau) < 2) stop(paste("The patch size with tau", tau, "is too small."))  
    if (!any(process[1] == c("average", "median", "envelope"))) stop("The patch transform must be 'average', 'median', or 'envelope'.")
    if (!any(process[2] == c("average", "median", "envelope"))) stop("The ensemnle transform must be 'average', 'median', or 'envelope'.")
    if (any(process[1] == c("average", "median")) & process[2] == "envelope") stop("If the patch transform is 'average' or 'median', the ensemble transform must be 'average' or 'median'.")
    if (pquantile[1] < 0 | pquantile[1] > 1) stop("Quantile for lower envelope of patch transform should be [0, 1].")
    if (pquantile[2] < 0 | pquantile[2] > 1) stop("Quantile for upper envelope of patch transform should be [0, 1].")
    if (pquantile[1] > pquantile[2]) stop("Quantile for lower envelope should be smaller than quantile for upper envelop of patch transform.")
    if (equantile[1] < 0 | equantile[1] > 1) stop("Quantile for lower envelope of ensemble transform should be [0, 1].")
    if (equantile[2] < 0 | equantile[2] > 1) stop("Quantile for upper envelope of ensemble transform should be [0, 1].")
  
    if (gamma < 0) stop("Scale parameter of patch transform should be positive.")
    if (!(any(boundary == c("none", "periodic", "symmetric")))) stop("The boundary condition should be 'none', 'periodic' or 'symmetric'." )
    if (maxiter < 0) stop("The maximum iteration of sifting should be positive.")  
  
    ndata <- length(signal)

    if (check) 
        if (is.null(tindex)) {
            plotindex <- 1:ndata 
        } else 
            plotindex <- tindex    

    if (any(process[1] == c("average", "median"))) {
        eptransform <- "Epstat"
    } else if (process[1] == "envelope") {
        eptransform <- "EpM"
    }
    
    #lower <- upper <- eptcomp <- NULL
    eptcomp <- NULL
    residue <- rep(0, ndata)
    
    iter <- 1
    input <- signal
    
    repeat {
        tmp <- eptransf(tindex=tindex, signal=input, type=type, tau=tau, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary)
        eptcomp <- cbind(eptcomp, tmp[[eptransform]]) 
        residue <- residue + tmp[[eptransform]] 
            
        #if (process[1] == "envelope") {
        #    lower <- cbind(lower, tmp$EpL)
        #    upper <- cbind(upper, tmp$EpU)
        #}

        if (check) {
            #if (process[1] == "envelope") ylimit <- range(c(tmp$EpL, tmp$EpU))
            ylimit <- range(c(input, residue), na.rm=TRUE)
            plot(plotindex, input, type = "l", xlab = "", ylab = "", ylim=ylimit, main=paste0("Frequency component at iter=", iter))
            #if (process[1] == "envelope") {
                #lines(plotindex, tmp$EpL, col = 3)
                #lines(plotindex, tmp$EpU, col = 3)
                #lines(plotindex, tmp$EpM, col = 2)
                #locator(1)
            #}
            
            plot(plotindex, residue, type = "l", xlab = "", ylab = "", ylim=ylimit, main=paste0("Residue at iter=", iter))
            locator(1)    
        }   
        
        if (stoprule == "type1" && (all(abs(eptcomp[, iter]) < tol) || iter >= maxiter)) {
            FC <- signal - residue                
            break
        } else if (stoprule == "type2" && iter >= 2) {
            if (sum((eptcomp[, iter-1] - eptcomp[, iter])^2) < tol * sum(eptcomp[, iter-1]^2) || iter >= maxiter) {
                FC <- signal - residue                
                break
            }
        }
        
        input <- input - tmp[[eptransform]] 
        
        iter <- iter + 1
    }
    
    if (any(process[1] == c("average", "median"))) {
        pquantile <- equantile <- gamma <- NULL 
    } else if (process[1] == "envelope") {
        if (process[2] != "envelope") equantile <- NULL
    }
    
    list(eptcomp=eptcomp, FC = FC, residue = residue, 
             parameters=list(type=type, tau=floor(tau), process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary, niter = iter))
}



eptplot <- function(eptransf, taus=eptransf$parameters$tau)
{
    taus <- sort(taus)
    nlevel <- length(taus)
    
    if (nlevel > eptransf$nlevel) stop ("The number of taus is larger than the number of multiscale levels.")
    if (!all(is.element(taus, eptransf$parameters$tau))) stop("The patch transform is not implemented for given tau's.")
    
    index <- match(taus, eptransf$parameters$tau)
    #colgray <- gray(seq(0.9, 0.45, length=nlevel))
    colgray <- gray(seq(1, 0, length=nlevel))
    
    if (is.null(eptransf$tindex)) {
        plotindex <- 1:length(eptransf$signal)
    } else
        plotindex <- eptransf$tindex
 
    if (any(eptransf$parameters$process[1] == c("average", "median"))) {
        eptcomp <- eptransf$Epstat
        
        if (eptransf$nlevel == 1) {
            plot(plotindex, eptransf$signal, ylim=range(c(eptransf$signal, eptcomp)), type="l", lty=1, lwd=0.5, main="", xlab="", ylab="")
            lines(plotindex, eptcomp, lwd=1.5, col=2)
            mtext(paste("tau = ", taus))
        } else {
            plot(plotindex, eptransf$signal, ylim=range(c(eptransf$signal, eptcomp[, index])), type="l", lty=1, lwd=0.5, main="", xlab="", ylab="")  #ylim=range(eptransf$signal), 
            for (j in nlevel:1) 
                lines(plotindex, eptcomp[, index[j]], lwd=1.5, lty=j+1, col=j+1) #col=colgray[j]) 
            mtext(paste("tau = ", paste(taus, collapse=",")))
        }        
    } else if (eptransf$parameters$process[1] == "envelope") {
        if (eptransf$nlevel== 1) {
            plot(plotindex, eptransf$signal, ylim=range(c(eptransf$EpL, eptransf$EpU)), type="l", lty=1, lwd=0.5, main="", xlab="", ylab="")
            #polygon(c(plotindex, rev(plotindex)), c(eptransf$EpL, rev(eptransf$EpU)), col="gray90", border="gray90")
            lines(plotindex, eptransf$EpM, lwd=1.5, col=2)
            lines(plotindex, eptransf$EpU, lwd=0.5, col=3)
            lines(plotindex, eptransf$EpL, lwd=0.5, col=3)
            mtext(paste("tau = ", taus))
        } else {
            if (nlevel > 1) {
                plot(plotindex, eptransf$signal, ylim=range(c(eptransf$EpL[, index], eptransf$EpU[, index])), type="n", main="", xlab="", ylab="") #ylim=range(eptransf$signal), 
                for (j in nlevel:1) 
                    polygon(c(plotindex, rev(plotindex)), c(eptransf$EpL[, index[j]], rev(eptransf$EpU[, index[j]])), col=colgray[j], border=colgray[j])
                lines(plotindex, eptransf$signal, lty=1, lwd=0.5)
                mtext(paste("tau = ", paste(taus, collapse=",")))
                
                for (j in nlevel:1)  
                    lines(plotindex, eptransf$EpM[,index[j]], lwd=1.5, lty=j+1, col=j+1)       
            } else {
                plot(plotindex, eptransf$signal, ylim=range(c(eptransf$EpL[, index], eptransf$EpU[, index])), 
                     type="l", lty=1, lwd=0.5, main="", xlab="", ylab="")
                #polygon(c(plotindex, rev(plotindex)), c(eptransf$EpL[, index], rev(eptransf$EpU[, index])), col="gray90", border="gray90")
                lines(plotindex, eptransf$EpM[,index], lwd=1.5, col=2)
                lines(plotindex, eptransf$EpU[,index], lwd=0.5, col=3)
                lines(plotindex, eptransf$EpL[,index], lwd=0.5, col=3)
                mtext(paste("tau =", taus))
            }
        }        
    } 
}



   
    
eptmap <- function(eptransf, taus=eptransf$parameters$tau, maptype=c("C", "D", "DC", "DD"), stat=c("pstat", "Epstat", "pM", "EpM", "psd", "Epsd"), der=c("time", "tau"), ncolor=100, ...)
{
    if (!(any(maptype == c("C", "D", "DC", "DD")))) stop("Map type is not properly specified." )
    if (!(any(stat == c("pstat", "Epstat", "pM", "EpM", "psd", "Epsd")))) stop("Summary stat is not properly specified." )
    if (!(any(der == c("time", "tau")))) stop("Derivative should be with respect to time or tau." )
    
    if (maptype == "C" || maptype =="DC") {
        if (any(eptransf$parameters$process[1] == c("average", "median")) & !any(stat == c("pstat", "Epstat"))) 
            stop("Centrality statistics for average or median of patch transform is not properly specified.")    
        if (eptransf$parameters$process[1] == "envelope" & !any(stat == c("pM", "EpM"))) 
            stop("Centrality statistics for mean envelope of patch transform is not properly specified.")     
    }
    if (maptype == "D" || maptype =="DD") {
        if (any(eptransf$parameters$process[1] == c("average", "median")) & !any(stat == c("psd", "Epsd"))) 
            stop("Dispersion statistics for average or median of patch transform is not properly specified.")
        if (eptransf$parameters$process[1] == "envelope" & !any(stat == c("psd", "Epsd"))) 
            stop("Dispersion statistics for mean envelope of patch transform is not properly specified.")     
    }

    taus <- sort(taus)
    nlevel <- length(taus)
    
    if (is.null(eptransf$tindex)) {
        plotindex <- 1:length(eptransf$signal)
    } else
        plotindex <- eptransf$tindex
    ndata <- length(eptransf$signal)

    if (nlevel > eptransf$nlevel) stop ("The number of taus is larger than the number of multiscale levels.")
    if (!all(is.element(taus, eptransf$parameters$tau))) stop("The patch transform is not implemented for given tau's.")
    
    index <- match(taus, eptransf$parameters$tau)
    
    #colgray <- gray(seq(0.9, 0.45, length=nlevel))
    colgray <- gray(seq(1, 0, length=nlevel))
  
    statindex <- match(stat, names(eptransf))
    if (maptype == "C" || maptype =="D") {
        image(plotindex, taus, eptransf[[statindex]][, index], xlab="", ylab="", col=heat.colors(ncolor), ...)      
    } else if (maptype == "DC" || maptype =="DD") {
        if (der == "time") {
            image(plotindex, taus, sweep(eptransf[[statindex]][, index][-1,] - eptransf[[statindex]][, index][-ndata,], 1, diff(plotindex), "/"), xlab="", ylab="", col=heat.colors(ncolor), ...)#, xlim=range(plotindex), ylim=range(taus)) 
        } else if (der == "tau") {
            image(plotindex, taus, sweep(eptransf[[statindex]][, index][,-1] - eptransf[[statindex]][, index][, -nlevel], 2, diff(taus), "/"), xlab="", ylab="", col=heat.colors(ncolor), ...)#, xlim=range(plotindex), ylim=range(taus))            
        }
    }
}

eptransf <- function(tindex=NULL, signal, type="rectangle", tau, process=c("average", "average"), pquantile=c(0, 1), equantile=c(0, 1), gamma=1, boundary="symmetric")
{
    if (!(any(type == c("rectangle", "oval")))) stop("Patch type should be 'rectangle', or 'oval'." )
    if (!any(process[1] == c("average", "median", "envelope"))) stop("The patch transform must be 'average', 'median', or 'envelope'.")
    if (!any(process[2] == c("average", "median", "envelope"))) stop("The ensemnle transform must be 'average', 'median', or 'envelope'.")
    if (any(process[1] == c("average", "median")) & process[2] == "envelope") stop("If the patch transform is 'average' or 'median', the ensemble transform must be 'average' or 'median'.")
    if (pquantile[1] < 0 | pquantile[1] > 1) stop("Quantile for lower envelope of patch transform should be [0, 1].")
    if (pquantile[2] < 0 | pquantile[2] > 1) stop("Quantile for upper envelope of patch transform should be [0, 1].")
    if (pquantile[1] > pquantile[2]) stop("Quantile for lower envelope should be smaller than quantile for upper envelop of patch transform.")
    if (equantile[1] < 0 | equantile[1] > 1) stop("Quantile for lower envelope of ensemble transform should be [0, 1].")
    if (equantile[2] < 0 | equantile[2] > 1) stop("Quantile for upper envelope of ensemble transform should be [0, 1].")
  
    if (gamma < 0) stop("Scale parameter of patch transform should be positive.")
    if (!(any(boundary == c("none", "periodic", "symmetric")))) stop("The boundary condition should be 'none', 'periodic' or 'symmetric'." )
    
    orindata <- length(signal)
    #if (is.null(tindex)) tindex <- 1:orindata    
    
    if (boundary == "symmetric") {
        signal <- c(rev(signal[-1]), signal, rev(signal)[-1])   
        if (!is.null(tindex)) 
            tindex <- c(tindex[1] - rev(cumsum(diff(tindex))), tindex, tindex[orindata] + cumsum(rev(diff(tindex)))) 
        obsindex <- orindata:(2*orindata-1)
    } else if (boundary == "periodic") {
        signal <- c(signal[-ndata], signal, signal[-1])   
        if (!is.null(tindex)) 
            tindex <- c(rev(tindex[1] - cumsum(rev(diff(tindex)))), tindex, tindex[orindata] + cumsum(diff(tindex)))     
        obsindex <- orindata:(2*orindata-1)
    } else {
        signal <- signal
        obsindex <- 1:orindata
    } 
    ndata <- length(signal)   

    pstat <- psd <- pL <- pU <- pM <- pR <- NULL

    tau0 <- floor(tau)
    if (tau0 < 2) stop(paste("The patch size with tau", tau, "is too small."))
    
    movetoleft <- floor((tau0+1)/2) # move the patch at (t, t+1, ..., t+tau0) to the left by movetoleft
    k <- 0:tau0 - movetoleft
    
    startindex <- 1 + (-k[1]); endindex <- ndata - k[tau0+1]
    envindex <- 1:(tau0+1)
    
    if (type == "rectangle") {
        envelopeweight <- gamma * tau0/2 * rep(1, length(k))
    } else if (type == "oval") {
        envelopeweight <- gamma * sqrt(tau0^2/4 - (k + (tau0 %% 2) / 2)^2)  
        #sqrt(tau0^2/4 - c(seq(-tau0/2, 0, length=-k[1]+1), seq(0, tau0/2, length=k[tau0+1]+1)[-1])^2) #gamma * sqrt(tau0^2/4 - k^2)  
    }        
    
    if (startindex > 1) {
        for (i in 1:(startindex-1)) {
            index <- (i + k)[(i + k) > 0]; envindex2 <- envindex[(i + k) > 0]
            if (process[1] == "average") {
                pstat[i] <- mean(signal[index], na.rm=TRUE)        
            } else if (process[1] == "median") {
                pstat[i] <- median(signal[index], na.rm=TRUE)             
            } else if (process[1] == "envelope") {
                 pL[i] <- quantile(signal[index] - envelopeweight[envindex2], pquantile[1], na.rm=TRUE) 
                 pU[i] <- quantile(signal[index] + envelopeweight[envindex2], pquantile[2], na.rm=TRUE)                    
            }
            psd[i] <- sd(signal[index], na.rm=TRUE)
        }           
    }
    if (endindex < ndata) {
        for (i in (endindex+1):ndata) { 
            index <- (i + k)[(i + k) <= ndata]; envindex2 <- envindex[(i + k) <= ndata]
            if (process[1] == "average") {
                 pstat[i] <- mean(signal[index], na.rm=TRUE)      
            } else if (process[1] == "median") {
                 pstat[i] <- median(signal[index], na.rm=TRUE)             
            } else if (process[1] == "envelope") {
                 pL[i] <- quantile(signal[index] - envelopeweight[envindex2], pquantile[1], na.rm=TRUE)   
                 pU[i] <- quantile(signal[index] + envelopeweight[envindex2], pquantile[2], na.rm=TRUE)                  
            }
            psd[i] <- sd(signal[index], na.rm=TRUE)
        }              
    }
    for (i in startindex:endindex) {
        index <- i + k    
        if (process[1] == "average") {
             pstat[i] <- mean(signal[index], na.rm=TRUE)
        } else if (process[1] == "median") {
             pstat[i] <- median(signal[index], na.rm=TRUE)        
        } else if (process[1] == "envelope") {
             pL[i] <- quantile(signal[index] - envelopeweight, pquantile[1], na.rm=TRUE)    
             pU[i] <- quantile(signal[index] + envelopeweight, pquantile[2], na.rm=TRUE)             
        }
        psd[i] <- sd(signal[index], na.rm=TRUE)
    }   
    
    #### For ensemble filter, the number of a patch touching a given point is "tau".
   
    Epstat <- Epsd <- EpL <- EpU <- EpM <- EpR <- NULL
    
    k <- -rev(k)
    startindex <- 1 + (-k[1]); endindex <- ndata - k[tau0+1]
    
    if (startindex > 1) {
        for (i in 1:(startindex-1)) {
            index <- (i + k)[(i + k) > 0]    
            if (any(process[1] == c("average", "median"))) {
                if (process[2] == "average") {
                     Epstat[i] <- mean(pstat[index])
                } else if (process[2] == "median") {
                     Epstat[i] <- median(pstat[index])               
                }
            } else if (process[1] == "envelope") {
                 if (process[2] == "average") {
                     EpL[i] <- mean(pL[index])
                     EpU[i] <- mean(pU[index])                     
                 } else if (process[2] == "median") {
                     EpL[i] <- median(pL[index])
                     EpU[i] <- median(pU[index])              
                 } else if (process[2] == "envelope") {
                     EpL[i] <- quantile(pL[index], equantile[1])    
                     EpU[i] <- quantile(pU[index], equantile[2])                       
                 }          
            }
            Epsd[i] <- mean(psd[index])
        }
    }   
    if (endindex < ndata) {
        for (i in (endindex+1):ndata) { 
            index <- (i + k)[(i + k) <= ndata]     
            if (any(process[1] == c("average", "median"))) {
                  if (process[2] == "average") {
                      Epstat[i] <- mean(pstat[index])
                  } else if (process[2] == "median") {
                      Epstat[i] <- median(pstat[index])               
                  }
            } else if (process[1] == "envelope") {
                  if (process[2] == "average") {
                      EpL[i] <- mean(pL[index])
                      EpU[i] <- mean(pU[index])                     
                  } else if (process[2] == "median") {
                      EpL[i] <- median(pL[index])
                      EpU[i] <- median(pU[index])              
                  } else if (process[2] == "envelope") {
                      EpL[i] <- quantile(pL[index], equantile[1])    
                      EpU[i] <- quantile(pU[index], equantile[2])                       
                  }          
            }
            Epsd[i] <- mean(psd[index])
        }
    }
    for (i in startindex:endindex) {
        index <- i + k 
        if (any(process[1] == c("average", "median"))) {
             if (process[2] == "average") {
                 Epstat[i] <- mean(pstat[index])
             } else if (process[2] == "median") {
                Epstat[i] <- median(pstat[index])               
          }
        } else if (process[1] == "envelope") {
            if (process[2] == "average") {
                EpL[i] <- mean(pL[index])
                EpU[i] <- mean(pU[index])                     
            } else if (process[2] == "median") {
                EpL[i] <- median(pL[index])
                EpU[i] <- median(pU[index])              
            } else if (process[2] == "envelope") {
                EpL[i] <- quantile(pL[index], equantile[1])    
                EpU[i] <- quantile(pU[index], equantile[2])                       
            }          
        }
        Epsd[i] <- mean(psd[index])
    }

    if (boundary == "symmetric" || boundary == "periodic") {
        if (!is.null(tindex)) tindex <- tindex[obsindex]
        signal <- signal[obsindex]
        pL <- pL[obsindex]
        pU <- pU[obsindex]
        pstat <- pstat[obsindex]
        psd <- psd[obsindex]
        
        EpL <- EpL[obsindex]
        EpU <- EpU[obsindex]
        Epstat <- Epstat[obsindex]
        Epsd <- Epsd[obsindex]
    } 

    if (process[1] == "envelope") {
        pM <- (pL + pU) / 2
        pR <- pU - pL          
        EpM <- (EpL + EpU) / 2
        EpR <- EpU - EpL    
    }
    
    if (any(process[1] == c("average", "median"))) {
        pquantile <- equantile <- gamma <- NULL 
    } else if (process[1] == "envelope") {
        if (process[2] != "envelope") equantile <- NULL
    }
    
    list(tindex=tindex, signal=signal, pstat=pstat, psd=psd, pL=pL, pU=pU, pM=pM, pR=pR, Epstat=Epstat, EpL=EpL, EpU=EpU, EpM=EpM,
         Epsd=Epsd, EpR=EpR, 
         parameters=list(type=type, tau=tau0, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary), nlevel=1)
}

meptransf <- function(tindex=NULL, signal, type="rectangle", taus, process=c("average", "average"), pquantile=c(0, 1), equantile=c(0, 1), gamma=1, boundary="symmetric") {
  
    if (any(floor(taus) < 2)) stop(paste("Some patch size with tau", paste(taus, collapse=","), "is too small."))
    if (!(any(type == c("rectangle", "oval")))) stop("Patch type should be 'rectangle', or 'oval'." )
    if (!any(process[1] == c("average", "median", "envelope"))) stop("The patch transform must be 'average', 'median', or 'envelope'.")
    if (!any(process[2] == c("average", "median", "envelope"))) stop("The ensemnle transform must be 'average', 'median', or 'envelope'.")
    if (any(process[1] == c("average", "median")) & process[2] == "envelope") stop("If the patch transform is 'average' or 'median', the ensemble transform must be 'average' or 'median'.")
    if (pquantile[1] < 0 | pquantile[1] > 1) stop("Quantile for lower envelope of patch transform should be [0, 1].")
    if (pquantile[2] < 0 | pquantile[2] > 1) stop("Quantile for upper envelope of patch transform should be [0, 1].")
    if (pquantile[1] > pquantile[2]) stop("Quantile for lower envelope should be smaller than quantile for upper envelop of patch transform.")
    if (equantile[1] < 0 | equantile[1] > 1) stop("Quantile for lower envelope of ensemble transform should be [0, 1].")
    if (equantile[2] < 0 | equantile[2] > 1) stop("Quantile for upper envelope of ensemble transform should be [0, 1].")
  
    if (gamma < 0) stop("Scale parameter of patch transform should be positive.")
    if (!(any(boundary == c("none", "periodic", "symmetric")))) stop("The boundary condition should be 'none', 'periodic' or 'symmetric'." )
  
    #ndata <- length(signal)
    #if (is.null(tindex)) tindex <- 1:ndata    
    taus <- sort(taus)
  
    pstat <- psd <- pL <- pU <- pM <- pR <- Epstat <- Epsd <- EpL <- EpU <- EpM  <- EpR <- NULL  
  
    taus <- sort(taus)
    n <- length(taus)
  
    for (i in 1:n) {
        tmp <- eptransf(tindex=tindex, signal=signal, type=type, tau=taus[i], process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary)  
        if (any(process[1] == c("average", "median"))) {
            pstat <- cbind(pstat, tmp$pstat) 
            Epstat <- cbind(Epstat, tmp$Epstat)  
        } else if (process[1] == "envelope") {
            pL <- cbind(pL, tmp$pL)
            pU <- cbind(pU, tmp$pU)  
            EpL <- cbind(EpL, tmp$EpL)
            EpU <- cbind(EpU, tmp$EpU)  
        } 
        psd <- cbind(psd, tmp$psd)  
        Epsd <- cbind(Epsd, tmp$Epsd) 
    }
    
    if (process[1] == "envelope") {
        pM <- (pL + pU) / 2
        pR <- pU - pL          
        EpM <- (EpL + EpU) / 2
        EpR <- EpU - EpL    
    }

    if (any(process[1] == c("average", "median"))) {
        rho <- diag(cor(signal-Epstat, Epstat))  
        pquantile <- equantile <- gamma <- NULL 
    } else if (process[1] == "envelope") {
        rho <- diag(cor(signal-EpM, EpM)) 
        if (process[2] != "envelope") equantile <- NULL
    }
  
    list(tindex=tindex, signal=signal, pstat=pstat, psd=psd, pL=pL, pU=pU, pM=pM, pR=pR, Epstat=Epstat, Epsd=Epsd, EpL=EpL, EpU=EpU, EpM=EpM, EpR=EpR, rho=rho, 
        parameters=list(type=type, tau=floor(taus), process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary), nlevel=n)
}

localextrema <- function(y) 
{
    ndata = length(y)
    minindex <- maxindex <- NULL
    nextreme <- 0
    cross <- NULL
    ncross <- 0
    z1 <- sign(diff(y))
    index1 <- seq(1, ndata - 1)[z1 != 0]
    z1 <- z1[z1 != 0]
    if (!(is.null(index1) || all(z1 == 1) || all(z1 == -1))) {
        index1 <- index1[c(z1[-length(z1)] != z1[-1], FALSE)] + 1
        z1 <- z1[c(z1[-length(z1)] != z1[-1], FALSE)]
        nextreme <- length(index1)
        if (nextreme >= 2) 
            for (i in 1:(nextreme - 1)) {
                tmpindex <- index1[i]:(index1[i + 1] - 1)
                if (z1[i] > 0) {
                    tmpindex <- tmpindex[y[index1[i]] == y[tmpindex]]
                    maxindex <- rbind(maxindex, c(min(tmpindex), max(tmpindex)))
                }
                else {
                    tmpindex <- tmpindex[y[index1[i]] == y[tmpindex]]
                    minindex <- rbind(minindex, c(min(tmpindex), max(tmpindex)))
                }
            }
        tmpindex <- index1[nextreme]:(ndata - 1)
        if (z1[nextreme] > 0) {
            tmpindex <- tmpindex[y[index1[nextreme]] == y[tmpindex]]
            maxindex <- rbind(maxindex, c(min(tmpindex), max(tmpindex)))
        }
        else {
            tmpindex <- tmpindex[y[index1[nextreme]] == y[tmpindex]]
            minindex <- rbind(minindex, c(min(tmpindex), max(tmpindex)))
        }
        if (!(all(sign(y) >= 0) || all(sign(y) <= 0) || all(sign(y) == 0))) {
            index1 <- c(1, index1)
            for (i in 1:nextreme) {
                if (y[index1[i]] == 0) {
                    tmp <- c(index1[i]:index1[i + 1])[y[index1[i]:index1[i + 1]] == 0]
                    cross <- rbind(cross, c(min(tmp), max(tmp)))
                }
                else if (y[index1[i]] * y[index1[i + 1]] < 0) {
                    tmp <- min(c(index1[i]:index1[i + 1])[y[index1[i]] * y[index1[i]:index1[i + 1]] <= 0])
                    if (y[tmp] == 0) {
                        tmp <- c(tmp:index1[i + 1])[y[tmp:index1[i + 1]] == 0]
                        cross <- rbind(cross, c(min(tmp), max(tmp)))
                    }
                    else cross <- rbind(cross, c(tmp - 1, tmp))
                }
            }
            if (any(y[index1[nextreme + 1]] * y[index1[nextreme + 1]:ndata] <= 0)) {
                tmp <- min(c(index1[nextreme + 1]:ndata)[y[index1[nextreme + 1]] * y[index1[nextreme + 1]:ndata] <= 0])
                if (y[tmp] == 0) {
                    tmp <- c(tmp:ndata)[y[tmp:ndata] == 0]
                    cross <- rbind(cross, c(min(tmp), max(tmp)))
                }
                else cross <- rbind(cross, c(tmp - 1, tmp))
            }
            ncross <- nrow(cross)
        }
    }
    list(minindex = minindex, maxindex = maxindex, nextreme = nextreme, cross = cross, ncross = ncross)
}  

