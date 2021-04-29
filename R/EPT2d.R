#### The second version
  
eptdecomp2d <- function(x=NULL, y=NULL, z, type="rectangle", tau, theta=0, process=c("average", "average"), pquantile=c(0, 1), equantile=c(0, 1), gamma=1, boundary="reflexive",   
                          stoprule = "type2", tol=0.1^2, maxiter=10, check=FALSE)
{   
    if (!(any(type == c("rectangle", "oval")))) stop("Patch type should be 'rectangle', or 'oval'." )
    if (any(floor(tau) < 2)) stop(paste("The patch size with tau", tau, "is too small."))
    if (theta < 0 || theta > 180) stop("Rotation of Patch should be clockwise by degree of (0, 180)." )   
    if (!any(process[1] == c("average", "median", "envelope"))) stop("The patch transform must be 'average', 'median', or 'envelope'.")
    if (!any(process[2] == c("average", "median", "envelope"))) stop("The ensemnle transform must be 'average', 'median', or 'envelope'.")
    if (any(process[1] == c("average", "median")) & process[2] == "envelope") stop("If the patch transform is 'average' or 'median', the ensemble transform must be 'average' or 'median'.")  
    if (pquantile[1] < 0 | pquantile[1] > 1) stop("Quantile for lower envelope of patch transform should be [0, 1].")
    if (pquantile[2] < 0 | pquantile[2] > 1) stop("Quantile for upper envelope of patch transform should be [0, 1].")
    if (pquantile[1] > pquantile[2]) stop("Quantile for lower envelope should be smaller than quantile for upper envelop of patch transform.")  
    if (equantile[1] < 0 | equantile[1] > 1) stop("Quantile for lower envelope of ensemble transform should be [0, 1].")
    if (equantile[2] < 0 | equantile[2] > 1) stop("Quantile for upper envelope of ensemble transform should be [0, 1].")  

    if (gamma < 0) stop("Scale parameter of patch transform should be positive.")
    if (!(any(boundary == c("none", "periodic", "reflexive")))) stop("The boundary condition should be 'none', 'periodic' or 'reflexive'." )
  
    if (maxiter < 0) stop("The maximum iteration of sifting should be positive.")

    if (any(process[1] == c("average", "median"))) {
        eptransform <- "Epstat"
    } else if (process[1] == "envelope") {
        eptransform <- "EpM"
    }
    
    #lower <- upper <- eptcomp <- NULL
    eptcomp <- list(NULL); residue <- matrix(0, dim(z)[1], dim(z)[2])
    
    iter <- 1
    input <- z
    repeat { 
        tmp <- eptransf2d(x=x, y=x, z=input, type=type, tau=tau, theta=theta, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary)
        
        eptcomp[[iter]] <- tmp[[eptransform]] 
        residue <- residue + tmp[[eptransform]] 
        
        if (check) {
            if (iter == 1) rangez <- range(z)
            image(input, xlab="", ylab="", col=gray(0:100/100), axes=F, zlim=rangez, main=paste0("Frequency component at iter=", iter))
            image(tmp[[eptransform]], xlab="", ylab="", col=gray(0:100/100), axes=F, zlim=rangez, main=paste0("Additional frequency component at iter=", iter))
            image(residue, xlab="", ylab="", col=gray(0:100/100), axes=F, zlim=rangez, main=paste0("Residue at iter=", iter))
            locator(1)          
        }   
        
        if (stoprule == "type1" && (all(abs(eptcomp[[iter]]) < tol) || iter >= maxiter)) {
            FC <- z - residue  
            break
        } else if (stoprule == "type2" && iter >= 2) {
            if (sum((eptcomp[[iter-1]] - eptcomp[[iter]])^2) < tol * sum(eptcomp[[iter-1]]^2) || iter >= maxiter) {
                FC <- z - residue  
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
    
    if (check) 
        list(eptcomp=eptcomp, FC=FC, residue=residue, 
             parameters=list(type=type, tau=floor(tau), theta=theta, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary, niter=iter))
    else list(FC=FC, residue=residue, parameters=list(type=type, tau=floor(tau), theta=theta, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary, niter=iter))
    
}

meptransf2d <- function(x=NULL, y=NULL, z, type="rectangle", taus, theta=0, process=c("average", "average"), pquantile=c(0, 1), equantile=c(0, 1), gamma=1, boundary="reflexive") {
    
    if (!(any(type == c("rectangle", "oval")))) stop("Patch type should be 'rectangle', or 'oval'." )
    if (any(floor(taus) < 2)) stop(paste("Some patch size with tau", paste(taus, collapse=","), "is too small."))
    if (theta < 0 || theta > 180) stop("Rotation of Patch should be clockwise by degree of (0, 180)." )   
    if (!any(process[1] == c("average", "median", "envelope"))) stop("The patch transform must be 'average', 'median', or 'envelope'.")
    if (!any(process[2] == c("average", "median", "envelope"))) stop("The ensemnle transform must be 'average', 'median', or 'envelope'.")
    if (any(process[1] == c("average", "median")) & process[2] == "envelope") stop("If the patch transform is 'average' or 'median', the ensemble transform must be 'average' or 'median'.")  
    if (pquantile[1] < 0 | pquantile[1] > 1) stop("Quantile for lower envelope of patch transform should be [0, 1].")
    if (pquantile[2] < 0 | pquantile[2] > 1) stop("Quantile for upper envelope of patch transform should be [0, 1].")
    if (pquantile[1] > pquantile[2]) stop("Quantile for lower envelope should be smaller than quantile for upper envelop of patch transform.")
    if (equantile[1] < 0 | equantile[1] > 1) stop("Quantile for lower envelope of ensemble transform should be [0, 1].")
    if (equantile[2] < 0 | equantile[2] > 1) stop("Quantile for upper envelope of ensemble transform should be [0, 1].")  
  
    if (gamma < 0) stop("Scale parameter of patch transform should be positive.")
    if (!(any(boundary == c("none", "periodic", "reflexive")))) stop("The boundary condition should be 'none', 'periodic' or 'reflexive'." )
    
    if (is.vector(taus)) {
        taus <- sort(taus)
        taus <- cbind(taus, taus)
        n <-  nrow(taus)         
    } else if (is.matrix(taus)) {
        taus <- taus[order(taus[,1]),]
        n <- nrow(taus)        
    }
 
    if (process[1] != "envelope") {
        pstat <- psd <- Epstat <- Epsd <- list(NULL); pL <- pU <- pM <- pR <- EpL <- EpU <- EpM <- EpR <- NULL
    } else {
        psd <- Epsd <- pL <- pU <- pM <- pR <- EpL <- EpU <- EpM <- EpR <- list(NULL); pstat <- Epstat <- NULL            
    }
  
    rho <- NULL

    for (i in 1:n) {
        tmp <- eptransf2d(x=x, y=y, z=z, type=type, tau=taus[i,], theta=theta, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary)  
        
        if (any(process[1] == c("average", "median"))) {
            pstat[[i]] <- tmp$pstat
            Epstat[[i]] <- tmp$Epstat
            rho <- c(rho, cor(c(z - Epstat[[i]]), c(Epstat[[i]])))
        } else if (process[1] == "envelope") {
            pL[[i]] <- tmp$pL
            pU[[i]] <- tmp$pU
            pM[[i]] <- tmp$pM
            pR[[i]] <- tmp$pR              
            EpL[[i]] <- tmp$EpL
            EpU[[i]] <- tmp$EpU
            EpM[[i]] <- tmp$EpM
            EpR[[i]] <- tmp$EpR
            rho <- c(rho, cor(c(z - EpM[[i]]), c(EpM[[i]])))
        } 
        psd[[i]] <- tmp$psd
        Epsd[[i]] <- tmp$Epsd
    }

   if (any(process[1] == c("average", "median"))) {
        pquantile <- equantile <- gamma <- NULL 
    } else if (process[1] == "envelope") {
        if (process[2] != "envelope") equantile <- NULL
    }
    
    list(x=x, y=y, z=z, pstat=pstat, psd=psd, pL=pL, pU=pU, pM=pM, pR=pR, Epstat=Epstat, Epsd=Epsd, EpL=EpL, EpU=EpU, EpM=EpM, EpR=EpR, rho=rho,
         parameters=list(type=type, taus=floor(taus), theta=theta, process=process, pquantile=pquantile, equantile=equantile, gamma=gamma, boundary=boundary, corr=rho), nlevel=n)
}

eptransf2d <- function(x=NULL, y=NULL, z, type="rectangle", tau, theta=0, process=c("average", "average"), pquantile=c(0, 1), equantile=c(0, 1), gamma=1, boundary="reflexive") {
    
    if (!(any(type == c("rectangle", "oval")))) stop("Patch type should be 'rectangle', or 'oval'." )
    if (theta < 0 || theta > 180) stop("Rotation of Patch should be clockwise by degree of (0, 180)." )   
    if (!any(process[1] == c("average", "median", "envelope"))) stop("The patch transform must be 'average', 'median', or 'envelope'.")
    if (!any(process[2] == c("average", "median", "envelope"))) stop("The ensemnle transform must be 'average', 'median', or 'envelope'.")
    if (any(process[1] == c("average", "median")) & process[2] == "envelope") stop("If the patch transform is 'average' or 'median', the ensemble transform must be 'average' or 'median'.")
    if (pquantile[1] < 0 | pquantile[1] > 1) stop("Quantile for lower envelope of patch transform should be [0, 1].")
    if (pquantile[2] < 0 | pquantile[2] > 1) stop("Quantile for upper envelope of patch transform should be [0, 1].")
    if (pquantile[1] > pquantile[2]) stop("Quantile for lower envelope should be smaller than quantile for upper envelop of patch transform.")
    if (equantile[1] < 0 | equantile[1] > 1) stop("Quantile for lower envelope of ensemble transform should be [0, 1].")
    if (equantile[2] < 0 | equantile[2] > 1) stop("Quantile for upper envelope of ensemble transform should be [0, 1].")
    if (gamma < 0) stop("Scale parameter of patch transform should be positive.")
    if (!(any(boundary == c("none", "periodic", "reflexive")))) stop("The boundary condition should be 'none', 'periodic' or 'reflexive'." )

    nr <- nrow(z) 
    nc <- ncol(z)
    
    boundperc <- 0.3
    
    if (boundary == "periodic" || boundary == "reflexive") {
        rowext <- round(nr * boundperc)
        colext <- round(nc * boundperc)
        nrext <- 2 * rowext + nr
        ncext <- 2 * colext + nc
        
        if (boundary == "periodic") {
            if (!is.null(x)) 
                x <- c((rev(x[1] - cumsum(rev(diff(x)))))[(nr-rowext):(nr-1)], x, (x[nr] + cumsum(diff(x)))[1:rowext])   
            if (!is.null(y))
                y <- c((rev(y[1] - cumsum(rev(diff(y)))))[(nc-colext):(nc-1)], y, (y[nc] + cumsum(diff(y)))[1:colext])            
            z <- cbind(z[, (nc - colext + 1):nc], z, z[, 1:colext])
            z <- rbind(z[(nr - rowext + 1):nr, ], z, z[1:rowext, ])
        } else if (boundary == "reflexive") {
            if (!is.null(x)) 
                x <- c((x[1] - rev(cumsum(diff(x))))[(nr-rowext):(nr-1)], x, (x[nr] + cumsum(rev(diff(x))))[1:rowext])
            if (!is.null(y)) 
                y <- c((y[1] - rev(cumsum(diff(y))))[(nc-colext):(nc-1)], y, (y[nc] + cumsum(rev(diff(y))))[1:colext])
            
            z <- cbind(z[, colext:1], z, z[, nc:(nc - colext + 1)])
            z <- rbind(z[rowext:1, ], z, z[nr:(nr - rowext + 1), ])
        }
    } else if (boundary == "none") {
        rowext <- 0; nrext <- nr
        colext <- 0; ncext <- nc
    } 
    
    if (process[1] != "envelope") {
        pstat <- psd <- matrix(0, nrext, ncext); pL <- pU <- pM <- pR <- EpL <- EpU <- EpM <- EpR <- NULL
    } else {
        psd <- pL <- pU <- matrix(0, nrext, ncext); pstat <- NULL            
    }

    
    if (length(tau) == 1) tau <- c(tau, tau)
    tau0 <- floor(tau)
    if (any(tau0 < 2)) stop(paste("The patch size with tau", tau, "is too small."))

    index1 <- 0:tau0[1] - ceiling(tau0[1]/2); index2 <- 0:tau0[2] - ceiling(tau0[2]/2)     
  
    if (theta == 0) {
        index12 <- as.matrix(expand.grid(index1, index2))
        n <- range(index1); m <- range(index2)
        if (type == "rectangle") {
            weight2d <- rep(gamma * sqrt(sum(tau0^2)), nrow(index12)) 
       } else if (type == "oval") {
         #weight2d <- outer(index1, index2, function(x, y) gamma * sqrt(sum(tau^2)/8 - 0.5 * (x^2 + y^2)))   
           weight2d <- apply(t(t(index12) +  tau0 %% 2/2), 1, function(tmp) gamma * sqrt(sum(tau0^2)/4 - (sum(tmp^2))))
       }
    } else {
        oripatch <- matrix(c(index1[tau0[1]+1], index2[tau0[2]+1], index1[tau0[1]+1], index2[1], 
                             index1[1], index2[1], index1[1], index2[tau0[2]+1], index1[tau0[1]+1], index2[tau0[2]+1]), nrow=2)
        theta <- pi/180 * theta # degree to radian
        rotmatrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2) 
        rotpatch <- rotmatrix %*% oripatch
        n <- ceiling(range(rotpatch[1,])+c(-1,0))
        m <- ceiling(range(rotpatch[2,])+c(-1,0))
        candidate <- expand.grid(n[1]:n[2], m[1]:m[2])
      
        rotpatchindex <- NULL
        for (i in 1:nrow(candidate)) {
            result <- (rotpatch[2, -1] - rotpatch[2, -5]) * (candidate[i, 1]-rotpatch[1,-5]) - (rotpatch[1, -1] - rotpatch[1, -5]) * (candidate[i, 2] - rotpatch[2,-5])
            if (all(result >= -1.0e-10)) 
            rotpatchindex <- rbind(rotpatchindex, candidate[i, ])
        }
        index12 <- as.matrix(rotpatchindex)
       
        if (type == "rectangle") {
            weight2d <- rep(gamma * sqrt(sum(tau0^2)), nrow(index12)) 
        } else if (type == "oval") {
            oriindex12 <- t(t(rotmatrix) %*% t(index12) +  (tau0 %% 2)/2)
            weight2d <- apply(oriindex12, 1, function(tmp) sqrt(sum(tau0^2)/4 - (sum(tmp^2))))
        }
    }

    n <- abs(n); m <- abs(m)
    for (i in (n[1]+1):(nrext-n[2])) {
        for (j in (m[1]+1):(ncext-m[2])) {
            findex <- list(index12=cbind(i + index12[, 1], j + index12[, 2]))
            if (process[1] == "average") {
                pstat[i,j] <- mean(z[findex$index12], na.rm=TRUE)  
            } else if (process[1] == "median") {
                pstat[i,j] <- median(z[findex$index12], na.rm=TRUE)          
            } else if (process[1] == "envelope") {
                pL[i,j] <- quantile(z[findex$index12] + weight2d, pquantile[1], na.rm=TRUE)   
                pU[i,j] <- quantile(z[findex$index12] + weight2d, pquantile[2], na.rm=TRUE)
            }
            psd[i,j] <- sd(z[findex$index12], na.rm=TRUE)
        }
    }
    
    for (i in c(1:n[1], (nrext-n[2]+1):nrext)) {
        for (j in 1:ncext) {
            findex <- chooseindex2d(i, j, index12, c(nrext, ncext)) 
            if (process[1] == "average") {
                pstat[i,j] <- mean(z[findex$index12], na.rm=TRUE)  
            } else if (process[1] == "median") {
                pstat[i,j] <- median(z[findex$index12], na.rm=TRUE)          
            } else if (process[1] == "envelope") {
                pL[i,j] <- quantile(z[findex$index12] + weight2d[findex$windex12], pquantile[1], na.rm=TRUE)   
                pU[i,j] <- quantile(z[findex$index12] + weight2d[findex$windex12], pquantile[2], na.rm=TRUE)
            }
            psd[i,j] <- sd(z[findex$index12], na.rm=TRUE)
        }
    }
    
    for (i in (n[1]+1):(nrext-n[2])) {
        for (j in c(1:m[1], (ncext-m[2]+1):ncext)) {
            findex <- chooseindex2d(i, j, index12, c(nrext, ncext)) 
            if (process[1] == "average") {
                pstat[i,j] <- mean(z[findex$index12], na.rm=TRUE)  
            } else if (process[1] == "median") {
                pstat[i,j] <- median(z[findex$index12], na.rm=TRUE)          
            } else if (process[1] == "envelope") {
              pL[i,j] <- quantile(z[findex$index12] + weight2d[findex$windex12], pquantile[1], na.rm=TRUE)   
              pU[i,j] <- quantile(z[findex$index12] + weight2d[findex$windex12], pquantile[2], na.rm=TRUE)
            }
            psd[i,j] <- sd(z[findex$index12], na.rm=TRUE)
        }
    }
 
    if (process[1] != "envelope") {
        Epstat <- Epsd <- matrix(0, nrext, ncext); EpL <- EpU <- NULL
    } else {
        Epsd <- EpL <- EpU <- matrix(0, nrext, ncext); Epstat <- NULL            
    }
            
    for (i in (n[1]+1):(nrext-n[2])) {
        for (j in (m[1]+1):(ncext-m[2])) {
            findex <- list(index12=cbind(i + index12[, 1], j + index12[, 2]))
            if (any(process[1] == c("average", "median"))) {
                if (process[2] == "average") {
                    Epstat[i, j] <- mean(pstat[findex$index12])
                } else if (process[2] == "median") {
                    Epstat[i, j] <- median(pstat[findex$index12])               
                }
            } else if (process[1] == "envelope") {
                if (process[2] == "average") {
                    EpL[i, j] <- mean(pL[findex$index12])
                    EpU[i, j] <- mean(pU[findex$index12])                     
                } else if (process[2] == "median") {
                    EpL[i, j] <- median(pL[findex$index12])
                    EpU[i, j] <- median(pU[findex$index12])              
                } else if (process[2] == "envelope") {
                    EpL[i, j] <- quantile(pL[findex$index12], equantile[1])    
                    EpU[i, j] <- quantile(pU[findex$index12], equantile[2])                       
                }          
            }
            Epsd[i] <- mean(psd[findex$index12])        
        }
    }
    
    for (i in c(1:n[1], (nrext-n[2]+1):nrext)) {
        for (j in 1:ncext) {
            findex <- chooseindex2d(i, j, index12, c(nrext, ncext)) 
            if (any(process[1] == c("average", "median"))) {
                if (process[2] == "average") {
                    Epstat[i, j] <- mean(pstat[findex$index12])
                } else if (process[2] == "median") {
                    Epstat[i, j] <- median(pstat[findex$index12])               
                }
            } else if (process[1] == "envelope") {
                if (process[2] == "average") {
                     EpL[i, j] <- mean(pL[findex$index12])
                     EpU[i, j] <- mean(pU[findex$index12])                     
                 } else if (process[2] == "median") {
                     EpL[i, j] <- median(pL[findex$index12])
                     EpU[i, j] <- median(pU[findex$index12])              
                 } else if (process[2] == "envelope") {
                     EpL[i, j] <- quantile(pL[findex$index12], equantile[1])    
                     EpU[i, j] <- quantile(pU[findex$index12], equantile[2])                       
                 }          
            }
            Epsd[i] <- mean(psd[findex$index12])  
        }
    }
    
    for (i in (n[1]+1):(nrext-n[2])) {
        for (j in c(1:m[1], (ncext-m[2]+1):ncext)) {
        findex <- chooseindex2d(i, j, index12, c(nrext, ncext)) 
            if (any(process[1] == c("average", "median"))) {
                if (process[2] == "average") {
                     Epstat[i, j] <- mean(pstat[findex$index12])
                } else if (process[2] == "median") {
                     Epstat[i, j] <- median(pstat[findex$index12])               
                }
            } else if (process[1] == "envelope") {
                 if (process[2] == "average") {
                     EpL[i, j] <- mean(pL[findex$index12])
                     EpU[i, j] <- mean(pU[findex$index12])                     
                 } else if (process[2] == "median") {
                     EpL[i, j] <- median(pL[findex$index12])
                     EpU[i, j] <- median(pU[findex$index12])              
                 } else if (process[2] == "envelope") {
                     EpL[i, j] <- quantile(pL[findex$index12], equantile[1])    
                     EpU[i, j] <- quantile(pU[findex$index12], equantile[2])                       
                 }          
            }
            Epsd[i] <- mean(psd[findex$index12])  
        }
    }
    
    obsindex1 <- (rowext+1):(rowext+nr); obsindex2 <- (colext+1):(colext+nc)
    if (!is.null(x)) x <- x[obsindex1]
    if (!is.null(y)) y <- y[obsindex2]
    z <- z[obsindex1, obsindex2]

    psd <- psd[obsindex1, obsindex2]
    Epsd <- Epsd[obsindex1, obsindex2]

    if (process[1] != "envelope") {
        pstat <- pstat[obsindex1, obsindex2]
        Epstat <- Epstat[obsindex1, obsindex2]
    } else {
        pL <- pL[obsindex1, obsindex2]
        pU <- pU[obsindex1, obsindex2]
        EpL <- EpL[obsindex1, obsindex2]
        EpU <- EpU[obsindex1, obsindex2]     
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
    
    list(x=x, y=y, z=z, pstat=pstat, psd=psd, pL=pL, pU=pU, pM=pM, pR=pR, Epstat=Epstat, EpL=EpL, EpU=EpU, EpM=EpM,
         Epsd=Epsd, EpR=EpR, 
         parameters=list(type=type, tau=tau0, theta=theta*180/pi, process=process, pquantile=pquantile, equantile=equantile, 
                         gamma=gamma, boundary=boundary), nlevel=1)
}

chooseindex2d <- function(i, j, index12, ndim) {
  index12 <-  cbind(i + index12[, 1], j + index12[, 2])
  windex12 <- (index12[, 1] >= 1 & index12[, 1] <= ndim[1] & index12[, 2] >= 1 & index12[, 2] <= ndim[2])
  index12 <- index12[windex12, ]
  
  list(index12=index12, windex12=windex12)
}
