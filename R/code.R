##### DISTANCE-BASED GENERAL (REGIONALISED) SENSITIVITY ANALYSIS #####

# ADAPTED FROM THE MATLAB CODE DEVELOPPED BY CELINE SCHEIDT
# (STANFORD UNIVERSITY)

# IDEA:
#		- use quantile fx from RCPP

require(RColorBrewer)


# Kolmogorov-Smirnov test (KS)
# Earth Mover’s Distance (EMD), also known as the first Wasserstein distance
# The Cramér-von Mises Distance
distCDF <- function(x, y, method = c("EMD", "KS", "CMD")){
  method <- match.arg(method, c("EMD", "KS", "CMD"))
  x_cdf <- ecdf(x)
  y_cdf <- ecdf(y)

  xy_knots <- sort(c(knots(x_cdf), knots(y_cdf)))
  y_val <- y_cdf(xy_knots)
  x_val <- x_cdf(xy_knots)
  n <- length(x_val)
  if(method == "EMD"){
    dCDF <- sum(abs(x_val[-n] - y_val[-n]) * diff(xy_knots))
  }else if(method == "KS"){
    dCDF <- max(abs(x_val - y_val))
  }else if(method == "CMD"){
    dCDF <- sqrt(um(abs(x_val[-n] - y_val[-n])^2 * diff(xy_knots)))
  }
  return(dCDF)
}

resampleTest <- function(x,y, nx, n = 1, ...){
  xPerm <- sample(x, size = nx, replace = FALSE)
  qxPerm <- quantile(xPerm, ...)
  return(sum(abs(qxPerm - y)^n)^(1/n))
}

# #  measure of distance between unconditional and conditional CDFs
# # + if nPerm > 0, compute p-values
# distCDF <- function(x, cl = NULL, nPerm = 1000, n = 100, type = 1){
#   probs <- seq(0, 1, length.out = n)
#   qP <- quantile(x, probs = probs, type = type)
#   cla <- unique(cl) # classes
#   ncl <- length(cla)
#   dCDF <- numeric(ncl) # initialisation distance between CDFs
#   names(dCDF) <- cla
#   pValue <- dCDF       # initialisation p-values
#   for(i in seq_len(ncl)){
#     inCli <- cl == cla[i]
#     nci <- sum(inCli)
#     qPc <- quantile(x[inCli], probs = probs, type = type)
#     dCDF[i] <- sum( abs(qPc - qP) )
#     dCDFPerm <- replicate(nPerm, resampleTest(x, qP, nx = nci,
#                                               probs = probs, type = type))
#     # dCDFPerm <- numeric(nPerm)
#     # for(j in seq_len(nPerm)){
#     #   xPerm <- sample(x, size = nci, replace = FALSE)
#     #   qPcPerm <- quantile(xPerm, probs = probs, type = type)
#     #   dCDFPerm[j] <- sum(abs(qPcPerm - qP))
#     # }
#     # pValue[i] <- (sum(dCDFPerm > dCDF[i])  + 1)/(nPerm + 1)
#     pValue[i] <- (sum(dCDFPerm > dCDF[i])  )/(nPerm )
#   }
#      return(list(dist = dCDF, pvalue = pValue))
# }


#  measure of distance between "CDF conditioned to cluster and
# parameter 'i' (bin)"
# and  CDF conditioned to cluster
# @ param pCond (optional) parameter to which to compute the conditional
#               interaction
distCDFint <- function(P, i, cl = NULL, nBins = 3, pCond = NULL,
                       n = 100, type = 1, nPerm = 100){
  x <- P[,i]  # the contioning parameter p
  probs <- seq(0, 1, length.out = n)
  cla <- unique(cl) # classes
  ncl <- length(cla)
  nP <- ncol(P)
  if(length(nBins) == 1){
    nb <- nBins[1]
  }else{
    nb <- nBins[i]
  }
  # definitions of the level bins for the conditional parameter
  if(length(uP <- unique(x)) <= nb){
    lvls <- sort(uP)
  }else{
    lvls <- quantile(x, probs = (1:(nb-1))/nb)
  }

  if(is.null(pCond)){
    intSens <- array(0,dim=c(nP, ncl, nb), dimnames = list(colnames(P),
                                                          paste0("cl",1:ncl),
                                                          paste0("bin",1:nb)))
  }else{
    intSens <- array(0,dim=c(ncl, nb))
  }

  pValue <- intSens

  # over the cluster
  for(j in seq_len(ncl)){
    idx_cl <- cl == cla[j]

    # over the bins
    for(k in seq_len(nb)){
      if(k ==1){
        idx_cl_bin <- x[idx_cl] <= lvls[1]
      }else if( k == nb){
        idx_cl_bin <- x[idx_cl] > lvls[nb-1]
      }else{
        idx_cl_bin <- x[idx_cl] > lvls[k-1] & x[idx_cl] < lvls[k]
      }
      if( (nci <- sum(idx_cl_bin)) <= 3){
        intSensPerm[,j,k,] <- NA
      }else{
        nci2 <- sum(idx_cl_bin)
        if(nci != nci2) stop("pborlme\n")
        # compute the L1-pseudo-norm
        # over the other parameters
        if(is.null(pCond)){
          for(l in seq_len(nP)){
            if(l != i){
            # distribution of parameter Pl conditioned to the cluster: F(p_i|c)
            qPi_cl <- quantile(P[idx_cl,l],probs=probs,type=type)
            # distribution of parameter Pl conditioned to Pi in bin l: F(p|i
            qPi_cl_Bin <- quantile(P[idx_cl_bin,l],probs=probs,type=type)
            intSens[l,j,k] <- sum(abs(qPi_cl_Bin - qPi_cl))
            dCDFPerm <- numeric(nPerm)
            for(u in seq_len(nPerm)){
#               Pi_cl_bin_perm <- sample(P[idx_cl,l],size=nci,replace=FALSE)
#               qPi_cl_bin_perm <- quantile(Pi_cl_bin_perm,probs=probs,type=type
# )
#               l1NormPerm[u] <- sum(abs(qPi_cl_bin_perm - qPi_cl))
              xPerm <- sample(P[idx_cl,l], size=nci, replace=FALSE)
              qPcPerm <- quantile(xPerm, probs=probs, type=type)
              dCDFPerm[u] <- sum(abs(qPcPerm - qPi_cl))
            }
        #       pValue[i] <- (sum(dCDFPerm > intSens[l,j,k])  + 1)/(nPerm + 1)
            pValue[l,j,k]  <- (sum(dCDFPerm > intSens[l,j,k])  )/(nPerm )
              #       cat(sum(abs(qPi_cl_Bin - qPi_cl)), " - ")
            }
          }
        }else{
          # distribution of parameter Pl conditioned to the cluster: F(p_i|c)
          qPi_cl <- quantile(P[idx_cl,l],probs=probs,type=type)
          # distribution of parameter Pl conditioned to Pi in bin l: F(p|i
          qPi_cl_Bin <- quantile(P[idx_cl_bin,l],probs=probs,type=type)
          intSens[j,k] <- sum(abs(qPi_cl_Bin - qPi_cl))
        }
      }
    }
  }
#   return(intSens[-i,,])
  return(list(dist = intSens, pvalue=pValue))
}




# colPal = c(rgb(0.95,0.95,0.95),
#                         brewer.pal(3,"Reds"))(4))


dGSABarplot <- function(x, colPal = brewer.pal(4,"Reds"), main="",
                        xlab="p-value"){
  xpval <- lapply(x, function(x) x$pvalue)
  xmax <- 1 - sapply(xpval, min)

  # sensitivity ranking
  rk <- order(xmax, decreasing = FALSE)
  rkx <- xmax[rk]
  XLIM <- c(0, 1)

  colVar <- rkx
  colVar[] <- 1
  colVar[rkx >= 0.90] <- 2
  colVar[rkx >= 0.95] <- 3
  colVar[rkx >= 1] <- 4

  op <- par(no.readonly = TRUE)
  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,1), heights=c(1))
  # layout.show(4)par(mai=c(1.1,1.1,0.2,0.2))
  par(mai=c(1.2,1.2,0.5,0.2))
  plot(0,0,type="n",xlim=XLIM, ylim=c(0,length(rkx)-0.6),ylab="",
      bty="n",yaxt='n', xaxt = "n",xlab=xlab)
  for(i in seq_along(rkx)){
    z <- i - 1
    rect(xleft=0,ybottom=z-0.4,xright=rkx[i],ytop=z+0.4,
         col=colPal[colVar[i]], lwd=1)
  }
  segments(x0=1,y0=-0.5,x1=1,y1=length(rkx)-0.6,col="dodgerblue",lwd=2)
  segments(x0=0.95,y0=-0.5,x1=0.95,y1=length(rkx)-0.6,col="dodgerblue",
        lwd=2, lty=3)
  title(main)
  axis(1, at= c(0,0.25,0.5,0.75,0.95,1), labels =TRUE)
  par(xpd=TRUE)
    text(x= XLIM[1], y=seq_along(rkx)-1,names(rkx),pos=2)
  par(xpd=NA)
  par(mai=c(1.2,0.2,0.5,0.7), xaxs="i")
  w <- 2
  val <- c(0, 0.9, 0.95, 1, 1.01)
  plot(c(0,w), range(val), type='n', bty='n', xaxt='n', xlab="",
      yaxt='n', ylab='')#, main=title)
  for (i in 1:(length(val)-1)) {
    ybot = val[i]
    ytop = val[i+1]
    rect(0,ybot,w,ytop, col=colPal[i], border=NA)
  }
  rect(0,val[1],w,tail(val,1), border="black")
  axis(4, at = val[-4], labels = val[-4], las=2)
  par(xpd = FALSE)
  par(op)
}

GSABarplotClust <- function(x, colPal =  brewer.pal(4,"Reds"), main="",
                        xlab="p-value", FUN=min, spc=1.15){
  op <- par(no.readonly = TRUE)
  ### H-errorbar
  xdist <- t(sapply(x, function(x) x$dist))
  xpval <- t(sapply(x, function(x) x$pvalue))
  x_max <- 1-apply(xpval,1,min)
  # sensitivity ranking
  rk <- order(x_max,decreasing = FALSE)
  rkx <- 1-xpval[rk,]

  ncl <- ncol(xpval)
  nP <- nrow(xpval)

  XLIM <- c(0, 1)
  YLIM <- c(0, spc*(nP- 1)*ncl+ncl)

#   colvar <- (1 - rkx) * (rkx >= deltam)*(rkx < 1)
#   maxcloser <- ((1 - rkx)[rkx < deltam & rkx < 1])
#   maxcloser <- min(maxcloser[maxcloser > max(colvar)])
#   colvar <- 100 * (maxcloser - colvar)/maxcloser*(rkx < 1)*(rkx >= deltam)
  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,1), heights=c(1))
  par(mai=c(1.2,1.2,0.5,0.2))
  plot(0,0,type="n",xlim=XLIM, ylim=YLIM, ylab="", xlab=xlab,
      bty="n",yaxt='n')
  for(i in seq_len(nP)){
    for(j in seq_len(ncl)){
    z <- spc * (i - 1)*ncl + (1+ncl - j) -1
#     if(rkx[j,i] < deltam[j,i]){
#       myCol <- head(colPal,1)
#   #     }else if(rkx[i] > deltap[i]){
#     }else if(rkx[j,i] > 1){
#       myCol <- tail(colPal,1)
#     }else{
#       myCol <- colPal[colvar[j,i]]
#     }
    rect(xleft=0,ybottom=z-0.45,xright=rkx[i,j],ytop=z+0.45, col="red", lwd=1)
#     arrows(deltam[j,i], z, deltap[j,i], z, length=0.05, angle=90, code=3)
    }
  }
  if(!is.null(colnames(rkx))){
    clName <- colnames(rkx)
    par(xpd=TRUE)
    for(j in seq_len(ncl)){
    zlab <- spc * (0 - 1)*ncl + 2.5*(ncl-j+1) -1
    text(XLIM[2],5+zlab,clName[j],pos=2)
    }
    text(XLIM[2],5+spc * (0 - 1)*ncl + 2.5*(ncl+1) -1,"clusters:",pos=2)
  par(xpd=NA)

  }
   segments(x0=1,y0=-0.5,x1=1,y1=YLIM[2] ,col="dodgerblue",lwd=2)
  segments(x0=0.95,y0=-0.5,x1=0.95,y1=YLIM[2],col="dodgerblue",
          lwd=2, lty=3)
  title(main)
  par(xpd=TRUE)
    text(x= XLIM[1], y=spc * (1:nP - 1)*ncl + ncl/2-0.5,rownames(rkx),pos=2)
  par(xpd=NA)
  # image.scale(z=rkx, col = heat.colors(12))

  par(mai=c(1.2,0.2,0.5,0.7), xaxs="i")

  w <- 2
  val <- c(0, 0.9, 0.95,0.99,1)
#   pos <- c(
#   scale = (length(lut)-1)/(max-min)
  plot(c(0,w), range(val), type='n', bty='n', xaxt='n', xlab="",
      yaxt='n', ylab='')#, main=title)
#   axis(4, round(ticks,n), las=1)
  for (i in 1:(length(val)-1)) {
    ybot = val[i]
    ytop = val[i+1]
    rect(0,ybot,w,ytop, col=colPal[i], border=NA)
  }
  rect(0,val[1],w,tail(val,1), border="black")
  axis(4, at = val[-4], labels = TRUE, las=2)

  par(op)
}

dGSABarplotClustOld <- function(x, colPal = c(rgb(0.95,0.95,0.95),
                  colorRampPalette(brewer.pal(3,"Reds"))(101)),
                  main="", FUN=max, spc=1.15){
	### H-errorbar
	x_max <- apply(x,2,FUN)
	# sensitivity ranking
	rk <- order(x_max,decreasing = FALSE)
	rkx <- x[,rk]

	ncl <- nrow(x)
	nP <- ncol(x)

# 	isSens <- isSens[rk]

	delta <- (abs(xmin[,rk] - xmax[,rk]))
	deltap <- 1 + delta/2
	deltam <- 1 - delta/2
	XLIM <- c(min(0,deltam), max(rkx, deltap))

	colvar <- (1 - rkx) * (rkx >= deltam)*(rkx < 1)
	maxcloser <- ((1 - rkx)[rkx < deltam & rkx < 1])
	maxcloser <- min(maxcloser[maxcloser > max(colvar)])
	colvar <- 100 * (maxcloser - colvar)/maxcloser*(rkx < 1)*(rkx >= deltam)
	layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,1.5), heights=c(1))
	# layout.show(4)par(mai=c(1.1,1.1,0.2,0.2))
	plot(0,0,type="n",xlim=XLIM, ylim=c(0, spc*(nP- 1)*ncl+ncl),ann=FALSE,
      bty="n",yaxt='n')
	for(i in seq_len(nP)){
	  for(j in seq_len(ncl)){
		z <- spc * (i - 1)*ncl + (1+ncl - j) -1
		if(rkx[j,i] < deltam[j,i]){
			myCol <- head(colPal,1)
  # 	  }else if(rkx[i] > deltap[i]){
		}else if(rkx[j,i] > 1){
			myCol <- tail(colPal,1)
		}else{
			myCol <- colPal[colvar[j,i]]
		}
		rect(xleft=0,ybottom=z-0.45,xright=rkx[j,i],ytop=z+0.45, col=myCol, lwd=1)
		arrows(deltam[j,i], z, deltap[j,i], z, length=0.05, angle=90, code=3)
	  }
	}
	if(!is.null(rownames(x))){
	  clName <- rownames(x)
	  par(xpd=TRUE)
	  for(j in seq_len(ncl)){
		zlab <- spc * (0 - 1)*ncl + 2.5*(ncl-j+1) -1
		text(XLIM[2],5+zlab,clName[j],pos=2)
	  }
	  text(XLIM[2],5+spc * (0 - 1)*ncl + 2.5*(ncl+1) -1,"clusters:",pos=2)
	par(xpd=NA)

	}
	segments(x0=1,y0=-1,x1=1,y1=z+2,col="dodgerblue",lwd=2)
	title(main)
	par(xpd=TRUE)
# 	text(x= XLIM[1], y=(ncl+0.70)*(seq_along(rkx[1,])-1)+2.5,colnames(rkx),pos=2)
	text(x= XLIM[1], y=spc * (1:nP - 1)*ncl + ncl/2-0.5,colnames(rkx),pos=2)
	par(xpd=NA)
	# image.scale(z=rkx, col = heat.colors(12))
	color.bar(lut=colPal,min=min(rkx), max=max(rkx),n=2)
}


color.bar <- function(lut, min, max=-min, nticks=11,
                      ticks=seq(min, max, len=nticks),
                      title='', n = 3, w = 10) {
    scale = (length(lut)-1)/(max-min)
#     dev.new(width=1.75, height=5)
#     par(mai=c(0.2,0.2,0.2,0.2))
    plot(c(0,w), c(min,max), type='n', bty='n', xaxt='n', xlab='',
        yaxt='n', ylab='', main=title)
    axis(4, round(ticks,n), las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,w,y+1/scale, col=lut[i], border=NA)
    }
}


#
# GSABarplot <- function(x,xmin,xmax,isSens, colPal = c(rgb(0.95,0.95,0.95),
# colorRampPalette(brewer.pal(3,"Reds"))(101)),
#             main=""){
#   ### H-errorbar
#   # sensitivity ranking
#   rk <- order(x,decreasing = FALSE)
#   isSens <- isSens[rk]
#   rkx <- x[rk]
#
#   delta <- as.vector(abs(xmin[rk] - xmax[rk]))
#   deltap <- 1 + delta/2
#   deltam <- 1 - delta/2
#   XLIM <- c(min(0,deltam), max(rkx, deltap))
#
# #   maxcloser <- min((1 - rkx)[rkx < deltam & rkx < 1])
# #   colvar <- (1 - rkx) * (rkx >= deltam)*(rkx < 1)
# #   colvar <- 100 * (max(colvar) - (colvar))/maxcloser*(rkx < 1)
#   colvar <- (1 - rkx) * (rkx >= deltam)*(rkx < 1)
#   maxcloser <- ((1 - rkx)[rkx < deltam & rkx < 1])
#   maxcloser <- min(maxcloser[maxcloser > max(colvar)])
#   colvar <- 100 * (maxcloser - colvar)/maxcloser*(rkx < 1)*(rkx >= deltam)
#   layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,1.5), heights=c(1))
#   # layout.show(4)par(mai=c(1.1,1.1,0.2,0.2))
#   plot(0,0,type="n",xlim=XLIM, ylim=c(0,length(rkx)+1),ann=FALSE,
#         bty="n",yaxt='n')
#   for(i in seq_along(rkx)){
#     z <- i - 1
#     if(rkx[i] < deltam[i]){
#       myCol <- head(colPal,1)
# #     }else if(rkx[i] > deltap[i]){
#     }else if(rkx[i] > 1){
#       myCol <- tail(colPal,1)
#     }else{
#       myCol <- colPal[colvar[i]]
#     }
#     rect(xleft=0,ybottom=z-0.4,xright=rkx[i],ytop=z+0.4, col=myCol,
#           lwd=1+isSens[i]*3)
#     arrows(deltam[i], z, deltap[i], z, length=0.05, angle=90, code=3)
#   }
#   segments(x0=1,y0=-1,x1=1,y1=length(rkx)-0.6,col="dodgerblue",lwd=2)
#   title(main)
#   par(xpd=TRUE)
#   text(x= XLIM[1], y=seq_along(rkx)-1,names(rkx),pos=2)
#   par(xpd=NA)
#   # image.scale(z=rkx, col = heat.colors(12))
#   color.bar(lut=colPal,min=min(rkx), max=max(rkx),n=2)
# }
#
# GSABarplotClust <- function(x,xmin,xmax, colPal = c(rgb(0.95,0.95,0.95),
#                   colorRampPalette(brewer.pal(3,"Reds"))(101)),
#                   main="", FUN=max, spc=1.15){
#   ### H-errorbar
#   x_max <- apply(x,2,FUN)
#   # sensitivity ranking
#   rk <- order(x_max,decreasing = FALSE)
#   rkx <- x[,rk]
#
#   ncl <- nrow(x)
#   nP <- ncol(x)
#
# #   isSens <- isSens[rk]
#
#   delta <- (abs(xmin[,rk] - xmax[,rk]))
#   deltap <- 1 + delta/2
#   deltam <- 1 - delta/2
#   XLIM <- c(min(0,deltam), max(rkx, deltap))
#
#   colvar <- (1 - rkx) * (rkx >= deltam)*(rkx < 1)
#   maxcloser <- ((1 - rkx)[rkx < deltam & rkx < 1])
#   maxcloser <- min(maxcloser[maxcloser > max(colvar)])
#   colvar <- 100 * (maxcloser - colvar)/maxcloser*(rkx < 1)*(rkx >= deltam)
#   layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,1.5), heights=c(1))
#   # layout.show(4)par(mai=c(1.1,1.1,0.2,0.2))
#   plot(0,0,type="n",xlim=XLIM, ylim=c(0, spc*(nP- 1)*ncl+ncl),ann=FALSE,
#       bty="n",yaxt='n')
#   for(i in seq_len(nP)){
#     for(j in seq_len(ncl)){
#     z <- spc * (i - 1)*ncl + (1+ncl - j) -1
#     if(rkx[j,i] < deltam[j,i]){
#       myCol <- head(colPal,1)
#   #     }else if(rkx[i] > deltap[i]){
#     }else if(rkx[j,i] > 1){
#       myCol <- tail(colPal,1)
#     }else{
#       myCol <- colPal[colvar[j,i]]
#     }
#     rect(xleft=0,ybottom=z-0.45,xright=rkx[j,i],ytop=z+0.45, col=myCol, lwd=1)
#     arrows(deltam[j,i], z, deltap[j,i], z, length=0.05, angle=90, code=3)
#     }
#   }
#   if(!is.null(rownames(x))){
#     clName <- rownames(x)
#     par(xpd=TRUE)
#     for(j in seq_len(ncl)){
#     zlab <- spc * (0 - 1)*ncl + 2.5*(ncl-j+1) -1
#     text(XLIM[2],5+zlab,clName[j],pos=2)
#     }
#     text(XLIM[2],5+spc * (0 - 1)*ncl + 2.5*(ncl+1) -1,"clusters:",pos=2)
#   par(xpd=NA)
#
#   }
#   segments(x0=1,y0=-1,x1=1,y1=z+2,col="dodgerblue",lwd=2)
#   title(main)
#   par(xpd=TRUE)
# #   text(x= XLIM[1],
# y=(ncl+0.70)*(seq_along(rkx[1,])-1)+2.5,colnames(rkx),pos=2)
#   text(x= XLIM[1], y=spc * (1:nP - 1)*ncl + ncl/2-0.5,colnames(rkx),pos=2)
#   par(xpd=NA)
#   # image.scale(z=rkx, col = heat.colors(12))
#   color.bar(lut=colPal,min=min(rkx), max=max(rkx),n=2)
# }
#
