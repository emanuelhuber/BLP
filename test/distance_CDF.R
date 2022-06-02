
#---------------------- DISTANCE BETWEEN TWO CDF ------------------------------#
# study the statistical variability.
x <- rnorm(100)
y <- rnorm(100, mean = 0, sd = 1)
plot(y)
x_cdf <- ecdf(x)
y_cdf <- ecdf(y)

xy_knots <- sort(c(knots(x_cdf), knots(y_cdf)))
y_val <- y_cdf(xy_knots)
x_val <- x_cdf(xy_knots)

plot(xy_knots, x_val, type = "s")
lines(xy_knots, y_val, col = "red", type = "s")
lines(xy_knots, abs(x_val - y_val), type = "s", col ="blue")
n <- length(x_val)

dCDF <- sum(abs(x_val[-n] - y_val[-n]) * diff(xy_knots))
dCDF



nit <- 100
vmean <- seq(0, by = 0.1, to = 15)
vdist <- matrix(nrow = length(vmean), ncol = nit)
for(i in seq_along(vmean)){
  for(j in seq_len(nit)){
    y <- rnorm(100, mean = vmean[i], sd = 1)
    vdist[i, j] <- distCDF(x, y)
  }
}

matplot(vmean, vdist, type = "l", pch = 20, col = rgb(0.2, 0.2, 0.2, 0.2))
boxplot(t(vdist), names = vmean)

#-------------------- DISTANCE BETWEEN TWO sub-CDF ----------------------------#

ns <- 50
x <- c(rnorm(ns), rnorm(ns, mean = 2), rnorm(ns, mean = 4))
cl <- rep(1:3, each = ns)

distCDFCl <- function(x, cl, method = c("EMD", "KS", "CMD")){
  cla <- unique(cl) # classes
  ncl <- length(cla)
  x_cdf <- ecdf(x)

  xsub <- vector(length = ncl, mode = "list")
  for(i in seq_len(ncl)){
    xsub[[i]] <- ecdf(x[cl == cla[i]])
  }
  x_knots <- knots(x_cdf)
  # x_knots <- sort(as.vector(sapply(xsub, knots)))
  x_val <- x_cdf(x_knots)
  xi_val <- sapply(xsub, function(x, xknts){ x(xknts)}, x_knots)
  n <- length(x)

  test <- abs(sweep(xi_val[-n, ], 1, x_val[-n], "-"))
  if(method == "EMD"){
    test2 <- sweep(test, 1, diff(x_knots), "*")
    dCDF <- colSums(test2)
  }else if(method == "KS"){
    dCDF <- max(test)
  }else if(method == "CMD"){
    test2 <- sweep(test^2, 1, diff(x_knots), "*")
    dCDF <- sqrt(colSums(test2))
  }
  return(dCDF)
}

cla <- unique(cl) # classes
ncl <- length(cla)
x_cdf <- ecdf(x)

xsub <- vector(length = ncl, mode = "list")
for(i in seq_len(ncl)){
  xsub[[i]] <- ecdf(x[cl == cla[i]])
}
x_knots <- knots(x_cdf)
# x_knots <- sort(as.vector(sapply(xsub, knots)))
xi_val <- sapply(xsub, function(x, xknts){ x(xknts)}, x_knots)


plot(x_knots, x_cdf(x_knots), type = "s")
lines(x_knots, xi_val[,1], type = "s", col = "blue")
lines(x_knots, xi_val[,2], type = "s", col = "red")
lines(x_knots, xi_val[,3], type = "s", col = "green")

dd <- n
dd <- distCDF(x, x[cl == cla[i]])

distCDFCl(x, cl, method = "EMD")
