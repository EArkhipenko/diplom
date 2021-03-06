require(nleqslv)
require(xlsx)
library(rgl)

# Aaia?e?oai aoiaiua aaiiua
generate <- function(n, tet, sigma1=0.05, sigma2=0.5, sigma3=0.05, sigma4=0.8) {
  m <- length(tet)
  c(runif(n, -1, 1)) -> ksi
  eta <- apply(sapply(1:m, function(i) ksi^(i-1)), 1, function(r) sum(r * tet))

  rbinom(n, 1, 0.05) -> rb1
  rbinom(n, 1, 0.05) -> rb2

  sapply(rb1, function(x) rnorm(1, 0, switch(x + 1, sigma1, sigma2))) -> delta
  sapply(rb2, function(x) rnorm(1, 0, switch(x + 1, sigma3, sigma4))) -> eps
  X <- cbind(rep(1,n), apply(sapply(2:m, function(i) ksi^(i-1)), 2, function(r) r + delta))
  list(X = matrix(X, nrow = n, ncol = m), y = eta + eps, outliersX = rb1, outliersY = rb2)
}

mnk <- function (X, y) {
  c(solve(t(X) %*% (X)) %*% t(X) %*% matrix(y)) 
}

# Au?eneyai aaooa
rcr <- function(X, y, tet.init, g0=0.001, h=0.001, eps=1.e-8, g=NULL) {
  b0 <- tet.init[2:length(tet.init)]
  p <- length(tet.init) - 1

  r <- sqrt(apply(sapply(1:p, function(i) (X[,i+1] - mean(X[,i+1]))^2), 1, sum) + (y - mean(y))^2)
  fSxx <- function (i, j) sum((X[,i+1]-mean(X[,i+1])) * (X[,j+1]-mean(X[,j+1])) / r^2)
  fSxy <- function (j) sum((X[,j+1]-mean(X[,j+1])) * (y - mean(y)) / r^2)
  Syy <- sum((y-mean(y))^2 / r^2)
  
  Sxx <- matrix(0, nrow = p, ncol = p)
  Sxy <- numeric(p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sxx[i,j] <- fSxx(i, j)
    }
    Sxy[i] <- fSxy(i)
  }


  SS <- function(b, gamma) {
    #gamma0 <- 1 - sum(gamma)
    gamma0 <- gamma[1]
    gamma.other <- gamma[2:length(gamma)]
    Sxx.sum <- 0
    Sxy.sum <- 0
    for (i in 1:p) {
      for (j in 1:p) {
        Sxx.sum <- Sxx.sum + b[i] * b[j] * Sxx[i, j]
      }
      Sxy.sum <- Sxy.sum + b[i] * Sxy[i]
    }
    #sum(sapply(1:5, function (j) sum(sapply(1:5, Sxx, j=j))))
    (gamma0 + sum(gamma.other / b^2)) * (Syy + Sxx.sum - 2 * Sxy.sum)
  }
# f <- function(v) SS(v[1:p], v[(p+1):(2*p)])
  fg <- function(b) function(g) SS(b, g)
  fb <- function(g) function(b) SS(b, g)
  df <- function(f) {
    function (v) {
      l <- length(v)
      grad <- numeric(l)
      for (i in 1:l) {
        v.plus.h <- v
        v.plus.h[i] <- v.plus.h[i] + h 
        grad[i] <- (f(v.plus.h) - f(v)) / h
      }
      grad
    }
  }
  norm <- function(x) sqrt(sum(x^2))
  basis <- function(n, i) {
  	a <- numeric(n)
  	a[i] <- 1
  	a
  }

  g.init <- rep(g0, p)
  b.init <- b0
  calc_b <- function (g, b.init=b0) {
  	res <- nleqslv(b.init, df(fb(g)), method = "Newton")
  	res$x
  }
  
  F <- function(g.p) {
  	g <- c((1-sum(g.p)), g.p)
    b <- calc_b(g, b.init)
    b.init <- b
    #for (i in 1:(p+1)) if (g[i] < 0 || g[i] > 1) return(1000000.0)
    -sum(sapply(1:(p+1), function(i) SS(calc_b(basis(p+1, i), b.init), basis(p+1, i)) / SS(b, basis(p+1, i))))
  }

  theta <- function(gamma) {
  	b <- calc_b(gamma)
  	alpha <- mean(y) - sum(b * sapply(1:p, function(i) (mean(X[,i+1]))))
    c(alpha, b)
  }

 # plot ((1:100) / 100, sapply((1:100) / 100, F))
 #  x <- (1:100) / 100
 #  xn = length(x)
 #  z <- matrix(0, nrow = xn, ncol = xn)
 #  for (i in 1:xn){
 #  for (j in 1:xn){
 #    z[i,j] <- f(x[i], x[j])
 #  	}
	# }
 #  persp3d(z)
  if (!is.null(g)) return(theta(g))
  
  res <- optim(par=g.init, F)
  #res <- optimize(F, interval=c(0,1))
  g <- c((1-sum(res$par)), res$par)
  #g <- c((1-res$minimum), res$minimum)
  print(c(g, sum(g)))
  #print(c(-F(res$minimum),F(numeric(p)), theta(g)))
  print(c(-F(res$par), -F(numeric(p))))
  theta(g)
  


  # repeat {
  #   f <- fg(b)
  #   plot (sapply((0:100) / 100, f), (0:100) / 100)

  #   res <- nleqslv(g, df(fg(b)), method = "Newton")
  #   g <- res$x
  #   print (g)
  #   res <- nleqslv(b, df(fb(g)), method = "Newton")
  #   b_new <- res$x
  #   if (norm(b_new - b) / norm(b) <= eps) {
  #     break
  #   }
  #   b <- b_new
  # }


  # psi <- function (gamma) {
  #   t <- theta(gamma)
  #   sum((t - te)^2 /  te^2)
  # }


  # for (g in (1:100) / 100) {
  #   print(c(g, psi(c(1-g,g)), -F(g), theta(c(1-g,g))))
  # }

	# for (g1 in (0:100) / 100) {
	# 	for (g2 in (0:10) / 10) {
	# 		g <- c(1-g1-g2, g1, g2)
	# 		print(c(g, psi(g), -F(c( g1, g2)), theta(g)))
	# 	}
	# }
  



  #plot(psi, 0, 1)

  # g <- 
  # p <- sapply(g, psi)
  # a <- sapply(g, 

  # x <- data.frame(g = g, psi = p)
  # write.table(x, file = "rcr.csv", sep = ",",  qmethod = "double")

  # for (g1 in (0:10) / 50) {
  #   for (g2 in (0:10) / 50) {
  #     g <- c(g1, g2)
  #     print(c(g, psi(g), theta(g)))
  #   }
  # }



# F <- function(gamma) {
# t <- theta(gamma)
# b <- t[2]
# ((Syy - b * Sxy) / (Syy - b * Sxy + b^4 * Sxx - b^3 * Sxy) - gamma)^2
# }

# o = optimize(F, interval=c(-1, 100000))
# o$minimum

# repeat {

# gamma_new <- (Syy - b * Sxy) / (Syy - b * Sxy + b^4 * Sxx - b^3 * Sxy)


# f <- function(b) (gamma0 + gamma_new / b^2) * (Syy + b^2 * Sxx - 2 * b * Sxy)
# df <- function(b) (f(b + h) - f(b)) / h
# res <- nleqslv(b0, df, method = "Newton")
# b_new <- res$x

# print(gamma)
# if (abs(b_new - b) / abs(b) <= eps || k > 100) {
# break 
# }
# #print(k)
# #print(abs(b_new - b) / abs(b))
# b <- b_new
# gamma <- gamma_new
# #k = k+1

# }
}


rcr_err <- function(X, y, tet.init, g0=0.001, h=0.001, eps=1.e-8, g=NULL) {
  b0 <- tet.init[2:length(tet.init)]
  p <- length(tet.init) - 1
  n <- length(y)
  
  r <- sqrt(apply(sapply(1:p, function(i) (X[,i+1] - mean(X[,i+1]))^2), 1, sum) + (y - mean(y))^2)
  fSxx <- function (i, j) sum((X[,i+1]-mean(X[,i+1])) * (X[,j+1]-mean(X[,j+1])) / r^2)
  fSxy <- function (j) sum((X[,j+1]-mean(X[,j+1])) * (y - mean(y)) / r^2)
  Syy <- sum((y-mean(y))^2 / r^2)
  
  Sxx <- matrix(0, nrow = p, ncol = p)
  Sxy <- numeric(p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sxx[i,j] <- fSxx(i, j)
    }
    Sxy[i] <- fSxy(i)
  }
  
  
  SS <- function(b, gamma) {
    #gamma0 <- 1 - sum(gamma)
    gamma0 <- gamma[1]
    gamma.other <- gamma[2:length(gamma)]
    Sxx.sum <- 0
    Sxy.sum <- 0
    for (i in 1:p) {
      for (j in 1:p) {
        Sxx.sum <- Sxx.sum + b[i] * b[j] * Sxx[i, j]
      }
      Sxy.sum <- Sxy.sum + b[i] * Sxy[i]
    }
    (gamma0 + sum(gamma.other / b^2)) * (Syy + Sxx.sum - 2 * Sxy.sum)
  }

  fg <- function(b) function(g) SS(b, g)
  fb <- function(g) function(b) SS(b, g)
  df <- function(f) {
    function (v) {
      l <- length(v)
      grad <- numeric(l)
      for (i in 1:l) {
        v.plus.h <- v
        v.plus.h[i] <- v.plus.h[i] + h 
        grad[i] <- (f(v.plus.h) - f(v)) / h
      }
      grad
    }
  }
  norm <- function(x) sqrt(sum(x^2))
  basis <- function(n, i) {
    a <- numeric(n)
    a[i] <- 1
    a
  }
  
  g.init <- rep(g0, p)
  b.init <- b0
  calc_b <- function (g, b.init=b0) {
    res <- nleqslv(b.init, df(fb(g)), method = "Newton")
    res$x
  }
  
  
  theta <- function(gamma) {
    b <- calc_b(gamma)
    alpha <- mean(y) - sum(b * sapply(1:p, function(i) (mean(X[,i+1]))))
    c(alpha, b)
  }
  
  F <- function(g.p) {
    g <- c((1-sum(g.p)), g.p)
    b <- calc_b(g, b.init)
    y.est <- X %*% matrix(theta(g))
    for (i in 1:(p+1)) if (g[i] < 0 || g[i] > 1) return(1000000.0)
    sum(abs(y - y.est)) / n
  }

  res <- optim(par=g.init, F)
  #res <- optimize(F, interval=c(0,1))
  g <- c((1-sum(res$par)), res$par)
  #g <- c((1-res$minimum), res$minimum)
  print(c(g, sum(g)))
  #print(c(-F(res$minimum),F(numeric(p)), theta(g)))
  print(c(F(res$par), F(numeric(p))))
  theta(g)
}

rcr_g0 <- function(X, y, tet.init) {
  m <- length(tet.init)
  g <- numeric(m)
  g[1] <- 1
  rcr(X, y, tet.init, g=g)
}

rcr2 <- function(x, y, b0=0.001, h=0.001, eps=1.e-8, k=0) {
  r <- sqrt((x-mean(x))^2 + (y - mean(y))^2)
  Sxx <- sum((x-mean(x))^2 / r^2)
  Sxy <- sum((x-mean(x)) * (y - mean(y)) / r^2)
  Syy <- sum((y-mean(y))^2 / r^2)
  my <- mean(y)
  mx <- mean(x)

  theta <- function() {
    f <- function(b) ((1 - (Syy - b * Sxy) / (Syy - b * Sxy + b^4 * Sxx - b^3 * Sxy)) + ((Syy - b * Sxy) / (Syy - b * Sxy + b^4 * Sxx - b^3 * Sxy)) / b^2) * (Syy + b^2 * Sxx - 2 * b * Sxy)
    df <- function(b) (f(b + h) - f(b)) / h
    res <- nleqslv(b0, df, method = "Newton")
    b <- res$x
    c(my - b * mx, b)
  }

# for (g in (1:100) / 100) {
# print(c(g, theta(g)))
# }

  theta()
# repeat {

# gamma_new <- (Syy - b * Sxy) / (Syy - b * Sxy + b^4 * Sxx - b^3 * Sxy)


# f <- function(b) (gamma0 + gamma_new / b^2) * (Syy + b^2 * Sxx - 2 * b * Sxy)
# df <- function(b) (f(b + h) - f(b)) / h
# res <- nleqslv(b0, df, method = "Newton")
# b_new <- res$x

# print(gamma)
# if (abs(b_new - b) / abs(b) <= eps || k > 100) {
# break 
# }
# #print(k)
# #print(abs(b_new - b) / abs(b))
# b <- b_new
# gamma <- gamma_new
# #k = k+1

# }
}

rcr3 <- function(x, y, gamma_init=0.1, b0=0.001, h=0.001, eps=1.e-8) {

  r <- sqrt((x-mean(x))^2 + (y - mean(y))^2)
  Sxx <- sum((x-mean(x))^2 / r^2)
  Sxy <- sum((x-mean(x)) * (y - mean(y)) / r^2)
  Syy <- sum((y-mean(y))^2 / r^2)

  gamma <- gamma_init
  gamma0 <- 1 - gamma
  f  <- function(b) (gamma0 + gamma / b^2) * (Syy + b^2 * Sxx - 2 * b * Sxy)
  df <- function(b) (f(b + h) - f(b)) / h
  res <- nleqslv(b0, df, method = "Newton")
  b <- res$x

  repeat {
  
    gamma_new <- (Syy - b * Sxy) / (Syy - b * Sxy + b^4 * Sxx - b^3 * Sxy)

    gamma0 <- 1 - gamma_new
    f  <- function(b) (gamma0 + gamma_new / b^2) * (Syy + b^2 * Sxx - 2 * b * Sxy)
    df <- function(b) (f(b + h) - f(b)) / h
    res <- nleqslv(b0, df, method = "Newton")
    b_new <- res$x

    #print(gamma)
    if (abs(b_new - b) / abs(b) <= eps) {
      break
    }
    b <- b_new
    gamma <- gamma_new
  }
  c(mean(y) - b * mean(x), b)
}

als <- function(X, y, tet.init, sigma_init=0.01, eps=1.e-8) {
  m <- length(tet.init)
  n <- length(y)
  matrix(1, nrow=n, ncol=2*m) -> t
  matrix(0, nrow=m, ncol=m) -> P
  rep(1,m) -> R

  t[,2] <- X[,2]

  theta <- function(sd2) {
    for (i in 3:(2*m)) { t[,i] <- X[,2] * t[,i-1] - (i - 2) * sd2 * t[,i-2] }
    for (i in 1:m) { R[i] <- sum(t[,i] * y) }
    for (i in 1:m) { for (j in 1:m) { P[i,j] <- sum(t[,i+j-1]) } }
    c(solve(P) %*% matrix(R))
  }
  # F <- function(sd2) {
  #   tet <- theta(sd2)
  #   #((n * sd2)/(sum(y^2) - sum(R * tet)) - (var(x) / var(y)))^2
  #   #y.est <- apply(sapply(1:m, function(i) x^(i-1) * tet[i]), 1, sum)
  #   y.est <- apply(sapply(1:m, function(i) x^(i-1)), 1, function(r) sum(r * tet))
  #   (((n-m) * sd2)/sum((y - y.est)^2) - var(x)/var(y))^2
  #   #(((n-m) * sd2)/sum((y - sum(R * tet))^2) - 1)^2
  # }

  tet <- tet.init

  norm <- function(x) sqrt(sum(x^2))

  theta(0.01^2)

  #g <- var(x) / var(y)
  # g <- 1

  # print (c(var(x) / var(y), 1))
  # repeat {
    
  #   y.est <- apply(sapply(1:m, function(i) x^(i-1)), 1, function(r) sum(r * tet))

  #   # print (c(y- y.est))
  #   # break
  #   sd2 <- g * sum((y - y.est)^2) / (n - m)
  #   tet_new <- theta(sd2)
  #   nev <- norm(tet_new - tet) / norm(tet)
  #   print(c(nev, sd2))
  #   if (nev < eps) {
  #     break
  #   }
  #   tet <- tet_new
  # }
  # tet

  # o = optimize(F, interval=c(0, 10))
  # #o = BBoptim(par=0, fn=F)
  # sd2 <- o$minimum
  # print(sd2)
  # print(F(sd2))
  # #plot(F, 0, 10)
  # theta(0.0001)
}

err <- function (data, tet) {
  m <- length(tet)
  n <- length(data$y)
  y.est <- apply(sapply(1:m, function(i) data$X[,i]), 1, function(r) sum(r * tet))
  c(sum((data$y - y.est)^2) / n, sum(abs(data$y - y.est)) / n)
}

report <- function (N, te) {
  psi <- function (tet) sum((tet - te)^2 /  te^2)
  m <- length(te)
  out <- matrix(0, nrow=N, ncol=(m+2)*4)
  for (i in 1:N) {
    data <- generate(500, te)
    tet.mnk <- mnk(data$X, data$y)
    tet.als <- als(data$X, data$y, tet.init=tet.mnk)
    tet.rcr <- rcr(data$X, data$y, tet.init=tet.mnk)
    tet.rcr_g0 <- rcr_g0(data$X, data$y, tet.init=tet.mnk)
    out[i,] <- c(tet.mnk, err(data, tet.mnk), tet.als, err(data, tet.als), tet.rcr, err(data, tet.rcr), tet.rcr_g0, err(data, tet.rcr_g0))
  }
  x <- data.frame(out)
  write.xlsx(x, file = "report.xlsx")
}

model <- function(tet) {
  m <- length(tet)
  Vectorize(function (x) sum(sapply(1:m, function(i) x^(i-1)) * tet))
}

research <- function(f, m=3) {
  dx <- read.table(sprintf("data/x_%s.txt", f))
  dy <- read.table(sprintf("data/y_%s.txt", f))
  data <- list(X = sapply(1:m, function(i) dx$V1^(i-1)), y = dy$V1)

  tet.mnk <- mnk(data$X, data$y)
  print(sprintf("tet.mnk = (%s), err = %f", paste(tet.mnk, collapse=", "), err(data, tet.mnk)[2]))
  plot(model(tet.mnk), 0, 1)
  
  tet.als <- als(data$X, data$y, tet.init=tet.mnk)
  print(sprintf("tet.als = (%s), err = %f", paste(tet.als, collapse=", "), err(data, tet.als)[2]))
  plot(model(tet.als), 0, 1)
  
  tet.rcr <- rcr(data$X, data$y, tet.init=tet.mnk)
  print(sprintf("tet.rcr = (%s), err = %f", paste(tet.rcr, collapse=", "), err(data, tet.rcr)[2]))
  plot(model(tet.rcr), 0, 1)
  
  tet.rcr_err <- rcr_err(data$X, data$y, tet.init=tet.mnk)
  print(sprintf("tet.rcr_err = (%s), err = %f", paste(tet.rcr_err, collapse=", "), err(data, tet.rcr_err)[2]))
  plot(model(tet.rcr_err), 0, 1)
  
  tet.rcr_g0 <- rcr_g0(data$X, data$y, tet.init=tet.mnk)
  print(sprintf("tet.rcr_g0 = (%s), err = %f", paste(tet.rcr_g0, collapse=", "), err(data, tet.rcr_g0)[2]))
  plot(model(tet.rcr_g0), 0, 1)
}


outliersPlot <- function(te) {
  data <- generate(50, te)
  X.clean = data$X[data$outliersX == 0 & data$outliersY == 0,]
  y.clean = data$y[data$outliersX == 0 & data$outliersY == 0]
  X.outliers = data$X[data$outliersX == 1 | data$outliersY == 1,]
  y.outliers = data$y[data$outliersX == 1 | data$outliersY == 1]
  
  tet.mnk <- mnk(data$X, data$y)
  tet.als <- als(data$X, data$y, tet.init=tet.mnk)
  tet.als.clean <- als(X.clean, y.clean, tet.init=tet.mnk)
  
  plot(X.clean[,2], y.clean, col="black")
  plot(model(tet.als), -1, 1, col="red", lwd=2, add=TRUE)
  plot(model(tet.als.clean), -1, 1, col="black", lwd=2, add=TRUE)
  points(X.outliers[,2], y.outliers, col="red")
}

te <- c(1.3, 2.1, 0.4)
#report(1, te)

#research("2_year_421")
outliersPlot(te)

#data <- generate(500, te)
# plot(data$X[,2], data$y)
#b0 <- mnk(data$x, data$y, te=te)
#tet <- rcr(data$x, data$y, te=te, tet.init=b0)



