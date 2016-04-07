require(nleqslv)
require(BB)

# Aaia?e?oai aoiaiua aaiiua
generate <- function(n, tet, sigma1=0.01,sigma2=0.01, sigma3=0.01, sigma4=0.5) {
	c(runif(n, -1, 1)) -> ksi
	eta <- tet[1] + ksi * tet[2]
	#rnorm(n, 0, sigma) -> delta

	rbinom(n, 1, 1) -> rb1
	rbinom(n, 1, 1) -> rb2
	sapply(rb1, function(x) rnorm(1, 0, switch(x + 1, sigma1, sigma2))) -> delta
	sapply(rb2, function(x) rnorm(1, 0, switch(x + 1, sigma3, sigma4))) -> eps
	list(x = ksi + delta, y = eta + eps)
}

mnk <- function (x, y, te = c(1.8, 0.4)) {
	p <- length(te) - 1
	n <- length(y)
	x <- cbind(rep(1,n), sapply(1:p, function(i) x^i))
	c(solve(t(x) %*% (x)) %*% t(x) %*% matrix(y)) 
}

# Au?eneyai aaooa
rcr <- function(x, y, te = c(1.8, 0.4), b0=0.001, g0=0.001, h=0.0001, eps=1.e-8, k=0) {
	p <- length(te) - 1

	r <- sqrt(apply(sapply(1:p, function(i) (x^i - mean(x^i))^2), 1, sum) + (y - mean(y))^2)
	Sxx <- function (i, j) sum((x^i-mean(x^i)) * (x^j-mean(x^j)) / r^2)
	Sxy <- function (j) sum((x^j-mean(x^j)) * (y - mean(y)) / r^2)
	Syy <- sum((y-mean(y))^2 / r^2)


	SS <- function(b, gamma) {
		gamma0 <- 1 - sum(gamma)
		#gamma0 <- gamma[1]
		#gamma.other <- gamma[2:length(gamma)]
		Sxx.sum <- 0
		Sxy.sum <- 0
		for (i in 1:p) {
			for (j in 1:p) {
				Sxx.sum <- Sxx.sum + b[i] * b[j] * Sxx(i, j)
			}
			Sxy.sum <- Sxy.sum + b[i] * Sxy(i)
		}
		#sum(sapply(1:5, function (j) sum(sapply(1:5, Sxx, j=j))))
		(gamma0 + sum(gamma / b^2)) * (Syy + Sxx.sum - 2 * Sxy.sum)
	}
#	f <- function(v) SS(v[1:p], v[(p+1):(2*p)])
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
	b <- b0
	g <- rep(g0, p)
	repeat {
		#res <- nleqslv(rep(g0, p), df(fg(b)), method = "Newton")
		res <- BBsolve(par=g, fn=df(fg(b)))
		g <- res$x
		print(g)
		# res <- nleqslv(b, df(fb(g)), method = "Newton")
		res <- BBsolve(par=b, fn=df(fb(g)))
		b_new <- res$x
		if (norm(b_new - b) / norm(b) <= eps) {
			break
		}
		b <- b_new
	}

	alpha <- mean(y) - sum(b * sapply(1:p, function(i) (mean(x^i))))
	c(alpha, b)


	# psi <- function (gamma) {
	# 	t <- theta(gamma)
	# 	sum((t - te)^2 /  te^2)
	# }


	#plot(psi, 0, 1)

	# g <- 
	# p <- sapply(g, psi)
	# a <- sapply(g, 

	# x <- data.frame(g = g, psi = p)
	# write.table(x, file = "rcr.csv", sep = ",", col.names = NA, qmethod = "double")

	# for (g1 in (0:10) / 50) {
	# 	for (g2 in (0:10) / 50) {
	# 		g <- c(g1, g2)
	# 		print(c(g, psi(g), theta(g)))
	# 	}
	# }

	# for (g in (0:100) / 100) {
	# 	print(c(g, psi(g), theta(g)))
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

# Вычисляем бетта
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

als <- function(x, y, te, sigma_init=0.01, eps=1.e-8) {
	m <- length(te)
	n <- length(y)
	matrix(1, nrow=n, ncol=2*m) -> t
	matrix(0, nrow=m, ncol=m) -> P
	rep(1,m) -> R

	t[,2] <- x

	theta <- function(sd2) {
		for (i in 3:(2*m)) { t[,i] <- x * t[,i-1] - (i - 2) * sd2 * t[,i-2] }
		for (i in 1:m) { R[i] <- sum(t[,i] * y) }
		for (i in 1:m) { for (j in 1:m) { P[i,j] <- sum(t[,i+j-1]) } }
		c(solve(P) %*% matrix(R))
	}
	F <- function(sd2) {
		tet <- theta(sd2)
		#((n * sd2)/(sum(y^2) - sum(R * tet)) - (var(x) / var(y)))^2
		#y.est <- apply(sapply(1:m, function(i) x^(i-1) * tet[i]), 1, sum)
		#		(((n-m) * sd2)/sum((y - y.est)^2) - 1)^2
		(((n-m) * sd2)/sum((y - sum(R * tet))^2) - 1)^2
	}

	o = optimize(F, interval=c(0, 10))
	#o = BBoptim(par=0, fn=F)
	sd2 <- o$minimum
	print(sd2)
	print(F(sd2))
	#plot(F, 0, 10)
	theta(0.0001)
}

te <- c(1.8, 0.4, 0.8)

data <- generate(500, te)
b0 <- mnk(data$x, data$y, te=te)
tet <- rcr(data$x, data$y, te=te, b0=b0[2:length(b0)])
#tet <- mnk(data$x, data$y)

# n <- 100
# for (i in 1:n) {
# 	data <- generate(500, c(1.8, 0.4))
# 	t <- als(data$x, data$y)
# 	tet <- tet+t
# }

# tet <- tet / n