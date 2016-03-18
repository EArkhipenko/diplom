require(nleqslv)

# Aaia?e?oai aoiaiua aaiiua
generate <- function(n, tet, sigma1=0.01,sigma2=0.01, sigma3=0.01, sigma4=0.5) {
	c(runif(n, -1, 1)) -> ksi
	eta <- tet[1] + ksi * tet[2]
	#rnorm(n, 0, sigma) -> delta

	rbinom(n, 1, 0.95) -> rb1
	rbinom(n, 1, 0.95) -> rb2
	sapply(rb1, function(x) rnorm(1, 0, switch(x + 1, sigma1, sigma2))) -> delta
	sapply(rb2, function(x) rnorm(1, 0, switch(x + 1, sigma3, sigma4))) -> eps
	list(x = ksi + delta, y = eta + eps)
}

mnk <- function (x, y) { 
	n <- length(y)
	x <- cbind(rep(1,n), x)
	c(solve(t(x) %*% (x)) %*% t(x) %*% matrix(y)) 
}

# Au?eneyai aaooa
rcr <- function(x, y, te = c(1.8, 0.4), b0=0.001, h=0.001, eps=1.e-8, k=0) {
	r <- sqrt((x-mean(x))^2 + (y - mean(y))^2)
	Sxx <- sum((x-mean(x))^2 / r^2)
	Sxy <- sum((x-mean(x)) * (y - mean(y)) / r^2)
	Syy <- sum((y-mean(y))^2 / r^2)
	my <- mean(y)
	mx <- mean(x)

	theta <- function(gamma) {
		gamma0 <- 1 - gamma
		f <- function(b) (gamma0 + gamma / b^2) * (Syy + b^2 * Sxx - 2 * b * Sxy)
		df <- function(b) (f(b + h) - f(b)) / h
		res <- nleqslv(b0, df, method = "Newton")
		b <- res$x
		c(my - b * mx, b)
	}

	psi <- function (gamma) {
		t <- theta(gamma)
		sum((t - te)^2 /  te^2)
	}


	#plot(psi, 0, 1)

	# g <- 
	# p <- sapply(g, psi)
	# a <- sapply(g, 

	# x <- data.frame(g = g, psi = p)
	# write.table(x, file = "rcr.csv", sep = ",", col.names = NA, qmethod = "double")

	for (g in (0:100) / 100) {
		print(c(g, psi(g), theta(g)))
	}

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

# ��������� �����
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

als <- function(x, y, sigma_init=0.01, eps=1.e-8) {
	m <- 2
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
		((n * sd2)/(sum(y^2 - R * tet)) - (var(x) / var(y)))^2
	}

	o = optimize(F, interval=c(0, 10))
	print(o$minimum)
	plot(F, 0, 10)
	theta(o$minimum)
}


data <- generate(500, c(1.8, 0.37))
tet <- rcr(data$x, data$y)
tet <- mnk(data$x, data$y)

# n <- 100
# for (i in 1:n) {
# 	data <- generate(500, c(1.8, 0.4))
# 	t <- als(data$x, data$y)
# 	tet <- tet+t
# }

# tet <- tet / n