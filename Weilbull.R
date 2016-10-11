## Weibull mixture models

## 1. Newton-Raphson algorithm 
##    (n can be set to 1 for raw data, but be a freq vector for grouped data)
newton_weibull <- function(n, x, t, r_init = 1) {
	g <- function(r) {
		A <- sum(n * t * x^r * log(x))
		B <- sum(n * t * x^r)
		C <- sum(n * t * log(x))
		D <- sum(n * t)
		A / B - 1 / r - C / D
	}
	g_diff <- function(r) {
		A <- sum(n * t * x^r * log(x)^2)
		B <- sum(n * t * x^r)
		C <- sum(n * t * x^r * log(x))
		D <- sum(n * t * x^r)	# B and D are the same
		A / B - (C / D)^2 + 1 / r^2	
	}
	r <- r_init
	 repeat { 	
		r_new <- r - g(r) / g_diff(r)
		if(is.nan(r_new)) break
		if(abs(r_new - r) < 1e-4) break
		r <- r_new
	 }

	theta <- sum(n * t * x^r) / sum(n * t)
	
	k <- r
	lambda <- theta^{1 / r}
	res <- c(k, lambda)
	# names(res) <- c("k", "lambda")
	return(res)
}

# 2. parameter transformation
to.mu.sd <- function(k, lambda) {
	mu <- lambda * gamma(1 + 1/k)
	var <- lambda^2 * (gamma(1 + 2/k) - gamma(1 + 1/k)^2)
	sd <- sqrt(var)
	list(mu  = mu, sd = sd)
}
to.k.lambda <- function(mu, sd) {
	k <- (sd / mu)^(-1.086)
	lambda <- mu / gamma(1 + 1/k)
	list(k = k, lambda = lambda)
}


## 3. function for e-step (this function is not for users)
expZweib <- function(x, pi, k, lambda) {
	n <- length(x)
	pi <- pi / sum(pi)
	ncomp <- length(pi)
	res <- matrix(NA, nrow = n, ncol = ncomp)
	for(i in 1:n) {
		tmp <- pi * dweibull(x[i], k, lambda)
		res[i, ] <- tmp / sum(tmp)
	}
	res
}

## 4. a function for generating random numbers from weibull mixture
rmixweib <- function(n, pi, k, lambda, mu, sd) {
	if(missing(k) | missing(lambda)) {
		k <- to.k.lambda(mu = mu, sd = sd)$k
		lambda <- to.k.lambda(mu = mu, sd = sd)$lambda
	}
	nvec <- length(pi)
	rn <- as.vector(rmultinom(1, n, pi))
	res <- numeric()
	for(i in 1:nvec) {
		res <- c(res, rweibull(rn[i], k[i], lambda[i]))
	}
	res
}

## 5. log-likelihood for weibull mixture (raw data)
loglik_weib <- function(x, pi, k, lambda) {
	n <- length(x)
	p <- length(pi)
	tmp <- matrix(NA, nrow = n, ncol = p)
	for(i in 1:p) {
		tmp[, i] <- pi[i] * dweibull(x, k[i], lambda[i])
	}
	res <- apply(tmp, 1, sum)
	sum(log(res))
}


######################################################################
## Weibull mixture model fitting (raw data)
######################################################################

## 1. define generic function
weibEM <- function(x, ...) UseMethod("weibEM", x)

## 2. model fitting
weibEM.default <- function(x, pi, mu, sd, k, lambda, ncomp = NULL, tol = 1e-4, max_iter = 1000) {
	
	## initial values by k-means (mu, sd) will overwrite (k, lamdba)
	## if none of (pi, mu, sd) or (pi, k, lambda) is provided, ncomp should be provided
	cond1 <- missing(pi) | missing(mu) | missing(sd)
	cond2 <- missing(pi) | missing(k) | missing(lambda)
	if(!cond1) {
		k <- to.k.lambda(mu, sd)$k
		lambda <- to.k.lambda(mu, sd)$lambda
	} else if(!cond2) {
		# nothing
	} else {
		init <- initz(x, ncomp)
		pi <- init$pi
		mu <- init$mu
		sd <- init$sd
		k <- to.k.lambda(mu, sd)$k
		lambda <- to.k.lambda(mu, sd)$lambda
	}
	
	ncomp <- length(pi)
	iter <- 1
	repeat {
		## e-step
		tx <- expZweib(x, pi, k, lambda)
		
		## m-step
		# pi
		pi_new <- apply(tx, 2, mean)

		# Newton method to get estimate of k_j and lambda_j
		par <- matrix(NA, nrow = ncomp, ncol = 2)
		for(i in 1:ncomp) {
			par[i, ] <- newton_weibull(1, x, tx[, i], k[i])
		}
		
		k_new <- par[, 1]; lambda_new <- par[, 2]

		diff <- c(pi_new, k_new, lambda_new) - c(pi, k, lambda)
		if(sum(is.nan(diff)) > 0) break
		if(max(abs(diff)) < tol | iter > max_iter) break
		sd <- to.mu.sd(k_new, lambda_new)$sd
		if(max(sd) / min(sd) > 10) break			
		pi <- pi_new
		k <- k_new
		lambda <- lambda_new	
		iter <- iter + 1		
	}
	mu <- to.mu.sd(k_new, lambda_new)$mu
	sd <- to.mu.sd(k_new, lambda_new)$sd
	
	loglik <- loglik_weib(x, pi_new, k_new, lambda_new)
	aic <- -2* loglik + 2 * (3 * ncomp - 1)
	bic <- -2* loglik + log(n) * (3 * ncomp - 1)
	
	res <- list(pi = pi_new, mu = mu, sd = sd, k = k_new, lambda = lambda_new, 
				iter = iter, loglik = loglik, aic = aic, bic = bic, data = x)
	class(res) <- "mixweib"
	res			   	
}

## 2. Plotting (base R plotting system)
# function to find highest point in a histogram
mhist <- function(data) {
	count <- data[, 3]
	n <- nrow(data)
	height <- numeric(n)
	for(i in 1:n) {
		area <- count[i] / sum(count)
		height[i] <- area / (data[i, 2] - data[i, 1])
	}
	max(height)
}

plot.mixweib <- function(object, detail = FALSE, smooth = 300, col = "black", 
						 breaks, xlim, ylim, lty = 1, lwd, ...) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	k <- object$k
	lambda <- object$lambda
		
	data <- object$data
	xlow <- max(min(mu - 3.5 * sd), 0)
	xupp <- max(mu + 3.5 * sd)
	xseq <- seq(xlow, xupp, length = smooth)
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dweibull(xseq, k[i], lambda[i])
	}
	yt <- apply(res, 1, sum)
	
	# plot parameters checking
	if(missing(breaks)) {
		breaks <- 30
	}
	if(missing(xlim)) {
		xlim <- c(xlow, xupp)
	}	
	if(missing(ylim)) {
		if(is.matrix(data)) {
			ylim <- c(0, max(c(yt, mhist(data))))
		} else {
			brks <- seq(min(data), max(data), length = breaks + 1)
			tmp <- bin(data, brks = brks)
			count <- tmp[, 3]
			max_freq <- max(count) / (sum(count) * (brks[2] - brks[1]))
			ylim <- c(0, max(c(yt, max_freq)))
		}		
	}	
	if(is.matrix(data)) {
		breaks <- sort(unique(c(data[, 1], data[, 2])))
		data <- reinstate(data)	
	}
	if(missing(lwd)) {
		lwd <- 2
	}
	# plot	
	hist(data, freq = FALSE, breaks = breaks, xlim = xlim, ylim = ylim, ...)
	if(detail) {
		for(i in 1:length(pi)) {
			lines(xseq, res[, i], col = i + 1, lwd = 1.5)
		}
	}
	lines(xseq, yt, lty = lty, lwd = lwd, col = col)
}


## 3. lines (base R plotting system)
lines.mixweib <- function(object, detail = FALSE, smooth = 300, lwd = 2, ...) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	k <- object$k
	lambda <- object$lambda

	xlow <- max(min(mu - 3.5 * sd), 0)
	xupp <- max(mu + 3.5 * sd)
	xseq <- seq(xlow, xupp, length = smooth)
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dweibull(xseq, k[i], lambda[i])
	}
	yt <- apply(res, 1, sum)
	
	if(detail) {
		for(i in 1:length(pi)) {
			lines(xseq, res[, i], lty = i + 1)
		}
	}
	lines(xseq, yt, lwd = lwd, ...)
}

## 4. Plotting (base R plotting system)
gplot.mixweib <- function(object, detail = FALSE, smooth = 300, title = NULL, xlim, ylim, 
                 xlab, ylab, breaks, ..., theme = c("grey", "bw"), alpha = 0.5) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	k <- object$k
	lambda <- object$lambda
	data <- object$data
	ncomp <- length(pi)
		
	if(missing(xlab)) {
		xlab <- "Data"
	}
	if(missing(ylab)) {
		ylab <- "Density"
	}
	if(missing(breaks)) {
		breaks <- 30
	}
	if(missing(xlim)) {
		xlim <- c(max(min(mu - 3.5 * sd), 0), max(mu + 3.5 * sd))
	}
	
	# binwidth = (xlim[2] - xlim[1]) / breaks
	xseq <- seq(xlim[1], xlim[2], length = smooth)	
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dweibull(xseq, k[i], lambda[i])
	}
	yt <- apply(res, 1, sum)
	
	if(missing(ylim)) {
		if(is.matrix(data)) {
			count <- data[, 3]
			max_freq <- max(count) / (sum(count) * max(diff(data[, 1])))
			ylim <- c(0, max(c(yt, max_freq)))
		} else {
			brks <- seq(min(data), max(data), length = breaks + 1)
			tmp <- bin(data, brks = brks)
			count <- tmp[, 3]
			max_freq <- max(count) / (sum(count) * (brks[2] - brks[1]))
			ylim <- c(0, max(c(yt, max_freq)))
		}		
	}	
	if(is.matrix(data)) {
		breaks <- sort(unique(c(data[, 1], data[, 2])))
		data <- reinstate(data)	
	} else {
		breaks <- brks
	}
	
	# prepare data frames
	df1 <- data.frame(x = rep(xseq, ncomp), comp = rep(1:ncomp, each = smooth), 
					  y = as.vector(res))
	df2 <- data.frame(x = xseq, y = yt)
	
	# plot
	if(is.null(title)) title = "Weibull Mixture Density"
	if(detail) {
		add <- geom_polygon(data = df1, aes(x, y, fill = as.factor(comp)), alpha = alpha)
	} else {
		add <- NULL
	}
	if(theme[1] == "bw") {
		theme <- theme_bw()
	} else {
		theme <- theme_grey()
	}
	ggplot(as.data.frame(data)) + 
		geom_histogram(aes(x = data, y = ..density..),breaks = breaks, color = "black", 
				   fill = "white", size = 0.3) + add + theme + 
		geom_path(data = df2, aes(x, y), ...) +
	scale_fill_discrete(guide = guide_legend(title = "Comp")) +
	labs(title = title, x = xlab, y = ylab)
}

## 5. glines function
glines.mixweib <- function(object, smooth = 300, ...) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	k <- object$k
	lambda <- object$lambda

	xlow <- max(min(mu - 3.5 * sd), 0)
	xupp <- max(mu + 3.5 * sd)
	xseq <- seq(xlow, xupp, length = smooth)
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dweibull(xseq, k[i], lambda[i])
	}
	yt <- apply(res, 1, sum)
	
	tmp <- data.frame(xseq, yt)
	list(
	geom_line(data = tmp, aes(xseq, yt), ...)
	)
}

######################################################################
## Weibull mixture model fitting (grouped data)
######################################################################

## 1. expected value of X
EXweib <- function(data, k, lambda) {
	n <- nrow(data)
	ncomp <- length(k)	
	a = data[, 1]
	b = data[, 2]
	res <- matrix(NA, nrow = n, ncol = ncomp)	
	for(i in 1:n) {
		for(j in 1:ncomp) {
			tryCatch({
				res[i, j] <- extrunc(spec = "weibull", a[i], b[i], k[j], lambda[j])
			}, error = function(e) a[i]
			)			
		}
	}	
	res
}

## 2. expected value of Z
TXweib <- function(data, k, lambda, pi, ex) {
	n <- nrow(data)
	ncomp <- length(k)
	res <- matrix(NA, nrow = n, ncol = ncomp)
	for(i in 1:n) {
		A = pi * dweibull(ex[i, ], k, lambda)
		res[i, ] <- A / sum(A)
	}
	res
}

## 3. loglikelihood(for grouped data)
loglik_weib_g <- function(data, pi, k, lambda) {
	n <- nrow(data)
	res <- numeric(n)
	for(i in 1:n) {
		low <- pweibull(data[i, 1], k, lambda)
		upp <- pweibull(data[i, 2], k, lambda)
		tmp <- sum((upp - low) * pi)
		res[i] <- data[i,3] * log(tmp)
	}
	sum(res)	
}


## 4. model fitting
weibEM.matrix <- function(data, pi, mu, sd, k, lambda, ncomp = NULL, 
						  tol = 1e-4, max_iter = 1000) {	

	## initial values by k-means (mu, sd) will overwrite (k, lamdba)
	## if none of (pi, mu, sd) or (pi, k, lambda) is provided, ncomp should be provided
	cond1 <- missing(pi) | missing(mu) | missing(sd)
	cond2 <- missing(pi) | missing(k) | missing(lambda)
	if(!cond1) {
		k <- to.k.lambda(mu, sd)$k
		lambda <- to.k.lambda(mu, sd)$lambda
	} else if(!cond2) {
		# nothing
	} else {
		init <- initz(x, ncomp)
		pi <- init$pi
		mu <- init$mu
		sd <- init$sd
		k <- to.k.lambda(mu, sd)$k
		lambda <- to.k.lambda(mu, sd)$lambda
	}

	count = data[, 3]
	ncomp <- length(pi)
	
	## EM algorithm
	iter <- 1
	repeat {
		## e-step
		ex <- EXweib(data, k, lambda)
		tx <- TXweib(data, k, lambda, pi, ex)


		## m-step
		# pi
		pi_new <- numeric(ncomp)
		for(j in 1:ncomp) {
			pi_new[j] <- sum(tx[, j] * count) / sum(apply(tx, 1, sum) * count)
		}

		# Newton method to get estimate of k1, k2, ..., lambda1, lambda2,...
		par <- matrix(NA, nrow = ncomp, ncol = 2)
		for(j in 1:ncomp) {
			par[j, ] <- newton_weibull(n = count, ex[, j], tx[, j], k[j])
		}

		k_new <- par[, 1]
		lambda_new <- par[, 2]

		diff <- c(pi_new, k_new, lambda_new) - c(pi, k, lambda)
		# diff <- loglik_weib_g(data, pi_new, k_new, lambda_new) - loglik_weib_g(data, pi, k, lambda)

		if(sum(is.nan(diff)) > 0) break
		if(max(abs(diff)) < tol | iter > max_iter) break
		sd <- to.mu.sd(k_new, lambda_new)$sd
		if(max(sd) / min(sd) > 10) break			
		pi <- pi_new
		k <- k_new
		lambda <- lambda_new	
		iter <- iter + 1
	}
	
	mu <- to.mu.sd(k_new, lambda_new)$mu
	sd <- to.mu.sd(k_new, lambda_new)$sd
	
	loglik <- loglik_weib_g(data, pi_new, k_new, lambda_new)
	aic <- -2 * loglik + 2 * (3 * ncomp - 1)
	bic <- -2 * loglik + log(sum(count)) * (3 * ncomp - 1)
	
	res <- list(pi = pi_new, mu = mu, sd = sd, k = k_new, lambda = lambda_new, 
				iter = iter, loglik = loglik, aic = aic, bic = bic, data = data)
	class(res) <- "mixweib"
	res						  	
}

## 5. plotting (see plot.mixweib, we use the same plotting function as for raw data)








