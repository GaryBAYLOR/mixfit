## 1. function for generating random data from normal mixture model
rmixnorm <- function(n, pi, mu, sd) {
	res <- numeric(0)
	for(i in 1:length(pi)) {
		tmp <- rnorm(n*pi[i], mu[i], sd[i])
		res <- c(res, tmp)
	}
	res
}


## 2. use k-means to set inital values (for both raw data and grouped data)
initz <- function(x, ncomp) {
	# check if 'x' is a matrix (from grouped data)
	if(is.matrix(x)) {
		x <- reinstate(x)
	} 
	a <- kmeans(x, centers = ncomp)$cluster
	res <- list()
	for(i in 1:ncomp) {
		res[[i]] <- x[a == i]
	}
	count <- sapply(res, length)
	pi <- count / sum(count)
	mu <- sapply(res, mean)
	sd <- sapply(res, sd)
	order <- order(mu)
	
	pi <- pi[order]
	mu <- mu[order]
	sd <- sd[order]
	list(pi = pi, mu = mu, sd = sd)
}


## 3. function for e-step (this function is not for users)
expz <- function(x, pi, mu, sd) {
	n <- length(x)
	pi <- pi / sum(pi)
	ncomp <- length(pi)
	res <- matrix(NA, nrow = n, ncol = ncomp)
	for(i in 1:n) {
		tmp <- pi * dnorm(x[i], mu, sd)
		res[i, ] <- tmp / sum(tmp)
	}
	res
}

## 4. loglikelihood calucation (not for users)
loglik <- function(x, pi, mu, sd) {
	n <- length(x)
	p <- length(pi)
	tmp <- matrix(NA, nrow = n, ncol = p)
	for(i in 1:p) {
		tmp[, i] <- pi[i] * dnorm(x, mu[i], sd[i])
	}
	res <- apply(tmp, 1, sum)
	sum(log(res))
}

## 5. reinstate pseudo original data for grouped data
## (used for plotting grouped data)
reinstate <- function(data) {
	rdata <- numeric()
	for(i in 1:nrow(data)) {
		rmin <- data[i, 1]
		rmax <- data[i, 2]
		r_num <- data[i, 3]
		rdata <- c(rdata, runif(r_num, rmin, rmax))
	}
	rdata
}	

######################################################################
## normal mixture model fitting (Using original data)
######################################################################

## 1. define generic function
normEM <- function(x, ...) UseMethod("normEM", x)

## 2. default method (use original data)
normEM.default <- function(x, pi, mu, sd, ncomp = NULL, ev = FALSE,
				   tol = 1e-4, max_iter = 1000) {	
	# check if initial values are provided
	if(missing(pi) | missing(mu) | missing(sd)) {
		init <- initz(x, ncomp = ncomp)
		pi <- init$pi
		mu <- init$mu
		sd <- init$sd
	}
	
	n <- length(x)
	ncomp <- length(pi)
	iter <- 1
	
	if(!ev) {
		repeat {
			# e-step
			Z <- expz(x, pi, mu, sd)
		
			# m-step
			pi_new <- apply(Z, 2, mean)
			mu_new <- as.vector(x %*% Z / apply(Z, 2, sum))
			var_new <- numeric(ncomp)
			for(j in 1:ncomp) {
				var_new[j] <- sum(Z[, j] * (x - mu_new[j])^2) / sum(Z[, j])
			}
			sd_new <- sqrt(var_new)
		
			# break (two ways)
			# diff <- c(pi_new, mu_new, sd_new) - c(pi, mu, sd)
			# if(max(abs(diff)) < tol | iter > max_iter) break
			diff2 <- abs(loglik(x, pi_new, mu_new, sd_new) - loglik(x, pi, mu, sd))
			if(max(abs(diff2)) < tol | iter > max_iter) break
			pi <- pi_new
			mu <- mu_new
			sd <- sd_new
			iter <- iter + 1
		}
	} else {
		repeat {
			# e-step
			Z <- expz(x, pi, mu, sd)
		
			# m-step
			pi_new <- apply(Z, 2, mean)
			mu_new <- as.vector(x %*% Z / apply(Z, 2, sum))
			tmp <- numeric(ncomp)
			for(j in 1:ncomp) {
				tmp[j] <- sum(Z[, j] * (x - mu_new[j])^2)
			}
			var_new <- sum(tmp) / sum(Z)
			var_new <- rep(var_new, ncomp)			
			sd_new <- sqrt(var_new)
		
			# break (two ways)
			# diff <- c(pi_new, mu_new, sd_new) - c(pi, mu, sd)
			# if(max(abs(diff)) < tol | iter > max_iter) break
			diff2 <- abs(loglik(x, pi_new, mu_new, sd_new) - loglik(x, pi, mu, sd))
			if(max(abs(diff2)) < tol | iter > max_iter) break
			pi <- pi_new
			mu <- mu_new
			sd <- sd_new
			iter <- iter + 1
		}
	}
		
	loglik <- loglik(x, pi_new, mu_new, sd_new)
	aic <- -2* loglik + 2 * ifelse(ev, 2 * ncomp, 3 * ncomp - 1)
	bic <- -2* loglik + log(n) * ifelse(ev, 2 * ncomp, 3 * ncomp - 1)
	
	res <- list(pi = pi_new, mu = mu_new, sd = sd_new, iter = iter, loglik = loglik, 
			    aic = aic, bic = bic, data = x)
	class(res) <- "mixnorm"
	return(res)
}

## 3. plotting (base R plotting system)
plot.mixnorm <- function(object, detail = FALSE, smooth = 300, col = "black", 
						 breaks, xlim, ylim, lwd, lty, ...) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	data <- object$data
	xlow <- min(mu - 3.5 * sd)
	xupp <- max(mu + 3.5 * sd)
	xseq <- seq(xlow, xupp, length = smooth)
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dnorm(xseq, mu[i], sd[i])
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

## 4. lines (base R plotting system)
lines.mixnorm <- function(object, detail = FALSE, smooth = 300, lwd = 2, ...) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	# data <- object$data
	xlow <- min(mu - 3.5 * sd)
	xupp <- max(mu + 3.5 * sd)
	xseq <- seq(xlow, xupp, length = smooth)
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dnorm(xseq, mu[i], sd[i])
	}
	yt <- apply(res, 1, sum)
	
	if(detail) {
		for(i in 1:length(pi)) {
			lines(xseq, res[, i], lty = i + 1)
		}
	}
	lines(xseq, yt, lwd = lwd, ...)
}

## 5. plotting (ggplot system)

# define generic function
gplot <- function(x, ...) UseMethod("gplot", x)

# ggplot (we use a name 'gplot')
gplot.mixnorm <- function(object, detail = FALSE, smooth = 300, title = NULL, xlim, ylim, 
                 xlab, ylab, breaks, ..., theme = c("grey", "bw"), alpha = 0.5) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
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
		xlim <- c(min(mu - 3.5 * sd), max(mu + 3.5 * sd))
	}
	
	# binwidth = (xlim[2] - xlim[1]) / breaks
	xseq <- seq(xlim[1], xlim[2], length = smooth)	
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dnorm(xseq, mu[i], sd[i])
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
	if(is.null(title)) title = "Mixture density"
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

# glines (used to add fitted density curve for another model)
glines <- function(x, ...) UseMethod("glines", x)

glines.mixnorm <- function(object, smooth = 300, ...) {
	pi <- object$pi
	mu <- object$mu
	sd <- object$sd
	
	xlow <- min(mu - 3.5 * sd)
	xupp <- max(mu + 3.5 * sd)
	xseq <- seq(xlow, xupp, length = smooth)
	res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
	for(i in 1:length(pi)) {
		res[ ,i] <-  pi[i] * dnorm(xseq, mu[i], sd[i])
	}
	yt <- apply(res, 1, sum)
	
	tmp <- data.frame(xseq, yt)
	list(
	geom_line(data = tmp, aes(xseq, yt), ...)
	)
}

######################################################################
## normal mixture model fitting (Using grouped data)
######################################################################
# library(truncdist)

## -1. initial value for grouped data (should be deprecated)
# initz_g <- function(data, ncomp) {
# 	rdata <- numeric()
# 	for(i in 1:nrow(data)) {
# 		rmin <- data[i, 1]
# 		rmax <- data[i, 2]
# 		r_num <- data[i, 3]
# 		set.seed(101)
# 		gen <- runif(r_num, rmin, rmax)
# 		# rm(.Random.seed, envir=.GlobalEnv)
# 		rdata <- c(rdata, gen)
# 	}
# 	a <- kmeans(rdata, centers = ncomp)$cluster
# 	res <- list()
# 	for(i in 1:ncomp) {
# 		res[[i]] <- rdata[a == i]
# 	}
# 	count <- sapply(res, length)
# 	pi <- count / sum(count)
# 	mu <- sapply(res, mean)
# 	sd <- sapply(res, sd)
# 	order <- order(mu)
# 	
# 	pi <- pi[order]
# 	mu <- mu[order]
# 	sd <- sd[order]
# 	list(pi = pi, mu = mu, sd = sd)
# }


## 0. binning data
bin <- function(x, brks) {
	k <- length(brks)
	res <- matrix(NA, nrow = k-1, ncol = 3)
	colnames(res) <- c("a", "b", "freq")
	res[,1] <- brks[-k]
	res[,2] <- brks[-1]
	temp <- .bincode(x, brks)
	temp <- temp[!is.na(temp)]
	
	for(i in 1:(k-1)) {
		res[i, 3] <- sum(temp == i)
	}	
	res[res[, 3] != 0, ]
}


## 1. conditional expectation of X (actually, this function is the same as that in weibull)
EXnorm <- function(data, mu, sd) {
	n <- nrow(data)
	ncomp <- length(mu)	
	a <- data[, 1]
	b <- data[, 2]
	res <- matrix(NA, nrow = n, ncol = ncomp)	
	for(i in 1:n) {
		for(j in 1:ncomp) {
			res[i, j] <- tryCatch({
				extrunc(spec = "norm", a[i], b[i], mu[j], sd[j])
			}, error = function(e) a[i]
			)			
		}
	}	
	res
}

## 2. conditional expectation of Z
TXnorm <- function(pi, mu, sd, ex) {
	n <- nrow(ex)	# 
	ncomp <- length(mu)
	res <- matrix(NA, nrow = n, ncol = ncomp)
	for(i in 1:n) {
		A = pi * dnorm(ex[i, ], mu, sd)
		res[i, ] <- A / sum(A)
	}
	res
}

## 3. log-likelihood
loglik_norm_g <- function(data, pi, mu, sd) {
	ndata <- nrow(data)
	res <- numeric(ndata)
	for(i in 1:ndata) {
		low <- pnorm(data[i, 1], mu, sd)
		upp <- pnorm(data[i, 2], mu, sd)
		tmp <- sum((upp - low) * pi)
		res[i] <- data[i,3] * log(tmp)
	}
	sum(res)	
}

## 4. model fitting (use grouped data)
normEM.matrix <- function(data, pi, mu, sd, ncomp = NULL, ev = FALSE,
				   tol = 1e-4, max_iter = 1000) {
	# check if initial values are missing
	if(missing(pi) | missing(mu) | missing(sd)) {
		init <- initz(data, ncomp = ncomp)
		pi <- init$pi
		mu <- init$mu
		sd <- init$sd
	}
	
	count <- data[, 3]
	ncomp <- length(pi)
	
	iter <- 1
	if(!ev) {
		repeat {
			## e-step
			ex <- EXnorm(data, mu, sd)
			tx <- TXnorm(pi, mu, sd, ex)

			## m-step
			# pi
			pi_new <- numeric(ncomp)
			for(j in 1:ncomp) {
				pi_new[j] <- sum(tx[ ,j] * count) / sum(apply(tx, 1, sum) * count)
			}

			# mu
			mu_new <- numeric(ncomp)
			for(j in 1:ncomp) {
				mu_new[j] <- mean(count * tx[, j] * ex[, j]) / mean(count * tx[, j])
			}
		
			# sd
			var_new <- numeric(ncomp)
			for(j in 1:ncomp) {
				var_new[j] <- mean(count * tx[, j] * (ex[, j] - mu_new[j])^2) / mean(count * tx[, j])
			}
			sd_new <- sqrt(var_new)

			# diff <- c(pi_new, mu_new, sd_new) - c(pi, mu, sd)
			diff2 <- loglik_norm_g(data, pi_new, mu_new, sd_new) - loglik_norm_g(data, pi, mu, sd)
				 
			if(abs(diff2) < tol | iter > max_iter) break		
			pi <- pi_new
			mu <- mu_new
			sd <-sd_new	
			iter <- iter + 1
		}		
	} else {
		repeat {
			## e-step
			ex <- EXnorm(data, mu, sd)
			tx <- TXnorm(pi, mu, sd, ex)

			## m-step
			# pi
			pi_new <- numeric(ncomp)
			for(j in 1:ncomp) {
				pi_new[j] <- sum(tx[ ,j] * count) / sum(apply(tx, 1, sum) * count)
			}

			# mu
			mu_new <- numeric(ncomp)
			for(j in 1:ncomp) {
				mu_new[j] <- mean(count * tx[, j] * ex[, j]) / mean(count * tx[, j])
			}
		
			# sd
			tmp <- numeric(ncomp)
			for(j in 1:ncomp) {
				tmp[j] <- sum(count * tx[, j] * (ex[, j] - mu_new[j])^2) #/ mean(count * tx[, j])
			
			}
			var_new <- sum(tmp) / (count %*% tx %*% matrix(1, nrow = ncomp))
			var_new <- rep(var_new, ncomp)
			sd_new <- sqrt(var_new)

			# diff <- c(pi_new, mu_new, sd_new) - c(pi, mu, sd)
			diff2 <- loglik_norm_g(data, pi_new, mu_new, sd_new) - loglik_norm_g(data, pi, mu, sd)
				 
			if(abs(diff2) < tol | iter > max_iter) break		
			pi <- pi_new
			mu <- mu_new
			sd <-sd_new	
			iter <- iter + 1
		}
	}
		
	loglik <- loglik_norm_g(data, pi_new, mu_new, sd_new)
	aic <- -2* loglik + 2 * ifelse(ev, 2 * ncomp, 3 * ncomp - 1)
	bic <- -2* loglik + log(sum(count)) * ifelse(ev, 2 * ncomp, 3 * ncomp - 1)
	
	res <- list(pi = pi_new, mu = mu_new, sd = sd_new, iter = iter, loglik = loglik, 
			    aic = aic, bic = bic, data = data)
	class(res) <- "mixnorm"
	return(res)
}

## we will contiue from here next time (Sep 24)





