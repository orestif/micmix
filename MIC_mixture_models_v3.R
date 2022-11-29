# ==============================================================
#
#             Mixture model analysis for MIC data
#
# ==============================================================

library(parallel)
options(mc.cores=detectCores())


# ----------------------- MAIN FUNCTION -----------------------------------

# Fit a series of mixture models ranging from 2 to K components.
# Arguments: MIC = vector of recorded MIC values, range = vector of concentrations tested.
# MIC can contain the value Inf, meaning MIC > max(range).
# We assume that the true MIC value is lower than or equal to the measured MIC. For example, if measurements are made at concentrations 16 and 32, a reported value of 32 means that the bacteria grew with concentrations up to 16, and were inhibited at concentrations from 32 and above. The true MIC could be anywhere in the range 16 < x <= 32.
# Hence MIC should not contain any values equal to 0 (they should be set to min(range))
# The concentrations are then transformed using trans(), which is log10() by default.

fit.multi.mix <- function(mic,range, K=4, trans=log10, inv.trans=function(x){10^x}, trace=F, n.runs=c(8,3)){
	
	if(K<2) stop("The number of mixture models must be greater than 1.")

	# Format the data for the MLE: transforms a vector of MIC values into a dataframe with three columns: Low, High, N.
	count <- data.frame(Low=c(0,range), High=c(range,Inf))
	count$N <- apply(count,1,function(x){
		length(mic[mic>x[1] & mic<=x[2]])
	})
	# Process data
	trans.count <- count %>% transmute(Low=trans(Low),High=trans(High),N=N)
	
	min.v <- min(trans.count$High)
	max.v <- max(trim(trans.count$High))
	mean.v <- with(trans.count,sum(trim(High)*trim(N))/sum(trim(N)))
	sd.v <- with(trans.count,sqrt(sum((trim(High)-mean.v)^2*trim(N)))/sum(trim(N)))
	
	# Caculate MLE with 1:K normal distributions
	fit <- lapply(1:K, function(k){
		init <- c(
			m = min.v + (1:k)*(max.v-min.v)/(k+1),
			sd = rep(sd.v/k,k),
			w = rep(1/k,k-1)
		)
		cat("\n---------------------------\nk=",k,"\n","init: ",init,"\n")
		# MLE:
		fit.k <- mix.norm.mle(k,init,trans.count,trace,n.runs)
		if(k==1) fit.k$par["w"] <- 1 else fit.k$par[paste("w",k,sep="")] <- 1-sum(fit.k$par[(2*k+1):(3*k-1)])
		print(fit.k$par,digits=3)
		return(fit.k)
	})
	# AIC values
	AIC <- sapply(fit,function(fit.k) {2*(fit.k$value + length(fit.k$par))})
	cat("\nBest model: K=",which.min(AIC),"\n")
	
	# Calculate predicted frequencies
	values <- trans.count$High
	pred <- lapply(1:K, function(k){
		X <- with(fit[[k]],sapply(1:nrow(trans.count), function(i){pmix(trans.count$High[i],par,k) - pmix(trans.count$Low[i],par,k)}))
		if (k==1) t(X) else X
	})
	# List of K dataframes, each showing the Predicted number of isolates, and the Classification of existing isolates based on the relative probabilities
	pred.n <- lapply(pred, function(pred.k){
		sum.p <- apply(pred.k,2,sum)
		all <- count
		for(i in 1:nrow(pred.k)){
			all <- all %>% bind_cols(data.frame(p=sum(count$N,na.rm=T)*pred.k[i,],
									n=ifelse(sum.p>0,count$N*pred.k[i,]/sum.p,0)))
			names(all)[2*i+2:3] <- paste(c("Pred","Class"),i,sep=".")
		}
		all
	})
	# Single cutoff (separating two types), only for k>1
	CO.2 <- do.call(rbind,lapply(2:K, function(k){ 
		# Compute a matrix of cumulative probabilities with each component in a column
		X <- with(fit[[k]],t(sapply(1:nrow(count), function(i){pmix(trans.count$High[i],par,k)})))
		as.data.frame(t(sapply(1:(k-1), function(i){
				# Identify the optimum position
				j <- which.min(apply(as.matrix(X[,(i+1):k]),1,sum) - apply(as.matrix(X[,1:i]),1,sum))
				# Pick a cutoff just above the optimum position
				co <- if(j<nrow(trans.count)) inv.trans((trans.count$High[j]+trans.count$High[j+1])/2) else inv.trans(trans.count$High[j])
				# Number of strains falling below or above
				N.wt <- sum(count$N[count$High<co])
				N.nwt <- sum(count$N[count$High>co])
				# Number of isolates mis-classified by the cutoff
				N.mis <- pred.n[[k]] %>% filter(High>co) %>% select(starts_with("Class")) %>% select(1:i) %>% sum(na.rm=T) +
					pred.n[[k]] %>% filter(High<co) %>% select(starts_with("Class")) %>% select((i+1):k) %>% sum(na.rm=T)
				c(k=k,split=i,cut.off=co,WT=N.wt,NWT=N.nwt,N.mis=N.mis)
		})))
	}))

	return(list(
		par = lapply(fit,function(fit.k) fit.k$par),
		LL = sapply(fit,function(fit.k) {-fit.k$value}),
		AIC = AIC,
		Best = which.min(AIC),
		convergence=sapply(fit,function(fit.k) {(fit.k$convergence==0)}),
		hessian=sapply(fit,function(fit.k) {fit.k$hessian}),
		eigen.hessian=lapply(fit,function(fit.k) {fit.k$eigen.hessian}),
		pred.n = pred.n,
		cut.offs = CO.2
	))
}


# =========================================== HELPER FUNCTIONS - NOT TO BE CALLED BY USERS ========================================

# Return a vector with NA removed
rm.na <- function(x) {x[which(is.finite(x))]}

# Remove the last element(s) of a vector
trim <- function(x,n=1) {
	k <- length(x)
	x[-seq(k-n+1,k)]
}


# Helper function to calculate the cumulative probability of a mixed model
# Return the vector of weighted CDF
pmix <- function(x,par,k){
	par[3*k] <- if(k==1) 1 else 1-sum(par[(2*k+1):(3*k-1)])
	sapply(1:k,function(i){
		par[2*k+i]*pnorm(x,par[i],par[k+i])
	})
}


# ----------------------------------- MLE sub-routine -------------------------------------------------
# Fits a mixture of k normal distributions to the transformed MIC count data
# This version uses the R function optim()
# Try random variations of init values and return the best. They are all run in parallel using mclapply.
mix.norm.mle <- function(k,init,count, trace, n.runs=c(8,3)){
	fit <- mclapply(1:(n.runs[1]), function(zz) {
		for(xx in 1:(n.runs[2])){
			sub <- optim(init*rnorm(length(init),1,0.2/xx), function(par){
				if(any(par[(k+1):(2*k)]<0) | any(par[1:(k-1)]>par[min(2,k):k])) return(1E10)
				if(k>1) {if(any(par[(2*k+1):(3*k-1)]<0) | sum(par[(2*k+1):(3*k-1)])>1) return(1E10)}
				
				probs <- sapply(1:nrow(count), function(i){
					with(count,sum(pmix(High[i],par,k) - pmix(Low[i],par,k)))
				})
				if(!is.finite(sum(probs))) return(1E10)
				ll <- -dmultinom(count$N,prob=probs,log = T)
				if(!is.finite(ll)) 1E10 else ll
			},hessian = T)
			init <- sub$par
		}
		sub
	})
	if(trace) cat("Attempts: ",round(sapply(fit,function(X) X$value),3),"\n")
	best <- which.min(sapply(fit,function(X) X$value))
	if(fit[[best]]$value>1E9){
		if(trace) cat("** Warning: MLE failed **\n")
		warning("MLE failed")
	} 
	return(fit[[best]])
}
