library(matrixStats)
library(xtable)

BinaryESS2 <- function(V,M){
###Function to calculate ESS for binary data
#M: posterior mean
#V: posterior variance

	u <- (M*(1-M))/V - 1
	a <- M*u
	b <- (1-M)*u
	return(a+b)

}


calc.MEM.betabin <- function(xvec, nvec, avec, bvec, prior, constraint=1){
###Function to calculate the MEM model weights for binomial data with beta(alpha,beta) prior
#xvec: vector with counts of those with event (primary source first, supplemental afterwards)
#nvec: vector of total number of subjects
#avec: vector of alpha parameters for beta priors for each source
#bvec: vector of beta parameters for beta priors for each source
#prior: prior to use on MEM source inclusion probability: equal (pi_e), pool (naive pooling), opt.source (pi_EB), opt.source_constrain (pi_EBc)
#constraint: value to limit maximum value of empirical Bayes prior to for pi_EBc

	mod.mat <- as.matrix(expand.grid( rep(list(c(0,1)), length(xvec)-1) )) #create matrix of all models with potential combinations of supplemental sources
	mod.mat <- mod.mat[order(rowSums(mod.mat)),] #group by number of supplemental sources
	mod.mat <- cbind(1, mod.mat) #add column before 1st column for primary source indicator
	colnames(mod.mat) <- c('p',paste0('s',seq(1,length(xvec)-1)))

	prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
	p.vec <- apply( t(sapply(1:dim(mod.mat)[1], function(x) prod.vec^(1-mod.mat[x,]))), 1, prod) #calculate the product portion of the integrated marginal likelihood corresponding to each model by taking power w/indicator for each source in a model and then using apply over rows of resulting matrix

	###Calculate the integrated marginal likelihood given data:
	marg.vec <- (beta(avec[1] + mod.mat%*%xvec, bvec[1] + mod.mat%*%(nvec-xvec)) / beta(avec[1],bvec[1]) ) * p.vec

	###Calculate prior:
	if(prior=='equal'){
		prior1 <- rep(.5, length(xvec)-1)
		prior0 <- rep(.5, length(xvec)-1)
	}else if(prior=='pool'){
		prior1 <- rep(1, length(xvec)-1)
		prior0 <- rep(0, length(xvec)-1)
	}else if(prior=='opt.source'){
		#this prior identifies the optimal MEM and gives it a weight of 1

		s.mat_opt.source <- mod.mat[,-1]
		min.vec <- -log( marg.vec ) 
		if(length(xvec)>2){
			if( length(which(min.vec == min(min.vec))) > 1){
				prior1.sum <- colSums(s.mat_opt.source[ which(min.vec == min(min.vec)) ,])
				prior1.sum[which(prior1.sum != 0)] <- 1
				prior1 <- prior1.sum
			}else{
				prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ,]
			}
			prior0 <- 1 - prior1
		}else{
			prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ]
			prior0 <- 1 - prior1
		}
	}else if(prior=='opt.source_constrain'){
		#this prior identifies the optimal MEM and gives it a weight of 1

		s.mat_opt.source <- mod.mat[,-1]
		min.vec <- -log( marg.vec ) 
		if(length(xvec)>2){
			if( length(which(min.vec == min(min.vec))) > 1){
				prior1.sum <- colSums(s.mat_opt.source[ which(min.vec == min(min.vec)) ,])
				prior1.sum[which(prior1.sum != 0)] <- 1
				prior1 <- prior1.sum*constraint
			}else{
				prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ,]*constraint
			}
			prior0 <- 1 - prior1
		}else{
			prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ]*constraint
			prior0 <- 1 - prior1
		}
	}else{print('Prior not valid, please enter a valid prior option.')}

	###Calculate model weights given priors

	if(length(xvec)==2){mps <- matrix( rbind(prior0,prior1)[ paste0('prior',(mod.mat[,2:length(xvec)])), 1], ncol=1)} #extract priors for sources for each MEM
	if(length(xvec)>2){mps <- sapply(1:(length(xvec)-1), function(x) rbind(prior0,prior1)[ paste0('prior',(mod.mat[,2:length(xvec)])[,x]), x])} #extract priors for sources for each MEM

	q.vec <- marg.vec * ( rowProds(mps)/sum(rowProds(mps)) ) / sum(marg.vec * ( rowProds(mps)/sum(rowProds(mps)) ) ) #model weights, note: need to include sum(rowProds(mps)) otherwise some cases result in NaN because values are so small

	ret <- list(q = q.vec, mod.mat = mod.mat, prior = prior )

	return(ret)
}

summary_calc.betabin <- function(xvec,nvec,avec,bvec,prior,constraint=1){
###Function to calculate summary outcomes with MEM binary data with beta prior
#p: probability of the outcome for the primary cohort
#xvec: vector with counts of those with event (primary source first, then supplemental)
#nvec: vector of total number of subjects (primary source first, supplemental afterwards)
#avec: vector of alpha parameters for beta priors for each source
#bvec: vector of beta parameters for beta priors for each source
#prior: prior to use on MEM source inclusion probability: equal
#constraint: value to limit maximum value of empirical Bayes prior to

	est.mu <- est.q <- est.ehss <- NULL

	res <- calc.MEM.betabin(xvec=xvec, nvec=nvec, avec=avec, bvec=bvec, prior=prior, constraint=constraint)

	###Calculate MEM mean
	a <- avec[1] + res$mod.mat %*% xvec #Beta distribution alpha term for each model
	b <- bvec[1] + res$mod.mat %*% (nvec - xvec) #Beta distribution beta term for each model
	model.mean <- a/(a+b) #Posterior mean from Beta distribution
	mem.mean <- sum(res$q * model.mean)
	est.q <- res$q #estimated model weights, res$mod.mat for which sources are included in each model
	
	###Calculate model variances
	est.var <- (a*b)/((a+b)^2 * (a+b+1)) #Posterior variance from Beta distribution
	est.esss <- sum((BinaryESS2(M = model.mean, V = est.var) - nvec[1]) * res$q) 

	###Calculate MEM variance
	var.est <- NA #need to think on this more
	
	###Note, ESSS provided twice in results (both in $est and $esss)
	ret <- list(est=c(mem.mean, var.est, est.esss), mod.weight = est.q, mod.mat=res$mod.mat, esss=est.esss)

	return(ret)

}
 
ebola.posterior.prob <- function(na, nb, xa, xb, a, b, delta=0, mod.mat=NULL, mod.weight=NULL){
###Function to calculate the posterior probability that arm A has lower mortality than arm B
#na: vector of sample sizes for arm A
#nb: vector of sample sizes for arm B
#xa: vector with number experiencing event in arm A
#xb: vector with number experiencing event in arm B
#a: beta dist. parameter alpha for primary source
#b: beta dist. parameter beta for primary source
#delta: clinically meaningful effect size to detect, default is 0
#mod.mat: matrix which specifies which supplemental sources are included in each model (only needed for MEM, not barely Bayesian)
#mod.weight: vector of model weights calculated from MEM

	if(length(nb)>1){
		mem.post.prob <- NULL
		n.mem <- mod.mat%*%nb
		x.mem <- mod.mat%*%xb
		for(u in 1:dim(mod.mat)[1]){
			xb.f <- x.mem[u]
			nb.f <- n.mem[u]
			f.mem <- function(x) { dbeta(x, xb.f+a, nb.f+b-xb.f)*pbeta(x-delta, xa+a, na+b-xa) } #function used to calculate posterior probability by integration
			post.prob.part <- integrate(f.mem, delta, 1, rel.tol=1e-4)$value
			mem.post.prob <- c(mem.post.prob, post.prob.part)
		}
		post.prob <- sum(mem.post.prob * mod.weight) #multiple each model posterior probability by MEM weight and sum to come up with MEM probability
	}else{
		f <- function(x) { dbeta(x, xb+a, nb+b-xb)*pbeta(x-delta, xa+a, na+b-xa) } #function used to calculate posterior probability by integration
		post.prob <- integrate(f, delta, 1, rel.tol=1e-4)$value
	}

	return(post.prob)

}

ebola.sim <- function(p, arm.RR, nvec, avec=NA, bvec=NA, pp, prior=NA, seed, MEM, AR=NULL, n.burn=NULL, A.allocation=1, block.num=NULL, pi.min=1, constraint=.67, int.mont=F, fb=.01){
###Function to simulate ebola trial scenarios
###NOTE: currently, for interim analyses with AR for MEM case you need to choose combination of n.burn and block.num such that resulting block sizes are all even numbers
###NOTE: cannot use A.allocation/=1 if using interim monitoring
#p: vector of 28-day probability of mortality for the SOC group
#arm.RR: relative risk for each treatment arm with respect to SOC only
#nvec: sample size for each arm (note, with adaptive randomization the sum of the values of nvec for the arms represents the total sample size over that segment with sample size adaptively allocated based on ESSS)
#avec: beta prior alpha parameter (if none given, 1's assumed for all arms)
#avec: beta prior beta parameter (if none given, 1's assumed for all arms)
#pp: the threshold to use for the posterior probability calculations to determine if arm A is better than arm B
#prior: prior to be used for MEM estimation
#seed: number to use for simulation seed
#MEM: indicator (TRUE or FALSE) if simulation is for MEM approach or for barely Bayesian approach
#AR: indicator (TRUE) if we should implement adaptive randomization for our MEM scenario
#n.burn: number of subjects to enroll (total for both arms) before adaptively randomizing 
#A.allocation: ratio to assign treatment arm A to 1 subject in arm B for period of n.burn with default of 1:1 (i.e., A.allocation:1 for arm A:arm B)
#block.num: number of blocks to use after n.burn for patients after n.burn period (groups will be divided into approximately equal groups of subjects for a given block.num)
#pi.min: set the maximum value for proportion of patients assigned to arm A in new block with default of 1 (i.e., all patients in that block would be assigned to arm A), lower values force at least some subjects to be randomized to the control arm B
#constraint: value to limit maximum value of empirical Bayes prior ('opt.source_constrain'), default of 0.67
#int.mont: if TRUE, interim monitoring implemented (every 2 subjects for BB, end of adaptive randomization windows for AR cases, every 40 subjects without AR)
#fb: futility boundary specified on the posterior probability scale such that values below the setting indicate early termination due to futility with default of 0 (i.e., no futility boundary)

	set.seed(515+seed)

	if(is.null(AR)==TRUE){AR <- FALSE}

	if(is.na(sum(avec))){ avec <- rep(1, length(nvec))}
	if(is.na(sum(bvec))){ bvec <- rep(1, length(nvec))}

	if(length(p)==1){ p <- rep(p, length(nvec)-1)} #if only one value given for p, create vector assuming constant SOC only mortality across all segments
	if(length(pp)==1){ pp <- rep(pp, length(nvec)-1)} #if 100 subjects enrolled in arm, use this value to determine if arm A improves mortality compared to arm B

	if(AR==TRUE & MEM==FALSE){stop('Cannot adaptively randomize with Barely Bayesian approach. Either set RA=FALSE for the BB case or MEM=TRUE and RA=TRUE for adaptive randomization of MEMs.')}

	n.tracker <- NULL
	if(MEM==F){int.analysis <- c(0,seq(6,20,1),seq(40,100,20))} #vector of points at which to complete interim analyses for BB
	if(MEM==T & AR==F){int.analysis <- seq(0,100,20)}
	if(MEM==T & AR==T){
		N.segment <- nvec[2] + nvec[1] #assuming nvec is constant for each segment
		R.subjects <- N.segment - n.burn #remaining subjects left to group into blocks after burn-in period
				
		block.size_vec <- diff(round(c(0,seq(from=R.subjects/block.num,to=R.subjects,length.out=block.num)))) #vector of block sizes

		n.burn_perarm <- n.burn/2
		int_perarm <- seq(0,n.burn_perarm,20) #number of times to check interim analysis per arm during burn-in
				
		if(n.burn_perarm %in% seq(0,100,20)){int_blocks <- cumsum( c(n.burn_perarm, (block.size_vec/2)))[2:(length(block.size_vec)+1)]}else{ int_blocks <- cumsum( c(n.burn_perarm, (block.size_vec/2))) }

		int.analysis <- c(int_perarm, int_blocks)
	} 

	#Initialize objects to record if termination early due to futility or efficacy when it occurs
	term.fut <- rep(0, length(nvec)-1 )
	term.eff <- rep(0, length(nvec)-1 )

	############################
	### Code for MEM case (both adaptive randomization (AR=TRUE) and 1:1 randomization cases (AR=FALSE/NULL))
	############################
	if(MEM==TRUE){

		###Initialize arms and vector to save past segment information:
		arm.a <- 2
		arm.b <- 1
		segment.num <- 1
		p.b <- p[segment.num] #initialize the probability for the first group
		p.a <- p[segment.num]*arm.RR[segment.num] #initialize the probability for the first group
		best.trt <- 'SOC'

		sup.n <- NULL #initialize supplementary segment sample size vector
		sup.x <- NULL #initialize supplementary segment number of events

		###Create NULL vector to keep track if new drug in segment was accepted(better) or rejected(not considered an improvement)
		seg.res <- NULL #1=drug improves survival (post.prob>0.975), 0=drug doesn't show significant improvement

		###Create vector to keep track of ESSS at each segment
		esss.res <- NULL
		ar.res <- NULL
		na.res <- NULL
		nb.res <- NULL
		xa.res <- NULL
		xb.res <- NULL

		###Create matrices to keep track of treatment in which arm
		arm.mat <- matrix(nrow=length(nvec), ncol=length(nvec)-1)
		colnames(arm.mat) <- paste0('Segment ',1:(length(nvec)-1))
		rownames(arm.mat) <- c('SOC alone',paste0('SOC+drug ',1:(length(nvec)-1)))
		arm.mat[c(arm.a,arm.b), 'Segment 1'] <- c('a','b')

		for(i in 1:(length(nvec)-1)){
			if(AR==FALSE){
				if(int.mont==FALSE){
					###Generate indicator whether individual survived to 28 days and sum to calculate total number of events
					#NOTE: arm A should be the new arm being compared to the current SOC (arm B)
					n.a <- nvec[arm.a]
					n.b <- nvec[arm.b] 
					x.a <- sum(rbinom(n=n.a, size=1, prob=p.a))
					x.b <- sum(rbinom(n=n.b, size=1, prob=p.b))
				}
				if(int.mont==TRUE){
					#note, nvec not needed here since sample size determined by interim analysis counts (currently set to 100 patients per arm, but can be changed as desired in int.analysis above)
					n.a <- n.b <- x.a <- x.b <- NULL #intialize objects to store sample size and number of events for each arm

					for(im in 1:(length(int.analysis)-1) ){
 						n.a_im <- n.b_im <- int.analysis[im+1] - int.analysis[im] #calculate number to be randomized to each arm
						 
						x.a_im <- sum(rbinom(n=n.a_im, size=1, prob=p.a))
						x.b_im <- sum(rbinom(n=n.b_im, size=1, prob=p.b))

						x.a <- sum(c(x.a, x.a_im)) #add new observations to previous interim analyses
						x.b <- sum(c(x.b, x.b_im))

						n.a <- sum(c(n.a, n.a_im))
						n.b <- sum(c(n.b, n.b_im)) 

						if(i==1){
							post.prob_im <- ebola.posterior.prob(na=n.a, nb=n.b, xa=x.a, xb=x.b, a=avec[1], b=bvec[1])
						}else{
							mem.est_im <- summary_calc.betabin(xvec = c(x.b, sup.x), nvec = c(n.b, sup.n), avec=rep(avec[arm.b],length(sup.x)+1), bvec=rep(bvec[arm.b],length(sup.x)+1), prior=prior, constraint=constraint)
							post.prob_im <- ebola.posterior.prob(nb=c(n.b,sup.n), xb=c(x.b,sup.x), na=n.a, xa=x.a, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=avec[1], b=bvec[1])
						}

						if(post.prob_im >= 0.999){
							term.eff[i] <- 1
							break
						}
						if(post.prob_im < fb){
							term.fut[i] <- 1
							break
						}

					}
				}

				ar.res <- NA
			}

			if(AR==TRUE){
				#NOTE: if i=1 there is no supplemental segment information so all n.a and n.b are assigned to their respective arms
				if(i==1){ #segment 1 (no supplemental information)
					if(int.mont==F){
						n.a <- nvec[arm.a]
						n.b <- nvec[arm.b] 
						x.a <- sum(rbinom(n=n.a, size=1, prob=p.a))
						x.b <- sum(rbinom(n=n.b, size=1, prob=p.b))
						ar.res <- c(n.a/n.b)
					}
					if(int.mont==T){
						#note, nvec not needed here since sample size determined by interim analysis counts (currently set to 100 patients per arm, but can be changed as desired in int.analysis above)
						n.a <- n.b <- x.a <- x.b <- NULL #intialize objects to store sample size and number of events for each arm
	
						for(im in 1:(length(int.analysis)-1) ){
 							n.a_im <- n.b_im <- int.analysis[im+1] - int.analysis[im] #calculate number to be randomized to each arm
						 
							x.a_im <- sum(rbinom(n=n.a_im, size=1, prob=p.a))
							x.b_im <- sum(rbinom(n=n.b_im, size=1, prob=p.b))

							x.a <- sum(c(x.a, x.a_im)) #add new observations to previous interim analyses
							x.b <- sum(c(x.b, x.b_im))

							n.a <- sum(c(n.a, n.a_im))
							n.b <- sum(c(n.b, n.b_im)) 

							post.prob_im <- ebola.posterior.prob(na=n.a, nb=n.b, xa=x.a, xb=x.b, a=avec[1], b=bvec[1])
							if(post.prob_im >= 0.999){
								term.eff[i] <- 1
								break
							}
							if(post.prob_im < fb){
								term.fut[i] <- 1
								break
							}

						}
					}
				}else{ #segments 2 onward below
					N.segment <- nvec[arm.a] + nvec[arm.b] #total sample size for segment
					R.subjects <- N.segment - n.burn #remaining subjects left to group into blocks after burn-in period
					R.blck <- R.subjects #remaining subjects at start of block
				
					block.size_vec <- diff(round(c(0,seq(from=R.subjects/block.num,to=R.subjects,length.out=block.num)))) #vector of block sizes

					if(int.mont==F){
						#Set sample sizes for each arm based on allocation during burn-in period
						n.a <- round(n.burn * A.allocation/(A.allocation+1)) 
						n.b <- n.burn - n.a
	
						x.a <- sum(rbinom(n=n.a, size=1, prob=p.a))
						x.b <- sum(rbinom(n=n.b, size=1, prob=p.b))

						for(blck in 1:length(block.size_vec)){

							esss.blck <- summary_calc.betabin(xvec = c(x.b, sup.x), nvec = c(n.b, sup.n), avec=rep(avec[arm.b],length(sup.x)+1), bvec=rep(bvec[arm.b],length(sup.x)+1), prior=prior, constraint=constraint)$esss
							pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
							pi.a <- max(0, min(c(pi.est, pi.min)))

							n.a_blck <- round(block.size_vec[blck]*pi.a)
							n.b_blck <- block.size_vec[blck] - n.a_blck
							R.blck <- R.blck - block.size_vec[blck]

							if(n.a_blck != 0){x.a <- x.a + sum(rbinom(n=n.a_blck, size=1, prob=p.a))}else{x.a <- x.a}
							if(n.b_blck != 0){x.b <- x.b + sum(rbinom(n=n.b_blck, size=1, prob=p.b))}else{x.b <- x.b}

							n.a <- n.a + n.a_blck
							n.b <- n.b + n.b_blck

						}
					}

					if(int.mont==T){
						#note, nvec not needed here since sample size determined by interim analysis counts (currently set to 100 patients per arm, but can be changed as desired in int.analysis above)
						n.a <- n.b <- x.a <- x.b <- NULL #intialize objects to store sample size and number of events for each arm
						ia_ar.burnin <- int.analysis[which(int.analysis <= n.burn/2)]
	
						for(im in 1:(length(ia_ar.burnin)-1) ){
 							n.a_im <- n.b_im <- int.analysis[im+1] - int.analysis[im] #calculate number to be randomized to each arm
						 
							x.a_im <- sum(rbinom(n=n.a_im, size=1, prob=p.a))
							x.b_im <- sum(rbinom(n=n.b_im, size=1, prob=p.b))

							x.a <- sum(c(x.a, x.a_im)) #add new observations to previous interim analyses
							x.b <- sum(c(x.b, x.b_im))

							n.a <- sum(c(n.a, n.a_im))
							n.b <- sum(c(n.b, n.b_im)) 

							post.prob_im <- ebola.posterior.prob(na=n.a, nb=n.b, xa=x.a, xb=x.b, a=avec[1], b=bvec[1])
							if(post.prob_im >= 0.999){
								term.eff[i] <- 1
								break
							}
							if(post.prob_im < fb){
								term.fut[i] <- 1
								break
							}
						}

						if(post.prob_im < 0.999 & post.prob_im >= fb){
							for(blck in 1:length(block.size_vec)){

								esss.blck <- summary_calc.betabin(xvec = c(x.b, sup.x), nvec = c(n.b, sup.n), avec=rep(avec[arm.b],length(sup.x)+1), bvec=rep(bvec[arm.b],length(sup.x)+1), prior=prior, constraint=constraint)$esss
								pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
								pi.a <- max(0, min(c(pi.est, pi.min)))

								n.a_blck <- round(block.size_vec[blck]*pi.a)
								n.b_blck <- block.size_vec[blck] - n.a_blck
								R.blck <- R.blck - block.size_vec[blck]

								if(n.a_blck != 0){x.a <- x.a + sum(rbinom(n=n.a_blck, size=1, prob=p.a))}else{x.a <- x.a}
								if(n.b_blck != 0){x.b <- x.b + sum(rbinom(n=n.b_blck, size=1, prob=p.b))}else{x.b <- x.b}

								n.a <- n.a + n.a_blck
								n.b <- n.b + n.b_blck

								mem.est_im <- summary_calc.betabin(xvec = c(x.b, sup.x), nvec = c(n.b, sup.n), avec=rep(avec[arm.b],length(sup.x)+1), bvec=rep(bvec[arm.b],length(sup.x)+1), prior=prior, constraint=constraint)
								post.prob_im <- ebola.posterior.prob(nb=c(n.b,sup.n), xb=c(x.b,sup.x), na=n.a, xa=x.a, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=avec[1], b=bvec[1])

								if(post.prob_im >= 0.999){
									term.eff[i] <- 1
									break
								}
								if(post.prob_im < fb){
									term.fut[i] <- 1
									break
								}

							}
						}
					}

					ar.res <- c(ar.res, n.a/n.b)

				}
			}

			###Calculate MEM estimate if at least one supplemental source and posterior probability
			if(length(sup.n) >= 1){
				mem.est <- summary_calc.betabin(xvec = c(x.b, sup.x), nvec = c(n.b, sup.n), avec=rep(avec[arm.b],length(sup.x)+1), bvec=rep(bvec[arm.b],length(sup.x)+1), prior=prior, constraint=constraint)
				esss.res <- c(esss.res, mem.est$esss)
				na.res <- c(na.res, n.a)
				nb.res <- c(nb.res, n.b)
				xa.res <- c(xa.res, x.a)
				xb.res <- c(xb.res, x.b)
				post.prob <- ebola.posterior.prob(nb=c(n.b,sup.n), xb=c(x.b,sup.x), na=n.a, xa=x.a, mod.mat=mem.est$mod.mat, mod.weight=mem.est$mod.weight, a=avec[1], b=bvec[1])
			}else{
				esss.res <- c(esss.res, 0) #no supplemental sources to borrow from, ESSS must be 0
				na.res <- c(na.res, n.a)
				nb.res <- c(nb.res, n.b)
				xa.res <- c(xa.res, x.a)
				xb.res <- c(xb.res, x.b)
				post.prob <- ebola.posterior.prob(na=n.a, nb=n.b, xa=x.a, xb=x.b, a=avec[1], b=bvec[1]) #calculate posterior probability with new arm A and arm B with potential supplementary information
			}

			###Update the supplementary segment vectors, arm specification, etc. if not at last segment
			if(i < (length(nvec)-1)){
				if(post.prob < pp[i]){
					sup.n <- c(sup.n, n.b); sup.x <- c(sup.x, x.b); arm.a <- arm.a+1; arm.mat[c(arm.a,arm.b),(i+1)] <- c('a','b')
					seg.res <- c(seg.res, 0)
					segment.num <- segment.num + 1
					p.b <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,0))
					p.a <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,1)) #updated mortality for new drug on arm A
				}else{
					best.trt <- paste0(best.trt, '+Drug ', arm.a-1)
					#NOTE: arm.b_RR updated to multiply current p by RR of new group (each subsequent group alters this RR), works because trial design stacks treatments as it proceeds
					segment.num <- segment.num + 1
					seg.res <- c(seg.res, 1)
					sup.n <- c(n.a); sup.x <- c(x.a); arm.b <- arm.a; arm.a <- arm.a+1; arm.mat[c(arm.a,arm.b),(i+1)] <- c('a','b')
					p.b <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,0)) #arm B's mortality becomes the estimate for arm A's mortality in next segment
					p.a <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,1)) #arm A's new mortality is arm B's mortality times new drug's RR
				}	
			}
		}

		if(post.prob >= pp[i]){	
			best.trt <- paste0(best.trt, '+Drug ', arm.a-1) 
			seg.res <- c(seg.res, 1)
			p.final <- p.a
		}else{
			seg.res <- c(seg.res, 0)
			p.final <- p.b
		}
	}

	############################
	### Code for PREVAIL II design (no incorporation of supplemental information)
	############################
	if(MEM==FALSE){

		###Initialize arms and vector to save past segment information:
		arm.a <- 2
		arm.b <- 1
		segment.num <- 1
		p.b <- p[segment.num] #initialize the probability for the first group
		p.a <- p[segment.num]*arm.RR[segment.num] #initialize the probability for the first group
		best.trt <- 'SOC'

		###Create NULL vector to keep track if new drug in segment was accepted(better) or rejected(not considered an improvement)
		seg.res <- NULL #let 1=drug improves survival (post.prob>0.975), 0=drug doesn't show significant improvement
		na.res <- NULL
		nb.res <- NULL
		xa.res <- NULL
		xb.res <- NULL
		esss.res <- NULL

		###Create matrices to keep track of treatment in which arm
		arm.mat <- matrix(nrow=length(nvec), ncol=length(nvec)-1)
		colnames(arm.mat) <- paste0('Segment ',1:(length(nvec)-1))
		rownames(arm.mat) <- c('SOC alone',paste0('SOC+drug ',1:(length(nvec)-1)))
		arm.mat[c(arm.a,arm.b), 'Segment 1'] <- c('a','b')

		for(i in 1:(length(nvec)-1)){
			if(int.mont==FALSE){
				###Generate indicator whether individual survived to 28 days and sum to calculate total number of events
				#NOTE: arm A should be the new arm being compared to the current SOC (arm B)
				n.a <- nvec[arm.a]
				n.b <- nvec[arm.b] 
				x.a <- sum(rbinom(n=n.a, size=1, prob=p.a))
				x.b <- sum(rbinom(n=n.b, size=1, prob=p.b))
			}
			if(int.mont==TRUE){
				#note, nvec not needed here since sample size determined by interim analysis counts (currently set to 100 patients per arm, but can be changed as desired in int.analysis above)
				n.a <- n.b <- x.a <- x.b <- NULL #intialize objects to store sample size and number of events for each arm

				for(im in 1:(length(int.analysis)-1) ){
					n.a_im <- n.b_im <- int.analysis[im+1] - int.analysis[im] #calculate number to be randomized to each arm
						 
					x.a_im <- sum(rbinom(n=n.a_im, size=1, prob=p.a))
					x.b_im <- sum(rbinom(n=n.b_im, size=1, prob=p.b))

					x.a <- sum(c(x.a, x.a_im)) #add new observations to previous interim analyses
					x.b <- sum(c(x.b, x.b_im))

					n.a <- sum(c(n.a, n.a_im))
					n.b <- sum(c(n.b, n.b_im)) 

					post.prob_im <- ebola.posterior.prob(na=n.a, nb=n.b, xa=x.a, xb=x.b, a=avec[1], b=bvec[1])
					if(post.prob_im >= 0.999){
						term.eff[i] <- 1
						break
					}
					if(post.prob_im < fb){
						term.fut[i] <- 1
						break
					}
				}
			}

			###Calculate the posterior probability
			na.res <- c(na.res, n.a)
			nb.res <- c(nb.res, n.b)
			xa.res <- c(xa.res, x.a)
			xb.res <- c(xb.res, x.b)
			post.prob <- ebola.posterior.prob(na=n.a, nb=n.b, xa=x.a, xb=x.b, a=avec[1], b=bvec[1]) #calculate posterior probability with new arm A and arm B with potential supplementary information

			###Update the supplementary segment vectors, arm specification, etc. if not at last segment
			if(i < (length(nvec)-1)){
				if(post.prob < pp[i]){
					arm.a <- arm.a+1; arm.mat[c(arm.a,arm.b),(i+1)] <- c('a','b')
					seg.res <- c(seg.res, 0)
					segment.num <- segment.num + 1
					p.b <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,0))
					p.a <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,1)) #updated mortality for new drug on arm A
				}else{
					best.trt <- paste0(best.trt, '+Drug ', arm.a-1)
					#NOTE: arm.b_RR updated to multiply current p by RR of new group (each subsequent group alters this RR), works because trial design stacks treatments as it proceeds
					segment.num <- segment.num + 1
					seg.res <- c(seg.res, 1)
					arm.b <- arm.a; arm.a <- arm.a+1; arm.mat[c(arm.a,arm.b),(i+1)] <- c('a','b')
					p.b <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,0)) #arm B's mortality becomes the estimate for arm A's mortality in next segment
					p.a <- p[segment.num] * prod(arm.RR[1:segment.num]^c(seg.res,1)) #arm A's new mortality is arm B's mortality times new drug's RR
				}	
			}
		}

		if(post.prob >= pp[i]){	
			best.trt <- paste0(best.trt, '+Drug ', arm.a-1) 
			seg.res <- c(seg.res, 1)
			p.final <- p.a
		}else{
			seg.res <- c(seg.res, 0)
			p.final <- p.b
		}

		esss.res <- rep(0, length(nvec)-1) #BB case has ESSS of 0 for all segments
		ar.res <- NA

	}

	###Values to return
	trt.n <- sum(na.res) #number of people over entire simulation who are randomized to treatment arm
	overall.n <- sum(na.res) + sum(nb.res) #total sample size across all segments
	overall.n25 <- sum(na.res[2:5]) + sum(nb.res[2:5]) #total sample size for segments 2-5
	overall.death <- sum(xa.res) + sum(xb.res) #total number of deaths observed in both arms
	seg.surv_a <- (na.res - xa.res)/na.res #percent surviving by segment in treatment arm
	seg.surv_b <- (nb.res - xb.res)/nb.res #percent surviving by segment in control arm
	seg.surv <- (na.res+nb.res - (xa.res+xb.res)) / (na.res + nb.res) #percent surviving overall by segment
	seg.surv_25a <- (sum(na.res[2:5]) - sum(xa.res[2:5]) )/sum(na.res[2:5]) #proportion surviving in treatment from segments 2-5
	seg.surv_25b <- (sum(nb.res[2:5]) - sum(xb.res[2:5]) )/sum(nb.res[2:5]) #proportion surviving in control from segments 2-5
	seg.surv_25 <- ( sum(na.res[2:5])+sum(nb.res[2:5]) - ( sum(xa.res[2:5])+sum(xb.res[2:5]) ) ) / overall.n25 #percent surviving overall fom segments 2-5

	seg_N_props <- c(seg.res, overall.n, sum(na.res[2:5])/overall.n25, seg.surv_a, seg.surv_25a, seg.surv_b, seg.surv_25b, seg.surv, seg.surv_25, (na.res+nb.res), term.eff, term.fut )
	ret <- list(mem = MEM, best.trt = best.trt, arm.mat = arm.mat, final.trt.p = p.final, seg.res = seg.res, seg_ntrt_res = c(seg.res, trt.n), seg_n_res = c(seg.res, overall.n), esss.res = esss.res, ar.res = ar.res, nb.res = nb.res, na.res = na.res, xb.res = xb.res, xa.res = xa.res, seg_N_props = seg_N_props, term.eff = term.eff, term.fut = term.fut)

}

tab.segN <- function(vn,v2,v3,v4,v5,cn,c2,c3,c4,c5){
###Function to take matrices of results for 5 scenarios (varying and constant), then create summary table
#Note: need power/type-1 error in first 5 columns, overall N in 6th column, proportion randomized to treatment arm in 7th column, proportion surviving in non-null segment (or 2-5 for null case) in 8th column
	asdf <- rbind( c( format(round(colMeans(vn)[1:5],3), nsmall=3), paste0(round( mean(vn[,6]),0)," (",round(sd(vn[,6]),2),")"), paste0(round( mean(vn[,7]), 3)," (",round(sd(vn[,7]),3),")"), paste0(round( mean(vn[,25]), 3)," (",round(sd(vn[,25]),3),")")),
		c( format(round(colMeans(v2)[1:5],3), nsmall=3), paste0(round( mean(v2[,6]),0)," (",round(sd(v2[,6]),2),")"), paste0(round( mean(v2[,7]), 3)," (",round(sd(v2[,7]),3),")"), paste0(round( mean(v2[,21]), 3)," (",round(sd(v2[,21]),3),")")),
		c( format(round(colMeans(v3)[1:5],3), nsmall=3), paste0(round( mean(v3[,6]),0)," (",round(sd(v3[,6]),2),")"), paste0(round( mean(v3[,7]), 3)," (",round(sd(v3[,7]),3),")"), paste0(round( mean(v3[,22]), 3)," (",round(sd(v3[,22]),3),")")),
		c( format(round(colMeans(v4)[1:5],3), nsmall=3), paste0(round( mean(v4[,6]),0)," (",round(sd(v4[,6]),2),")"), paste0(round( mean(v4[,7]), 3)," (",round(sd(v4[,7]),3),")"), paste0(round( mean(v4[,23]), 3)," (",round(sd(v4[,23]),3),")")),
		c( format(round(colMeans(v5)[1:5],3), nsmall=3), paste0(round( mean(v5[,6]),0)," (",round(sd(v5[,6]),2),")"), paste0(round( mean(v5[,7]), 3)," (",round(sd(v5[,7]),3),")"), paste0(round( mean(v5[,24]), 3)," (",round(sd(v5[,24]),3),")")),
		c( format(round(colMeans(cn)[1:5],3), nsmall=3), paste0(round( mean(cn[,6]),0)," (",round(sd(cn[,6]),2),")"), paste0(round( mean(cn[,7]), 3)," (",round(sd(cn[,7]),3),")"), paste0(round( mean(cn[,25]), 3)," (",round(sd(cn[,25]),3),")")), 
		c( format(round(colMeans(c2)[1:5],3), nsmall=3), paste0(round( mean(c2[,6]),0)," (",round(sd(c2[,6]),2),")"), paste0(round( mean(c2[,7]), 3)," (",round(sd(c2[,7]),3),")"), paste0(round( mean(c2[,21]), 3)," (",round(sd(c2[,21]),3),")")),
		c( format(round(colMeans(c3)[1:5],3), nsmall=3), paste0(round( mean(c3[,6]),0)," (",round(sd(c3[,6]),2),")"), paste0(round( mean(c3[,7]), 3)," (",round(sd(c3[,7]),3),")"), paste0(round( mean(c3[,22]), 3)," (",round(sd(c3[,22]),3),")")),
		c( format(round(colMeans(c4)[1:5],3), nsmall=3), paste0(round( mean(c4[,6]),0)," (",round(sd(c4[,6]),2),")"), paste0(round( mean(c4[,7]), 3)," (",round(sd(c4[,7]),3),")"), paste0(round( mean(c4[,23]), 3)," (",round(sd(c4[,23]),3),")")),
		c( format(round(colMeans(c5)[1:5],3), nsmall=3), paste0(round( mean(c5[,6]),0)," (",round(sd(c5[,6]),2),")"), paste0(round( mean(c5[,7]), 3)," (",round(sd(c5[,7]),3),")"), paste0(round( mean(c5[,24]), 3)," (",round(sd(c5[,24]),3),")")))
	return(asdf)
}

tab.segN_twotrt <- function(vn,v2,v3,v4,v5,v6,v7,cn,c2,c3,c4,c5,c6,c7){
###Function to take matrices of results for 5 scenarios (varying and constant), then create summary table
#Note: need power/type-1 error in first 5 columns, overall N in 6th column, proportion randomized to treatment arm in 7th column, proportion surviving in non-null segment (or 2-5 for null case) in 8th column
	###Calculate weighted average of proportion surviving the non-null segments
	v2s <- (v2[,21]*v2[,27] + v2[,22]*v2[,28]) / (v2[,27] + v2[,28]) #S2/S3
	v3s <- (v2[,21]*v2[,27] + v2[,23]*v2[,29]) / (v2[,27] + v2[,29]) #S2/S4
	v4s <- (v2[,21]*v2[,27] + v2[,24]*v2[,30]) / (v2[,27] + v2[,30]) #S2/S5
	v5s <- (v2[,22]*v2[,28] + v2[,23]*v2[,29]) / (v2[,28] + v2[,29]) #S3/S4
	v6s <- (v2[,22]*v2[,28] + v2[,24]*v2[,30]) / (v2[,28] + v2[,30]) #S3/S5
	v7s <- (v2[,23]*v2[,29] + v2[,24]*v2[,30]) / (v2[,29] + v2[,30]) #S4/S5

	c2s <- (c2[,21]*c2[,27] + c2[,22]*c2[,28]) / (c2[,27] + c2[,28]) #S2/S3
	c3s <- (c2[,21]*c2[,27] + c2[,23]*c2[,29]) / (c2[,27] + c2[,29]) #S2/S4
	c4s <- (c2[,21]*c2[,27] + c2[,24]*c2[,30]) / (c2[,27] + c2[,30]) #S2/S5
	c5s <- (c2[,22]*c2[,28] + c2[,23]*c2[,29]) / (c2[,28] + c2[,29]) #S3/S4
	c6s <- (c2[,22]*c2[,28] + c2[,24]*c2[,30]) / (c2[,28] + c2[,30]) #S3/S5
	c7s <- (c2[,23]*c2[,29] + c2[,24]*c2[,30]) / (c2[,29] + c2[,30]) #S4/S5

	asdf <- rbind( c( format(round(colMeans(vn)[1:5],3), nsmall=3), paste0(round( mean(vn[,6]),0)," (",round(sd(vn[,6]),2),")"), paste0(round( mean(vn[,7]), 3)," (",round(sd(vn[,7]),3),")"), paste0(round( mean(vn[,25]), 3)," (",round(sd(vn[,25]),3),")")),
		c( format(round(colMeans(v2)[1:5],3), nsmall=3), paste0(round( mean(v2[,6]),0)," (",round(sd(v2[,6]),2),")"), paste0(round( mean(v2[,7]), 3)," (",round(sd(v2[,7]),3),")"), paste0(round( mean(v2s), 3)," (",round(sd(v2s),3),")")),
		c( format(round(colMeans(v3)[1:5],3), nsmall=3), paste0(round( mean(v3[,6]),0)," (",round(sd(v3[,6]),2),")"), paste0(round( mean(v3[,7]), 3)," (",round(sd(v3[,7]),3),")"), paste0(round( mean(v3s), 3)," (",round(sd(v3s),3),")")),
		c( format(round(colMeans(v4)[1:5],3), nsmall=3), paste0(round( mean(v4[,6]),0)," (",round(sd(v4[,6]),2),")"), paste0(round( mean(v4[,7]), 3)," (",round(sd(v4[,7]),3),")"), paste0(round( mean(v4s), 3)," (",round(sd(v4s),3),")")),
		c( format(round(colMeans(v5)[1:5],3), nsmall=3), paste0(round( mean(v5[,6]),0)," (",round(sd(v5[,6]),2),")"), paste0(round( mean(v5[,7]), 3)," (",round(sd(v5[,7]),3),")"), paste0(round( mean(v5s), 3)," (",round(sd(v5s),3),")")),
		c( format(round(colMeans(v6)[1:5],3), nsmall=3), paste0(round( mean(v6[,6]),0)," (",round(sd(v6[,6]),2),")"), paste0(round( mean(v6[,7]), 3)," (",round(sd(v6[,7]),3),")"), paste0(round( mean(v6s), 3)," (",round(sd(v6s),3),")")),
		c( format(round(colMeans(v7)[1:5],3), nsmall=3), paste0(round( mean(v7[,6]),0)," (",round(sd(v7[,6]),2),")"), paste0(round( mean(v7[,7]), 3)," (",round(sd(v7[,7]),3),")"), paste0(round( mean(v7s), 3)," (",round(sd(v7s),3),")")),
		c( format(round(colMeans(cn)[1:5],3), nsmall=3), paste0(round( mean(cn[,6]),0)," (",round(sd(cn[,6]),2),")"), paste0(round( mean(cn[,7]), 3)," (",round(sd(cn[,7]),3),")"), paste0(round( mean(cn[,25]), 3)," (",round(sd(cn[,25]),3),")")), 
		c( format(round(colMeans(c2)[1:5],3), nsmall=3), paste0(round( mean(c2[,6]),0)," (",round(sd(c2[,6]),2),")"), paste0(round( mean(c2[,7]), 3)," (",round(sd(c2[,7]),3),")"), paste0(round( mean(c2s), 3)," (",round(sd(c2s),3),")")),
		c( format(round(colMeans(c3)[1:5],3), nsmall=3), paste0(round( mean(c3[,6]),0)," (",round(sd(c3[,6]),2),")"), paste0(round( mean(c3[,7]), 3)," (",round(sd(c3[,7]),3),")"), paste0(round( mean(c3s), 3)," (",round(sd(c3s),3),")")),
		c( format(round(colMeans(c4)[1:5],3), nsmall=3), paste0(round( mean(c4[,6]),0)," (",round(sd(c4[,6]),2),")"), paste0(round( mean(c4[,7]), 3)," (",round(sd(c4[,7]),3),")"), paste0(round( mean(c4s), 3)," (",round(sd(c4s),3),")")),
		c( format(round(colMeans(c5)[1:5],3), nsmall=3), paste0(round( mean(c5[,6]),0)," (",round(sd(c5[,6]),2),")"), paste0(round( mean(c5[,7]), 3)," (",round(sd(c5[,7]),3),")"), paste0(round( mean(c5s), 3)," (",round(sd(c5s),3),")")),
		c( format(round(colMeans(c6)[1:5],3), nsmall=3), paste0(round( mean(c6[,6]),0)," (",round(sd(c6[,6]),2),")"), paste0(round( mean(c6[,7]), 3)," (",round(sd(c6[,7]),3),")"), paste0(round( mean(c6s), 3)," (",round(sd(c6s),3),")")),
		c( format(round(colMeans(c7)[1:5],3), nsmall=3), paste0(round( mean(c7[,6]),0)," (",round(sd(c7[,6]),2),")"), paste0(round( mean(c7[,7]), 3)," (",round(sd(c7[,7]),3),")"), paste0(round( mean(c7s), 3)," (",round(sd(c7s),3),")")))
	return(asdf)
}

meansd_fnx <- function(val){
###Function to take matrix of 6 values and return mean (sd) for table
	ma <- format(round(colMeans(val),3),nsmall=3)
	sa <- format(round(colSds(val),3),nsmall=3)

	ret <- c( paste0(ma[1], " (",sa[1],")"), paste0(ma[2], " (",sa[2],")"), paste0(ma[3], " (",sa[3],")"), paste0(ma[4], " (",sa[4],")"), paste0(ma[5], " (",sa[5],")"), paste0(ma[6], " (",sa[6],")"))
	return(ret)
}

tab.segN_propsurv <- function(vn,v2,v3,v4,v5,cn,c2,c3,c4,c5){
###Function to take matrices of results for 5 scenarios (varying and constant), then create summary table
#Note: columns 1-5 represent each segment, column 6 represents the mean (sd) from segments 2-5
#Note: rows 1-10 represent arm A, 11-20 represent arm B (control), 21-30 represent survivors for both arms combined
	asdf <- rbind( meansd_fnx(val=vn[,8:13]),meansd_fnx(val=v2[,8:13]),meansd_fnx(val=v3[,8:13]),meansd_fnx(val=v4[,8:13]),meansd_fnx(val=v5[,8:13]),
		meansd_fnx(val=cn[,8:13]),meansd_fnx(val=c2[,8:13]),meansd_fnx(val=c3[,8:13]),meansd_fnx(val=c4[,8:13]),meansd_fnx(val=c5[,8:13]),
		meansd_fnx(val=vn[,14:19]),meansd_fnx(val=v2[,14:19]),meansd_fnx(val=v3[,14:19]),meansd_fnx(val=v4[,14:19]),meansd_fnx(val=v5[,14:19]),
		meansd_fnx(val=cn[,14:19]),meansd_fnx(val=c2[,14:19]),meansd_fnx(val=c3[,14:19]),meansd_fnx(val=c4[,14:19]),meansd_fnx(val=c5[,14:19]),
		meansd_fnx(val=vn[,20:25]),meansd_fnx(val=v2[,20:25]),meansd_fnx(val=v3[,20:25]),meansd_fnx(val=v4[,20:25]),meansd_fnx(val=v5[,20:25]),
		meansd_fnx(val=cn[,20:25]),meansd_fnx(val=c2[,20:25]),meansd_fnx(val=c3[,20:25]),meansd_fnx(val=c4[,20:25]),meansd_fnx(val=c5[,20:25]))
	return(asdf)
}

mat.segN_nonnull <- function(sim.num, vn,v2,v3,v4,v5,cn,c2,c3,c4,c5){
###Function to take matrices of results for 5 scenarios (varying and constant), take relevant output for use in making figures (boxplots, etc.)
#columns 1-5 are N per scenario (null/2/3/4/5), 6-10 are proportion assigned to treatment in segments N/2/3/4/5, 11-15 are proportion surviving segments N/2/3/4/5 in given scenario

	vary.mat <- cbind(sim.num, vn[,6], v2[,6], v3[,6], v4[,6], v5[,6], vn[,7], v2[,7], v3[,7], v4[,7], v5[,7], vn[,25], v2[,21], v3[,22], v4[,23], v5[,24])
	constant.mat <- cbind(sim.num, cn[,6], c2[,6], c3[,6], c4[,6], c5[,6], cn[,7], c2[,7], c3[,7], c4[,7], c5[,7], cn[,25], c2[,21], c3[,22], c4[,23], c5[,24])

	write.table(vary.mat, file ='paper2_simresults_varying_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
	write.table(constant.mat, file ='paper2_simresults_constant_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

}

mat.segN_nonnull_futility <- function(sim.num, vn,v2,v3,v4,v5,cn,c2,c3,c4,c5){
###Function to take matrices of results for 5 scenarios (varying and constant), take relevant output for use in making figures (boxplots, etc.)
#columns 1-5 are N per scenario (null/2/3/4/5), 6-10 are proportion assigned to treatment in segments N/2/3/4/5, 11-15 are proportion surviving segments N/2/3/4/5 in given scenario

	vary.mat <- cbind(sim.num, vn[,6], v2[,6], v3[,6], v4[,6], v5[,6], vn[,7], v2[,7], v3[,7], v4[,7], v5[,7], vn[,25], v2[,21], v3[,22], v4[,23], v5[,24])
	constant.mat <- cbind(sim.num, cn[,6], c2[,6], c3[,6], c4[,6], c5[,6], cn[,7], c2[,7], c3[,7], c4[,7], c5[,7], cn[,25], c2[,21], c3[,22], c4[,23], c5[,24])

	write.table(vary.mat, file ='paper2_simresults_futility_varying_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
	write.table(constant.mat, file ='paper2_simresults_futility_constant_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

}

mat.segN_propearly <- function(sim.num, vn,v2,v3,v4,v5,cn,c2,c3,c4,c5){
###Function to take matrices of results for 5 scenarios (varying and constant), take relevant output for use in making figures (boxplots, etc.)
#columns 1-5 are N per scenario (null/2/3/4/5), 6-10 are proportion assigned to treatment in segments N/2/3/4/5, 11-15 are proportion surviving segments N/2/3/4/5 in given scenario

	vary.mat <- cbind(sim.num, vn[,31:40], v2[,31:40], v3[,31:40], v4[,31:40], v5[,31:40])
	constant.mat <- cbind(sim.num, cn[,31:40], c2[,31:40], c3[,31:40], c4[,31:40], c5[,31:40])

	write.table(vary.mat, file ='paper2_simresults_futility_props_varying_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
	write.table(constant.mat, file ='paper2_simresults_futility_props_constant_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

}

server.sim <- function(sim.num){
###Function to implement simulations on server with parallel programming
#sim.num: combination of elements to simulate (ranges from 1 to 16)

	p.vary <- c(.74,.61,.48,.36,.23)

	type <- c( rep('EB10', 2), rep('pie',2), 'BB', rep('Pooled',2), rep('EB10', 3), 'EB20', rep('EB10',2), 'BB', rep('EB10',2))[sim.num]
	pp.vec <- list( c(.975,.97125,.96625,.95875,.9575), c(.975,.975,.9775,.975,.975),
		c(.975,.96375,.95875,.94375,.9325), c(.975,.9815,.985,.9875,.9875), rep(.975,5), c(.975,.97375,.96375,.95375,.94),
		c(.975,.995,.999,.9975,.999), c(.975,.97125,.96625,.95875,.9575), c(.975,.97125,.96625,.95875,.9575), 
		c(.975,.97125,.96625,.95875,.9575), c(.975,.96375,.96125,.95625,.95), c(.975,.97125,.96625,.95875,.9575), 
		c(.975,.975,.9775,.975,.975), rep(.975,5),
		c(.975,.97125,.96625,.95875,.9575), c(.975,.975,.9775,.975,.975) )[[sim.num]]

	rr <- c( rep(.7, 11), rep(.5,3), rep(.7, 2) )[sim.num]
	mem.ind <- c(rep(TRUE,4), FALSE, rep(TRUE,8), FALSE, rep(TRUE,2) )[sim.num]
	prior <- c( rep('opt.source_constrain',2), rep('equal',2), NA, rep('pool',2), rep('opt.source_constrain', 6), NA, rep('opt.source_constrain', 2) )[sim.num]
	constraint.val <- c(.1,.1,rep(NA,5), rep(.1,3), .2, .1, .1, NA, .1, .1)[sim.num]
	AR.ind <- c( rep(TRUE, 4),FALSE,rep(TRUE,4),  FALSE, rep(TRUE, 3), FALSE, rep(TRUE,2) )[sim.num]
	n.burn <- c( rep(60, 4), NA, rep(60,2), 30, 90, NA, rep(60,3), NA, rep(60,2) )[sim.num]
	block.num <- c( rep(5, 4), NA, rep(5,2), 6, 4, NA, rep(5, 3), NA, rep(5,2) )[sim.num]
	int.ind <- c( rep(TRUE,14), rep(FALSE,2) )[sim.num] #indicator for if interim monitoring should be implemented

	vn <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=rep(1,5), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v2 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,rr,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v3 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,rr,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v4 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v5 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))

	cn <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c2 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,rr,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c3 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,rr,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c4 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c5 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))

	cap.use <- paste0(type,': Int Analyses=',int.ind,', AR=',AR.ind,', nburn=', n.burn,', blocks=',block.num,', RR=',rr,', PP thresholds: (',pp.vec[1],',',pp.vec[2],',',pp.vec[3],',',pp.vec[4],',',pp.vec[5],')')

	asd <- print(xtable(tab.segN(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5), caption=cap.use))
	asd_prop <- print(xtable(tab.segN_propsurv(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5), caption=cap.use))
	asd_temp <- mat.segN_nonnull(sim.num, vn,v2,v3,v4,v5,cn,c2,c3,c4,c5)

	write.table(asd, file ='paper2_simresults_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
	write.table(asd_prop, file ='paper2_propsurv_results_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

}

server.sim_futility <- function(sim.num){
###Function to implement simulations on server with parallel programming which incorporate futility boundaries
#sim.num: combination of elements to simulate (ranges from 1 to 16 as of now)

	p.vary <- c(.74,.61,.48,.36,.23)

	type <- rep(c( 'EB10','pie','BB','Pooled'),4)[sim.num]
	pp.vec <- rep(list( c(.975,.97125,.96625,.95875,.9575), c(.975,.96375,.95875,.94375,.9325), rep(.975,5), c(.975,.97375,.96375,.95375,.94) ),4)[[sim.num]]

	rr <- .7
	mem.ind <- rep(c(T,T,F,T),4)[sim.num]
	prior <- rep(c( rep('opt.source_constrain',1), rep('equal',1), NA, rep('pool',1) ),4)[sim.num]
	constraint.val <- rep(c(.1,rep(NA,3)),4)[sim.num]
	AR.ind <- rep(c( rep(TRUE, 2),FALSE, TRUE ),4)[sim.num]
	n.burn <- rep(c( rep(60, 2), NA, rep(60,1) ),4)[sim.num]
	block.num <- rep(c( rep(5, 2), NA, rep(5,1) ),4)[sim.num]
	fb <- c( rep(.1,4), rep(.01,4), rep(.001,4), rep(.2,4) )[sim.num]
	int.ind <- TRUE #indicator for if interim monitoring should be implemented

	vn <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=rep(1,5), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v2 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,rr,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v3 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,rr,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v4 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v5 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))

	cn <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c2 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,rr,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c3 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,rr,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c4 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c5 <- t(sapply(1:25000, function(x) ebola.sim(fb=fb, int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))

	cap.use <- paste0(type,': Int Analyses=',int.ind,', AR=',AR.ind,', nburn=', n.burn,', blocks=',block.num,', RR=',rr,', PP thresholds: (',pp.vec[1],',',pp.vec[2],',',pp.vec[3],',',pp.vec[4],',',pp.vec[5],'), Futility Bound=', fb)

	asd_temp <- mat.segN_nonnull_futility(sim.num, vn,v2,v3,v4,v5,cn,c2,c3,c4,c5)
	asd_temp2 <- mat.segN_propearly(sim.num, vn,v2,v3,v4,v5,cn,c2,c3,c4,c5)

	asd <- print(xtable(tab.segN(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5), caption=cap.use))
	asd_prop <- print(xtable(tab.segN_propsurv(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5), caption=cap.use))

	write.table(asd, file ='paper2_simresults_futility_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
	write.table(asd_prop, file ='paper2_propsurv_results_futility_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

}


server.sim_2nonnull <- function(sim.num){
###Function to implement simulations on server with parallel programming for cases where 2 drugs are effective
#sim.num: combination of elements to simulate (ranges from 1 to 4 as of now)

	p.vary <- c(.74,.61,.48,.36,.23)

	type <- c( 'EB10','pie','BB','Pooled')[sim.num]
	pp.vec <- list( c(.975,.97125,.96625,.95875,.9575), c(.975,.96375,.95875,.94375,.9325), rep(.975,5), c(.975,.97375,.96375,.95375,.94) )[[sim.num]]

	rr <- .7
	mem.ind <- c(rep(TRUE,2), FALSE, TRUE )[sim.num]
	prior <- c( rep('opt.source_constrain',1), rep('equal',1), NA, rep('pool',1) )[sim.num]
	constraint.val <- c(.1,rep(NA,3))[sim.num]
	AR.ind <- c( rep(TRUE, 2),FALSE, TRUE )[sim.num]
	n.burn <- c( rep(60, 2), NA, rep(60,1) )[sim.num]
	block.num <- c( rep(5, 2), NA, rep(5,1) )[sim.num]
	int.ind <- TRUE #indicator for if interim monitoring should be implemented

	vn <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=rep(1,5), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v2 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,rr,rr,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v3 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,rr,1,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v4 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,rr,1,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v5 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,rr,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v6 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,rr,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	v7 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,rr,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))

	cn <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c2 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,rr,rr,1,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c3 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,rr,1,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c4 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,rr,1,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c5 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,rr,rr,1), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c6 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,rr,1,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))
	c7 <- t(sapply(1:25000, function(x) ebola.sim(int.mont=int.ind, pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,rr,rr), nvec=rep(100,6), MEM=mem.ind, seed=x, prior=prior,constraint=constraint.val,AR=AR.ind,n.burn=n.burn,block.num=block.num)$seg_N_props))

	cap.use <- paste0(type,': Two Effective Treatments, Int Analyses=',int.ind,', AR=',AR.ind,', nburn=', n.burn,', blocks=',block.num,', RR=',rr,', PP thresholds: (',pp.vec[1],',',pp.vec[2],',',pp.vec[3],',',pp.vec[4],',',pp.vec[5],')')

	asd <- print(xtable(tab.segN_twotrt(vn, v2, v3, v4, v5,v6,v7,cn, c2, c3, c4, c5,c6,c7), caption=cap.use))

	write.table(asd, file ='paper2_simresults_mmddyy.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

}
