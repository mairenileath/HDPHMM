backwards_message_vec <- function(likelihood, blockEnd, pi_z, pi_s){
	# Allocate storage space
	Kz = NCOL(pi_z)
	Ks = NCOL(pi_s)
	T  = NCOL(blockEnd)
	#print(paste("T in bmv=",T,sep=""))
	#print(paste("bmv blockEnd=",blockEnd))

	bwds_msg     = array(1,c(Kz,T))
	partial_marg = array(0,c(Kz,T))

	# Compute marginalized likelihoods for all times, integrating s_t
	if (Kz==1 & Ks==1){
		marg_like = likelihood[1,1,]
	}else{
		#print(dim(likelihood))
		#print(dim(array(pi_s,c(Kz,Ks,T))))
		#pi_s_times_likelihood = likelihood * array(pi_s,c(Kz,Ks,T))
		#use log-space to avoid values going to zero
		log_pi_s_times_likelihood = log(likelihood+1e-22) + array(log(pi_s+1e-22),c(Kz,Ks,T))
		norm_factor = sapply(1:T, function(i,l) max(l[,,i],na.rm=T), l=log_pi_s_times_likelihood)
		pi_s_times_likelihood = exp(log_pi_s_times_likelihood - array(norm_factor,c(Kz,Ks,T)))
		
		if(Ks==1){
			#marg_like = sapply(1:T,function(i,p) p[,1,i],p=pi_s_times_likelihood)
			#050312 - should be faster
			marg_like = pi_s_times_likelihood[,1,]
		}else{
			#marg_like = sapply(1:T,function(i,p) colSums(p[,,i]),p=pi_s_times_likelihood)
			marg_like = sapply(1:T,function(i,p) rowSums(p[,,i]),p=pi_s_times_likelihood)
		}
	}
	#print(paste("dim(marg_like=",dim(marg_like),sep=""))
	marg_like[is.na(marg_like) | is.nan(marg_like)] = 0
	
	# If necessary, combine likelihoods within blocks, avoiding underflow
	if (T < blockEnd[T]){
		marg_like = log(marg_like+1e-22)

		block_like = array(0,c(Kz,T))
		block_like[,1] = rowSums(marg_like[,1:blockEnd[1]])
	  	#for (tt in 2:T){
		#	block_like[,tt] = rowSums(marg_like[,(blockEnd[tt-1]+1):blockEnd[tt]])
		#}
		#050312 - should be faster
		block_like[,2:T] = t(sapply(2:T, function(i,m,b) rowSums(m[,(bb[i-1]+1):bb[i]]),m=marg_like,b=blockEnd))

	  	#block_norm = sapply(1:Kz, function(i,b) max(b[,i],na.rm=T), b=block_like)
	  	block_norm = sapply(1:T, function(i,b) max(b[,i],na.rm=T), b=block_like)
	  	block_like = exp(block_like - array(block_norm,c(Kz,T)))
	}else{
	  	block_like = marg_like
	}
	#print(paste("min/max block_like=",min(block_like),",",max(block_like),sep=""))
	#print(paste("dims pi_z=",dim(pi_z)))
	#print(paste("dims partial_marg=",dim(partial_marg)))
	#print(paste("dims bwds_msg=",dim(bwds_msg)))
	#print(paste("first block_like=",block_like[,1]))
	#print(paste("first bwds_msg=",bwds_msg[,1]))

	bwds_msg[,T] = bwds_msg[,T] / sum(c(bwds_msg[,T],1e-22),na.rm=T)
	bwds_msg[is.na(bwds_msg[,T]) | is.nan(bwds_msg[,T]),T] = 0
	
	#050312 do multiplication in one go - should be quicker. THIS WON'T WORK
	#partial_marg[,] = block_like[,] * bwds_msg[,]
	#partial_marg[is.na(partial_marg[,]) | is.nan(partial_marg[,]),]=0
	#norm_factor = sapply(1:T, function(i,p) max(c(p[,,i],1e-22),na.rm=T), p=partial_marg)
	#partial_marg[,] = partial_marg[,] / t(array(norm_factor,c(T,Kz)))
	
	# Compute messages backwards in time
	for (tt in seq(T-1,1,-1)){
	  # Multiply likelihood by incoming message:
	  partial_marg[,tt+1] = block_like[,tt+1] * bwds_msg[,tt+1]
	  #288212 remove NAs (underflow)
	  partial_marg[is.na(partial_marg[,tt+1]) | is.nan(partial_marg[,tt+1]),tt+1]=0
	  #290212 scale partial_marg
	  partial_marg[,tt+1] = partial_marg[,tt+1] / max(c(partial_marg[,tt+1],1e-22),na.rm=T)#maybe we should sum rather than max??

	  # Integrate out z_t:
	  bwds_msg[,tt] = pi_z %*% partial_marg[,tt+1]
	  bwds_msg[,tt] = bwds_msg[,tt] / sum(c(bwds_msg[,tt],1e-22),na.rm=T)
	  bwds_msg[is.na(bwds_msg[,tt]) | is.nan(bwds_msg[,tt]),tt] = 0
	}

	# Compute marginal for first time point
	partial_marg[,1] = block_like[,1] * bwds_msg[,1]	
	#288212 remove NAs (underflow)
	partial_marg[is.na(partial_marg[,1]) | is.nan(partial_marg[,1]),1]=0
	#290212 scale partial_marg
	partial_marg[,1] = partial_marg[,1] / max(c(partial_marg[,1],1e-22),na.rm=T)#maybe we should sum rather than max??
	  
	return(list(bwds_msg = bwds_msg,partial_marg=partial_marg))
}