sample_zs_init <- function(data_struct,dist_struct,obsModelType)
{

	####################################
	# Define and initialize parameters #
	####################################

	# Define parameters:
	pi_z = dist_struct$pi_z
	pi_s = dist_struct$pi_s
	pi_init = dist_struct$pi_init

	Kz = NCOL(pi_z)
	Ks = NCOL(pi_s)

	# Initialize state count matrices:
	N = array(0,c(Kz+1,Kz))
	Ns = array(0,c(Kz,Ks))
	
	if(!("test_cases" %in% names(data_struct[[1]]))){
		data_struct[[1]]$test_cases = 1:length(data_struct)
	}

	# Preallocate INDS
	INDS=vector("list",length(data_struct))
	stateSeq=vector("list",length(data_struct))
	for (ii in 1:length(data_struct)){
	  #T = NROW(data_struct[[ii]]$blockSize)
	  T = NCOL(data_struct[[ii]]$blockSize)
	  #print(paste("T in sample_zs_init=",T,sep=""))
	  INDS[[ii]]$obsIndzs = array(list(NULL),c(Kz,Ks))
	  #INDS[[ii]]$obsIndzs[1:Kz,1:Ks] = list(inds=array(0,c(1,T)),tot=0)
	  for(jj in 1:Kz){
		  for(kk in 1:Ks){
			  #INDS[[ii]]$obsIndzs[[jj,kk]] = list(inds=array(0,c(1,T)),tot=0)
			  #simplify inds
			  INDS[[ii]]$obsIndzs[[jj,kk]] = list(inds=array(0,T),tot=0)
		  }
	  }
	}

	#print("OK1")

	for (ii in data_struct[[1]]$test_cases){
		if("z_init" %in% names(data_struct[[1]])){
			tempList = setZtoFixedSeq(data_struct[[ii]],dist_struct,N,Ns,data_struct[[ii]]$z_init,1)
		}else{
			tempList = sampleZfromPrior(data_struct[[ii]],dist_struct,N,Ns)
		}
		stateSeq[[ii]]$z = tempList$z
		stateSeq[[ii]]$s = tempList$s
		totSeq = tempList$totSeq
		indSeq = tempList$indSeq
		N = tempList$N
		Ns = tempList$Ns
		
		#print("OK1.1")
		#print(paste("dims(totseq=",dim(totSeq),")",sep=""))
		#print(paste("dims(indseq=",dim(indSeq),")",sep=""))
		#print(paste("dims(tot=",dim(INDS[[ii]]$obsIndzs[[1,1]]$tot),")",sep=""))
		#print(paste("dims(inds=",dim(INDS[[ii]]$obsIndzs[[1,1]]$inds),")",sep=""))
		#print(paste("dim(INDS[[ii]]$obsIndzs)=",dim(INDS[[ii]]$obsIndzs),")",sep=""))
		for (jj in 1:Kz){
			for (kk in 1:Ks){
				INDS[[ii]]$obsIndzs[[jj,kk]]$tot  = totSeq[jj,kk]
				INDS[[ii]]$obsIndzs[[jj,kk]]$inds = indSeq[,jj,kk]
			}
		}

	}
	
#print("OK2")

	for (ii in setdiff(1:length(data_struct),data_struct[[1]]$test_cases)){
		tempList = setZtoFixedSeq(data_struct[[ii]],dist_struct,N,Ns,data_struct[[ii]]$true_labels,0)
		stateSeq[[ii]]$z = tempList$z
		stateSeq[[ii]]$s = tempList$s
		totSeq = tempList$totSeq
		indSeq = tempList$indSeq
		N = tempList$N
		Ns = tempList$Ns

		for (jj in 1:Kz){
			for (kk in 1:Ks){
				INDS[[ii]]$obsIndzs[[jj,kk]]$tot  = totSeq[jj,kk]
				INDS[[ii]]$obsIndzs[[jj,kk]]$inds = indSeq[,jj,kk]
			}
		}

	}
	
#print("OK3")

	binNs = array(0,dim(Ns))
	binNs[which(Ns>0)] = 1
	uniqueS = rowSums(binNs)

	stateCounts = list(uniqueS = uniqueS, N = N, Ns = Ns)
	return(list(stateSeq=stateSeq, INDS=INDS, stateCounts=stateCounts))
}


sampleZfromPrior <- function(data_struct,dist_struct,N,Ns){

	# Define parameters:
	pi_z = dist_struct$pi_z
	pi_s = dist_struct$pi_s
	pi_init = dist_struct$pi_init

	Kz = NCOL(pi_z)
	Ks = NCOL(pi_s)

	T = NCOL(data_struct$blockSize)
	blockSize = data_struct$blockSize
	blockEnd = data_struct$blockEnd

	# Initialize state and sub-state sequences:
	z = array(0,T)
	s = array(0,sum(blockSize))

	############################################
	# Sample the state and sub-state sequences #
	############################################

	# Sample (z(1),{s(1,1)...s(1,N1)}).  We first sample z(1) given the
	# observations u(1,1)...u(1,N1) having marginalized over the associated s's
	# and then sample s(1,1)...s(1,N1) given z(1) and the observations.

	totSeq = array(0,c(Kz,Ks))
	indSeq = array(0,c(T,Kz,Ks))

	for (t in 1:T){
		# Sample z(t):
		if (t == 1){
			Pz = pi_init
			obsInd = 1:blockEnd[1]
		}else{
			Pz = pi_z[z[t-1],]
			obsInd = (blockEnd[t-1]+1):blockEnd[t]
		}
		Pz   = cumsum(Pz)
		z[t] = 1 + sum(Pz[length(Pz)]*runif(1) > Pz);

		# Add state to counts matrix:
		if (t > 1){
			N[z[t-1],z[t]] = N[z[t-1],z[t]] + 1
		}else{
			N[Kz+1,z[t]] = N[Kz+1,z[t]] + 1  # Store initial point in "root" restaurant Kz+1
		}

		# Sample s(t,1)...s(t,Nt) and store sufficient stats:
		#This could be done quicker without a loop!!
		for (k in 1:blockSize[t]){
			# Sample s(t,k):
			if (Ks > 1){
				Ps = pi_s[z[t],]
				Ps = cumsum(Ps)
				s[obsInd[k]] = 1 + sum(Ps[length(Ps)]*runif(1) > Ps)
			}else{
				s[obsInd[k]] = 1
			}

			# Add s(t,k) to count matrix and observation statistics:
			Ns[z[t],s[obsInd[k]]] = Ns[z[t],s[obsInd[k]]] + 1
			totSeq[z[t],s[obsInd[k]]] = totSeq[z[t],s[obsInd[k]]] + 1
			indSeq[totSeq[z[t],s[obsInd[k]]],z[t],s[obsInd[k]]] = obsInd[k]
		}
	}
	return(list(z=z, s=s, totSeq=totSeq, indSeq=indSeq, N=N, Ns=Ns))
}

setZtoFixedSeq <- function(data_struct,dist_struct,N,Ns,z_fixed,sampleS){        
	# Define parameters:
	pi_z = dist_struct$pi_z
	pi_s = dist_struct$pi_s
	pi_init = dist_struct$pi_init

	Kz = NCOL(pi_z)
	Ks = NCOL(pi_s)

	T = NCOL(data_struct$blockSize)
	blockSize = data_struct$blockSize
	blockEnd = data_struct$blockEnd

	totSeq = array(0,c(Kz,Ks))
	indSeq = array(0,c(T,Kz,Ks))

	# Initialize state and sub-state sequences:
	z = z_fixed
	if (sampleS == 1){
		for (t in 1:T){
			# Sample z(t):
			if (t == 1){
				obsInd = 1:blockEnd[1]
			}else{
				obsInd = (blockEnd[t-1]+1):blockEnd[t]
			}

			# Sample s(t,1)...s(t,Nt) and store sufficient stats:
			for (k in 1:blockSize[t]){
				# Sample s(t,k):
				if (Ks > 1){
					Ps = pi_s[z[t],]
					Ps = cumsum(Ps)
					s[obsInd[k]] = 1 + sum(Ps[length(Ps)]*runif(1) > Ps)
				}else{
					s[obsInd[k]] = 1
				}
			}
		}
	}else{
		s = array(1,sum(blockSize))
	}


	for (t in 1:T){
		# Sample z(t):
		if (t == 1){
			obsInd = 1:blockEnd[1]
		}else{
			obsInd = (blockEnd[t-1]+1):blockEnd[t]
		}

		# Add state to counts matrix:
		if (t > 1){
			N[z[t-1],z[t]] = N[z[t-1],z[t]] + 1
		}else{
			N[Kz+1,z[t]] = N[Kz+1,z[t]] + 1  # Store initial point in "root" restaurant Kz+1
		}

		# Sample s(t,1)...s(t,Nt) and store sufficient stats:
		for (k in 1:blockSize[t]){
			# Add s(t,k) to count matrix and observation statistics:
			Ns[z(t),s[obsInd[k]]] = Ns[z(t),s[obsInd[k]]] + 1
			totSeq[z[t],s[obsInd[k]]] = totSeq[z[t],s[obsInd[k]]] + 1
			indSeq[totSeq[z[t],s[obsInd[k]]],z[t],s[obsInd[k]]] = obsInd[k]
		}
	}

	return(list(z=z,s=s,totSeq=totSeq,indSeq=indSeq,N=N,Ns=Ns))
}
