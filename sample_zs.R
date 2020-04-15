# function function [stateSeq INDS stateCounts] = sample_zs(data_struct,dist_struct,theta,obsModelType)
# Sample the mode and sub-mode sequence given the observations, transition
# distributions, and emission parameters. If SLDS model, the "observations"
# are the sampled state sequence.

sample_zs <- function(data_struct,dist_struct,theta,obsModelType){
	####################################
	# Define and initialize parameters #
	####################################

	# Define parameters:
	pi_z = dist_struct$pi_z  # transition distributions with pi_z(i,j) the probability of going from i->j
	pi_s = dist_struct$pi_s  # mixture weights with pi_s(i,j) the probability of s_t=j when z_t=i
	pi_init = dist_struct$pi_init  # initial distribution on z_1

	Kz = NCOL(pi_z) # truncation level for transition distributions
	Ks = NCOL(pi_s) # truncation level for MoG emissions

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
	  #print(paste("T in sample_zs=",T,sep=""))
	  INDS[[ii]]$obsIndzs = array(list(NULL),c(Kz,Ks))
	  #INDS[[ii]]$obsIndzs[1:Kz,1:Ks] = list(inds=array(0,c(1,T)),tot=0)
	  for(jj in 1:Kz){
		  for(kk in 1:Ks){
			  #INDS[[ii]]$obsIndzs[[jj,kk]] = list(inds=array(0,c(1,T)),tot=0)
			  #simplify inds
			  INDS[[ii]]$obsIndzs[[jj,kk]] = list(inds=array(0,T),tot=0)
		  }
	  }
	  
	  #print(paste("init tot=",INDS[[ii]]$obsIndzs[[1,1]]$tot,sep=""))
	  
	  # Initialize state sequence structure:
	  #stateSeq[[ii]] = list(z=array(0,c(1,T)),s=array(0,c(1,data_struct[ii]$blockEnd[length(data_struct[ii]$blockEnd)])))
	  #stateSeq has been simplified
	  stateSeq[[ii]] = list(z=array(0,T),s=array(0,data_struct[ii]$blockEnd[length(data_struct[ii]$blockEnd)]))
	}

	#for (ii in 1:length(data_struct){
	for (ii in data_struct[[1]]$test_cases){
		#print(paste("ii=",ii))
		T = NCOL(data_struct[[ii]]$blockSize)
		#print(paste("T=",T,sep=""))
		blockSize = data_struct[[ii]]$blockSize
		blockEnd = data_struct[[ii]]$blockEnd

		# Initialize state and sub-state sequences:
		#z = array(0,c(1,T))
		#s = array(0,c(1,sum(blockSize)))
		#simplify z and s
		z = array(0,T)
		s = array(0,1,sum(blockSize))
		####################################
		# Compute likelihoods and messages #
		####################################

		# Compute likelihood(kz,ks,u_i) of each observation u_i under each
		# parameter theta(kz,ks):
		likelihood = compute_likelihood(data_struct[[ii]],theta,obsModelType,Kz,Ks)[[1]]

		# Compute backwards messages:
		tempList = backwards_message_vec(likelihood, blockEnd, pi_z, pi_s)
		bwds_msg = tempList$bwds_msg
		partial_marg = tempList$partial_marg
		############################################
		# Sample the state and sub-state sequences #
		############################################

		# Sample (z(1),{s(1,1)...s(1,N1)}).  We first sample z(1) given the
		# observations u(1,1)...u(1,N1) having marginalized over the associated s's
		# and then sample s(1,1)...s(1,N1) given z(1) and the observations.

		totSeq = array(0,c(Kz,Ks))
		indSeq = array(0,c(T,Kz,Ks))

		#print(paste("first partial_marg=",partial_marg[,1],sep=""))
		#print(paste("last partial_marg=",partial_marg[,T],sep=""))
		#print(paste("pi_init=",pi_init,sep=""))

		for (t in 1:T){
			# Sample z(t):
			if (t == 1){
				Pz = pi_init * partial_marg[,1]
				obsInd = 1:blockEnd[1]
			}else{
				Pz = pi_z[z[t-1],] * partial_marg[,t]
				obsInd = (blockEnd[t-1]+1):blockEnd[t]
			}
			Pz   = cumsum(Pz)
			z[t] = 1 + sum(Pz[length(Pz)]*runif(1) > Pz)

			# Add state to counts matrix:
			if (t > 1){
				N[z[t-1],z[t]] = N[z[t-1],z[t]] + 1
			}else{
				N[Kz+1,z[t]] = N[Kz+1,z[t]] + 1  # Store initial point in "root" restaurant Kz+1
			}

			# Sample s(t,1)...s(t,Nt) and store sufficient stats:
			for (k in 1:blockSize[t]){
				# Sample s(t,k):
				if (Ks > 1){
					Ps = pi_s[z[t],] * likelihood(z[t],,obsInd[k])
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

		#print(paste("first 100 zs from sample_zs:",z[1:100],sep=""))

		stateSeq[[ii]]$z = z
		stateSeq[[ii]]$s = s

		#print(paste("length of totals:",length(INDS[[ii]]$obsIndzs[[jj,kk]]$tot),",",length(totSeq[jj,kk]),sep=""))

		for (jj in 1:Kz){
			for (kk in 1:Ks){
				#print(paste("first indSeq=",indSeq[1]))
				INDS[[ii]]$obsIndzs[[jj,kk]]$tot  = totSeq[jj,kk]
				#INDS[[ii]]$obsIndzs[[jj,kk]]$inds = array(0,indSeq[,jj,kk])
				#270212 not sure if this is correct
				INDS[[ii]]$obsIndzs[[jj,kk]]$inds = indSeq[,jj,kk]
			}
		}
	}
		
	#We don't normally have any known states, so this is not used
	for(ii in setdiff(1:length(data_struct),data_struct[[1]]$test_cases)){ # for sequences ii with fixed z_{1:T}
		print("KNOWN STATES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		T = NROW(data_struct[[ii]]$blockSize)
		blockSize = data_struct[[ii]]$blockSize
		blockEnd = data_struct[[ii]]$blockEnd
		#INDS(ii).obsIndzs(1:Kz,1:Ks) = struct('inds',sparse(1,T),'tot',0);

		# Initialize state and sub-state sequences:
		z = data_struct[[ii]]$true_labels
		s = array(1,c(1,sum(blockSize)))

		# Add s(1,1)...s(1,N1) counts and store sufficient stats:
		for (i in 1:blockSize[1]){
			# Add s(t,i) to counts matrix:
			Ns[z[1],s[i]] = Ns[z[1],s[i]] + 1
			#ERROR here 240212
			INDS[[ii]]$obsIndzs[[z[1],s[i]]]$tot = INDS[[ii]]$obsIndzs[[z[1],s[i]]]$tot + 1
			INDS[[ii]]$obsIndzs[[z[1],s[i]]]$inds[INDS[[ii]]$obsIndzs[[z[1],s[i]]]$tot] = i
		}
		# Add z(1) count:
		N[Kz+1,z[1]] = N[Kz+1,z[1]] + 1

		# Sample (z(t),{s(t,1)...s(t,Nt)}).  We first sample z(t) given the
		# observations u(t,1)...u(t,Nt) having marginalized over the associated s's
		# and then sample s(t,1)...s(t,Nt) given z(t) and the observations.
		for (t in 2:T){
			# Add state to counts matrix:
			N[z[t-1],z[t]] = N[z[t-1],z[t]]+1

			# Sample s(t,1)...s(t,Nt) and store sufficient stats:
			for (i in 1:blockSize[t]){
				obsInd = blockEnd[t-1] + i

				# Add s[t,i] to counts matrix:
				Ns[z[t],s[obsInd]] = Ns[z[t],s[obsInd]] + 1

				INDS[[ii]]$obsIndzs[z[t],s[obsInd]]$tot = INDS[[ii]]$obsIndzs[z[t],s[obsInd]]$tot + 1
				INDS[[ii]]$obsIndzs[z(t),s[obsInd]]$inds[INDS[[ii]]$obsIndzs[z[t],s[obsInd]]$tot] = obsInd
			}
		}

		stateSeq[[ii]]$z = z
		stateSeq[[ii]]$s = s

	}

	binNs = array(0,dim(Ns))
	binNs[which(Ns>0)] = 1
	uniqueS = rowSums(binNs)

	stateCounts = list()
	stateCounts$uniqueS = uniqueS
	stateCounts$N = N
	stateCounts$Ns = Ns

	return(list(stateSeq=stateSeq, INDS=INDS, stateCounts=stateCounts))
}