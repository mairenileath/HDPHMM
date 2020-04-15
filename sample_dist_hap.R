sample_dist <-function(stateCounts,hyperparams,model){
	source("randdirichlet.R")
	Kz = NROW(stateCounts$Ns) # truncation level for transition distributions
	Ks = NCOL(stateCounts$Ns) # truncation level for mode-specific MoG emission distributions

	#print(hyperparams)
	
	# Define alpha0 and kappa0 in terms of alpha0+kappa0 and rho0:
	alpha0 = hyperparams$alpha0_p_kappa0*(1-hyperparams$rho0-hyperparams$rho0_hap)
	kappa0 = hyperparams$alpha0_p_kappa0*hyperparams$rho0
	kappa0_hap = hyperparams$alpha0_p_kappa0*hyperparams$rho0_hap
	sigma0 = hyperparams$sigma0

	#print(paste("rho0_hap in sample_dist=",hyperparams$rho0_hap,sep=""))

	N = stateCounts$N  # N(i,j) = # z_t = i to z_{t+1}=j transitions. N(Kz+1,i) = 1 for i=z_1.
	Ns = stateCounts$Ns  # Ns(i,j) = # s_t = j given z_t=i
	barM = stateCounts$barM  # barM(i,j) = number of tables in restaurant i that considered dish j

	#print(paste(NROW(barM),NCOL(barM),sep=","))
	
	if(model$HMMmodel$type == "HDP"){
			# Sample beta, the global menu, given new barM:
			gamma0 = hyperparams$gamma0
			beta_vec = randdirichlet(colSums(barM) + gamma0/Kz)
	}else if(model$HMMmodel$type == "finite"){
			# A finite HMM model with a sparse Dirichlet prior is exactly
			# equivalent to the truncated HDP-HMM model with a uniform global
			# menu:
			beta_vec = matrix(1/Kz,c(1,Kz))
	}

	pi_z = array(0,c(Kz,Kz))
	pi_s = array(0,c(Kz,Ks))
	for (j in 1:Kz){
		kappa_vec = array(0,c(1,Kz))
		# Add an amount \kappa to Dirichlet parameter corresponding to a
		# self-transition:
		kappa_vec[j] = kappa0
		#050312 give preference to haplotype block switch over CNV segment switch using a hyperparameter
		if (j %% 2 == 0){
			kappa_vec[j-1] = kappa0_hap
		}else{
			kappa_vec[j+1] = kappa0_hap
		}
		# Sample \pi_j's given sampled \beta_vec and counts N, where
		# DP(\alpha+\kappa,(\alpha\beta+\kappa\delta(j))/(\alpha+\kappa)) is
		# Dirichlet distributed over the finite partition defined by beta_vec:
		pi_z[j,] = randdirichlet(alpha0*beta_vec + kappa_vec + N[j,])
		# Sample HMM-state-specific mixture weights \psi_j's with truncation
		# level Ks given sampled s stats Ns:
		pi_s[j,] = randdirichlet(Ns[j,] + sigma0/Ks)
	}
	alpha_vec = alpha0*beta_vec + N[Kz+1,]
	sum_alpha = sum(alpha_vec)
	#very low values of alpha cause rdirichlet to crash, so scale up alpha
	if(sum_alpha<0.01){
		alpha_vec = alpha_vec * 0.01/sum_alpha
	}
	pi_init = randdirichlet(alpha_vec)
	#print(paste("beta_vec=",beta_vec,sep=""))
	#print(paste("pi_init=",pi_init,sep=""))
	#print(paste("alpha0=",alpha0,sep=""))
	#print(paste("alpha0*beta_vec=",alpha0*beta_vec,sep=""))
	#print(paste("alpha_vec=",alpha_vec,sep=""))
	#print(paste("N[Kz+1,]=",N[Kz+1,],sep=""))
	dist_struct = list(pi_z = pi_z,pi_init = pi_init,pi_s = pi_s,beta_vec = beta_vec)
	return(dist_struct)
}