sample_hyperparams_init <-function(stateCounts,hyperparams,hyperhyperparams,HMMmodelType,resample_kappa){

	# Hyperparams for gamma dist over \alpha+\kappa, where transition distributions
	# \pi_j \sim DP(\alpha+\kappa,(\alpha\beta+\kappa\delta(j))/(\alpha+\kappa))
	# = DP(alpha_p_kappa, (1-\rho)*\beta + \rho*\delta(j)):
	a_alpha=hyperhyperparams$a_alpha
	b_alpha=hyperhyperparams$b_alpha

	# Hyperparams for beta dist over \rho, where \rho relates \alpha+\kappa to
	# \alpha and \kappa individually.
	c=hyperhyperparams$c
	c_hap=hyperhyperparams$c_hap
	d=hyperhyperparams$d

	Ns = stateCounts$Ns

	if(HMMmodelType=="HDP"){

			# Hyperparams for gamma dist over \gamma, where avg transition distribution
			# \beta \sim stick(\gamma):
			a_gamma=hyperhyperparams$a_gamma
			b_gamma=hyperhyperparams$b_gamma

			# Resample concentration parameters:
			alpha0_p_kappa0 = a_alpha / b_alpha    # Gj concentration parameter
			gamma0 = a_gamma / b_gamma    # G0 concentration parameter

			hyperparams$gamma0 = gamma0

	}else if(HMMmodelType=="finite"){
			#         alpha0_p_kappa0 = alpha0_p_kappa0;
			# Resample concentration parameters:
			alpha0_p_kappa0 = a_alpha / b_alpha    # Gj concentration parameter
	}
	
	# MAYBE CHANGE THIS WHEN HDP-HMM with FINITE EMISSION STUFF IS ADDED
	if (NROW(Ns)>1){ 
		# Hyperparams for gamma dist over \sigma, where HMM-state-specific mixture
		# weights \psi_j \sim stick(\sigma):
		a_sigma=hyperhyperparams$a_sigma
		b_sigma=hyperhyperparams$b_sigma

		sigma0 = a_sigma / b_sigma
	}else{
		sigma0 = 1
	}

	if (resample_kappa){
		rho0 = c/(c+c_hap+d)
		rho0_hap = c_hap/(c+c_hap+d)
	}else{
		rho0 = 0 # betarnd(0,1);
	}

	hyperparams$alpha0_p_kappa0 = alpha0_p_kappa0
	hyperparams$sigma0 = sigma0
	hyperparams$rho0 = rho0
	hyperparams$rho0_hap = rho0_hap

	#print(paste("rho0_hap in sample_hyperparams=",rho0_hap,sep=""))

	return(hyperparams)
}