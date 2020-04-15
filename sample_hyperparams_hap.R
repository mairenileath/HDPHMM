# function [hyperparams] = sample_hyperparams(stateCounts,hyperparams,hyperhyperparams,HMMmodelType,resample_kappa)
# Sample concentration parameters that define the distribution on transition
# distributions and mixture weights of the various model components.

sample_hyperparams <- function(stateCounts,hyperparams,hyperhyperparams,HMMmodelType,resample_kappa,resample_haplotypeSwitching=TRUE){

	print(paste("resample_haplotypeSwitching in sample_hyperparams=",resample_haplotypeSwitching,sep=""))

	# Hyperparams for Gamma dist over \alpha+\kappa, where transition distributions
	# \pi_j \sim DP(\alpha+\kappa,(\alpha\beta+\kappa\delta(j))/(\alpha+\kappa))
	# = DP(alpha_p_kappa, (1-\rho)*\beta + \rho*\delta(j)):
	a_alpha=hyperhyperparams$a_alpha
	b_alpha=hyperhyperparams$b_alpha

	# Hyperparams for Beta dist over \rho, where \rho relates \alpha+\kappa to
	# \alpha and \kappa individually.
	c=hyperhyperparams$c
	c_hap=hyperhyperparams$c_hap
	d=hyperhyperparams$d

	# Grab out last value of the hyperparameters:
	alpha0_p_kappa0 = hyperparams$alpha0_p_kappa0
	sigma0 = hyperparams$sigma0

	N = stateCounts$N # N(i,j) = # z_t = i to z_{t+1}=j transitions in z_{1:T}. N(Kz+1,i) = 1 for i=z_1.
	Ns = stateCounts$Ns # Ns(i,k) = # of obs assigned to mix component k in mode i (i.e., # s_t = k given z_t=i)
	uniqueS = stateCounts$uniqueS # uniqueS(i) = sum_j Ns(i,j) = # of mixture components for HMM-state i
	M = stateCounts$M # M(i,j) = # of tables in restaurant i serving dish k
	barM = stateCounts$barM # barM(i,j) = # of tables in restaurant i considering dish k
	sum_w = stateCounts$sum_w # sum_w(i) = # of overriden dish assignments in restaurant i
	sum_w_hap = stateCounts$sum_w_hap # sum_w(i) = # of overriden dish assignments in restaurant i

	Nkdot = rowSums(N)
	Mkdot = rowSums(M)
	Nskdot = rowSums(Ns)
	barK =sum((colSums(barM)>0))
	validindices = which(Nkdot>0)
	validindices2 = which(Nskdot>0)

	if(HMMmodelType=="HDP"){
			# Hyperparams for gamma dist over \gamma, where avg transition distribution
			# \beta \sim stick(\gamma):
			a_gamma=hyperhyperparams$a_gamma
			b_gamma=hyperhyperparams$b_gamma

			gamma0 = hyperparams$gamma0

			# Resample concentration parameters:
			if (length(validindices)==0){
				alpha0_p_kappa0 = rgamma(1,a_alpha) / b_alpha    # Gj concentration parameter
				gamma0 = rgamma(1,a_gamma) / b_gamma    # G0 concentration parameter
			}else{
				alpha0_p_kappa0  = gibbs_conparam(alpha0_p_kappa0, Nkdot[validindices],Mkdot[validindices],a_alpha,b_alpha,50)
				gamma0 = gibbs_conparam(gamma0,sum(barM),barK,a_gamma,b_gamma,50)
			}

			hyperparams$gamma0 = gamma0
			print(paste("gamma0 in sample_hyperparams=",hyperparams$gamma0,sep=""))
	}else if(HMMmodelType=="finite"){
			# Resample Dirichlet parameter for \pi_j \sim
			# Dir(\alpha/L,...,\alpha/L + \kappa,...,\alpha/L):
			if (length(validindices)==0){
				alpha0_p_kappa0 = rgamma(1,a_alpha) / b_alpha   
			}else{
				alpha0_p_kappa0  = gibbs_conparam(alpha0_p_kappa0, Nkdot[validindices],Mkdot[validindices],a_alpha,b_alpha,50)
			}
	}

	if (NCOL(Ns)>1){ # Only spend time resampling sigma0 if MoG mode-specific emission distribution

		# Hyperparams for Gamma dist over \sigma, where HMM-state-specific mixture
		# weights \psi_j \sim stick(\sigma):
		a_sigma=hyperhyperparams$a_sigma
		b_sigma=hyperhyperparams$b_sigma

		if (length(validindices2)==0){
			sigma0 = rgamma(1,a_sigma) / b_sigma
		}else{
			sigma0 = gibbs_conparam(sigma0,Nskdot(validindices2),uniqueS(validindices2),a_sigma,b_sigma,50)
		}

	}else{
		sigma0 = 1
	}


	if (resample_kappa==1){  # Only spend time resampling rho0 if sticky model

		# Resample self-transition proportion parameter:
		if(resample_haplotypeSwitching){		
			rho_vec = rdirichlet(1,c(c+sum(sum_w),c_hap+sum(sum_w_hap),d+(sum(M)-sum(sum_w)-sum(sum_w_hap))))
			rho0 = rho_vec[1]
			rho0_hap = rho_vec[2]
		}else{
			rho0_hap = c_hap/(c+c_hap+d)
			rho0 = rbeta(1,c+sum(sum_w),d+(sum(M)-sum(sum_w))) * (c+d)/(c+c_hap+d)
		}
		
	}else{
		rho0 = 0 #rbeta(1,0,1);
		rho0_hap=0
	}

	hyperparams$alpha0_p_kappa0 = alpha0_p_kappa0
	hyperparams$sigma0 = sigma0
	hyperparams$rho0 = rho0
	hyperparams$rho0_hap = rho0_hap
	
	#print(paste("rho0_hap in sample_hyperparams=",rho0_hap,sep=""))
	return(hyperparams)
 }
