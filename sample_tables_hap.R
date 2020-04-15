# function stateCounts = sample_tables(stateCounts,hyperparams,beta_vec,Kz)
# Sample the number of tables in restaurant i serving dish j given the mode
# sequence z_{1:T} and hyperparameters. Also sample the override variables.

sample_tables <- function(stateCounts,hyperparams,beta_vec,Kz){
	# Split \alpha and \kappa using \rho:
	rho0 = hyperparams$rho0
	rho0_hap = hyperparams$rho0_hap
	#print(paste("rho0_hap in sample_tables=",rho0_hap,sep=""))
	alpha0 = hyperparams$alpha0_p_kappa0*(1-rho0-rho0_hap)
	kappa0 = hyperparams$alpha0_p_kappa0*rho0
	kappa0_hap = hyperparams$alpha0_p_kappa0*rho0_hap

	N = stateCounts$N

	# Sample M, where M(i,j) = # of tables in restaurant i served dish j:
	#M = randnumtable(cbind(alpha0*array(beta_vec,c(Kz,Kz))+diag(kappa0,nrow=Kz,ncol=Kz),alpha0*t(beta_vec)),N)
	kappa_vec=diag(kappa0,nrow=Kz,ncol=Kz)
	kappa_vec[cbind(seq(1,Kz-1,2),seq(2,Kz,2))] = kappa0_hap
	kappa_vec[cbind(seq(2,Kz,2),seq(1,Kz-1,2))] = kappa0_hap
	M = randnumtable(cbind(alpha0*array(beta_vec,c(Kz,Kz))+kappa_vec,alpha0*t(beta_vec)),N)
	# Sample barM (the table counts for the underlying restaurant), where
	# barM(i,j) = # tables in restaurant i that considered dish j:
	tempList = sample_barM(M,beta_vec,rho0,rho0_hap)

	stateCounts$M = M
	stateCounts$barM = tempList$barM
	stateCounts$sum_w = tempList$sum_w
	stateCounts$sum_w_hap = tempList$sum_w_hap

	return(stateCounts)
}