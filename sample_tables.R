# function stateCounts = sample_tables(stateCounts,hyperparams,beta_vec,Kz)
# Sample the number of tables in restaurant i serving dish j given the mode
# sequence z_{1:T} and hyperparameters. Also sample the override variables.

sample_tables <- function(stateCounts,hyperparams,beta_vec,Kz){
	# Split \alpha and \kappa using \rho:
	rho0 = hyperparams$rho0
	alpha0 = hyperparams$alpha0_p_kappa0*(1-rho0)
	kappa0 = hyperparams$alpha0_p_kappa0*rho0

	N = stateCounts$N

	# Sample M, where M(i,j) = # of tables in restaurant i served dish j:
	M = randnumtable(cbind(alpha0*array(beta_vec,c(Kz,Kz))+diag(kappa0,nrow=Kz,ncol=Kz),alpha0*t(beta_vec)),N)
	# Sample barM (the table counts for the underlying restaurant), where
	# barM(i,j) = # tables in restaurant i that considered dish j:
	tempList = sample_barM(M,beta_vec,rho0)
	barM=tempList$barM
	sum_w=tempList$sum_w

	stateCounts$M = M
	stateCounts$barM = barM
	stateCounts$sum_w = sum_w

	return(stateCounts)
}