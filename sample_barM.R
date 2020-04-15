sample_barM <- function(M,beta_vec,rho0){

	barM = M
	sum_w = array(0,NCOL(M))

	for (j in 1:NCOL(M)){
		if (rho0>0){
			p = rho0/(beta_vec[j]*(1-rho0) + rho0)
		}else{
			p = 0
		}
		sum_w[j] = rbinom(1,M[j,j],p)
		barM[j,j] = M[j,j] - sum_w[j]
	}

	return(list(barM=barM,sum_w=sum_w))
}