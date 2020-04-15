sample_barM <- function(M,beta_vec,rho0,rho0_hap){

	barM = M
	sum_w = array(0,NCOL(M))
	sum_w_hap = array(0,NCOL(M))
	
	#for (j in 1:NCOL(M)){
	#	p1 = rho0/(beta_vec[j]*(1-rho0-rho0_hap) + rho0 + rho0_hap)
	#	if(j%%2 ==0){
	#		p2 = rho0_hap/(beta_vec[j-1]*(1-rho0-rho0_hap) + rho0 + rho0_hap)
	#	}else{
	#		p2 = rho0_hap/(beta_vec[j+1]*(1-rho0-rho0_hap) + rho0 + rho0_hap)
	#	}
	#	p3 = 1-(p1+p2)
	#	overridden = rmultinom(1,M[j,j],c(p1,p2,p3))
	#	sum_w[j] = overridden[1]
	#	sum_w_hap[j] = overridden[2]
	#	barM[j,j] = M[j,j] - sum_w[j] - sum_w_hap[j]
	#}
	
	for (j in 1:NCOL(M)){
		p = rho0/(beta_vec[j]*(1-rho0-rho0_hap) + rho0)
		sum_w[j] = rbinom(1,M[j,j],p)
		barM[j,j] = M[j,j] - sum_w[j]
		
		p = rho0_hap/(beta_vec[j]*(1-rho0-rho0_hap) + rho0_hap)
		
		if(j %% 2 == 0){
			sum_w_hap[j] = rbinom(1,M[j-1,j],p)
			barM[j-1,j] = M[j-1,j] - sum_w_hap[j]
		}else{
			sum_w_hap[j] = rbinom(1,M[j+1,j],p)
			barM[j+1,j] = M[j+1,j] - sum_w_hap[j]
		}
	}

	return(list(barM=barM,sum_w=sum_w,sum_w_hap=sum_w_hap))
}