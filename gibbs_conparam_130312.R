gibbs_conparam <- function(alpha, numdata, numclass, aa, bb, numiter=1){
	#gibbs_conparam     Auxiliary variable resampling of DP concentration parameter


	numgroup   = length(numdata)
	totalclass = sum(numclass)

	xx = array(0,c(1,numgroup))

	A=array(0,c(numgroup,2))
	A[,1]=alpha+1
	A[,2]=t(numdata)
	A=t(A)

	for (ii in 1:numiter){
	  # beta auxiliary variables
	  #for jj = 1:numgroup
	  #  xj     = dirichlet_rnd([alpha+1 numdata(jj)], 1);
	  #  xx(jj) = xj(1);
	  #end  
	  xj=randdirichlet(A)
	  xx=xj[1,]

	  # binomial auxiliary variables
	  zz = (runif(numgroup)*(alpha+numdata)) < numdata

	  # gamma resampling of concentration parameter
	  gammaa = aa + totalclass - sum(zz)
	  gammab = bb - sum(log(xx))
	  alpha  = rgamma(1,gammaa) / gammab
	}
	return(alpha)
}