randiwishart <-function(sigma,df,di){
	#RANDIWISHART Generate inverse Wishart random matrix
	#   W=RANDIWISHART(SIGMA,DF) generates a random matrix W from the inverse
	#   Wishart distribution with parameters SIGMA and DF.  The inverse of W
	#   has the Wishart distribution with covariance matrix inv(SIGMA) and DF
	#   degrees of freedom.
	#
	#   W=RANDIWISHART(SIGMA,DF,DI) expects DI to be the Cholesky factor of
	#   the inverse of SIGMA.
	#
	#   [W,DI]=RANDIWISHART(SIGMA,DF) returns DI so it can be used again in
	#   future calls to RANDIWISHART.

	n = NROW(sigma)
	if (df<n){ # require this to ensure invertibility
	   stop('randiwish:BadDf',paste('Degrees of freedom must be no smaller than the dimension of SIGMA. (',toString(df),',',toString(n),')',sep=""))
	}

	# Get Cholesky factor for inv(sigma) unless that's already done
	if (nargs()<3){
		#     [d,p] = chol(sigma,0);
		#     if p~=0
		#         error('stats:iwishrnd:BadCovariance',...
		#             'Covariance matrix must be symmetric and positive definite.');
		#     end
		#if(length(sigma)==1){
		#	d=sigma
		#}else{
			#d = Cholesky(sigma)
			d = chol(sigma)
		#}	
		di = solve(t(d),diag(nrow(d))) # either take inverse here and scale chol of
		#randwishart sample and then take inverse of sample, or take inverse of
		#sample and then scale after w/o the inverse.
	}

	a = rwish(df/2,n)
	
	print(paste("a=",a,sep=""))
	print(paste("n=",n,sep=""))
	print(paste("di=",di,sep=""))

	#sqrtinvx = sqrt(2)*a[1,1]*di
	sqrtinvx = sqrt(2)*a %*% di
	sqrtx = t(solve(sqrtinvx,diag(nrow(sqrtinvx))))

	# x = 2*(a'*a);
	# x = x\eye(size(x));
	# x = d'*(x*d);

	# sqrtx = sqrt(2)*a;
	# sqrtx = (sqrtx\eye(size(sqrtx)))';
	# sqrtx = sqrtx*d;


	return(list(sqrtx=sqrtx,sqrtinvx=sqrtinvx,di=di))
}