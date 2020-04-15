randdirichlet <-function(alpha){
	#library(MCMCpack)
	x=rdirichlet(1,alpha)
	return(x)
}