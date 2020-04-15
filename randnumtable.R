randnumtable <- function(alpha,numdata){
	numtable=array(0,dim(numdata))
	#print(paste("numdata len=",length(numtable),sep=""))
	for (ii in 1:length(numdata)){
		#if(ii %% 100 ==0){
		#	print(paste("numdata[",ii,"]=",numdata[ii],sep=""))
		#}
		if(numdata[ii]>0){
			numtable[ii]=1+sum(runif(numdata[ii]-1)<(array(1,numdata[ii]-1)*alpha[ii]/(alpha[ii]+(1:(numdata[ii]-1)))))
		}
	}
	numtable[which(numdata==0)]=0
	return(numtable)
}