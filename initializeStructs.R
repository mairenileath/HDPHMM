initializeStructs<-function(model,data_struct,settings){
	Kz = settings$Kz
	Ks = settings$Ks

	#print(paste("initializeStructs Kz=",Kz,sep=""))

	prior_params = model$obsModel$params
	#240212 - not necessary because blockSize is set further down
	#if(!("blockSize" %in% names(data_struct[1]))){
	#	data_struct[[1]]$blockSize = NA
	#}

	if(!("test_cases" %in% names(data_struct[1]))){
		data_struct[[1]]$test_cases = 1
	}

	if (is.null(data_struct[1]$test_cases)){
		data_struct[[1]]$test_cases = 1
	}


	if(model$obsModel$type == "Gaussian"){
		#dimu = NCOL(data_struct[[1]]$obs)
		dimu = NROW(data_struct[[1]]$obs)
		print(paste("dimu=",dimu,", Kz=",Kz,", Ks=",Ks,sep=""))
		#print(paste("length(data_struct)=",length(data_struct)))
		for (ii in 1:length(data_struct)){
			if (is.null(data_struct[[ii]]$blockSize)){
				data_struct[[ii]]$blockSize = matrix(1,1,NCOL(data_struct[[1]]$obs))
			}
			data_struct[[ii]]$blockEnd = matrix(cumsum(data_struct[[ii]]$blockSize),1,NCOL(data_struct[[1]]$obs))
		}
		#print(paste("init blockEnd dims=",dim(data_struct[[1]]$blockEnd)))
		#print(paste("datalen=",NCOL(data_struct[[1]]$obs)))

		
		theta=list(invSigma=array(0,c(dimu,dimu,Kz,Ks)),mu=array(0,c(dimu,Kz,Ks)))
		Ustats=list(card=array(0,c(Kz,Ks)),YY=array(0,c(dimu,dimu,Kz,Ks)),sumY=array(0,c(dimu,Kz,Ks)))	
	}
	
	stateCounts=list(
		N = array(0,c(Kz+1,Kz)),
		Ns = array(0,c(Kz,Ks)),
		uniqueS = array(0,c(Kz,1)),
		M = array(0,c(Kz+1,Kz)),
		barM = array(0,c(Kz+1,Kz)),
		sum_w = array(0,c(1,Kz)),
		sum_w_hap = array(0,c(1,Kz))
	)

	hyperparams=list(
		alpha0_p_kappa0 = 0,
		rho0 = 0,
		rho0_hap = 0
	)
	
	if(model$HMMmodel$type == "HDP"){
		hyperparams$gamma0 = 0
	}else if("alpha0" %in% names(model$obsModel$params)){
		hyperparams$alpha0_p_kappa0 = model$obsModel$params$alpha0
	}
	
	if(model$obsModel$mixtureType == "infinite"){
		 hyperparams$sigma0 = 0
	}else{
		hyperparams$sigma0 = array(0,c(1,Ks))
	}

	numSaves = settings$saveEvery %/% settings$storeEvery
	numStateSeqSaves = settings$saveEvery %/% settings$storeStateSeqEvery
	T = NROW(data_struct[[1]]$obs)
	S=list()
	#CHECK THESE 2 LINES
	stateSeq=array(list(),c(numStateSeqSaves,length(data_struct)))
	#stateSeq[1:numStateSeqSaves,length(data_struct)] = rep.int(list(z=array(0,c(1,T)),s=array(0,c(1,T))),numStateSeqSaves)
	for(n in 1:numStateSeqSaves){
		#stateSeq[[n,length(data_struct)]] = list(z=array(0,c(1,T)),s=array(0,c(1,T)))
		#simplify stateSeq
		stateSeq[[n,length(data_struct)]] = list(z=array(0,T),s=array(0,T))
	}
	
	S$stateSeq = stateSeq
	S$dist_struct[[1:numSaves]] = list(pi_z=array(0,c(Kz,Kz)),pi_init=array(0,c(1,Kz)),pi_s=array(0,c(Kz,Ks)),beta_vec=array(0,c(1,Kz)))

	#S$theta[[1:numSaves]] = theta
	#S$theta[[1:numSaves]] =list(invSigma=array(0,c(dimu,dimu,Kz,Ks)),mu=array(0,c(dimu,Kz,Ks)))
	S$theta[[1:numSaves]] = theta
	
	S$hyperparams[[1:numSaves]] = hyperparams
	S$m = 1
	S$n = 1

	return(list(theta=theta,Ustats=Ustats,stateCounts=stateCounts,hyperparams=hyperparams,data_struct=data_struct,model=model,S=S))
}