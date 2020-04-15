HDPHMMDPinference<-function(data_struct,model,settings,restart){
	library(MCMCpack)
	#library(MASS)
	library(Matrix)
	source("sample_hyperparams_init.R")
	#source("sample_dist.R")
	source("initializeStructs.R")
	source("sample_theta.R")
	source("sample_zs.R")
	source("sample_zs_init.R")
	source("compute_likelihood.R")
	source("backwards_message_vec.R")
	source("update_Ustats.R")
	#source("sample_tables.R")	
	source("randnumtable.R")
	#source("sample_barM.R")
	#source("sample_hyperparams.R")
	source("gibbs_conparam.R")
	source("store_stats.R")
	
	source("sample_barM_hap.R")
	source("sample_hyperparams_hap.R")
	source("sample_tables_hap.R")
	source("sample_dist_hap.R")
	
	trial = settings$trial
	if (!("saveMin" %in% names(settings))){
	    settings$saveMin = 1
	}
	resample_kappa = settings$resample_kappa
	resample_haplotypeSwitching = settings$resample_haplotypeSwitching
	Kz = settings$Kz
	Niter = settings$Niter
	
	print(paste("resample_haplotypeSwitching in HDPHMMInference=",resample_haplotypeSwitching,sep=""))
	print(paste('Trial: ',toString(trial),sep=""))
	
	if(restart==1){
		n = settings$lastSave
        
        # Build necessary structures and clear structures that exist as
        # part of the saved statistics:
        tempList = initializeStructs(model,data_struct,settings)
        data_struct = tempList[[5]]
        model = tempList[[6]]
        #rm(ls(theta,Ustats,stateCounts,hyperparams,S))
        
        #print(paste("saveDir=",settings$saveDir,sep=""))
        #print(paste("lastSave=",settings$lastSave,sep=""))
        # Load the last saved statistics structure S:
        lastSaveDir = gsub(paste("iters",Niter,sep=""),paste("iters",n,sep=""),settings$saveDir)
        if("filename" %in% names(settings)){
            filename = paste(lastSaveDir,'/',settings$filename,'iter',toString(n),'trial',toString(settings$trial),".Rdata",sep="")
        }else{
            filename = paste(lastSaveDir,'/HDPHMMDPstats','iter',toString(n),'trial',toString(settings$trial),".Rdata",sep="")
		}
        #print(paste("filename=",filename,sep=""))
        load(filename)
        #print("OK")
        
        obsModel = model$obsModel  # structure containing the observation model parameters
        obsModelType = obsModel$type   # type of emissions including Gaussian, multinomial, AR, and SLDS.
        HMMhyperparams = model$HMMmodel$params # hyperparameter structure for the HMM parameters
        HMMmodelType = model$HMMmodel$type # type of HMM including finite and HDP
        
        # Set new save counter variables to 1:
        S$m = 1
        S$n = 1
        
        # Grab out the last saved statistics from the S structure:
        numSaves = settings$saveEvery %/% settings$storeEvery
        numStateSeqSaves = settings$saveEvery %/% settings$storeStateSeqEvery
        theta = S$theta[[numSaves]]
        
        #print(paste("HDPHMMinference theta$invSigma=",theta$invSigma,sep=""))
        #print(paste("HDPHMMinference dim(theta$invSigma)=",dim(theta$invSigma),sep=""))
        
        dist_struct = S$dist_struct[[numSaves]]
        hyperparams = S$hyperparams[[numSaves]]
        stateSeq = S$stateSeq[[numStateSeqSaves]]
        
        # Set the new starting iteration to be lastSave + 1:
        n_start = n + 1		
	}else{
		# Set the starting iteration:
		n_start = 1

		# Build initial structures for parameters and sufficient statistics:
		tempList = initializeStructs(model,data_struct,settings)
		theta=tempList$theta
		Ustats=tempList$Ustats
		stateCounts=tempList$stateCounts
		hyperparams=tempList$hyperparams
		data_struct=tempList$data_struct
		#model=tempList$model
		S=tempList$S
		

		obsModel = model$obsModel  # structure containing the observation model parameters
		obsModelType = obsModel$type   # type of emissions including Gaussian, multinomial, AR, and SLDS.
		HMMhyperparams = model$HMMmodel$params # hyperparameter structure for the HMM parameters
		HMMmodelType = model$HMMmodel$type # type of HMM including finite and HDP

		# Resample concentration parameters:
		hyperparams = sample_hyperparams_init(stateCounts,hyperparams,HMMhyperparams,HMMmodelType,resample_kappa)

		# Sample the transition distributions pi_z, initial distribution
		# pi_init, emission weights pi_s, and global transition distribution beta
		# (only if HDP-HMM) from the priors on these distributions:
		dist_struct = sample_dist(stateCounts,hyperparams,model)

		# If the optional 'formZInit' option has been added to the settings
		# structure, then form an initial mode sequence in one of two ways.  If
		# 'z_init' is a field of data_struct, then the specified initial
		# sequence will be used. Otherwise, the sequence will be sampled from
		# the prior.
		if("formZInit" %in% names(settings)){
			if(settings$formZInit == 1){
				tempList = sample_zs_init(data_struct,dist_struct,obsModelType)
				stateSeq = tempList$stateSeq
				INDS = tempList$INDS
				stateCounts = tempList$stateCounts			
				Ustats = update_Ustats(data_struct,INDS,stateCounts,obsModelType)
				print('Forming initial z using specified z_init or sampling from the prior using whatever fixed data is available')
			}
		}else if (length(data_struct)>length(data_struct[1]$test_cases)){
			print('Do you want z_init set to truth for extra datasets?  If so, make setttings.formZinit =1 ')
		}

		# Sample emission params theta_{z,s}'s. If the above 'formZInit' option
		# was not utilized, the initial parameters will just be drawn from the
		# prior.
		theta = sample_theta(theta,Ustats,obsModel)

		# Create directory in which to save files if it does not currently exist
		if (!file.exists(settings$saveDir)){
			dir.create(settings$saveDir)
		}

		# Save initial statistics and settings for this trial:
		if ("filename" %in% names(settings)){
			settings_filename = paste(settings$saveDir,'/',settings$filename,'_info4trial',toString(trial),".Rdata",sep="")    # create filename for current iteration
			init_stats_filename = paste(settings$saveDir,'/',settings$filename,'initialStats_trial',toString(trial),".Rdata",sep="")    # create filename for current iteration
		}else{
			settings_filename = paste(settings$saveDir,'/info4trial',toString(trial),".Rdata",sep="")  # create filename for current iteration
			init_stats_filename = paste(settings$saveDir,'/initialStats_trial',toString(trial),".Rdata",sep="")   # create filename for current iteration
		}
		stats1<-list(data_struct=data_struct,settings=settings,model=model)
		save(stats1,file=settings_filename) # save current statistics
		#save(data_struct,settings,model,file=settings_filename) # save current statistics
		stats2=list(dist_struct=dist_struct,theta=theta,hyperparams=hyperparams)
		save(stats2,file=init_stats_filename) # save current statistics
	}
		
	########## Run Sampler ##########
	for (n in n_start:Niter){

		# Sample z and s sequences given data, transition distributions,
		# HMM-state-specific mixture weights, and emission parameters:

		# Block sample (z_{1:T},s_{1:T})|y_{1:T}
		temp_list = sample_zs(data_struct,dist_struct,theta,obsModelType)
		print("sample_zs OK")
		#stateSeq = temp_list[[1]]
		#INDS = temp_list[[2]]
		#stateCounts = temp_list[[3]]
		stateSeq = temp_list$stateSeq
		INDS = temp_list$INDS
		stateCounts = temp_list$stateCounts
		
		print(paste("unique states=",c(unique(stateSeq[[1]]$z))))
		print(paste("unique segmented BAF values=",c(theta$mu[cbind(1,unique(stateSeq[[1]]$z),1)]),sep=""))
		print(paste("unique segmented log(R) values=",c(theta$mu[cbind(2,unique(stateSeq[[1]]$z),1)]),sep=""))
		
		# Create sufficient statistics:
		#print("calling Ustats")
		Ustats = update_Ustats(data_struct,INDS,stateCounts,obsModelType)

		print("Ustats OK")
		# Based on mode sequence assignment, sample how many tables in each
		# restaurant are serving each of the selected dishes. Also sample the
		# dish override variables:
		stateCounts = sample_tables(stateCounts,hyperparams,dist_struct$beta_vec,Kz)

		print("stateCounts OK")

		# Sample the transition distributions pi_z, initial distribution
		# pi_init, emission weights pi_s, and avg transition distribution beta:
		dist_struct = sample_dist(stateCounts,hyperparams,model)

		print("sampleDist OK")

		# Sample theta_{z,s}'s conditioned on the z and s sequences and the
		# sufficient statistics structure Ustats:
		theta = sample_theta(theta,Ustats,obsModel)

		print("sampleTheta OK")

		# Resample concentration parameters:
		hyperparams = sample_hyperparams(stateCounts,hyperparams,HMMhyperparams,HMMmodelType,resample_kappa,resample_haplotypeSwitching)

		print("sample_hyperparams OK")

		# Build and save stats structure:
		S = store_stats(S,n,settings,stateSeq,dist_struct,theta,hyperparams)
		#print(paste("Unique segmented BAF values=",c(as.vector(S$theta[[1]]$mu[cbind(1,unique(S$stateSeq[1,1][[1]]$z),1)])),sep=""))

		print("store_stats OK")
	}
}