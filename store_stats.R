# function S = store_stats(S,n,settings,stateSeq_n,dist_struct_n,theta_n,hyperparams_n)   
# Store statistics into structure S and save if rem(n,saveEvery) = 0

store_stats <- function(S,n,settings,stateSeq_n,dist_struct_n,theta_n,hyperparams_n){

	# If we are at a storing iteration:
	if (n %% settings$storeEvery==0 & n>=settings$saveMin){
		# And if we are at a mode-sequence storing iteration:
		if (n %% settings$storeStateSeqEvery==0){
			# Store all sampled mode sequences:
			for (ii in 1:length(stateSeq_n)){
				S$stateSeq[S$n,ii] = stateSeq_n[ii]
			}
			# Increment counter for the mode-sequence store variable:
			S$n = S$n + 1;
		}
		# Store all sampled model parameters:
		S$dist_struct[[S$m]] = dist_struct_n
		S$theta[[S$m]] = theta_n
		S$hyperparams[[S$m]] = hyperparams_n
		# Increment counter for the regular store variable:
		S$m = S$m + 1

	}

	# If we are at a saving iteration:
	if (n %% settings$saveEvery==0){

		# Save stats to specified directory:
		if("filename" %in% names(settings)){
			filename = paste(settings$saveDir,'/',settings$filename,'iter',toString(n),'trial',toString(settings$trial),".Rdata",sep="")    # create filename for current iteration
		}else{
			filename = paste(settings$saveDir,'/HDPHMMDPstats','iter',toString(n),'trial',toString(settings$trial),".Rdata",sep="")    # create filename for current iteration
		}

		save(S,file=filename) # save current statistics

		# Reset S counter variables:
		S$m = 1
		S$n = 1

		#print(paste('Iteration: ',toString(n),sep=""))

	}
	print(paste('Iteration: ',toString(n),sep=""))
	
	return(S)
}