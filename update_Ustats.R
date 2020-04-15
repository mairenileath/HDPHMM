update_Ustats <-function (data_struct,INDS,stateCounts,obsModelType){

	Ns = stateCounts$Ns

	Kz = NROW(Ns)
	Ks = NCOL(Ns)

	if(obsModelType == "Gaussian"){

			unique_z = which(rowSums(Ns)>0)

			dimu = NROW(data_struct[[1]]$obs)
			
			#print(paste("Kz=",Kz,", Ks=",Ks,", dimu=",dimu,sep=""))

			store_YY = array(0,c(dimu,dimu,Kz,Ks))
			store_sumY = array(0,c(dimu,Kz,Ks))

			for (ii in 1:length(data_struct)){

				u = data_struct[[ii]]$obs

				for (kz in unique_z){
					unique_s_for_z = which(Ns[kz,]>0)
					for (ks in unique_s_for_z){
						#print(paste("INDS=",INDS[[ii]]$obsIndzs[[kz,ks]],sep=""))
						obsInd = INDS[[ii]]$obsIndzs[[kz,ks]]$inds[1:INDS[[ii]]$obsIndzs[[kz,ks]]$tot]
						#print(paste("dim(obsInd)=",dim(obsInd),sep=""))
						#print(paste("length(obsInd)=",length(obsInd),sep=""))
						#print(paste("dim(obsInd)=",dim(obsInd),sep=""))
						#if(is.null(dim(obsInd))){
						if(length(obsInd)==1){
							#store_YY[,,kz,ks] = store_YY[,,kz,ks] + sum(u[,obsInd]^2)
							#store_sumY[,kz,ks] = store_sumY[,kz,ks] + sum(u[,obsInd])
							store_sumY[,kz,ks] = store_sumY[,kz,ks] + u[,obsInd]
							store_YY[,,kz,ks] = store_YY[,,kz,ks] + array(u[,obsInd],c(dimu,1)) %*% u[,obsInd]
						}else if(dim(u)[1]>1){
							#store_YY[,,kz,ks] = store_YY[,,kz,ks] + rowSums(u[,obsInd]^2)
							store_sumY[,kz,ks] = store_sumY[,kz,ks] + rowSums(u[,obsInd])
							store_YY[,,kz,ks] = store_YY[,,kz,ks] + u[,obsInd] %*% t(u[,obsInd])
						}else{
							store_sumY[,kz,ks] = store_sumY[,kz,ks] + sum(u[,obsInd])
							store_YY[,,kz,ks] = store_YY[,,kz,ks] + array(u[,obsInd],c(1,length(obsInd))) %*% array(u[,obsInd],c(length(obsInd),1))						
						}
						
					}
				}
			}

			Ustats=list(card = Ns,YY = store_YY,sumY = store_sumY)
	}
        
	return(Ustats)
}