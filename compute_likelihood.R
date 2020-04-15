compute_likelihood <- function(data_struct,theta,obsModelType,Kz,Ks){
	if(obsModelType == "Gaussian"){        
        invSigma = theta$invSigma
        mu = theta$mu
        
        #T = NROW(data_struct$obs)
        #dimu = NCOL(data_struct$obs)
        T = NCOL(data_struct$obs)
        dimu = NROW(data_struct$obs)        
        
        log_likelihood = array(0,c(Kz,Ks,T))
        #print(paste("invSigma=",invSigma,sep=""))
        for (kz in 1:Kz){
            for (ks in 1:Ks){               
                cholinvSigma = chol(invSigma[,,kz,ks])
                #print(paste("cholinvsigma dims=",dim(cholinvSigma),sep=""))
                dcholinvSigma = diag(cholinvSigma)
                u = cholinvSigma %*% (data_struct$obs - array(mu[,kz,ks],c(dimu,T)))
                log_likelihood[kz,ks,] = -0.5*colSums(u^2) + sum(log(dcholinvSigma))
			}
		}
        normalizer = sapply(1:T,function(i,l) max(l[,,i],na.rm=T),l=log_likelihood)
        #log_likelihood = log_likelihood - array(rep(normalizer,Kz*Ks),c(Kz,Ks,T))
        log_likelihood = log_likelihood - as.vector(t(matrix(rep(normalizer,Kz*Ks),nrow=T,ncol=Kz*Ks)))
        likelihood = exp(log_likelihood)
        #print(paste("min/max likelihood=",min(likelihood),",",max(likelihood),sep=""))
        
        normalizer = normalizer - (dimu/2)*log(2*pi)
	}
       
    return(list(likelihood=likelihood,normalizer=normalizer))
}