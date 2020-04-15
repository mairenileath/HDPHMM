sample_theta <-function(theta,Ustats,obsModel){
	prior_params = obsModel$params
	source("randiwishart.R")

	if(obsModel$type == 'Gaussian'){
			theta = sample_theta_submodule(theta,Ustats,obsModel$priorType,prior_params)
			return(theta)
	}
}


sample_theta_submodule <-function(theta,Ustats,priorType,prior_params)
{
	nu = prior_params$nu
	nu_delta = prior_params$nu_delta
	store_card = Ustats$card

	#this may not be necessary
	if (NROW(store_card)==1){
		store_card = t(store_card)
	}
	Kz=NROW(store_card)
	Ks=NCOL(store_card)


	if(priorType =='IW-N'){

		invSigma = theta$invSigma
		mu = theta$mu

		store_YY = Ustats$YY
		store_sumY = Ustats$sumY
		
		if(!("numIter" %in% names(prior_params))){
			prior_params$numIter = 50
		}
		numIter = prior_params$numIter

		mu0 = prior_params$mu0
		cholSigma0 = prior_params$cholSigma0
		Lambda0 = ginv(t(prior_params$cholSigma0) %*% prior_params$cholSigma0)
		theta0 = Lambda0 %*% prior_params$mu0


		dimu = NROW(nu_delta)
		#020212 set all emissions to have the same invSigma
		for (n in 1:numIter){
			if (sum(store_card)>0){
				Syy=array(0,c(dim(store_YY)[1],dim(store_YY)[1]))
				for (kz in 1:Kz){                
					for (ks in 1:Ks){
						if (store_card[kz,ks]>0){
								# Given X, Y get sufficient statistics
								Syy = Syy + store_YY[,,kz,ks] - (array(mu[,kz,ks],c(dimu,1)) %*% store_sumY[,kz,ks]) - (array(store_sumY[,kz,ks],c(dimu,1)) %*% mu[,kz,ks]) + store_card[kz,ks] * (array(mu[,kz,ks],c(dimu,1)) %*% mu[,kz,ks])                             
						}
					}
				}
				# Sample Sigma given s.stats
				Sygx = (Syy + t(Syy))/2
				#print(paste("dim(Sygx=",dim(Sygx),sep=""))
				#tmpArr = randiwishart(Sygx + nu_delta,nu+sum(store_card))
				#sqrtSigma = tmpArr[[1]]
				#sqrtinvSigma = tmpArr[[2]]
				#print(paste("param1=",nu+sum(store_card),sep=""))
				#print(paste("param2=",Sygx + nu_delta,sep=""))
				TEST=riwish(nu+sum(store_card),Sygx + nu_delta)
				#TEST=riwish(nu+sum(store_card),ginv(Sygx + nu_delta))
			}else{
				#tmpArr = randiwishart(nu_delta,nu)
				#sqrtSigma = tmpArr[[1]]
				#sqrtinvSigma = tmpArr[[2]]    
				TEST=riwish(nu,nu_delta)
				#TEST=riwish(nu,ginv(nu_delta))
			}

			#print(class(sqrtinvSigma))
			#print(paste("sqrtinvSigma=",sqrtinvSigma,sep=""))
			#print(paste("sqrtSigma=",sqrtSigma,sep=""))
			#print(paste("TEST=",TEST,sep=""))
			#print(dim(invSigma))

			is = ginv(TEST)
			for (kz in 1:Kz){               
				for (ks in 1:Ks){          
					#invSigma[,,kz,ks] = t(sqrtinvSigma) %*% sqrtinvSigma
					#invSigma[,,kz,ks] = 1/TEST #this works for the 1-D case, but not 2-D
					#invSigma[,,kz,ks] = TEST
					invSigma[,,kz,ks] = is
					#invSigma[,,kz,ks] = 1/sqrt(TEST %*% t(TEST))
				}
			}
			
			#print(paste("dim(Lambda0)=",dim(Lambda0),sep=""))
			#print(paste("dim(invSigma)=",dim(invSigma),sep=""))
			#print(paste("dim(theta0)=",dim(theta0),sep=""))
			#print(paste("dim(store_card)=",dim(store_card),sep=""))
			
#             for kz=1:Kz                
#                 for ks=1:Ks
#                     if store_card(kz,ks)>0  #**
#                         # Sample mu given A and Sigma
#                         Sigma_n = inv(Lambda0 + store_card(kz,ks)*invSigma(:,:,kz,ks));
#                         mu_n = Sigma_n*(theta0 + invSigma(:,:,kz,ks)*store_sumY(:,kz,ks));
#                         mu(:,kz,ks) = mu_n + chol(Sigma_n)'*randn(dimu,1);
#                     else
#                         mu(:,kz,ks) = mu0 + cholSigma0'*randn(dimu,1);
#                     end
# 				end
#             end

			#force BAF to take values symmetrical around 0.5 and consider logRs for haplotype blocks together
			 for (kz in 1:(Kz %/% 2)){           
				for (ks in 1:Ks){
					if (store_card[2*kz-1,ks]>0 | store_card[2*kz,ks]>0){
						# Sample mu given A and Sigma
						Sigma_n = ginv(Lambda0 + (store_card[2*kz-1,ks] + store_card[2*kz,ks])*invSigma[,,2*kz,ks])
						#020312 treat BAF and logR differently
						if(dim(store_sumY)[1]==1){
							combined_store_sumY = array(store_sumY[,2*kz-1,ks]+store_card[2*kz,ks]-store_sumY[,2*kz,ks],c(1,1))
						}else{
							combined_store_sumY = array(0,c(2,1))
							combined_store_sumY[1,1] = store_sumY[1,2*kz-1,ks]+store_card[2*kz,ks]-store_sumY[1,2*kz,ks]
							combined_store_sumY[2,1] = store_sumY[2,2*kz-1,ks]+store_sumY[2,2*kz,ks]
						}
						#print(combined_store_sumY/(store_card[2*kz-1,ks]+store_card[2*kz,ks]))
						mu_n = Sigma_n %*% (theta0 + invSigma[,,2*kz,ks] %*% combined_store_sumY)
						#mu_n = Sigma_n %*% (theta0 + invSigma[,,2*kz,ks] %*% array(store_sumY[,2*kz-1,ks]+store_card[2*kz,ks]-store_sumY[,2*kz,ks],c(dimu,1)))
						
						#print(paste("theta0=",dim(theta0),sep=""))
						#print(paste("Sigma_n=",Sigma_n,sep=""))
						#print(paste("chol(Sigma_n)=",dim(chol(Sigma_n)),sep=""))
						#print(paste("mu_n=",dim(Sigma_n),sep=""))
						#print(paste("dimu=",dimu,sep=""))
						
						mu[,kz*2-1,ks] = mu_n + (t(chol(Sigma_n)) %*% array(rnorm(dimu),c(dimu,1)))
						#mu(:,kz*2,ks) = 1.0 - mu(:,kz*2-1,ks);
						mu[1,kz*2,ks] = 1.0 - mu[1,kz*2-1,ks] #BAFs should be 1-x for haplotype blocks
						if(dim(mu)[1]>1){
							mu[2,kz*2,ks] = mu[2,kz*2-1,ks] #logRs should be the same for haplotype blocks
						}
					}else{
						mu[,kz*2-1,ks] = mu0 + t(cholSigma0) %*% array(rnorm(dimu),c(dimu,1))
						#mu(:,kz*2,ks) = 1.0 - mu(:,kz*2-1,ks);
						mu[1,kz*2,ks] = 1.0 - mu[1,kz*2-1,ks] #BAFs should be 1-x for haplotype blocks
						if(dim(mu)[1]>1){
							mu[2,kz*2,ks] = mu[2,kz*2-1,ks] #logRs should be the same for haplotype blocks
						}
					}
				}
			}
		}

		theta$invSigma = invSigma
		theta$mu =  mu
	}
	return(theta)
}