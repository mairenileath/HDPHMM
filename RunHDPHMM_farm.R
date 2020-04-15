args=(commandArgs(TRUE)) 
samplename<-toString(args[1])
chr<-as.integer(args[2])
max_no_states<-as.integer(args[3])
no_iterations<-as.integer(args[4])
c<-as.integer(args[5])
c_hap<-as.integer(args[6])
conc_param_states<-as.integer(args[7])
conc_param_transitions<-as.integer(args[8])
resample_haplotypeSwitching = T
if(length(args)>=9){
	resample_haplotypeSwitching<-as.logical(args[9])
}

setwd("/lustre/scratch103/sanger/dw9/HDPHMM")
source("HDPHMMDPinference_farm.R")
source("MAD.R")

#settings=list(Kz=30,Ks=1,Niter = 2000,resample_kappa = 1,seqSampleEvery = 100,saveEvery = 100,storeEvery = 100,storeStateSeqEvery = 100,trial=1,formZInit =1)
#settings=list(Kz=30,Ks=1,Niter = 1000,resample_kappa = 1,seqSampleEvery = 100,saveEvery = 100,storeEvery = 100,storeStateSeqEvery = 100,trial=1,formZInit =1)
#settings=list(Kz=30,Ks=1,Niter = 5,resample_kappa = 1,seqSampleEvery = 1,saveEvery = 1,storeEvery = 1,storeStateSeqEvery = 1,trial=1,formZInit =1)
#settings=list(Kz=30,Ks=1,Niter = 20,resample_kappa = 1,seqSampleEvery = 5,saveEvery = 5,storeEvery = 5,storeStateSeqEvery = 5,trial=1,formZInit =1)
#settings=list(Kz=30,Ks=1,Niter = 100,resample_kappa = 1,seqSampleEvery = 10,saveEvery = 10,storeEvery = 10,storeStateSeqEvery = 10,trial=1,formZInit =1)
#settings=list(Kz=30,Ks=1,Niter = 20,resample_kappa = 1,saveEvery = 5,storeEvery = 5,storeStateSeqEvery = 5,formZInit =1)

saveDir = paste("/lustre/scratch103/sanger/dw9/HDPHMM/outputs_maxz",max_no_states,"_iters",no_iterations,"_c",c,"_chap",c_hap,"_concParamStates",conc_param_states,"_concParamTransitions",conc_param_transitions,"_resampleHap",resample_haplotypeSwitching,sep="")
if(!file.exists(saveDir)){
	dir.create(saveDir)
}
#settings=list(Kz=30,Ks=1,Niter = 25,resample_kappa = 1,saveEvery = 5,storeEvery = 5,storeStateSeqEvery = 5,formZInit =1,saveDir="C:/R_code/HDPHMM/TEST_OUTPUT")

settings=list(Kz=max_no_states,Ks=1,Niter = no_iterations,resample_kappa = 1,resample_haplotypeSwitching = resample_haplotypeSwitching, saveEvery = 5,storeEvery = 5,storeStateSeqEvery = 5,formZInit =1,saveDir=saveDir)
##moved down 060312
#hmmParams=list(a_alpha=1,b_alpha=0.01,a_gamma=1,b_gamma=1/conc_param,a_sigma=1,b_sigma=0.01,c=c,c_hap=c_hap,d=1)
#hmmModel=list(type="HDP",params=hmmParams)

##CREATE params
#params=list(mu0=array(c(0.5,0),c(2,1)),cholSigma0=0.5*diag(2),nu=4,nu_delta=2*diag(2))#2-D version
##params=list(mu0=0.5,cholSigma0=0.5,nu=3,nu_delta=1)
#obsModel=list(params=params,type="Gaussian",mixtureType="infinite",priorType="IW-N")
#model=list(obsModel=obsModel,HMMmodel=hmmModel)

	
chr_name=toString(chr)
if(chr==23){
	chr_name="X"
}
settings$trial=chr
BAFdata=read.table(paste("/nfs/cancer_archive02/dw9/haplotyped_data/",samplename,"_filtered_0109_ALL_noRefLoci/",samplename,"_chr",toString(chr),"_heterozygousMutBAFs_haplotyped.txt",sep=""),sep="\t",header=T)

BAFs=BAFdata[,3]
pos=BAFdata[,2]/1000000
#150212 - winsorise data
BAFs = madWins(BAFs,2.5,25)[[1]]

logRdata=read.table(paste("/nfs/cancer_archive02/dw9/haplotyped_data/",samplename,"_filtered_0109_ALL_noRefLoci/",samplename,"_chr",toString(chr),"_heterozygousLogRs.txt",sep=""),sep="\t",header=T)

logRs=logRdata[,3]
#150212 - winsorise data
logRs = madWins(logRs,2.5,25)[[1]]


if(is.na(c)){
	c = length(BAFs)
	c_hap = c / c_hap #if c is NA, c_hap is assumed to be the expected length of a haplotype block
}

#CREATE params
hmmParams=list(a_alpha=1,b_alpha=conc_param_transitions/(c+c_hap+1),a_gamma=1,b_gamma=conc_param_states/settings$Kz,a_sigma=1,b_sigma=0.01,c=c,c_hap=c_hap,d=1)
hmmModel=list(type="HDP",params=hmmParams)
params=list(mu0=array(c(0.5,0),c(2,1)),cholSigma0=0.5*diag(2),nu=4,nu_delta=2*diag(2))#2-D version
#params=list(mu0=0.5,cholSigma0=0.5,nu=3,nu_delta=1)
obsModel=list(params=params,type="Gaussian",mixtureType="infinite",priorType="IW-N")
model=list(obsModel=obsModel,HMMmodel=hmmModel)

#data1=list(obs=array(BAFs,c(1,length(BAFs))))#1-D
data1=list(obs=rbind(BAFs,logRs))
data_struct=list(data1)

if(length(args)>=10){
	settings$lastSave=as.integer(args[10])
	HDPHMMDPinference(data_struct,model,settings,1)
}else{
	HDPHMMDPinference(data_struct,model,settings,0)
}	

filename=paste(settings$saveDir,"/HDPHMMDPstatsiter",settings$Niter,"trial",settings$trial,".Rdata",sep="")
load(filename)
print(unique(S$stateSeq[1,1][[1]]$z))
print(S$theta[[1]]$mu)

png(filename=paste(settings$saveDir,"/chr",chr_name,"_segmented.png",sep=""),width=1500,height=1000)
par(mfrow=c(2,1))
segmentedBAFs = S$theta[[1]]$mu[cbind(1,S$stateSeq[1,1][[1]]$z,1)]
plot(pos,BAFs,type="p",col="red",pch=".",ylim=c(0,1),main=paste(samplename," chr ",chr_name,sep=""),xlab="position (Mb)",ylab="BAF")
points(pos,segmentedBAFs,col="green",pch=".",cex=2)

segmentedlogRs = S$theta[[1]]$mu[cbind(2,S$stateSeq[1,1][[1]]$z,1)]
plot(pos,logRs,type="p",col="red",pch=".",ylim=c(-1,1),main=paste(samplename," chr ",chr_name,sep=""),xlab="position (Mb)",ylab="log(R)")
points(pos,segmentedlogRs,col="green",pch=".",cex=2)
dev.off()

outdata=cbind(BAFdata[,1:2],segmentedBAFs)
names(outdata)[3]=samplename
write.table(outdata,file=paste(settings$saveDir,"/chr",chr_name,"_segmentedBAFs.txt",sep=""),sep="\t",row.names=F,quote=F)
outdata=cbind(logRdata[,1:2],segmentedlogRs)
names(outdata)[3]=samplename
write.table(outdata,file=paste(settings$saveDir,"/chr",chr_name,"_segmentedLogRs.txt",sep=""),sep="\t",row.names=F,quote=F)
write.table(data.frame(state=S$stateSeq[1,1][[1]]$z),file=paste(settings$saveDir,"/chr",chr_name,"_states.txt",sep=""),sep="\t",row.names=F,quote=F)

q(save="no")