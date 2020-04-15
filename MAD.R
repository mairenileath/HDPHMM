
#Get mad SD (based on KL code)
getMad <- function(x,k=25){

  if(length(x) > 2*k){
  
  	filtLengdeM1 <- 2*k
  	lengdeX <- length(x)
  	runMedian <- rep(0,lengdeX)
  	for(i in 1:(lengdeX-filtLengdeM1)){
  		runMedian[i+k] <- median(x[i:(i+filtLengdeM1)])
  	}
  	for(i in 1:k){
  		runMedian[i] <- median(x[1:(i+k)])
  		runMedian[lengdeX-i+1] <- median(x[(lengdeX-i-k+1):lengdeX])
  	}
  }else{
    runMedian <- median(x)
  }	
	dif <- x-runMedian
	SD <- mad(dif)
	
	return(SD)
}

#Perform MAD winsorization:
madWins <- function(x,tau,k){
	xhat <- medianFilter(x,k)
	d <- x-xhat
	SD <- mad(d)
	z <- tau*SD
	xwin <- xhat + psi(d, z)
	outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
	return(list(ywin=xwin,sdev=SD,outliers=outliers))
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)

}



psi <- function(x,z){
 xwin <- x
 xwin[x < -z] <- -z
 xwin[x > z] <- z
 return(xwin)
}