lpcresidual <- function(x,L,shift,order){
  # Function to derive the Linear Prediction residual signal
  #
  # Description
  #  Function to derive the Linear Prediction residual signal
  #
  # Inputs
  #  x               : [samples] [Nx1] Input signal
  #  L               : [samples] [1x1] window length (e.g., 25ms =>  25/1000*fs)
  #  shift           : [samples] [1x1] window shift (e.g., 5ms => 5/1000*fs)
  #  order           : [samples] [1x1] Order of Linear Prediction
  #
  # Outputs
  #  res             : [samples] [Mx1] Linear Prediction residual
  #  LPCcoeff        : [samples] [order+1xM] Linear Prediction coefficients
  #
  # LPC residual reference:
  # http://musicweb.ucsd.edu/~trsmyth/analysis/LPC_residual.html
  # http://www.otolith.com/otolith/olt/lpc.html
  # http://iitg.vlab.co.in/?sub=59&brch=164&sim=616&cnt=1108
  #
  # Theory
  #   The LP residual is the production error e(n).
  #   Obtained as the difference between the predicted speech sample 
  #   s^(n) and the current sample s(n).
  #   Alternatively, e(n) can be obtained by filtering the speech signal
  #   through the reciprocal of H(z), A(z).
  #   Thus, LP residual is obtained via inverse filtering of speech.
  
  # Implementation
 
  # Ensure x is in [Nx1] form
  if (!is.matrix(x)){
    x = as.matrix(x, ncol = 1)
  }
  
  # Initial settings for segmentation
  start = 1
  stop = start+L    
  
  # Allocate space for OUTPUT matricies
  res = zeros(1,length(x))  
  LPCcoeff = zeros((order+1),rnd(length(x)/shift))
 
  # Do processing
  n=1
  system.time({  
  while (stop < length(x)){
    
    # Signal frame  
    segment = x[start:stop]  # segmentation of input signal
   
    # LPC coefficients 
    # Perform LPC analysis. LPC function has in-built Hanning-windowed
    # function. No need to apply window on segment (unlike in MATLAB)
    A = lpc_modified(segment,order, preemph = FALSE, window = TRUE, mean = FALSE)  
    # LPC coefficients matrix
    LPCcoeff[,n] = matrix(A, ncol=1)    # assign predictor coeff values to column "n"
    
    # Hanning-windowed speech segment to be inversed filter
    segment = segment * hanningR(length(segment), type = "hanning")    # apply Hanning window
    
    # Inverse filtering of speech signal to produce LP residue signal
    # Residual signal r(n), obtained by estimating through Linear Prediction (LP) 
    # analysis the coefficients A of an Auto-Regressive model of the speech
    # signal s(n), and by removing the contribution of this spectral
    # envelope (vocal tract) by inverse filtering.
  
    # inv = filter(A,1,segment)                       
    inv=filter(A,1,segment) *sqrt(sum(segment^2)/sum((filter(A,1,segment)) ^2))
    
    # Compute residue signal
    res[start:stop] = res[start:stop]+inv # Overlap and add
    
    # Increment values in while-loop
    start = start+shift
    stop = stop+shift
    n = n+1
  }
  })   
  res = res/max(abs(res)) # Normalise amplitude 
                          # max(abs(res)) equals an integer, as res is a vector
  
  return(list("res" = res, "LPCcoeff" = LPCcoeff))  # Output list$(2 objects)
  # Produces same resutls as MATLAB
  #
  # lpcresidual.test to get some proper results. 
  # This is not how this file is used in pitch_srh.m
  # But this allows us to compare to results to MATLAB
  # y = loadsound("nat_speech.wav")
  # x = y$sound
  # fs = y$fs
  #
  # x = as.matrix(x, ncol = 1)
  #
  # L=400
  # shift = 80
  # order = 12
}