GetLPCresidual_WithTwoSignals <- function(wave, waveToComputeLPC, L, shift, order){
  # Used in polarity_reskew.R method
  #
  # INPUTS
  #   wave              original signal
  #   waveToCompute     elliptical high-passed signal
  #   L     = window length (samples) (typ.25ms)
  #   shift = window shift (samples) (typ.5ms)
  #   order = LPC order
  #   gci   = gci position (samples)
  #   type  = vector of voicing decisions (=0 if Unvoiced, =1 if Voiced)
  #   t0    = vector of period values (in samples)
  # 
  # Written by Thomas Drugman, TCTS Lab.
  # 
  
  # Initialize loop parameters
  start = 1
  stop = start+L
  
  # Initialize matricies
  res = zeros(1, length(wave))
  LPCcoeff = zeros((order+1),rnd(length(wave)/shift))
  
  # Do processing
  n=1
  while (stop<length(wave)){
    
    # Wave segment for LPC (no Hanning window)
    segment = waveToComputeLPC[start:stop]
    
    # LPC Coefficients computed from signal "waveToComputeLPC" (s_h)
    A = lpc_modified(segment,order, preemph = FALSE, window = TRUE, mean = FALSE)  
    LPCcoeff[,n] = matrix(A, ncol=1) 
    
    # Segmentation of signal "wave" (s)
    segment2 = wave[start:stop]
    segment2 = segment2 * hanningR(length(segment2), type = "hanning")    
    
    # Inverse filtering of speech signal to produce LP residue signal.
    # Uses coefficients computed from the "waveToCompute" signal
    inv = filter(A, 1, segment2)
    inv = inv * sqrt(sum(segment2^2)/sum(inv^2));
    
    #Compute residual
    res[start:stop]=res[start:stop]+matrix(inv, nrow = 1)
    
    # Increment 
    start=start+shift;
    stop=stop+shift;
    n=n+1;
  }
  
  # Normalise amplitude 
  res=res/max(abs(res))
  
  return(res)
  # Produces same results as MATLAB
  
}