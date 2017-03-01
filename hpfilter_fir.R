hpfilter_fir <- function (Fstop,Fpass,fs,N){
  
  # Function used in iaif.m in MATLAB
  # Funciton call:
  # function B = hpfilter_fir(Fstop,Fpass,fs,N)
  #
  # FIR least-squares Highpass filter design using the FIRLS function
  #
  # Tuomo Raitio (MATLAB)
  # 10.7.2012
  
  # Calculate the coefficients using the FIRLS function.
  frequencies = c(0, Fstop, Fpass, fs/2)/(fs/2)
  pass = c(0, 0, 1, 1)
  weight = matrix(c(1, 1), nrow = 1)
  
  b  = firls1(N, frequencies, pass, weight)
  return(b)
  
}