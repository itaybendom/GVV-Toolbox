fftn <- function (signal, N){
  #FFT Discrete Fourier transform.
  #   FFT(X) is the discrete Fourier transform (DFT) of vector X.  For
  #   matrices, the FFT operation is applied to each column. For N-D
  #   arrays, the FFT operation operates on the first non-singleton
  #   dimension.
  #
  #   When X is a vector, the value computed and returned by fft is the unnormalized univariate 
  #   DFT of the sequence of values in z. 
  #   Specifically, y <- fft(z) returns:
  #
  #                           y[h] = sum_{k=1}^n z[k]*exp(-2*pi*1i*(k-1)*(h-1)/n)
  #
  #   FFT(X,N) is the N-point FFT, padded with zeros if X has less than N points.
  
  #   Make sure signal is in vector/matrix form
  if (!is.matrix(signal)){
    signal = as.matrix(signal)
  }
  
  #   Zero padding to length N
  signalzp = c(signal, rep(0, (N-length(signal))))
  #   Perform FFT analysis on padded signal
  signal = fft(signalzp)
  
}