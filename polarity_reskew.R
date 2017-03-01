polarity_reskew <- function(s, fs){
  
  # RESKEW is an efficient technique to determine the polarity of the speech
  # signal.
  #
  # Octave compatible
  #
  # Description
  # The Residual Excitation Skewness (RESKEW) method is described in [1].
  # This algorithm determines the polarity of the speech signal by inspecting
  # the skewness of two residual excitation signals. Its advantages (shown in
  # [1]) are: i) its high performance; ii) its robustness to an additive
  # noise; iii) the fact that it does not depend any voicing or pitch
  # estimate.
  #
  #
  # Inputs
  #  s               : [samples] [Nx1] input signal (speech signal)
  #  fs              : [Hz]      [1x1] sampling frequency
  #
  # Outputs
  #  polarity              : [1x1] the speech polarity
  #
  # Example
  #  polarity = polarity_reskew(s,fs);
  #
  # References
  #  [1] T.Drugman, Residual Excitation Skewness for Automatic Speech Polarity
  #  Detection, IEEE Signal Processing Letters, vol. 20, issue 4, pp.
  #  387-390, 2013.
  #  Publication available at the following link:
  #  http://tcts.fpms.ac.be/~drugman/files/SPL-Polarity.pdf
  #
  # Copyright (c) 2013 University of Mons, FNRS
  #
  # License
  #  This code is a part of the GLOAT toolbox with the following
  #  licence:
  #  This program is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.
  #
  # This function is part of the Covarep project: http://covarep.github.io/covarep
  # 
  # Author 
  #  Thomas Drugman thomas.drugman@umons.ac.be
  
  ###########################################################################
  # Computation of high-passed version s_h of signal s using elliptic filter
  #
  # Filter parameters
  Wp = 480/(fs/2)   # Pass-band edge
  Ws = 500/(fs/2)   # Stop-band edge
  Rp = 3            # Allowable decibels of ripple in the pass band
  Rs = 60           # Minimum attenuation in the stop band in dB
  
  # Compute discrete elliptic filter order and cutoff
  filterlist = ellipordR(Wp,Ws,Rp,Rs)
  n = filterlist$n                     # Elliptic filter order
  Wn = filterlist$Wc                   # Elliptic filter cutoff frequency
  
  # Compute elliptic filter with new order and cutoff frequency
  elliplist =  ellip(n, Rp, Rs, Wn, type = "high")
  b_h = elliplist$b           # Moving average (MA) polynomial coefficients
  a_h = elliplist$a           # Autoregressive (AR) polynomial coefficients
  
  # Pass signal s through a high-pass elliptic filter to produce s_h
  # s_h spectral components are kept beyond Wn (cutoff frequency)
  s_h = filtfilt(b_h, a_h, s)   # Filtered signal, length(s_h) = length(s)
  ###########################################################################
  #
  # The next step is to extract the two residual signals
  # 1. residual signal obtained by performing LPC analysis on signal s,
  #    then removing the spectral envelope fo the vocal tract through inverse
  #    filtering, resulting in the residual signal [res2].
  #    Used function: lpcresidual
  #
  # 2. residual signal in the form of the glottal flow derivative signal.
  #    Obtained by Computing the LPC coefficients of the high-passed signal
  #    s_h and using those coefficients to remove the spectral envelope from 
  #    the original signal s (using inverse filtering). This results in the
  #    glottal flow derivative [res]
  #    Used function: GetLPCresidual_WithTwoSignals
  #
  ###########################################################################
  # Residual signals computation
  # Two excitation signals are considered in this work:
  #
  # First residual signal r(n)
  # Obtained by estimating thourgh Linear Prediction the coefficients 
  # of an Auto-Regressive model of speech signal s(n) (AR model has only poles ) 
  # and removing the contribution of the vocal tract spectral envelope
  # through inverse filtering (procedure of lpcresidual.R method)
  lpcResidualList = lpcresidual(s, rnd(25/1000*fs), rnd(5/1000*fs), (rnd(fs/1000)+2))
  res2 = lpcResidualList$res
  # View signal:
  # a = 1:length(res2)
  # plot(a, res2, type = "l")
  
  # Second residual signal
  # This signal is a rough approximation of the glottal flow derivative
  # s_h is a high passed version of signal s
  # This high-pass filtering is carried out with an elliptic filter (above)
  #
  # GetLPCresidual_WithTwoSignals:
  # Computes LPC coefficients from s_h signal
  # The estimate of the glottal flow derivative is then achieved by inverse
  # filtering the original speech signal s(n) using the LPC coefficients of
  # the high-passed filtered signal, s_h. 
  res = GetLPCresidual_WithTwoSignals(s, s_h, rnd(25/1000*fs), rnd(5/1000*fs), (rnd(fs/1000)+2));
  # View signal:
  # plot(a, res, type = "l")
  # This signal presents negative peaks at GCIs
  ###########################################################################
  
  if (FALSE){
    a = 1:length(res2)
    plot(a, res, type = "l")
    par(new=T)
    plot(a, res2, type = "l",col="red",
         main = "Red = res, Black = derivative")
  }
  
  # Compute skewness of residual signals
  # Skewness requires {e1071} package
  Val1 = skewness(res)
  Val2 = skewness(res2)
  
  # Return 1, 0, or -1 if number is +ve, zero, or -ve, respectively
  polarity = sign(Val2-Val1)
  
  return(polarity)
  # Produces same results as MATLAB (slight rounding difference)
  
}