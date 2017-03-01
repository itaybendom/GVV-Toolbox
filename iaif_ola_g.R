# IAIF Glottal Inverse Filtering
#  [g,dg,a,ag] = iaif_ola(x,fs,winLen,winShift,p_vt,p_gl,d,hpfilt)
#
# Description
#  This function estimates vocal tract linear prediction coefficients and
#  the glottal volume velocity waveform from a speech signal frame using
#  Iterative Adaptive Inverse Filtering (IAIF) method. Analysis is carried
#  out on a fixed frame basis and waveforms are generated using overlap and
#  add.
#
# Inputs
#  x       : Speech signal frame [samples]
#  fs      : Sampling frequency [Hz]
#  winLen  : Window Length [samples] 
#  winShift: Window Shift [samples]
#  p_vt    : Order of LPC analysis for vocal tract
#  p_gl    : Order of LPC analysis for glottal source
#  d       : Leaky integration coefficient (e.g. 0.99)
#  hpfilt  : High-pass filter flag (0: do not apply, 1...N: apply N times)
#
# Outputs
#  g       : Glottal volume velocity waveform
#
# Notes
#  This function does not perform pitch synchronous analysis. This ensures
#  the robustness of the method regardless of the GCI estimation
#  performance.
#
# Example
#  Simplest, type g = iaif(x,fs) to estimate the glottal flow of a speech
#  frame.
#  And see the HOWTO_glottalsource.m example file.
#
# References
#  [1] P. Alku, "Glottal wave analysis with pitch synchronous iterative
#      adaptive inverse filtering", Speech Communication, vol. 11, no. 2-3,
#      pp. 109-118, 1992.
#
# Copyright (c) 2013 Aalto University
#
# Author
#  Tuomo Raitio <tuomo.raitio@aalto.fi>

iaif_ola_g <- function(x,fs,winLen,winShift,p_vt,p_gl,d,hpfilt){
  
  # This is a faster implementation of iaif_ola.
  # The output is the Glottal volume velocity waveform.
  
  # The IAIF method was modified from iaif_modified.R to iaif_fast.R
  # The new iaif_fast.R method accepts the FIR values of B as an input
  # argument, thus cutting the computation time down.
  
  # Further optimisation was carried out by using sapply instead of for-loop
  # sapply outputs a matrix if all columns are the same length. otherwise it outputs a list.
  
  # Initial settings
  #if (nargs()<3){
    winLen = 25/1000*fs
    winShift = 5/1000*fs
  #}
  
  # Use default settings from iaif.m implementation
  #if (nargs() < 5){
    hpfilt = 1
    d = 0.99
    p_gl = 2*rnd(fs/4000)
    p_vt = 2*rnd(fs/2000)+4
  #}
  
  # Allocate space
  g = zeros(1,length(x))      # Vector
  wins = zeros(1,length(x))   # Vector
  #g = rep(0 ,length(wave))  
  #wins = rep(0 ,length(wave))  
  
  # Compute all the start and stop times for the segmentation of the speech wave.
  # This is similar to computing x_frame in iaif_ola.
  # Initialize vectors of start and stop index of windows
  startv = c()
  stopv = c()
  
  # Initialize parameters
  start = 1
  stop = start+winLen-1
  cnt = 1
  
  # Do segmentation processing
  while (stop <= length(x)){
    
    startv[cnt] = start
    stopv[cnt] = stop
    start = start+winShift-1
    stop = start+winLen-1
    cnt = cnt+1
  }
  
  # Compute High Pass FIR filter
  Fstop = 40
  Fpass = 70
  Nfir = round(300/16000*fs)
  B = hpfilter_fir(Fstop, Fpass, fs, Nfir)
  
  #get all available x_frames into a data frame. each column is an x_frame
  x_frames = data.frame(sapply(1:length(stopv), function(m) x[startv[m]:stopv[m]]))
  #user  system elapsed 
  #0       0       0 
  #same as previous implementation
  
  g_frames = data.frame(sapply(1:ncol(x_frames), function(m) iaif_g(x_frames[[m]],B)))
  #user  system elapsed 
  #0.20    0.00    0.21 
  #same as previous implementation
  
  # use [[]] in list to access values as vector, rather than one [] which accesses
  # values as list
  #sapply(1:cnt, function(i) g[startv[i]:stopv[i]] = g[startv[i]:stopv[i]]+ gval[[i]]*t(win))
  #sapply(1:cnt, function(i) wins[startv[i]:stopv[i]] = wins[startv[i]:stopv[i]]+t(win))
  
  # Returns the N-point symmetric Hanning window
  win = hanningR(winLen, "hanning")
  
  for (i in 1:ncol(g_frames)){
    g[startv[i]:stopv[i]] = g[startv[i]:stopv[i]]+g_frames[[i]]*t(win)
    wins[startv[i]:stopv[i]] = wins[startv[i]:stopv[i]]+t(win)
  }

  idx = which(wins>0)        # find non-zero values' index
  g[idx] = g[idx]/wins[idx]
 
  return(g)

}
  





