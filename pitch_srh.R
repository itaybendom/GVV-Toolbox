pitch_srh <- function(wave,fs,f0min,f0max,hopsize){
  # SRH is a robust pitch tracker.
  #
  # Octave compatible
  #
  # Description
  #  The Summation of the Residual Harmonics (SRH) method is described in [1].
  #  This algorithm exploits a criterion taking into the strength of the
  #  harmonics and subharmonics of the residual excitation signal in order to
  #  determine both voicing decision and F0 estimates. It is shown in [1] to
  #  be particularly interesting in adverse conditions (low SNRs with various
  #  types of additive noises).
  #
  #
  # Inputs
  #  wave            : [samples] [Nx1] input signal (speech signal)
  #  fs              : [Hz]      [1x1] sampling frequency
  #  f0min           : [Hz]      [1x1] minimum possible F0 value
  #  f0max           : [Hz]      [1x1] maximum possible F0 value
  #  hopsize         : [ms]      [1x1] time interval between two consecutive
  #                    frames (i.e. defines the rate of feature extraction).
  #                    hopsize = frame_shift
  #
  # Outputs
  #  F0s             : vector containing the F0 estimates (values are provided even in unvoiced parts).
  #  VUVDecisions    : vector containing the binary voicing decisions.
  #  SRHVal          : vector containing the SRH values (according the
  #                    harmonic criterion - voicing decision are derived from
  #                    these values by simple thresholding).
  #  time            : [s] Analysis instants of the features described above.
  #
  #
  # References
  #  [1] T.Drugman, A.Alwan, "Joint Robust Voicing Detection and Pitch Estimation
  #      Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011.
  #      Publication available at the following link:
  #      http://tcts.fpms.ac.be/~drugman/files/IS11-Pitch.pdf
  
  # Warnings
  if (length(wave)/fs < 0.1){
    stop("SRH error: the duration of your file should be at least 100ms long")
  }
  if (f0max<=f0min){
    stop("You look funny! Your f0min should be lower than f0max!!")
  }
  if (fs != 16000){
    print("Sample rate not equal to 16kHz. Audio is resampled.")
    wave = resample(wave, 16000, fs)
    fs = 16000
  }
  if (nargs() < 5){
    hopsize = 10
  }
  
  # Settings
  nHarmonics = 5
  SRHiterThresh = 0.1
  SRHstdThresh = 0.05
  VoicingThresh = 0.07
  VoicingThresh2 = 0.085
  LPCorder = rnd(3/4*fs/1000)
  Niter = 2
  
  # Compute LPC Residual
  list.res = lpcresidual(wave,rnd(25/1000*fs),rnd(5/1000*fs), LPCorder)
  res = list.res$res  # Extract res values from lpcresidual returned list
  
  # Create frame matrix
  waveLen = length(wave)
  #rm(wave)
  frameDuration = rnd(100/1000*fs)-2 # Minus 2 to make equivalent to original
  shift = rnd(hopsize/1000*fs)
  halfDur = rnd(frameDuration/2)
  time = seq( (halfDur+1), (waveLen-halfDur), shift )
  N = length(time)
  frameMat = zeros(frameDuration, N)
  for (n in 1:N){
    frameMat[,n] = res[ (time[n] - halfDur) : (time[n] + (halfDur-1)) ]
  }
  #rm(res)
  
  # Create window matrix and apply to frames
  win = blackman(frameDuration)
  winMat = repmat( matrix(win), 1 , N )
  frameMatWin = frameMat * winMat
  #rm(winMat)
  
  # Do mean subtraction
  frameMean =apply(frameMatWin,2, mean)   # MATLAB: mean(frameMatWin,1);
                                          # In MATLAB, mean(x,DIM) takes the mean along the dimension DIM of X
                                          # However, mean(x,1) takes the mean of the matrix per COLUMN
                                          # In R, apply(x,1,mean) takes the mean per ROW
                                          # Thus, use apply(x,2,mean) to compute mean per column
  frameMeanMat = repmat(frameMean,frameDuration,1)
  frameMatWinMean = frameMatWin - frameMeanMat
  #rm(frameMean, frameMeanMat, frameMatWin, frameMat)
  
  # Compute spectrogram matrix 
  specMat = zeros(fs, ncol(frameMatWinMean))  # [fs x ncol(frameMatWinMean)] dimansions
  
  for (i in 1:ncol(frameMatWinMean)){
    tempCol = abs( fftn(frameMatWinMean[,i],fs) )
    specMat[,i] = matrix(tempCol)  # Overall: specMat = abs( fft(frameMatWinMean,fs) )
  }
  # Debug: if specMat is a column vector, R turns it into a row vector when we try to remove
  # half of its rows, which leads to issues when using apply() function.
  # So first we check if it's a column vector. If so, we turn it into a matrix before indexing.
  if (ncol(specMat) == 1){
    specMat = matrix(specMat[ (1:(fs/2)), ])  # Remove half the matrix rows
  } else{ # if size(specMat) = [number>1 number>1] then we can do standard indexing
    specMat = (specMat[ (1:(fs/2)), ])  # Remove half the matrix rows
  }
  specDenom = matrix(sqrt( apply((specMat^2), 2, sum) ), nrow = 1)
  specDenomMat = repmat( specDenom, (fs/2), 1 )
  specMat = specMat / specDenomMat
  rm(specDenom, specDenomMat)
 
  ## Estimate the pitch track in 2 iterations 
  for (Iter in 1:Niter){
    
    list.SRH = SRH( specMat, nHarmonics, f0min, f0max )
    F0s = list.SRH$F0
    SRHVal = list.SRH$SRHVal
    
    if (max(SRHVal) > SRHiterThresh){
      F0medEst = median( F0s[ SRHVal > SRHiterThresh ] )
      
      # Only refine F0 limits if within the original limits
      if (rnd(0.5*F0medEst) > f0min){
        f0min = rnd(0.5*F0medEst)
      }
      if (rnd(2*F0medEst) < f0max){
        f0max = rnd(2*F0medEst)  
      } 
    } 
  }   

  return(F0s)
  
  #time = time/fs
  
  # Voiced-Unvoiced decisions are derived from the value of SRH (Summation of Residual Harmonics)
  #VUVDecisions = zeros(1,N)
  
  #if (std(SRHVal,1) > SRHstdThresh){
  #  VoicingThresh = VoicingThresh2;
  #}
  
  #VUVDecisions[SRHVal > VoicingThresh] = 1
  
  # Return Function Outputs as list
  #return(list("F0s" = F0s, "VUVDecisions" = VUVDecisions, "SRHVal" = SRHVal, "time" = time))
  
  
  # Same as MATLAB up to here
  # WORKS

}
