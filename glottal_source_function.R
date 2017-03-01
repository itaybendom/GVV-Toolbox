c# 30/11/16
# 
# The following method extract the Glottal Source, or Glottal Velocity Waveform, using Alku's Iterative Adaptive
# Inverse Filtering (IAIF) method.
# Additional filtering of the wave is carried out, outputting the locations of the minima and maxima points, to be
# used for the extraction of other time domain parameters:
#   Open Quotient
#   Speed Quotient
#   Normalised Peak Amplitude
#   Jitter
#   Shimmer
#
# INPUTS
# wave      A speech signal
#
# OUTPUTS
#   glottal     The glottal waveform derived from the input signal.
#   pitch       The pitch of the waveform.
#   peakIdx     A vector containing the sample index of the peaks
#   valIdx      A vector containing the sample index of the valleys
#   peakAmp     A vector containing the amplitude values of the peaks
#   valAmp      A vector containing the amplitude values of the valleys
#
# Author: Itay Ben-Dom @ University of Auckland 2016
#
###################################################################################################################

glottal_source_function <- function (wave, 
                                     fs, 
                                     fileName = "", 
                                     idx,
                                     polarity = 1,
                                     plotCycle = FALSE,
                                     plotGVV = TRUE) {
  
  # Polarity check
  #polarity = polarity_reskew(wave, fs)
  wave = polarity * wave
  print(polarity)
  # Using SRH to get robust pitch track and its mean value in the vowel period
  pitch = mean( pitch_srh(wave, fs, 50, 250, 60))
  
  #################
  # GVV Computation
  #################
  
  # Lowpass Butterworth filter with cutoff frequency dependent on the average pitch period.
  # The cutoff frequency is (2*pitch)/fs. We add an integer to increase the cutoff frequency so that we don't lose any valuable information.
  # Since some of the data has uncharacteristically low pitch, we need to make sure that the butterworth filter applied doesn't have a too
  # low of cutoff frequency (pitch = 100 for man, 200 for women, on average). Thus, we add an IF condition.
  if (pitch < 100) {
    cutoff = (2 * 6 * 100)
  } else{
    cutoff = 2 * 6 * pitch
  }
  butterlist = butter(16, cutoff / 16000, type = "low")
  
  # Erase the interference of the neighbouring consonants by truncating the wave
  wavin = wave[rnd(0.05 * length(wave)):rnd(0.95 * length(wave))]
  
  # Using Alku's Iterative Adaptive Inverse Filtering to obtain the GVV signal
  gvv = iaif_ola_g(wavin, fs)
  # Filter GVV to generate a more smooth wave
  gvv = filter(butterlist$b, butterlist$a, gvv)
  
  # Find the first instance of zero crossing (sgvv) and start analysis from there.
  for (i in 2:length(gvv)) {
    if (gvv[i - 1] > 0) {
      if (gvv[i] < 0) {
        sgvv = i - 1
        break
      }
    }
  }
  
  # gvv is a time-series objet. This results in a time shift when analysed. Thus, we convert it to a vector of type "double"
  gvvWave = c(gvv)
  
  # Define a window length
  winlen = round(length(gvvWave)/10)
  
  
  # start and stop smaple index for the window-shifting loop
  start = round(sgvv)
  stop =  round(length(gvvWave)-1-winlen)
  
  # Determine tolerance for peak/valley-detection method
  if (pitch >= 210) {
    tol = pitch / 5
  } else if (pitch < 210 & pitch >= 195) {
    tol = pitch / 4
  } else if (pitch < 195 & pitch >= 150) {
    tol = pitch / 3
  } else if (pitch < 150 & pitch >= 100) {
    tol = pitch / 2
  } else if (pitch < 100) {
    tol = pitch
  }
  
  # initialize peak/valley vectors 
  peakIdx = c()
  valIdx  = c()
  peakAmp = c()
  valAmp  = c()
  aa = 0
  
  # Peak detection processing
  for (i in seq(start, stop, winlen)){
    
    # wave segmentation
    if (aa > 0){
      if (i-sgvv > 0){
        sas = i-10
      }
    } else {
      sas = i
    }
    sae = i+winlen-1
    segmentGVV = gvvWave[sas:sae]

    # Extract peaks index & amplitude values
    peak_list = peak_detection(segmentGVV, "q", tol)
    peakIdxtemp = peak_list$index
    peakAmptemp = peak_list$amplitude
    peakIdx = c(peakIdx, peakIdxtemp+sas-1) # add new peaks to peak vector
    peakAmp = c(peakAmp, peakAmptemp) # add new peaks amplitude to peak vector

    # Extract valleys index & amplitude values
    valley_list = peak_detection(segmentGVV, "vq", tol)
    valIdxtemp = valley_list$index
    valAmptemp = valley_list$amplitude
    valIdx = c(valIdx, valIdxtemp+sas-1) # add new valleys to valley vector
    valAmp = c(valAmp, valAmptemp)
    
    aa = aa+1
   
    # plot window shift cycles
    if (plotCycle){
      plot(x = sas:sae,
           y = segmentGVV,
           col = "black")
      abline(v = valIdx, col = "blue", lty = 1)
      abline(v = peakIdx, col = "red", lty = 1)
    }
  
  }
  
  ######################################
  # Peak/Valley Detection - Improvements
  ######################################
  
  # Eliminate repeated values due to window overlap
  rmIdx = c()
  for (i in 1:(length(valIdx)-1)){
    bb = which(ceil(valIdx[i]) == ceil(valIdx))
    if (length(bb) > 1){
      rmIdx = c(unique(rmIdx), bb[-1])
    }
    rm(bb)
  }
  if (length(rmIdx) > 0){
    valIdx = valIdx[-rmIdx]
    valAmp = valAmp[-rmIdx]
  }
  
  rmIdx = c()
  for (i in 1:(length(peakIdx)-1)){
    cc = which(ceil(peakIdx[i]) == ceil(peakIdx))
    if (length(cc) > 1){
      rmIdx = c(unique(rmIdx), cc[-1])
    }
    rm(cc)
  }
  if (length(rmIdx) > 0){
    peakIdx = peakIdx[-rmIdx]
    peakAmp = peakAmp[-rmIdx]
  }
  
  rm(rmIdx)
  
  # Sort valley and peak index by ascending order
  valIdxSort = sort(valIdx, index.return = T)$ix
  peakIdxSort =  sort(peakIdx, index.return = T)$ix
  valIdx = valIdx[valIdxSort]
  valAmp= valAmp[valIdxSort]
  peakIdx = peakIdx[peakIdxSort]
  peakAmp= peakAmp[peakIdxSort]
  
  # AMPLITUDE THRESHOLD DETECTION
  
  rmIdx = c()
  
  # Eliminate valleys that are above the average amplitude level
  valAmpRange = mean(gvvWave[valIdx]) * 0.20       # Amplitude threshold
  rmIdx = which(gvvWave[valIdx] > valAmpRange)
  if (length(rmIdx) > 0) {
    valIdx = valIdx[-rmIdx]
    valAmp = valAmp[-rmIdx]
    if (any(rmIdx == 1)) {
      # Only remove first peak
      peakIdx= peakIdx[-1]
      peakAmp = peakAmp[-1]
    }
  }
  
  rm(rmIdx)
  rmIdx = c()
  
  peakAmpRange = mean(gvvWave[peakIdx]) * 0.20        # Amplitude threshold
  rmIdx = which(gvvWave[peakIdx] < peakAmpRange)
  if (length(rmIdx) > 0) {
    peakIdx = peakIdx[-rmIdx]
    peakAmp = peakAmp[-rmIdx]
  }
  rm(rmIdx)
  
  
  # INDEX THRESHOLD DETECTION
  
  ######################################
  # Eliminate Neighbouring Peaks/Valleys
  ######################################
  
  # Perfrom neighbouring peaks elimination twice, as the first iteration sometimes fails to detect all double 
  # occurances
  for (j in 1:2){
    rmIdx = c()
    
    # Remove remaining neighbouring peaks
    # If length(closeTest) == 0 that means there was no valley between 2 neighbouring peaks, and as a result we remove the peak with the
    # smaller amplitude between the two.
    for (i in 1:(length(peakIdx) - 1)) {
      closeTest = which(valIdx > peakIdx[i] & valIdx < peakIdx[i + 1])      # Find a valley between 2 neighbouring peaks
      if (length(closeTest) == 0) {
        # Identify index of smaller peak amplitude
        #rmIdx = c(rmIdx, which(peakAmp == min(peakAmp[i], peakAmp[i + 1])))  
        # Since the glottal closing duration should be shorter than the opening duration, we choose the peak that is 
        # further along the time axis, thus we eliminate neighbouring peaks based on their index, not amplitude.
        # this differs to the OQsub50.R method, which uses the above line.
        rmIdx = c(rmIdx, i)
      }
    }
    if (!is.null(rmIdx)) {
      # Remove neighbouring peaks
      peakIdx = peakIdx[-rmIdx]
      peakAmp = peakAmp[-rmIdx]
    }
    
    rm(rmIdx)
  }
 
  for (j in 1:2){
    rmIdx = c()
    # Remove remaining neighbouring valleys
    # If length(closeTest) == 0 that means there was no valley between 2 neighbouring valleys, and as a result we remove the valley with the
    # smaller amplitude between the two.
    for (i in 1:(length(valIdx) - 1)) {
      closeTest1 = which(peakIdx > valIdx[i] &
                           peakIdx < valIdx[i + 1])      # Find a peak between 2 neighbouring valleys
      if (length(closeTest1) == 0) {
        # Check if there is no peak between the 2 neighbouring valleys
        # Check for irregular spacing
        rmIdx = c(rmIdx, which(valAmp == max(valAmp[i], valAmp[i + 1])))    # Identify index of higher valley amplitude
        
      }
    }
    if (!is.null(rmIdx)) {
      # Remove neighbouring peaks
      valIdx = valIdx[-rmIdx]
      valAmp = valAmp[-rmIdx]
    }
    
    rm(rmIdx)
  }
  
  
  # Final step
  
  # We want to make sure both start and finish on a valley, as we define a cycle as occuring between two valleys
  i = 1
  while (peakIdx[i] < valIdx[1]) {
    # start on a valley
    peakIdx = peakIdx[-i]
    peakAmp = peakAmp[-i]
  }
  while (peakIdx[length(peakIdx)] > valIdx[length(valIdx)]) {
    # end on a valley
    peakIdx = peakIdx[-length(peakIdx)]
    peakAmp = peakAmp[-length(peakAmp)]
  }
  
  if (plotGVV){
    
    #pdf(
    #  paste("E:/Master/Prom_DB_Glottal_Plots/GVV_", fileName, "_", idx, ".pdf", sep = ""),
    #  height = 10.22,
    #  width = 16.23
    #)
    
    plot(gvvWave,
         main = paste(fileName, "_", idx, sep = ""),
         ylab = "Sample Index",
         xlab = "Amplitude")
    abline(v = valIdx, col = "blue")
    abline(v = peakIdx, col = "red")
    
    #dev.off()   # to avoide writing error. 
  }
  
  return(list("glottal" = gvvWave,
              "pitch"   = pitch,
              "peakIdx" = peakIdx, 
              "valIdx"  = valIdx, 
              "peakAmp" = peakAmp, 
              "valAmp"  = valAmp)
         )
}


