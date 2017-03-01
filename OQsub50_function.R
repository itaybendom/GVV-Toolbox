###############################################################################################################################################
# The following function is an implementation of Tian's OQ measurements code.
# This is a far more robust implementation of Amit's OQ_detection code.
# This file is a further optimization of OQpoints_original_implementation.R
# The window loop was removed as it was unnecessary for my purposes.

# 22/8/2016
# This current implementaion is different than OQpoints.R although some similarities are present i.e. the detection 
# of the "rough" GOI/GCI points location.

# 2/9/2016
# This is a modification of the OQpoints_tolerance.R function (also OQ50). In this version, all steps are the same 
# except for the GCI detection method which is based on the GVV waveform derivative and involves interpolation rather 
# than polynomial fitting.
#
# 5/12/2016
# This is an implementation of the OQsub50.R function, except that it accepts as an input the pre-computed glottal 
# signal, using the glottal_source_function.R method.
# 
#
# Created 5/12/2016
# By Itay Ben-Dom
###############################################################################################################################################
###############################################################################################################################################

OQsub50_function <- function(glottalList, threshold = 0.5){
  
  # First extract all the applicable information from the list generated via glottal_source_function.R
  winwave = glottalList[[1]]$glottal  # glottal waveform
  valIdx  = glottalList[[1]]$valIdx   # vector of valley instances
  valAmp  = glottalList[[1]]$valAmp   # vector of corresponding valley amplitudes
  peakIdx = glottalList[[1]]$peakIdx  # vector of peaks instances
  peakAmp = glottalList[[1]]$peakAmp  # vector of corresponding peak amplitudes
  
  # Next step is to set up a data-frame object containing all the information required for OQ calculations.
  
  ###################
  # Data Frame Set-Up
  ###################
  
  # Dataframe (with respect to Glottal velocity waveform)
  # start         end            peak             threshold value
  # valIdx[1]     valIdx[2]      peakIdx[1]              X
  # valIdx[2]     valIdx[3]      peakIdx[2]              Y
  
  # Create dataframe first 2 columns, indicating the start and end instances of each glottal cycle in the waveform.
  startIdx = list()
  endIdx = list()
  for (i in 1:(length(valIdx) - 1)) { # A cycle occurs between valIdx[i] and valIdx[i+1]
    startIdx[i] = valIdx[i]
    endIdx[i] = valIdx[i + 1]
  }
  
  # Create dataframe object
  results_df = data.frame(unlist(startIdx), unlist(endIdx), unlist(peakIdx))
  names(results_df)[1:3] = c("Start Idx", "End Idx", "Peak Idx")
  
  #######################################
  # Determine Threshold's Amplitude Value
  #######################################
  
  # The threshold value per cycle is determind by multiplying the threshold by the peak-to-peak amplitude value, and 
  # adding that value to the lowest valley.
  threVal = c()
  for (i in 1:length(peakAmp)) {
    threVal = c(threVal,  min(valAmp[i], valAmp[i + 1]) + (peakAmp[i] - min(valAmp[i], valAmp[i + 1])) * threshold)
  }
  
  # Add threshold amplitude column to data-frame
  results_df[, 4] = threVal
  names(results_df)[4] = "Threshold"
  
  # Now we can calculate the GOI and GCI instances.
  
  #######################
  # GOI and GCI Detection
  #######################
  
  # Note that GOI/GCI instances will not match GOI/GCI instances of original speech waveform, as that wave had been 
  # truncated and filtered.
  
  GOI_list = list()
  GCI_list = list()
  rmIdx = c()
  
  for (i in 1:nrow(results_df)) {
    flag = c()
    GOItmp  = GOI_instant(results_df[i, ], winwave,  i, showPlot = F)
    if (length(GOItmp > 0)) {
      GOI_list[i] = GOItmp
    } else{
      rmIdx = c(rmIdx, i)
      flag = 1
    }
    GCItmp = GCI_instant(results_df[i, ], winwave,  i, showPlot = T)
    if (length(GCItmp > 0)) {
      GCI_list[i] = GCItmp
    } else{
      if (flag != 1)
        rmIdx = c(rmIdx, i)
    }
    rm (GOItmp, GCItmp, flag)
  }
  
  # Remove unused cycles from dataframe and GOI/GCI list
  if (length(rmIdx) > 0) {
    results_df = results_df[-rmIdx, ]
    GOI_list = GOI_list[-rmIdx]
    GCI_list = GCI_list[-rmIdx]
  }
  # Add GOI and GCI to dataframe object, GOI and GCI correspond to sample index of GVV waveform
  results_df[, 5] = unlist(GOI_list)
  results_df[, 6] = unlist(GCI_list)
  names(results_df)[5:6] = c("GOI", "GCI")
  
  ###########################
  # OPEN QUOTIENT Calculation
  ###########################
  
  # Since we got GOI and GCI points we can now find the OQ.
  # It's OK the GOI and GCI points are not actualy values, we don't need to round them, since we got them via interpolation.
  # Compute OQ measurements via:
  # OQ = (time between GCI and GOI)/(time between consequtive GOIs)
  # q1 = GCI - GOI                     # time difference between GCI and GOI
  # cycleT = valIdx[i+1] - valIdx[i]   # find the time difference between two consequtive valleys (one period)
  # OQ = q1 / cycleT                  
  # Compute OQ and add measurements to dataframe:
  for (i in 1:nrow(results_df)) {
    OQ = (results_df$GCI[i] - results_df$GOI[i]) / (results_df$`End Idx`[i] - results_df$`Start Idx`[i])
    results_df[i, 7] = OQ
    rm(OQ)
  }
  names(results_df)[7] = "OQ"
  
  # Return OQsub50 values
  # Vaector of OQ values is outputted as list object. This allows the values to be inserted into another list
  # and then analyzed using lapply(). results_df$OQ is outputted as a list to mimic the output of trapply().
  return(list(results_df$OQ))
  
}




