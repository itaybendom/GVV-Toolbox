OQ50_function <- function(glottalList, threshold = 0.5){
  
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
  
  
  # The GVV's final cycle might end abruptly, with no GCI point as the GVV wave does not reach the threshold.
  # In that case, GOI_GCI_detection() returns an empty list, which means goi_gci_points will be length zero. If that is the case,
  # The corresponding entry in the dataframe object is removed, as that cycle is not valid for open quotient calculation.
  
  rmIdx = c()
  goi_gci_list = list()
  
  for (i in 1:nrow(results_df)) {
    goi_gci_points = GOI_GCI_detection(results_df[i,], winwave , i)  # Find GOI/GCI points for each cycle
    if (length(goi_gci_points) > 0) {
      # If GOI/GCI points are found, we add them to the list
      goi_gci_list[[i]] = goi_gci_points
    } else {
      rmIdx = c(rmIdx, i)        # If GOI/GCI points aren't found, we remove the cycle entry from the dataframe.
    }
    rm(goi_gci_points)
  }
  # Remove unused cycles from dataframe and GOI/GCI list
  if (length(rmIdx) > 0) {
    results_df = results_df[-rmIdx,]
    goi_gci_list = goi_gci_list[-rmIdx]
  }
  # Add GOI and GCI to dataframe object, GOI and GCI correspond to x-value instances i.e. samples
  for (i in 1:nrow(results_df)) {
    results_df[i, 5] = goi_gci_list[[i]]$GOI
    results_df[i, 6] = goi_gci_list[[i]]$GCI
  }
  names(results_df)[5:6] = c("GOI", "GCI")
  
  
  ###########################
  # OPEN QUOTIENT Calculation
  ###########################
  
  # Since we got GOI and GCI points we can now find the OQ.
  # It's OK the GOI and GCI points are not actualy values, we don't need to round them, since we got them via interpolation.
  # Compute OQ measurements via:
  # q1 = GCI - GOI                                 # time difference between GCI and GOI
  # cycleT = df$`Start Idx`[i] - df$`End Idx`[i]   # find the time difference between two consequtive valleys (one period)
  # OQ = q1 / cycleT                               # OQ = (time between GCI and GOI)/(time between consequtive GOIs)
  # Compute OQ and add measurements to dataframe:
  for (i in 1:nrow(results_df)) {
    OQ = (results_df$GCI[i] - results_df$GOI[i]) / (results_df$`End Idx`[i] - results_df$`Start Idx`[i])
    results_df[i, 7] = OQ
    rm(OQ)
  }
  names(results_df)[7] = "OQ"
  
  # Plot OQ vs. pitch period.
  # plot(16000/diff(results_df$`Start Idx`))
  # par(new=T)
  # plot(1:39,results_df$OQ[2:40],type="l",col="red")
  
  # Return OQ values
  return(list(results_df$OQ))
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  