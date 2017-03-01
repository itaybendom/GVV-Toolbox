# Created 2/9/2016
# By Itay Ben-Dom

GCI_instant <- function(instantsList, winwave, cycleNumber, showPlot = T){
  
  # A method to detect GCI occurance via initial interpolation and detecting the minima point of the interpolated GVV waveform's derivative.
  # A GCI instance is found between successive peak and valley.
  # The GVV is interpolated between the peak index and the following valley index. Upsampling is used to increase the resolution from
  # 16 kHz to 48 kHz. the interpolated waveform [interpWave] is differentiated using the gradient() function. 
  # Following the Liljencrants-Fant (LF) model, the minima point of the derivative GVV waveform corresponds to the exact point of glottal
  # closure time-parameter of the GVV waveform. 
  # Once the point is found, additional computations are carried out in order to detect the exact (approximately) point of GCI instance in
  # the GVV waveform. 
  #
  # INPUTS
  # instantsList    A row (list) of the dataframe object results_df, i.e. results_df[1,].
  # winwave         A GVV waveform.
  # cycleNumber     Corresponds to the row number of results_df. Provides an easy way to find the plot corresponding to the OQ value.
  #                 It also provides a visual explanation to why the OQ wasn't computed for the specified cycleNumber.
  # showPlots       If TRUE, the GVV, interpolated GVV, derivative GVV, and GCI point are all plotted.
  #
  # OUTPUT
  # GCI             The GCI instance [sample index]. 
  #
  # Related Functions:
  # OQpoints_tolerance
  # OQsub50
  #
  # Papers:
  # Glottal Airflow and Electroglottographic Measures of Vocal Function at Multiple Intensities (Dromey)
  #
  # Parameterisation Methods of the Glottal Flow Estimated by Inverse filtering (Alku)
  #
  # A Four-Parameter Model of the Glottal Flow (Fant)
  # http://www.speech.kth.se/prod/publications/files/qpsr/1985/1985_26_4_001-013.pdf
  
  # Extract information from dataframe
  startIdx = instantsList[[3]] # peakIdx
  endIdx   = instantsList[[2]] # valIdx
  
  # initialize variables
  GCI = c()
  
  # Interpolate wave. Upsample from 16 kHz to 48 kHz
  x = seq(floor(startIdx), floor(endIdx))
  y = winwave[x]
  interpWave = spline(x, y, n = 3 * length(x))
  # interpWave$x = index
  # interpWave$y = amplitude
  
  # Compute derivative of interpolated GVV waveform
  dgvv = gradient(interpWave$y)
  
  # Find location of valley which corresponds to GCI instance
  dgvv_val = peak_detection(dgvv, "qv", 10)
  dgvv_val_idx = dgvv_val$index[which.min(dgvv_val$amplitude)]
  
  # Find actual sample corresponding to the location of the valley in the GVV waveform
  # HOW TO:
  # Since we are dealing with discrete values, we cannot use (interpWave$x[dgvv_val_idx]) to find the sample index for the GCI occurance.
  # First, we find the sample index values the valley is located in-between:
  # tmp1 = interpWave$x[floor(dgvv_val_idx)]
  # tmp2 = interpWave$x[ceil(dgvv_val_idx)]
  # Second, we find the difference between them. 
  # tmp2 - tmp1
  # Third, we multiply that difference by the decimals of dgvv_val_idx.
  # Four, we add that value to tmp1 to find the exact sample index occurance of the GCI instance.
  
  diffVal = interpWave$x[ceil(dgvv_val_idx)] - interpWave$x[floor(dgvv_val_idx)] 
  fracVal = dgvv_val_idx - floor(dgvv_val_idx)
  GCI = interpWave$x[floor(dgvv_val_idx)] + diffVal * fracVal
  
  # Plots for reference 
  if (showPlot == TRUE){
    #plot(floor(startIdx):floor(endIdx),
    #     winwave[floor(startIdx):floor(endIdx)],
    #     col = "red",
    #     type = "s", 
    #     xlab = "",
    #     ylab = ""
    #     )
    #par(new=T)
    plot(interpWave$x,
         interpWave$y,
         xlab = "",
         ylab = ""
         )
    par(new=T)
    plot(interpWave$x,
         dgvv,
         type = "p",
         yaxt = "n",             # GVV derivative has different amplitude values that are not plotted for ease of viewing
         #main = paste("GVV Cycle", cycleNumber),
         main = "Glottal Flow and Glottal Derivative",
         xlab = "Sample Index",
         ylab = "Amplitude",
         col= "blue"
         )
    abline(v = GCI, lty = 2, col = "blue")
    legend("bottomleft", 
           #c("GVV", "Interp","dGVV", "GCI"),
           c("GVV", "dGVV", "GCI"),
           lwd = 1,
           lty = c(NA,NA,2),
           pch = c(1,1,NA),
           col = c("black","blue", "blue"))
  }
  
  return(GCI)
}




