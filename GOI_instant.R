# Created 2/9/2016
# By Itay Ben-Dom

GOI_instant <- function(instantsList, winwave, cycleNumber, showPlot = T){
  
  # A method to detect GOI occurance point via rough estimation followed by polynomial curve fitting. 
  # A GOI instance is found between two successive valleys.
  # First, an estimated instance for the location of GOI is computed by examining the GVV for the first sample index that is greater than 
  # the threshold specified for the cycle.
  # Second, a polynomial (up to 8th order) is fitted to 10 points around the initial GOI estimate. The polynomial is then solved to find the
  # x-instance (sample index) where the threshold and the fitted curve intersect. This is the polyGOI point.
  #
  # INPUTS
  # instantsList    A row (list) of the dataframe object results_df, i.e. results_df[1,].
  # winwave         A GVV waveform.
  # cycleNumber     Corresponds to the row number of results_df. Provides an easy way to find the plot corresponding to the OQ value.
  #                 It also provides a visual explanation to why the OQ wasn't computed for the specified cycleNumber.
  # showPlots       If TRUE, the GVV, rough GOI point and exact GOI point are all plotted.
  #
  # OUTPUT
  # polyGOI          The GOI instance. The output is the output of polynomial_fitting() for GOI.
  #
  # Derived from the following function:
  # polynomial_fitting
  #
  # Related Functions:
  # OQpoints_tolerance
  # OQ50
  # OQsub50
  
  # Extract information from dataframe
  startIdx = instantsList[[1]]
  endIdx   = instantsList[[3]]
  thre     = instantsList[[4]]
  
  # initialize variables
  GOI = c()
  
  # find GOI instant
  for (m in floor(startIdx):floor(endIdx)){
    if (winwave[m-1] < thre){
      if (winwave[m] > thre){
        GOI = m
        break
      }
    }
  } 
  
  if (length(GOI) == 0){
    return(GOI)
  } else {
    # Perfrom curve fitting to find exact location of GOI/GCI instances. 
    xRange = (GOI-6) : (GOI+6) 
    # polynomial order
    #n = 3
    #polynomial = polyfit(xRange, wave[xRange], n) # outputs the intercept as the right-most value i.e. A(x^2) B(x) C, opposite to lm() 
    polynomial <- lm(winwave[xRange]~xRange+I(xRange^2)+I(xRange^3)+I(xRange^4)+I(xRange^5)+I(xRange^6)+I(xRange^7)+I(xRange^8))
    polynomial = polynomial$coefficients
    polynomial[is.na(polynomial)] <- 0
    polynomial = rev(polynomial)        # Reverse polynomial since roots() required the coefficients in decreasing order 
    
    # To confirm via plot
    if(TRUE){
      plot(x = xRange, 
           y = winwave[xRange], 
           ylim=c(min(polyval(polynomial,xRange),winwave[xRange]), max(polyval(polynomial,xRange),winwave[xRange])),
           main = "GOI Polynomial Fitting",
           xlab="",
           ylab="",
           col = ifelse(winwave[xRange]==GOI,"blue","black")
      )
      par(new=T)
      plot(xRange,
           polyval(polynomial,xRange), 
           col="red", 
           type="l",
           ylim=c(min(polyval(polynomial,xRange),winwave[xRange]), max(polyval(polynomial,xRange),winwave[xRange])),
           xlab = "Sample Index",
           ylab= "Amplitude"
      )
      abline(h=thre, col="blue",lty=2)
      #abline(v = polyGOI, col = "blue", lty = 1)   # object polyGOI doesn't exist yet: run line after going through the entire function
      legend("bottomright", 
             c("GVV", "Threshold", "POlynomial Fit"), 
             lwd = 1,
             lty = c(NA,2,1),
             pch = c(1,NA,NA),
             col = c("black", "blue", "red"))
    }
 
    # STEP 2 + 3
    # Subtract threshold to get new polynomial who's roots are the intersections
    polynomial[length(polynomial)] = polynomial[length(polynomial)] - thre
    
    # Step 4
    # Find roots 
    polynomialRoots = roots(polynomial)
    
    # Step 5
    # Find roots within x-range
    polynomialRoots = abs(polynomialRoots)
    polynomialRoots = polynomialRoots[which(polynomialRoots > xRange[1] & polynomialRoots < xRange[length(xRange)])]
    
    # find polynomial closest to original rough GOI
    polyGOI = polynomialRoots[which.min(abs(GOI - polynomialRoots))]
    
    if (showPlot == TRUE){
      plot(floor(startIdx):floor(endIdx),
           winwave[floor(startIdx):floor(endIdx)],
           #main = paste("Glottal Pulse", cycleNumber),
           main = paste("Glottal Pulse Opening Phase"),
           xlab = "Sample Index",
           ylab = "Amplitude")
      abline(h = thre, lty = 2)
      abline(v = GOI, col = "blue", lty= 2)
      abline(v = polyGOI, col = "red", lty = 1)
      legend("bottomright", 
             c("GVV", "Threshold","Rough GOI", "Poly GOI"), 
             lwd = 1,
             lty = c(NA,2,3,1),
             pch = c(1,NA,NA,NA),
             col = c("black", "black","blue", "red"))
    }
    
    return(polyGOI)
  }
  

  
  
  
  
}