polynomial_fitting <- function(GOI = NULL,
                               GCI = NULL,
                               wave,
                               thre) {
  
  # The function is used with GOI_GCI_detection() and GOI_instant() functions.
  # The function fits an 8th order polynomial (maximum polynomial size) over the GVV samples in the xRange. Once the polynomial coefficients
  # are extracted, the equation is modified by subtracting the threshold value from the intercept value in order to find the possible GOI/GCI
  # instances in the xRange. If multiple values are found, the value closest to the rough GOI/GCI instance is used.
  # Polynomial fitting is used since we are dealing with discrete values but want to find the exact occurance of glottal opening/closing. 
  #
  # INPUTS
  # GOI           Rough GOI instance
  # GCI           Rough GCI instance
  # wave          GVV waveform
  # thre          The threshold of the cycle, where the GOI/GCI occur (theoretically)
  #
  # GOI and GCI points are going to be found via the following interpolation:
  # 1. Fit a curve to 5 points from each direction of the initially detected GOI/GCI point
  # 2. The fitted polynomial will have a funciton of the following nature (if second order):
  #           ax^2 + bx + c = y
  # 3. The threshold [thre] is therefore the y-value of this polynomial over the specified range
  #           ax^2 + bx + c = thre
  # Therefore, we can create a new polynomial:   ax^2 + bx + (c-thre) = 0
  # 4. After getting the new polynomial equation, we use the {roots} function from the {signal} package to find where
  # the threshold intersects with the polynomial i.e. the roots (can use either method ="polyroot" or "eigen").
  # 5. Finally, since the polynomial will have multiple intersections with thre, we only extract the roots value which is
  # applicable to our range.
  
  # There are two ways to generate polynomials.
  # 1. polyfit()
  # 2. lm()
  # Even though polyfit() is much more user friendly, we will use lm(), as polyfit() throws an exception:
  # Error in qr.solve(A, y) : singular matrix 'a' in solve
  # This is due to its inability to find the inverse.
  # Instead, we will use lm() up to order 10 and truncate it as needed
  
  #############
  # GOI FITTING
  #############
  GOIvec = c()
  
  # STEP 1
  # Define x-range
  xRange = ((rnd(GOI) - 5):(rnd(GOI) + 5))
  # polynomial order
  #n = 3
  #polynomial = polyfit(xRange, wave[xRange], n) # outputs the intercept as the right-most value i.e. A(x^2) B(x) C, opposite to lm()
  polynomial <-
    lm(
      wave[xRange] ~ xRange + I(xRange ^ 2) + I(xRange ^ 3) + I(xRange ^ 4) +
        I(xRange ^ 5) + I(xRange ^ 6) + I(xRange ^ 7) + I(xRange ^ 8)
    )
  polynomial = polynomial$coefficients
  polynomial[is.na(polynomial)] <- 0
  polynomial = rev(polynomial)        # Reverse polynomial since roots() required the coefficients in decreasing order
  
  
  # To confirm via plot
  #plot(xRange, wave[xRange])
  #par(new=T)
  #plot(polyval(polynomial,xRange), col="red", type="l")
  #abline(h=thre, col="green",lty=2)
  
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
  
  # GOI instance
  GOIvec = polynomialRoots
  
  rm(xRange, polynomial, polynomialRoots)
  
  
  #############
  # GCI FITTING
  #############
  GCIvec = c()
  
  # STEP 1
  # Define x-range
  xRange = ((rnd(GCI) - 5):(rnd(GCI) + 5))
  # polynomial order
  #n = 3
  #polynomial = polyfit(xRange, wave[xRange], n)
  polynomial <-
    lm(
      wave[xRange] ~ xRange + I(xRange ^ 2) + I(xRange ^ 3) + I(xRange ^ 4) +
        I(xRange ^ 5) + I(xRange ^ 6) + I(xRange ^ 7) + I(xRange ^ 8)
    )
  polynomial = polynomial$coefficients
  polynomial[is.na(polynomial)] <- 0
  polynomial = rev(polynomial)        # Reverse polynomial since roots() required the coefficients in decreasing order
  
  
  # To confirm via plot
  #plot(xRange, wave[xRange])
  #par(new=T)
  #plot(polyval(polynomial,xRange), col="red", type="l")
  #abline(h=thre, col="green",lty=2)
  
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
 
  # GCI instance
  GCIvec = polynomialRoots
  
  # polynomialRoots might contain several roots within the xRange.
  # Thus, we use initial, "rough" GOI/GCI points as reference for the desired polynomial root value
  # Only take value closest to the original "rough" GOI point
  if (length(GOIvec) > 1) {
    GOIvec = GOIvec[which.min(abs(GOIvec - GOI))]
  }
  # Only take value closest to the original "rough" GCI point
  if (length(GCIvec) > 1) {
    GCIvec = GCIvec[which.min(abs(GCIvec - GCI))]
  }
  
  return(list("GOI" = GOIvec, "GCI" = GCIvec))
  
}