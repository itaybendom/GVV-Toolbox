maxR <- function(x){
  #
  # Function emulates MATLAB max() function for matrix.
  #
  # MATLAB
  #   [Y,I] = MAX(X) returns the indices of the maximum values in vector I.
  #   If the values along the first non-singleton dimension contain more
  #   than one maximal element, the index of the first one is returned.
  #
  # INPUT
  #   x     Matrix
  #
  # OUTPUTS
  #   Y     returns the the maximum values of each column in vector form. 
  #   I     Returns the indices of the maximum values in vector form.
  #
  # Author: ITAY BEN-DOM
  # Date: 26/4/2016
  
  # Find maximum values in each column
  Y = apply(x, 2, max)
  
  # Index for each column maximum
  I = max.col(t(x))
  
  return(list("max" = Y, "index" = I))
  
}