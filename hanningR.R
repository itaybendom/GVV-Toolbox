hanningR <- function(n, type = ""){
  # The following function is equivalen to MATLAB's hanning() function
  # MATLAB has a hanning() function and a hann() funciton. 
  # The later is part of the Singal Processsing Toolbox
  # R has several window functions such as windowfunc in {phonTools package}
  #
  # The following 
  #
  # INPUTS:
  #
  #   n:        window length
  #   type:     determines if Hanning window is computed as hanning() funciton or hann() funciton  
  # 
  # OUTPUTS:
  #
  #   w:        Hanning window as column vector
  #
  # For reference: https://www.dsprelated.com/freebooks/sasp/Matlab_Hann_Window.html
  #
  # Author: Itay Ben-Dom @ 22 April, 2016
  
  # If a vector is given, the window will have the same length as the vector.
  if (length(n) > 1){
    n = length(n)
  }
  # Warning, window length must be an integer
  if (!is.numeric(n)){
    stop("Invalid number of points specified.")
  } 
  # Window of length n=1  
  if (n == 1){
    w = 1
  } 
  if (type == "hanning"){
    w = 0.5*(1 - cos(2 * pi * as.matrix((1:n), ncol = 1)/(n+1)))  
  }
  if (type == "hann"){
    w = 0.5*(1 - cos(2 * pi * as.matrix((0:(n-1)), ncol = 1)/(n-1)))  
  }
  if (nchar(type) == 0){
    stop("No window type provided.")
  } 
  if (!(type == "hanning" | type == "hann")){
    stop("Invalid window type provided.")
  }
  if (is.null(w)){
    stop("Invalid window type provided.")
  } 
    
  return(w)
  
}