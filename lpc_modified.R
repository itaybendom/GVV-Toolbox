#LPC function in {phonTools} package
#http://www.santiagobarreda.com/rstuff/lpc/lpc.html
#Code was adapted from MATLAB
#http://au.mathworks.com/help/signal/ref/lpc.html#f7-1107

#Linear predictive coding (autocorrelation) estimation of filter coefficients. 
#tested for nat_speech.wav file

lpc_modified <- function (sound, order = round(fs/1000) + 3, fs = 10000, show = FALSE, 
          add = FALSE, preemph = TRUE, mean = FALSE, window = TRUE) 
{
  # Sound: either a numeric vector representing a sequence of samples taken from a sound wave or a sound object
  # This if statement is not the case
  if (class(sound) == "ts")
    fs = frequency(sound)   #frequency returns the number of samples per unit time
  
  # Since we use loadsound() function, object is of class sound.
  # This if statement is executed
  if (class(sound) == "sound") 
  {
    fs = sound$fs       #Set sampleing frequency to that of Sound Object 
    sound = sound$sound #Set sound as sampled of Sound Object
  }
 
  # Debugging
  if (!is.numeric(sound)) 
    stop("Input must be numeric.")
  
  # Pre-emphasis
  if (preemph == TRUE) 
    sound = preemphasis(sound, fs = fs)
  
  # Object length
  n = length(sound)                 #n is length of vector
  
  # Normalization of sound samples around the mean
  if(mean == TRUE){
    sound = sound - mean(sound)    
  }
  
  # Hanning window function
  if (window == TRUE){
    sound = sound * hanningR(sound, type = "hanning") #apply a Hanning window of length n (length of sound)
  }
  
  # Zero Padding
  sound = c(sound, rep(0, order))   
  
  #sapply: apply a function over a list or vector
  #Create a sequence from 1 to length(sound) with increments of 1. Thus, created sequence 1,2,3...6186
  #(function (x) sound[(x):(x + order)]) creates a matrix of (LPC_order+1)=21 rows and length(sound)=6186
  #However, each row is missing one element, going from [1,6186] to [21,6167]
  #t() is transpose, resulting in a final 6186 rows by 21 columns matrix 
  #Thus, predictors is a 6186x21 matrix
  #However, it starst losing elements per row from [6166,21], [6167.20]... [6186,1] i.e. [6186,2] = NA
  #Overll, the X matrix in the MATLAB example will start on row p = order
  predictors = t(sapply(seq(1, n, 1), function(x) sound[(x):(x + order)]))
  
  #Computation of predictor coefficients
  y = sound[1:n]            # ignore padding
  r = y %*% predictors      # %*% denotes matrix multiplication. cross multiply with lagged copies
  tmp = c(rev(r), r[-1])    # mirror 
  #make the bottom matrix in the matlab example
  w = t(sapply(seq(order + 1, 2, -1), function(x) tmp[(x):(x + order - 1)]))
  
  #Final computation of linear predictor coefficients
  coeffs = -r[2:(order + 1)] %*% solve(w)   #Compute coefficients
  coeffs = c(1, coeffs)                     #Add 1 to coefficients vectos (a[0] = 1)
  
  #Add plots
  if (show == TRUE & add == TRUE) 
    freqresponse(1, coeffs, fs = fs, add = add)
  if (show == TRUE & add == FALSE) {
    spectralslice(sound, fs = fs, col = 4)
    freqresponse(1, coeffs, fs = fs, add = TRUE)
  }
  
  #Return linear predictor coefficients
  coeffs
}