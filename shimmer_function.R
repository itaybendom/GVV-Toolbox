#####################################################################################################################
#
# 5/12/2016
# This is an implementation of the shimmer_method.R function, except that it accepts as an input the pre-computed 
# glottal signal, using the glottal_source_function.R method
#
# The following function computes cycle-to-cucle variations in peak amplitude, also known as Shimmer.
#
# Shimmer Calculation
# http://www.jvoice.org/article/S0892-1997(13)00234-8/pdf
# This is the coputation of Jitter (relative): http://www.cs.upc.edu/~nlp/papers/far_jit_07.pdf
#
# Created 5/12/2016
# By Itay Ben-Dom
#####################################################################################################################
#####################################################################################################################

shimmer_function <- function(glottalList){
  
  # First extract all the applicable information from the list generated via glottal_source_function.R
  peakAmp  = glottalList[[1]]$peakAmp  # vector of peak amplitude values
  
  # Shimmer calculation
  
  # x = vocal fold amplitude vector
  x = (peakAmp)
  # N = total number of peak amplitudes
  N = length(x)
  # Average amplitude = sum(x)/N
  avgAmp = mean(x) 
  
  # Shimmer equation
  shimmerVal = (100 * ( 1 / (N-1) ) * sum(abs(diff(x))) ) / avgAmp 
  
  return(list(shimmerVal))
  
}