#####################################################################################################################
#
# 5/12/2016
# This is an implementation of the jitter_method.R function, except that it accepts as an input the pre-computed 
# glottal signal, using the glottal_source_function.R method
#
# The following function computes cycle-to-cucle variations in fundamental frequency, also known as Jitter.
#
# Jitter Calculation
# http://www.jvoice.org/article/S0892-1997(13)00234-8/pdf
# This is the coputation of Jitter (relative): http://www.cs.upc.edu/~nlp/papers/far_jit_07.pdf
#
# Created 5/12/2016
# By Itay Ben-Dom
#####################################################################################################################
#####################################################################################################################

jitter_function <- function(glottalList){
  
  # First extract all the applicable information from the list generated via glottal_source_function.R
  valIdx  = glottalList[[1]]$valIdx   # vector of valley instances
  
  # Jitter calculation

  
  # x = vocal fold period vector
  x = diff(valIdx)
  # N = total number of vocal fold periods
  N = length(x)
  # Average period = sum(x)/N
  avgT = mean(x) 
  
  # Jitter equation
  jitterVal = (100 * ( 1 / (N-1) ) * sum(abs(diff(x))) ) / avgT
  
  return(list(jitterVal))
  
}