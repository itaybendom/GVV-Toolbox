SRH <- function(specMat, nHarmonics, f0min, f0max) {
  # SRH (Summation of Residual Harmonics)
  # used in pitch_srh
  # Function to compute Summation of Residual harmonics function
  # on a spectrogram matrix, with each column corresponding to one
  # spectrum.
  #
  # Used in pitch_srh.m file. Part of gci_detection.m function.
  #
  # Settings
  #   Spectrogram Matrix (specMat) is computed in pitch_srh.m
  #   F0min = 80  # Minimum F0 set to 80 Hz
  #   F0max = 500
  #   nHarmonics = 5
  #
  # References
  #  [1] T.Drugman, A.Alwan, "Joint Robust Voicing Detection and Pitch Estimation
  #    Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011.
  #      Publication available at the following link:
  #      http://tcts.fpms.ac.be/~drugman/files/IS11-Pitch.pdf
  #
  # Copyright (c) 2011 University of Mons, FNRS
  #
  # License
  #  This code is a part of the GLOAT toolbox with the following
  #  licence:
  #  This program is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.
  #
  # This function is part of the Covarep project: http://covarep.github.io/covarep
  #
  # Author
  #  Thomas Drugman thomas.drugman@umons.ac.be
  #
  # Modified
  #  John Kane kanejo@tcd.ie September 27th 2014 - Bug fix and efficiency
  
  # Initial settings
  # MATLAB: [M,N] = size( specMat );
  #M = nrow(specMat)
  N = ncol(specMat)
  
  # Initialize matrix
  SRHmat = zeros(f0max, N)
  
  # Initialize frequency components
  fSeq = f0min:f0max      # F0 frequency sequence
  fLen = length(fSeq)     # length of sequence
  
  # Prepare harmonic indeces matrices
  # plusidx corresponds to nHarmonics size summation of E(k*f)
  plusIdx = repmat(matrix((1:nHarmonics), ncol = 1), 1, fLen) * repmat(fSeq, nHarmonics, 1)
  # subtridx corresponds to nHarmonics size summation of E((k-0.5)*f)
  subtrIdx = rnd(repmat( matrix((1:(nHarmonics - 1)) + 0.5, ncol = 1), 1, fLen) * repmat(fSeq, (nHarmonics - 1), 1))
  
  # The following two lines aren't included in MATLAB, but Required in R.
  # Indexing a matrix with a matrix doesn't automatically produce a matrix
  # Define row number for matrices.
  plusIdxRow = nrow(plusIdx)
  subtrIdxRow = nrow(subtrIdx)
  
  # Do harmonic summation (OLD)
  #for (n in 1:N){
  #  specMatCur = repmat((matrix(specMat[,n], ncol = 1)), 1, fLen)
  #  SRHmat[fSeq,n] = t(apply(matrix(specMatCur[plusIdx], nrow = plusIdxRow ), 2, sum)) #- apply(matrix(specMatCur[subtrIdx], nrow = subtrIdxRow), 2, sum))
  # Above - Needs to be converted into matrices to have same dimansions as MATLAB
  }
  
  # Do harmonic summation (NEW)
  # New implementation 14/9/16
  #aa = lapply( 1:ncol(specMat) , function(m) replicate(fLen, specMat[,m]))
  #for (n in 1:N){
  #  SRHmat[fSeq,n] = t( colSums(matrix(aa[[n]][plusIdx], nrow = plusIdxRow )) - colSums(matrix(aa[[n]][subtrIdx], nrow = subtrIdxRow)) )
  #}
  
  # Do harmonic summation (NEW - FASTEST) 14/9/16
  for (n in 1:ncol(specMat)) {
    SRHmat[fSeq, n] = t(unlist(lapply(1:ncol(plusIdx) , function (m)  sum(specMat[, n][plusIdx[, m]]))) - 
                          unlist(lapply(1:ncol(subtrIdx) , function (m) sum(specMat[, n][subtrIdx[, m]]))))
  }
  
  # Retrieve f0 and SRH value
  templist = maxR(SRHmat)
  # New implementation 14/9/16
  #SRHVal = templist$max
  #F0 = templist$index
  
  return(list("SRHVal" =  templist$max, "F0" = templist$index))
  # Same results as MATLAB
  # WORKS
}