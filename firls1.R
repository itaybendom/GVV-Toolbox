# R function for firls() function in MATLAB
# Code structure 
# The Octave script was verified to produce the same results as MATLAB
#http://octave.cvs.sourceforge.net/viewvc/octave/octave-forge/main/signal/inst/firls.m?content-type=text%2Fplain&revision=HEAD
#http://eeglab.googlecode.com/svn/!svn/bc/9956/eeglab/functions/octavefunc/signal/firls.m

# Theory for Linear-Phase FIR Filter Design by Least Squares
#http://cnx.org/contents/6x7LNQOp@7/Linear-Phase-Fir-Filter-Design

# MATLAB function call:
# b = firls(N, F, A);
# b = firls(N, F, A, W);

#  FIR filter design using least squares method. Returns a length N+1
#  linear phase filter such that the integral of the weighted mean
#  squared error in the specified bands is minimized.
#
#  F specifies the frequencies of the band edges, normalized so that
#  half the sample frequency is equal to 1.  Each band is specified by
#  two frequencies, to the vector must have an even length.
#
#  A specifies the amplitude of the desired response at each band edge.
#
#  W is an optional weighting function that contains one value for each
#  band that weights the mean squared error in that band. A must be the
#  same length as F, and W must be half the length of F.
#
# The least squares optimization algorithm for computing FIR filter
# coefficients is derived in detail in:
#
# I. Selesnick, 'Linear-Phase FIR Filter Design by Least Squares,'
# http://cnx.org/content/m10577

# F = frequencies
# A = pass
# W = weight
#
# For R make sure inputs are in matrix form.

firls1 <- function (N, frequencies, pass, weight){
  
  # Debugging
  #
  if (nargs() == 3){
    weight = ones(1, length(pass)/2);
  }
  if (length (frequencies) != length (pass)){
    stop("frequencies and pass must have equal lengths.")
  }
  if (2 * length (weight) != length (pass)){
    stop("weight must contain one weight per band.");
  }
 
  # Code start
  #
  M = N/2
 
  tmp = matrix(c(-1, 1), nrow=2)  # Y matrix for Kronecker product 
  w = kron(t(weight), tmp)        # Kronecker product 
  
  omega = matrix(frequencies * pi, nrow=1)
  
  i1 = seq(1, length(omega), 2)
  i2 = seq(2, length(omega), 2)
  
  # Generate the matrix Q
  # The matrix can be expressed as the sum of a Hankel and Toeplitz matrix.
  # A factor of 1/2 has been dropped and the final filter coefficients multiplied by 2 to compensate.
  #
  # The matrix Q1 is a symmetric Toeplitz matrix (constant along its diagonals)
  # The matrix Q2 is a Hankel matrix (constant along its anti-diagonals). 
  # Consequently, the matrices can be stored with less memory than arbitrary matrices.
  # There are fast algorithms to compute the solution to 'Toeplitz plus Hankel' systems with computational complexity O(M^2) instead of O(M^3). 
  
  cos_ints11 = omega
  cos_ints12 = sin(matrix(1:N) %*% omega)
  cos_ints1 = rbind(cos_ints11, cos_ints12)
  
  q1 = matrix(c(1, 1/(1:N)))
  q2 = cos_ints1 %*% w
  q =  q1 * q2
  
  # Compute the Q matrix  
  # Q1 = toeplitz(q[1:(M+1)])
  # Q2 = hankel(q[1:(M+1)], q[(M+1):length(q)])
  Q = toeplitz(q[1:(M+1)]) + hankel(q[1:(M+1)], q[(M+1):length(q)]);
  
  # The vector b is derived from solving the integral:
  #
  #           _ w
  #          /   2
  #  b  =   /       W(w) D(w) cos(kw) dw
  #   k    /    w
  #       -      1
  #
  # Since we assume that W(w) is constant over each band (if not, the
  # computation of Q above would be considerably more complex), but
  # D(w) is allowed to be a linear function, in general the function
  # W(w) D(w) is linear. The computations below are derived from the
  # fact that:
  #     _
  #    /                          a              ax + b
  #   /   (ax + b) cos(nx) dx =  --- cos (nx) +  ------ sin(nx)
  #  /                             2                n
  # -                             n
  #
  
  cos_ints21 = omega[i1]^2 - omega[i2]^2 
  cos_ints22 = cos(matrix(c(1:M), ncol = 1) %*% omega[i2]) - cos(matrix(c(1:M), ncol = 1) %*% omega[i1])
  cos_ints23 = (matrix(c(2,1:M), ncol=1) %*% (omega[i2] - omega[i1]))
  cos_ints24 = rbind(cos_ints21,cos_ints22)
  cos_ints2 = cos_ints24/cos_ints23
  
  d1 = -weight * pass[i1]
  d2 =  weight * pass[i2]
  dmat = rbind(d1, d2)
  d = matrix(dmat,ncol=1)
  
  b1 = matrix(c(1, 1./(1:M)))
  b2 = ((kron(cos_ints2, matrix(c(1, 1),nrow = 1)) + cos_ints1[1:(M+1),]) %*% d);
  b = b1 * b2
  
  # Having computed the components Q and b of the  matrix equation, 
  # solve for the filter coefficients.
  a = mldivide(Q,b)
  
  #break matrix rows:
  coef1 = matrix(rev(a[-1]), ncol = 1)
  coef2 = 2 * a[1]
  coef3 = matrix(a[-1], ncol = 1)
  
  #results in ->
  coef = rbind(coef1, coef2, coef3)   # return values
  
}