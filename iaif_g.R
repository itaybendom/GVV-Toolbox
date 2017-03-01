# IAIF Glottal Inverse Filtering
#   [G,DG,A,AG] = IAIF(X,FS,P_VT,P_GL,D,HPFILT) estimates the glottal
#   volume velocity waveform from a speech signal using Iterative Adaptive
#   Inverse Filtering (IAIF).
#
#   INPUTS: 
#   x      - Speech signal
#   B      - FIR least-squares Highpass filter coefficients
#
#
#   OUTPUTS: 
#   g      - Glottal volume velocity waveform
#
# Reference:
#
# P. Alku, "Glottal wave analysis with pitch synchronous iterative adaptive
# inverse filtering", Speech Commun., vol. 11, no. 2-3, pp. 109-118, 1992.
#
# Tuomo Raitio, 20.09.2011
# Revised 12.08.2012
# Revised 10.07.2013

iaif_g <- function (x, B){

  # Only outputs Glottal volume velocity waveform.
  # iaif_modified.R is a better implementation of iaif.R
  # However, a faster implementation is iaif_g.R
  # Since B becomes an input argument, we don't have to compute it everytime, which
  # shaves about .08 sec off the computation time.
  
  
  # Set default parameters
  # Debug number of function input arguments 
  fs = 16000
  hpfilt = 1
  d = 0.99
  p_gl = 2*round(fs/4000);
  p_vt = 2*round(fs/2000)+4;
  
  preflt = p_vt+1;

  Fstop = 40                 # Stopband Frequency
  Fpass = 70                 # Passband Frequency
  Nfir = round(300/16000*fs) # FIR numerator order
  if (mod(Nfir,2) == 1){
    Nfir = Nfir + 1
  }
  
  # ONly need to calculate B once
  x = filter(B, 1, append(x, zeros(rnd(length(B)/2)-1, 1)))
  x = x[rnd(length(B)/2):length(x)]

  # rbind(matrix(linspace(-x[1], x[1], preflt), ncol = 1) , x)
  x1 = append(linspace(-x[1], x[1], preflt) ,x)
  # Estimate the combined effect of the glottal flow and the lip radiation
  # (Hg1) and cancel it out through inverse filtering. Note that before
  # filtering, a mean-normalized pre-frame ramp is appended in order to
  # diminish ripple in the beginning of the frame. The ramp is removed after
  # filtering.
  
  Hg1 = lpc_modified(x, 1, preemph = FALSE)
  y = filter(Hg1, 1, x1)[-1:-(preflt)]  #at this point length(y) should equal length(x)
  
  # Estimate the effect of the vocal tract (Hvt1) and cancel it out through
  # inverse filtering. The effect of the lip radiation is canceled through
  # intergration. Signal g1 is the first estimate of the glottal flow.
  
  Hvt1 = lpc_modified(y, p_vt, preemph = FALSE)
  #g1 = filter(Hvt1, 1, x1)
  g1 = filter(1, c(1, -d), filter(Hvt1, 1, x1))[-1:-(preflt)]
  
  # Re-estimate the effect of the glottal flow (Hg2). Cancel the contribution
  # of the glottis and the lip radiation through inverse filtering and
  # integration, respectively.
  
  Hg2 = lpc_modified(g1, p_gl, preemph = FALSE)
  #y = filter(Hg2, 1, x1) 
  y = filter(1, c(1, -d), filter(Hg2, 1, x1))[-1:-(preflt)]
  
  # Estimate the model for the vocal tract (Hvt2) and cancel it out through
  # inverse filtering. The final estimate of the glottal flow is obtained
  # through canceling the effect of the lip radiation.
  
  Hvt2 = lpc_modified(y, p_vt, preemph = FALSE);
  dg = filter(Hvt2, 1, x1)
  g = filter(1, c(1, -d), dg)[-1:-(preflt)]    
  
 
  return(g)

}