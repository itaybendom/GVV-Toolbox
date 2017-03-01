peak_detection <- function(wave, type, tolerance){
  
  # Adapted from MATLAB's V_FINDPEAKS algorithm.
  #V_FINDPEAKS finds peaks with optional quadratic interpolation [K,V]=(X,M,W)
  #
  #  Inputs:  wave        is the input signal (does not work with UInt datatype)
  #           type        is mode:
  #                       'q' performs quadratic interpolation
  #                       'v' finds valleys instead of peaks
  #           toletance   is the width tolerance; a peak will be eliminated if there is
  #                       a higher peak within +-W samples
  #
  # Outputs:  K        are the peak locations in X (fractional if M='q')
  #           V        are the peak amplitudes: if M='q' the amplitudes will be interpolated
  #                    whereas if M~='q' then V=X(K). 
  
  # Outputs are column vectors regardless of whether X is row or column.
  # If there is a plateau rather than a sharp peak, the routine will place the
  # peak in the centre of the plateau. When the W input argument is specified,
  # the routine will eliminate the lower of any pair of peaks whose separation
  # is <=W; if the peaks have exactly the same height, the second one will be eliminated.
  # All peak locations satisfy 1<K<length(X).
  #
  # If no output arguments are specified, the results will be plotted.
  #
  
  #	   Copyright (C) Mike Brookes 2005
  #      Version: $Id: v_findpeaks.m 3601 2013-10-11 15:27:30Z dmb $
    #
  #   VOICEBOX is a MATLAB toolbox for speech processing.
  #   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
  #
  ###################################################################################
  #   This program is free software; you can redistribute it and/or modify
  #   it under the terms of the GNU General Public License as published by
  #   the Free Software Foundation; either version 2 of the License, or
  #   (at your option) any later version.
  #
  #   This program is distributed in the hope that it will be useful,
  #   but WITHOUT ANY WARRANTY; without even the implied warranty of
  #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #   GNU General Public License for more details.
  #
  #   You can obtain a copy of the GNU General Public License from
  #   http://www.gnu.org/copyleft/gpl.html or by writing to
  #   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
  ###################################################################################
  
  x = wave
  nx = length(x)
  
  # Check the Mode for required pattern
  # If looking for valleys, the wave is inverted in polarity. Now valleys are peaks
  if (grepl("v", type)){
    x = -x
  } else{
    x = x
  }
  
  dx = diff(x)
  r=which(dx>0)
  f=which(dx<0)
  
  lr = length(r)
  lf = length(f)
  
  if (lr>0 & lf>0){  # we must have at least one rise and one fall
    
    dr=r
    dr[-1]= diff(r)
    
    rc = repmat(1,nx,1)
    rc[r+1] = 1-dr 
    rc[1] = 0
    rs = cumsum(rc) # = time since the last rise
    
    df = f
    df[-1]= diff(f)
    
    fc = repmat(1,nx,1)
    fc[f+1] = 1-df
    fc[1] = 0
    fx = cumsum(fc) # = time since the last fall
    
    rp = repmat(-1,nx,1)
    rp[c(1,r+1)] = matrix(c(dr-1, nx-r[lr]-1), ncol=1)
    rq = cumsum(rp)  # = time to the next rise
    
    fp = repmat(-1,nx,1);
    fp[matrix(c(1,f+1),ncol=1)] = matrix(c(df-1, nx-f[lf]-1), ncol=1)
    fq = cumsum(fp) # = time to the next fall
    
    # k is index, v is amplitude
    k = which((rs<fx) & (fq<rq) & (floor((fq-rs)/2)==0));   # the final term centres peaks within a plateau
    v = x[k]
    
    
    if (grepl("q", type)){ # do quadratic interpolation
      
      b=0.5*(x[k+1]-x[k-1]);
      a=x[k]-b-x[k-1];
      j=(a>0);            # j=0 on a plateau
      v[j]=x[k[j]]+0.25*b[j]^2/a[j];
      k[j]=k[j]+0.5*b[j]/a[j];
      k[!j]=k[!j]+(fq[k[!j]]-rs[k[!j]])/2;    # add 0.5 to k if plateau has an even width
    }         
   
    
    w = tolerance
    if (length(k) > 1){
      # now purge nearby peaks
      j = which(diff(k)<=w)
      while (any(j)){
        j = j + (v[j] >= v[j+1])
        k = k[-j]
        v = v[-j]
        j = which(diff(k)<=w)
      }
    }
    
    
  } else{
    k=c()
    v=c()
  }
  
  if (grepl("v", type)){
    v = -v  # Invert value if searching for valleys
  }
  
  return(list("amplitude" = v, "index" = k))
  
  
  
  
}