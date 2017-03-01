ellipordR <-function (Wp, Ws, Rp, Rs) 
{
  if (!(length(Wp) <= 2 && length(Ws) <= 2)) 
    stop("ellipord: frequency must be given as w0 or [w0, w1]")
  if (!all(Wp >= 0 & Wp <= 1 & Ws >= 0 & Ws <= 1)) 
    stop("ellipord: critical frequencies must be in (0 1)")
  if (!(length(Wp) == 1 || length(Wp) == 2 & length(Ws) == 
        1 || length(Ws) == 2)) 
    stop("ellipord: only one filter band allowed")
  if (length(Wp) == 2 && !(Wp[1] < Wp[2])) 
    stop("ellipord: first band edge must be smaller than second")
  if (length(Ws) == 2 && !(length(Wp) == 2)) 
    stop("ellipord: you must specify band pass borders.")
  if (length(Wp) != 1 && length(Wp) != 2) 
    stop("ellipord: Wp,Ws must have length 1 or 2")
  if (length(Wp) == 2 && length(Ws) == 2 && !(Ws[1] < Wp[1] && 
                                              Ws[2] > Wp[2])) 
    stop("ellipord: ( Wp[1], Wp[2] ) must be inside of interval ( Ws[1], Ws[2] )")
  if (length(Wp) == 2 && length(Ws) == 1 && !(Ws < Wp[1] || 
                                              Ws > Wp[2])) 
    stop("ellipord: Ws must be out of interval ( Wp[1], Wp[2] )")
  T <- 2
  Wpw <- tan(pi * Wp/T)
  Wsw <- tan(pi * Ws/T)
  if (length(Wpw) == 2 && length(Wsw) == 2) {
    type <- "pass"
    wp <- 1
    w02 <- Wpw[1] * Wpw[2]
    w3 <- w02/Wsw[2]
    w4 <- w02/Wsw[1]
    if (w3 > Wsw[1]) {
      ws <- (Wsw[2] - w3)/(Wpw[2] - Wpw[1])
    }
    else if (w4 < Wsw[2]) {
      ws <- (w4 - Wsw[1])/(Wpw[2] - Wpw[1])
    }
    else {
      ws <- (Wsw[2] - Wsw[1])/(Wpw[2] - Wpw[1])
    }
  }
  else if (length(Wpw) == 2 && length(Wsw) == 1) {
    type <- "pass"
    wp <- 1
    w02 <- Wpw[1] * Wpw[2]
    if (Wsw > Wpw[2]) {
      w3 <- w02/Wsw
      ws <- (Wsw - w3)/(Wpw[2] - Wpw[1])
    }
    else {
      w4 <- w02/Wsw
      ws <- (w4 - Wsw)/(Wpw[2] - Wpw[1])
    }
  }
  else {
    type <- if (Wpw < Wsw) 
      "low"
    else "high"
    wp <- Wpw
    ws <- Wsw
  }
  k <- wp/ws
  if (type == "high")
    k <- ws/wp
  k1 <- sqrt(1 - k^2)
  q0 <- (1/2) * ((1 - sqrt(k1))/(1 + sqrt(k1)))
  q <- q0 + 2 * q0^5 + 15 * q0^9 + 150 * q0^13
  D <- (10^(0.1 * Rs) - 1)/(10^(0.1 * Rp) - 1)
  n <- ceiling(log10(16 * D)/log10(1/q))
  FilterOfOrder(n = n, Wc = Wp, type = type, Rp = Rp, Rs = Rs)
}