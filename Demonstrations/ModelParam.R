LHCSIZE=10
NLHC=3
NSCMS=30
NPARA=3
param.names=c("A_U","A_EPSILON","A_T")
param.lows=c(0.01,0.1,0.01)
param.highs=c(0.4,3.,1.)
param.defaults=c(0.126,0.85,0.14)
which.logs<-c()
  param.defaults <- param.defaults[1:NPARA]
  param.highs <- param.highs[1:NPARA]
  param.lows <- param.lows[1:NPARA]
  param.names <- param.names[1:NPARA]
