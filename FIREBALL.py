import math
from atm62 import *

class FIREBALL:

#  compute the diameter of a fireball using a chart of scaled radius vs scaled time
#  from AD0680842 Calculation of Retinal Burn....
#  used engague to digitize the points then used STEPIT for fit a curve (see FD)
#  this version output in meters
#  rise model taken from (30 yr old) SDLAMB runs

   def __init__(self,W,HOB):

#  W    yield (Kt)
#  HOB  height of burst (Km)

      self.fuzzyGraphFix = 10.0

      self.atm   = atm62()
      self.p1 = [  3.6046156E+04,  1.8216787E-01, -3.3106571E-02]
      self.p0 = [  4.4965491E+04,  3.5498392E-01, -1.1305412E-03]

      (temp,pres,self.rho0,sos,visc,visck,fmp,g,hs) = self.atm.returnConditions(0.0)
      (temp,pres,self.rhob,sos,visc,visck,fmp,g,hs) = self.atm.returnConditions(HOB)

      self.t2max = 0.037 * (self.rhob / self.rho0)**0.282 * W**0.47
      self.con   = 0.3504 * W**0.35 * (self.rho0 / self.rhob)**0.18
      self.con   = self.con * self.fuzzyGraphFix

#  output in meters

      self.con   = 1.0e-2 * self.con

      HOBm       = HOB * 1000.0
      self.HOBcm = HOBm * 100.0

#  computed average over 1 second intervals over yields and height of burst
#  now get the average growth rate of the fireball
#  we will decay the rate over time (just a little)

      self.growthDecay = -0.002
      a0               = 74.18 + 0.0056 * HOBm -2.873e-7 * HOBm * HOBm 
      self.growthRate  =  a0 * W**(0.1857 + 0.0161 * math.log(W))

#  computed average over 1 second intervals over yields and height of burst
#  now get the average rise rate of the fireball
#  we will decay the rate over time (just a little)

      self.riseDecay = -0.002
      a0             =  1440.87 + 0.0237 * HOBm -3.879e-6 * HOBm * HOBm 
      self.riseRate  =  a0 * W**(0.250 -3.0e-3 * math.log(W))

   def FD(self,T):
      tx    = max(0.0,min(T,10.0 * self.t2max))

      tstar = tx / self.t2max

#  broke the scale radius vs scaled time into 2 pieces.  

      if (tstar >= 0.03):
         fd = self.p1[0] * math.exp(self.p1[1] * math.log(tstar) + self.p1[2] * math.log(tstar)**2)
      else:
         fd = self.p0[0] * math.exp(self.p0[1] * math.log(tstar) + self.p0[2] * math.log(tstar)**2)

#  original fit include the burst height so include here then
#  subtract out to just return rise and convert to meters

      hfb  = (self.HOBcm + self.riseRate * T) * T**self.riseDecay
      rise = (hfb - self.HOBcm) / 100.0

      return (self.con * fd,rise)
