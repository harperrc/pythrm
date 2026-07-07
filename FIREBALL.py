import math
from atm62 import *

class FIREBALL:

#  compute the diameter of a fireball using a chart of scaled radius vs scaled time
#  from AD0680842 Calculation of Retinal Burn....
#  used engague to digitize the points then used STEPIT for fit a curve (see FD)
#  this version output in meters

#  rise model taken from (30 yr old) STLAMB runs

   def __init__(self,W,HOB):

#  W    yield (Kt)
#  HOB  height of burst (Km)

      self.fuzzyGraphFix = 10.0

      self.atm = atm62()

#  fits to Figure 9 compute DFB in ?? i think the plot is mislabeled
#  it appears the y axis is meters vs cm (compares favorabily with other data)

      self.p0 = [  1.0000000E+04,  1.2651851E-01, -6.3315969E-03]
      self.p1 = [  4.2439795E+04,  3.1559732E-01, -6.5007706E-03]
      self.p2 = [  3.6159920E+04,  1.8173914E-01, -3.4663785E-02]

      (temp,pres,self.rho0,sos,visc,visck,fmp,g,hs) = self.atm.returnConditions(0.0)
      (temp,pres,self.rhob,sos,visc,visck,fmp,g,hs) = self.atm.returnConditions(HOB)

      self.t2max = 0.037 * (self.rhob / self.rho0)**0.282 * W**0.47

#      self.con   = 0.3504 * W**0.35 * (self.rho0 / self.rhob)**0.18

#  removed the 0.3504 since a 1Kt at 0 hob should match figure 9 and not be
#  0.3504 smaller

      self.con   = W**0.35 * (self.rho0 / self.rhob)**0.18

#  computed average over 1 second intervals over yields and height of burst
#  now get the average growth rate of the fireball
#  we will decay the rate over time (just a little)

      HOBm  = HOB * 1000.0

      self.HOBcm = HOBm * 100.0

      self.growthDecay = -0.002
      a0               = 74.18 + 0.0056 * HOBm -2.873e-7 * HOBm * HOBm 
      self.growthRate  =  a0 * W**(0.1857 + 0.0161 * math.log(W))

#  computed average over 1 second intervals over yields and height of burst
#  now get the average rise rate of the fireball
#  we will decay the rate over time (just a little)

      self.riseDecay = -0.002
      a0             =  1440.87 + 0.0237 * HOBm -3.879e-6 * HOBm * HOBm 
      self.riseRate  =  a0 * W**(0.250 -3.0e-3 * math.log(W))
      
#  hack in a minimum fireball size using initial time point in fit

      tstar = 1.0e-4
      tmp   = self.p0[0] * math.exp(self.p0[1] * math.log(tstar) + self.p0[2] * math.log(tstar)**2)

      self.minDFB = tmp * self.con

   def FD(self,T):
      tx    = max(0.0,min(T,20.0 * self.t2max))
      tstar = tx / self.t2max

#  broke the scale radius vs scaled time into 3 pieces.  

      if (tstar < 0.0005):
         fd = self.p0[0] * math.exp(self.p0[1] * math.log(tstar) + self.p0[2] * math.log(tstar)**2)
      elif (tstar <= 0.1236):
         fd = self.p1[0] * math.exp(self.p1[1] * math.log(tstar) + self.p1[2] * math.log(tstar)**2)
      elif (tstar >= 0.0005):
         fd = self.p2[0] * math.exp(self.p2[1] * math.log(tstar) + self.p2[2] * math.log(tstar)**2)

#  original STLAMB fit include the burst height so include here then
#  subtract out to just return rise and convert to meters

      hfb  = (self.HOBcm + self.riseRate * T) * T**self.riseDecay
      rise = (hfb - self.HOBcm) / 100.0

      return (max(self.minDFB,self.con * fd),rise)
