#!/usr/bin/python3

import math
import numpy as np
import sys

class CRUDEFB:
   def __init__(self,W,H):

#  this is a crude fit to STLAMB data that was generated about 30+ years ago
#  using the McDonnell Douglas version of LAMB

#  somewhat representative of the hot region of a fireball

#  fits were over the following heights of burst (meters)
#      0 10 50 100 200 500 1000 2000 5000 10000
#  and yields (Kt)
#      1 5 10 20 50 100 200 500 1000 5000 10000 20000

#  N.B. return radius and rise in cm

#  W  yield (Kt)
#  H  height of burst (m)

      self.H    = H
      self.Hcm  = H * 100.0
      self.W    = W

#  start with the initial radius that was fit to height of burst (H)
#  and yield (W)

      a0      = 3421.37 + 0.05483 * H + 2.446e-6 * H * H

      self.r0 = a0 * W**0.22

#  computed average over 1 second intervals over yields and height of burst
#  now get the average growth rate of the fireball
#  we will decay the rate over time (just a little)

      self.growthDecay = -0.002
      a0               = 74.18 + 0.0056 * H -2.873e-7 * H * H
      self.growthRate  =  a0 * W**(0.1857 + 0.0161 * math.log(W))

#  computed average over 1 second intervals over yields and height of burst
#  now get the average rise rate of the fireball
#  we will decay the rate over time (just a little)

      self.riseDecay = -0.002
      a0             =  1440.87 + 0.0237 * H -3.879e-6 * H * H
      self.riseRate  =  a0 * W**(0.250 -3.0e-3 * math.log(W))

   def radiusHfb(self,age):

      if (age <= 0.0):
         return (self.r0,self.Hcm)
    
      radius = (self.r0  + self.growthRate * age) * age**self.growthDecay
      hfb    = (self.Hcm + self.riseRate * age) * age**self.riseDecay

      return (radius,hfb)
