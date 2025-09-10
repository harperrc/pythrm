#!/usr/bin/python3

import math
import numpy as np
import sys
from thermal import *

whichrun = 2

if (whichrun == 0):
   yld  =   10.0
   hob  =    0.0
   
   yld  =  20000.0
   hob  =      0.0
   
   gnd  =    0.0
   hdet =    0.0
   
   rdet = 10000.0
   
   th = THERMAL(yld,hob,gnd,True)
   
   #  get the suggested time array
   
   ta = th.ta
   
   fo = open('out.thermal','w')
   
   for t in ta:
      (a,b,c) = th.atATime(hdet,rdet,t)
      fo.write('%15.5e %15.5e %15.5e %15.5e %15.5f %15.5f\n' % (t,a,b,c,th.rfb,th.hfb))
   
   fo.close()

if (whichrun == 1):
   ylda  = [1,10,100,1000,10000]
   hoba  = [0]
   rdeta = [1.41,10]

   gnd   = 0.0
   hdet  = 0.0

   for yld in ylda:
      for hob in hoba:
         for rdetkm in rdeta:
            rdet = rdetkm * 1000.0

            th = THERMAL(yld,hob,gnd,True)
            ta = th.ta

            oname = '%d_%d_%.2f.out' % (yld,hob,rdetkm)
            fo = open(oname,'w')
       
            for t in ta:
               (a,b,c) = th.atATime(hdet,rdet,t)
               fo.write('%15.5e %15.5e %15.5e %15.5e %15.5f %15.5f\n' % (t,a,b,c,th.rfb,th.hfb))
            fo.close()

if (whichrun == 2):
   fo = open('t2max.out','w')
   ylda  = [1,5,10,20,100,200,500,1000,5000,10000]
   ylda  = np.logspace(0,4,256)
   hoba  = [0,3000]

   gnd   = 0.0
   hdet  = 0.0

   for hob in hoba:
      for yld in ylda:
         th = THERMAL(yld,hob,gnd,False)
         fo.write('%10.3f %10.3f %15.5e %s\n' % (yld,hob,th.pwr[1],th.pwrnames[1]))
      fo.write('\n')
