import numpy as np
from CRUDEFB import *
from atm62 import *

class THERMAL:
   def __init__(self,yld,hob,gnd,dump=False):

#  yld    yield (Kt)
#  hob    height of burst (m)
#  gnd    height of ground (m)
#  dump   true/false to write out thermal parameters

      self.fourpi  = 4.0 * np.pi

      self.atmos = atm62()

      self.fb    = CRUDEFB(yld,hob)

      self.w     = yld

      self.hob   = hob
      self.hobcm = hob * 100.0

      self.gnd   = gnd
      self.gndcm = gnd * 100.0

      self.rzero = self.atmos.returnConditions(self.gndcm)

      self.rho   = self.atmos.returnConditions(self.hobcm)

#  height to apply correction to thermal output (cm)

      self.hc = 5486.0 * self.w**0.40

      self.vis     = 25.7   # Km

#  convert visibility to cm

      self.vis     = self.vis * 1.0e5

      self.tolddt  = 0.0
      self.oldth   = 0.0
      self.intflux = 0.0

#  curvefits based on SPUTTER HYDRO runs
#  c.f. ADA093517
#       Estimates of Thermal Radiation Environments for Planning a Thermal Simulation on a HE Test
#       Burton S. Chambers
#       John A. Hasdal
#       SAIC 
#       30 March 1979

      self.pwrnames = ['null',\
            'time to principal maximum (sec)',\
            'time to principal minimum (sec)',\
            'power at principal minimum (watts)',\
            'curvefit coefficient c',\
            'time to shock exposure maximum (sec)',\
            'peak power at shock exposure or radiation max (watts)',\
            'shape smoothing factor',\
            'first pulse exponential scale factor (bnr)',\
            'first pulse exponential scale factor (bnl)',\
            'late time dip function scale factor (bmr)',\
            'late time dip function scale factor (bml)',\
            'peak amplitude of dip (watts)',\
            'time associated with peak amplitude of dip (sec)',\
            'power time scale factor prior to tmin (hintmn)',\
            'power time scale factor after tmin (hintmx)',\
            'curvefit coefficient s',\
            'power at 2nd max']
            
      self.pwr = [0.0] * 18
   
      sigma   = self.rho / self.rzero

      s            = 1.75995 * sigma**(-0.0432)
      self.pwr[6]  = 0.17e14 * sigma**(-0.6244) * self.w**0.53
      
      x       = np.log(sigma)
      self.pwr[5]  = (((-1.51933e-7 * x - 3.72132e-6) * x - 7.62326e-6) * x + 1.6828e-4) * self.w**0.3123
      
      self.pwr[10] = 2.0 * sigma**0.01008
      
      self.pwr[11] = sigma**0.1189
      
      tmax    = 0.03800 * self.w**0.44 * sigma**0.36
      tmin    = 0.00256 * self.w**0.39 * sigma**(-0.062)
   
      chk     = 0.99 * tmax
      if (tmin > chk):
         tmin = chk
      
      pmax    = 1.49e13 * self.w**0.59*sigma**(-0.45)
      pmin    = 6.82e11 * self.w**0.54*sigma**(-1.002)
      
      chk     = 0.99 * pmax
      if (pmin > chk):
         pmin = chk
      
      self.pwr[7]  = 0.0289 * self.w**0.3209 + 1.023 * self.w**(-0.0013979)
      
      self.pwr[8]  = 3.4071 * sigma**0.17827
      
      self.pwr[9]  = ((-0.007741 * x - 0.207) * x - 0.45) * x + 9.83382
      
      self.pwr[12] = 1.0 / (0.80 * sigma**0.3372) - 1.0
      
      self.pwr[13] = 1.523 * sigma**(-0.0752) * tmax
      
      tr      = tmin / tmax
      c       = np.log(tr) * tmax / (tmin - tmax)
      b       = s * np.exp(c) / c
      
      self.pwr[14] = pmin * (tr**s) * np.exp(b / np.exp(c * tr))
      
      self.pwr[15] = pmax * np.exp(b / np.exp(c))
      
      self.pwr[1]  = tmax
      
      self.pwr[2]  = tmin
      
      self.pwr[3]  = pmin
      
      self.pwr[4]  = c
      
      self.pwr[16] = s
   
      self.pwr[17] = pmax

#  build suggested time array based on above values

      self.ta = []

      t0 = math.log10(1.0e-5)
      t1 = math.log10(self.pwr[2])
      n  = 30
      dt = (t1 - t0) / (n - 1)
      t  = t0
      for i in range(n):
         self.ta.append(10.0**t)
         t = t + dt

      t0 = t1 + dt
      t1 = math.log10(2.0 * self.pwr[1])
      n  = 30
      dt = (t1 - t0) / (n - 1)
      t  = t0
      for i in range(n):
         self.ta.append(10.0**t)
         t = t + dt

      t0 = t1 + dt
      t1 = math.log10(100.0)
      n  = 50 
      dt = (t1 - t0) / (n - 1)
      t  = t0
      for i in range(n):
         self.ta.append(10.0**t)
         t = t + dt
      
      if (dump):
         fo = open('terms_thermal.out','w')
         fo.write('yield: %10.3f hob: %10.3f\n' % (yld,hob))
         fo.write('\n')

         for i in range(1,len(self.pwr)):
            fo.write('%5d %15.5e %s\n' % (i,self.pwr[i],self.pwrnames[i]))
         fo.close()

   def atATime(self,hdet,rdet,t):

#  hdet    height of detector (m)
#  rdet    distance from burst to detector (m)
#  t       time (sec)

      if (t <= 0.0):
         return (0.0,0.0,0.0)

      hdetcm = hdet * 100.0
      rdetcm = rdet * 100.0

#  get fireball at current time

      (rfb,rise) = self.fb.radiusHfb(t)

#  warn user if inside fireball

      if (rdetcm < rfb):
         print('at %15.5e (sec) detector inside fireball...results are meaningless!!!!!' % (t))

#  compute flux......

      r1     = self.pwr[13] / t
      bfac   = 2.0 * self.pwr[12] / (r1**self.pwr[9] + r1**(-self.pwr[8])) + 1.0

      bsem   = self.pwr[6] / self.pwr[3] - 1.0
      r2     = self.pwr[5] / t
      facm   = 2.0 * bsem / (r2**self.pwr[11] + r2**(-self.pwr[10])) + 1.0

      if (t < self.pwr[2]):
         thflux = self.pwr[3] * facm / bfac
      else:
         tsmax  = t / self.pwr[1]
         tr     = self.pwr[2] / self.pwr[1]
         hintmn = self.pwr[14]
         hintmx = self.pwr[15]

         if (tsmax >= 1.0):
            hint = hintmx
         else:
            cthet = np.cos(np.pi * (1.0 - np.log(tsmax) / np.log(tr)))
            hint  = hintmx + (hintmn - hintmx) * ((1.0 + cthet) / 2.0)**self.pwr[7]

         s = self.pwr[16]
         c = self.pwr[4]
 
         thflux = hint * facm / (bfac * (tsmax**s) * np.exp(s * np.exp(c * (1.0 - tsmax)) / c))

#  add in rise

      hfb  = self.hobcm + rise

      self.rfb  = rfb
      self.rise = rise
      self.hfb  = hfb

#  get density at current fireball height

      den  = self.atmos.returnConditions(hfb)

#  correction due to ground

      hbhp = hfb - self.gndcm + hdetcm - self.gndcm

      fcor = 1.0

      if (hbhp < self.hc):
         fcor = max(0.50,0.50 * (1.0 + hbhp / self.hc))

      dtovr2 = 0.0
      if (self.tolddt > 0.0):
         dtovr2 = 0.50 * (t - self.tolddt)

      self.tolddt = t

# compute slant range accounting for fireball rise (flat earth)

      srx = np.sqrt(hfb * hfb + rdetcm * rdetcm)

      sr   = max(rfb,srx)
      sr2  = sr * sr
      rfb2 = rfb * rfb

#  range term.... at large distance this looks just like 1 / (4 pi r2)

      rf  = 2.0 * (1.0 - np.sqrt(1.0 - rfb2 / sr2)) / (self.fourpi * rfb2)

#  account for transmission thru atmosphere

      rhor = sr * den
      tr   = max(0.01,1.0 - 875.0 * rhor / self.vis)

#  4.184 watts -> cal

      th   = fcor * tr * thflux / 4.184 * rf

      abc  = 0.0
      if (dtovr2 > 0.0):
         abc = dtovr2 * (self.oldth + th)

      self.oldth   = th

      self.intflux = self.intflux + abc

      return (th,thflux,self.intflux)
