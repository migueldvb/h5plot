#!/usr/bin/python
"""
Read namelist
"""

import sys
sys.path.append("/home/miguel/src/fortran-namelist")

import namelist

from plplot import *

from scipy import sqrt

plinit()

plenv(0,100,0,2.e-3,0,2)

for i in range(1,99):
   filename = "%s%02u" % (sys.argv[1], i)
   print filename
   mynmlfile = namelist.Namelist(filename)

#    mynmlfile = namelist.Namelist(sys.argv[1])

#    print mynmlfile

#    print mynmlfile["POINTMASS"]["par"]
#    print mynmlfile["POINTMASS1"]["par"][0]["angmom"][-1]

   angmom0 = float(mynmlfile["POINTMASS"]["par"][0]["angmom"][-1])
   m0 = float(mynmlfile["POINTMASS"]["par"][0]["m"][0])
   r0 = mynmlfile["POINTMASS"]["par"][0]["r"]
   ekin0 = float(mynmlfile["POINTMASS"]["par"][0]["ekin"][0])

   angmom1 = float(mynmlfile["POINTMASS1"]["par"][0]["angmom"][-1]) 
   m1 = float(mynmlfile["POINTMASS1"]["par"][0]["m"][0])
   r1 = mynmlfile["POINTMASS1"]["par"][0]["r"]
   ekin1 = float(mynmlfile["POINTMASS1"]["par"][0]["ekin"][0])

   angmom2 = float(mynmlfile["POINTMASS2"]["par"][0]["angmom"][-1]) 
   m2 = float(mynmlfile["POINTMASS2"]["par"][0]["m"][0])
   r2 = mynmlfile["POINTMASS2"]["par"][0]["r"]
   ekin2 = float(mynmlfile["POINTMASS2"]["par"][0]["ekin"][0])

   etot = ekin0 + ekin1 +ekin2 - \
   m0*m1/sqrt( (float(r0[0])-float(r1[0]))**2 + (float(r0[1])-float(r1[1]))**2 )
   print i, etot, angmom0+angmom1

   plpoin([i],[angmom0],2)
   plpoin([i],[angmom1],3)
   plpoin([i],[angmom2],3)
   plpoin([i],[angmom0+angmom1],4)

plend()
