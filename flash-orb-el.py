#!/usr/bin/python
"""
Plot orbital parameters of flash n-body simulations.
"""

import getopt, sys, os.path
from plplot import *
from numpy import *
from scipy.io import read_array

def usage():
    print "-aewt"

def semia(r,v2,mu):
    semia = 1./(2./r - v2/mu)
    return semia

def ecc(r,v,mu):
    ecc_arr = sum(square(v))*r/mu - dot(r,v)*v/mu - r/sqrt(sum(square(r)))
    return ecc_arr

def readfile(filename):
    time = read_array(file(filename),lines=(0,),columns=(0,))
    mass = read_array(file(filename),lines=(2,-9),columns=(1,),separator='=',atype='f')
    rarr = read_array(file(filename),lines=(3,-9),columns=(1,2,3),separator=('=',','),atype='f')
    varr = read_array(file(filename),lines=(4,-9),columns=(1,2,3),separator=('=',','),atype='f')
    mu = mass + mass[0]
    return time,mu,rarr,varr

def keplerian(r,v,mu):
    rxv = cross(r,v)
    hs = sum(square(rxv))
    h = sqrt(hs)
    vs = sum(square(v))
    rdotv = dot(r,v)
    rdot = rdotv/sqrt(sum(square(r)))
    vxh = cross(v,rxv)
    rhat = r/sqrt(sum(square(r)))
    e = vxh/mu - rhat
    ecc = sqrt(sum(square(e)))
    i = arccos(rxv[2]/h)
    if (rxv[0] !=0 or rxv[1] !=0):
        longnode = atan2(rxv[0],-rxv[1])
    else:
        longnode = 0.
    ecosargperi = e[0]*cos(longnode) + e[1]*sin(longnode)
    if (sin(i) != 0):
        esinargperi = e[2]/sin(i)
    else:
        esinargperi = 0.
    parameter = hs/mu
    ecostrueanom = parameter/sqrt(sum(square(r))) - 1.
    esintrueanom = rdot * h / mu
#     ecc = sqrt(square(ecostrueanom) + square(esintrueanom))
    if (esintrueanom !=0 or ecostrueanom !=0):
        trueanom = arctan2(esintrueanom, ecostrueanom)
    else:
        trueanom = 0.
    cosnode = cos(longnode)
    sinnode = sin(longnode)
    rcosu = r[0]*cosnode + r[1]*sinnode
    rsinu = (r[1]*cosnode - r[0]*sinnode)/cos(i)
    if (rsinu != 0 or rcosu != 0):
        u = arctan2(rsinu,rcosu)
    else:
        u = 0.
    argperi = u - trueanom
    a = 1./(2/sqrt(sum(square(r))) - vs/mu)
    eccanom = 2. * arctan(sqrt((1.-ecc)/(1.+ecc)) * tan(trueanom/2.))
    meananom = eccanom - ecc*sin(eccanom)
    return a,ecc,i,longnode,argperi,meananom

try:
    opts, args = getopt.getopt(sys.argv[1:], "haewt", ["help"])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err)
    usage()
    sys.exit(2)

tend = 2288
step = 10
nplanets = 4
w = zeros(nplanets)
argperi = zeros(nplanets)
meananom = zeros(nplanets)
p = 1
q = 1

plinit()
o, a = opts[0]
if o == "-a":
    plenv(0,tend/10,0.5,2.,0,2)
    pllab("Time [orbits]","a [AU]","")
elif o == "-e":
    plenv(0,tend/10,0.,.5,0,2)
    pllab("Time [orbits]","e","")
elif o == "-w":
    plenv(0,tend/10,-pi,pi,0,2)
    pllab("Time [orbits]","#gD #gw","")
elif o == "-t":
    plenv(0,tend/10,0.,2*pi,0,2)
    pllab("Time [orbits]","#gT","")

for i in range(1,tend,step):
    filename = "%s%04u" % (args[0], i)
    print "Reading %s" % filename
    if (os.path.isfile(filename) == False):
        break
    time,mu,rarr,varr = readfile(filename)
    for i in range(1,nplanets):
        a,e,inc,longnode,argperi[i],meananom[i] = keplerian(rarr[i,:],varr[i,:],mu[i])
        if o == "-a":
            plpoin(time/2/pi,[a],2)
        elif o == "-e":
            plpoin(time/2/pi,[e],2)
        elif o == "-w":
            w[i] = arctan2( ecc_arr[0], ecc_arr[1] )
    if o == "-w":
        plpoin(time/2/pi,[w[2]-w[1]],2)
        plpoin(time/2/pi,[argperi[2]-argperi[1]],2)
    if o == "-t":
        plpoin(time/2/pi,[(p+q)*(meananom[2]+argperi[2])-p*(meananom[1]+argperi[1])-q*argperi[1]],2)

plend()
