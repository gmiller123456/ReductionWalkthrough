import numpy as np
from datetime import datetime
from skyfield.api import Loader, utc, Topos

LOAD_FROM = '/home/mkbrewer/.skyfield'
APR_6_2017 = 1491436800.0

load = Loader(LOAD_FROM)
ts = load.timescale()
planets = load('de430_1850-2150.bsp')
jupiter = planets['jupiter barycenter']
moon = planets['moon']
earth = planets['earth']

# Polar Motion (arcseconds)
xpole = 0.010021 
ypole = 0.385588
# No Polar Motion
#xpole = 0.0
#ypole = 0.0

site = Topos('22.959748 S', '67.787260 W', elevation_m=5186.0, x=xpole, y=ypole)
place = earth + site


start = APR_6_2017
nsamples = 1441
delta = 60.0
times = []
for i in range(nsamples):
    ctime = start + i * delta
    times.append(datetime.utcfromtimestamp(ctime).replace(tzinfo=utc))

ut = ts.utc(times)
astrometric = place.at(ut).observe(moon)
ra_icrs, dec_icrs, distance = astrometric.radec()
dist = distance.m
apparent = astrometric.apparent()
ra, dec, distance = apparent.radec(epoch='date')
el, az, distance = apparent.altaz()
lst = ut.gast + site.longitude._hours
mask = (lst < 0.0)
lst[mask] += 24.0

print("UTC,LST,ICRS_RA,ICRS_DEC,RA,DEC,AZ,EL,DIST,DUT1,XP,YP")
for i in range(nsamples):
    times[i] = times[i].strftime("%Y-%m-%d %H:%M:%S.%f")
    print ('%s,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%.9lf,%.3lf,%.6lf,%.6lf,%.6lf' % (
           times[i], lst[i], 
           ra_icrs._degrees[i], dec_icrs._degrees[i], 
           ra._degrees[i], dec._degrees[i], 
           az._degrees[i], el._degrees[i], dist[i], 
           ut.dut1[i], xpole, ypole )) 



