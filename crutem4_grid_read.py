import numpy as np
import math
from datetime import date, time, datetime
import calendar
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from time_series import *

file = Dataset('Data/CRUTEM.4.3.0.0.anomalies.nc')
#file = Dataset('Data/CRUTEM3.nc')

allDimNames = file.dimensions.keys() 
variableNames = file.variables.keys()

print variableNames

temps = file.variables['temperature_anomaly'][:]
#temps = file.variables['temp'][:].squeeze()
lats = file.variables['latitude'][:]
lons = file.variables['longitude'][:]
time = file.variables['time'][:]
#time = file.variables['t'][:]

lons, lats = np.meshgrid(lons,lats)

weights = np.cos(lats*np.pi/180.)


years = []
months = []

for y in range(1850,2015):
    lmn  = 12
    if y == 2014:
        lmn = 9
    for m in range(1,lmn+1):
        years.append(y)
        months.append(m)


ts_nh = []
ts_sh = []
ts_europe = []
ts_globe = []
ts_samerica = []
ts_europe2 = []
ts_turkey = []
ts_africa = []
ts_africa2 = []

#ts = []

ts_eea_europe = []

nreg = 10
  #0 nh
   #1 sh
   #2 Europe
   #3 South America
   #4 EEA Europe
   #5 EEA Turkey
   #6 Africa
   #7 W Block of Africa
   #8 E Block of Africa
    #9 NE Block of Africa
regions = [ [-180,   0, 180,  90], \
            [-180, -90, 180,   0], \
            [-25,   35,  45,  75], \
            [-85,  -60,  30,  15], \
            [-25,   35,  30,  70], \
            [ 30,   35,  45,  40], \
            [-20,  -40,  50,  40], \
            [-20,  -40,  30,  40], \
            [ 30,  -40,  50,  15], \
            [ 30,   15,  40,  30] ]

for i in range(0,len(time)):
    
    sumdat = np.zeros(nreg)
    sumwei = np.zeros(nreg)

    for xx in range(0,72):
        for yy in range(0,36):
                if temps[i,yy,xx] != None:
                    for j in range(0,nreg):
                        if lats[yy,xx] > regions[j][1] and \
                           lats[yy,xx] < regions[j][3] and \
                           lons[yy,xx] > regions[j][0] and \
                           lons[yy,xx] < regions[j][2]:
                            sumdat[j] += temps[i,yy,xx]*weights[yy,xx]
                            sumwei[j] += weights[yy,xx]

    avs = []
    for j in range(0,nreg):     
        if sumwei[j] != 0.0:
            avs.append(sumdat[j]/sumwei[j])
        else:
            avs.append(-99)

    ts_nh.append(avs[0])
    ts_sh.append(avs[1])
    ts_europe.append(avs[2])
    ts_samerica.append(avs[3])
    ts_europe2.append(avs[4])
    ts_turkey.append(avs[5])
    ts_africa.append(avs[6])

    ts_africa2.append( (avs[7]*sumwei[7]+avs[8]*sumwei[8]+avs[9]*sumwei[9])/(sumwei[7]+sumwei[8]+sumwei[9]) )

    ts_eea_europe.append((ts_europe2[i]*sumwei[4] + ts_turkey[i]*sumwei[5])/(sumwei[4]+sumwei[5]))

    if ts_nh[i] != -99 and ts_sh[i] != -99:
        ts_globe.append((2*ts_nh[i]+ts_sh[i])/3.)
    else:
        ts_globe.append(-99)

#    print i,ts_nh[i],ts_sh[i],ts_globe[i]

#plt.plot(ts_europe)
#plt.show()

ts_africa_monthly = monthly_time_series(years,months,ts_africa)
ts_africa_annual = ts_africa_monthly.annualise()
ts_africa2_monthly = monthly_time_series(years,months,ts_africa2)
ts_africa2_annual = ts_africa2_monthly.annualise()

ts_africa_annual.plot_ts("Red")
ts_africa2_annual.plot_ts("Green")
plt.show()

ts_samerica_monthly = monthly_time_series(years,months,ts_samerica)
ts_samerica_annual = ts_samerica_monthly.annualise()

ts_europe_monthly = monthly_time_series(years, months, ts_europe)
ts_europe_annual = ts_europe_monthly.annualise()

ts_eea_europe_monthly = monthly_time_series(years, months, ts_eea_europe)
ts_eea_europe_annual = ts_eea_europe_monthly.annualise()

ts_eea_europe_annual.plot_ts("Blue")
ts_europe_annual.plot_ts("Red")
plt.show()

ts_nh_monthly     = monthly_time_series(years, months, ts_nh)     
ts_sh_monthly     = monthly_time_series(years, months, ts_sh)     
                
