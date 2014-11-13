import numpy as np
import math
from datetime import date, time, datetime
import calendar
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

class monthly_time_series:

    def __init__(self,years,months,data):
        self.years = years
        self.months = months
        self.data = data
        self.make_time_axis()

    def make_time_axis(self):
        self.time_axis = self.years[:]
        for i in range(0,len(self.time_axis)):
            self.time_axis[i] += (self.months[i]-1.)/12.

    def plot_ts(self):
        plt.plot(self.time_axis,self.data)
        plt.plot(self.time_axis,np.zeros(len(self.time_axis)),color="Black")
        plt.plot([2014,2014],[-1,1],color="Red")
        plt.show()

    def rebaseline(self,year1,year2):
        #choose new climatology period
        clim = np.zeros(12)
        climcounts = np.zeros(12)
        for i in range(0,len(self.years)):
            if self.years[i] >= year1 and self.years[i] <= year2:
                clim[self.months[i]-1] += self.data[i]
                climcounts[self.months[i]-1] += 1.0
        for i in range(0,12):
            clim[i] /= climcounts[i]
        for i in range(0,len(self.years)):
            self.data[i] -= clim[self.months[i]-1]
                       
    def annualise(self):
        #go from monthly averages to annual averages
        fyear = int(min(self.years))
        lyear = int(max(self.years))

        annual_years = []
        annual_data = []
        annual_data_ct = []

        for i in range(fyear,lyear+1):
            annual_years.append(0)
            annual_data.append(0)
            annual_data_ct.append(0)

        for i in range(0,len(self.data)):
            y = int(self.years[i])
            annual_data[y-fyear] += self.data[i]
            annual_data_ct[y-fyear] += 1.0
            annual_years[y-fyear] = y

        for i in range(0,len(annual_data)):
            annual_data[i] /= annual_data_ct[i]

        return time_series(annual_years,annual_data,annual_data,annual_data)

class time_series:

    def __init__(self,times,data,lounc,hiunc):
        self.times = times
        self.data = data
        self.lounc = lounc
        self.hiunc = hiunc
        self.name = ""

    def add_year(self,year,data,lounc,hiunc):
        self.times.append(year)
        self.data.append(data)
        self.lounc.append(lounc)
        self.hiunc.append(hiunc)

    def add_name(self,name):
        self.name = name

    def rebaseline(self,year1,year2):
        #change baseline for anomalies
        ind1 = self.times.index(year1)
        ind2 = self.times.index(year2)
        clim = np.mean(self.data[ind1:ind2+1])
        for i in range(0,len(self.data)):
            self.data[i] -= clim
            self.lounc[i] -= clim
            self.hiunc[i] -= clim

    def print_ordered_ts(self,topten):
        #print the warmest n years where n=topten
        print self.name+" Top "+str(topten)
        order = self.data[:]
        order.sort()
        for i in range(len(order)-topten,len(order)):
            print("%3d %4d %7.3f %7.3f %7.3f " % ( len(order)-i,self.times[self.data.index(order[i])], \
                               self.data[self.data.index(order[i])], \
                               self.lounc[self.data.index(order[i])], \
                               self.hiunc[self.data.index(order[i])]))

    def print_ts(self):
        print self.name+" Annual averages"
        for i in range(0,len(self.data)):
            print("%3d %4d %7.3f" % (i,self.times[i],self.data[i]))


    def print_running_mean(self,filter_width):
        print self.name+" Running mean"
        for i in range(filter_width-1,len(self.data)):
            print ("%3d %4d %7.3f" % (i,self.times[i],np.mean(self.data[i-filter_width+1:i])))

    def plot_ts(self, color):
        plt.plot(self.times, self.data, linewidth=2.0, color=color, label=self.name)

    def plot_ts_with_unc(self, colora, colorbk):
        plt.plot(self.times, self.data, linewidth=2.0, color=colora, label=self.name)
        plt.fill_between(self.times, self.lounc, self.hiunc,
                facecolor=colorbk,color=colorbk, alpha=0.5,
                label='1 sigma range')
        plt.plot(self.times, self.hiunc, linewidth=1.0, color=colora, alpha=0.5)
        plt.plot(self.times, self.lounc, linewidth=1.0, color=colora, alpha=0.5)

        mx = max(self.hiunc)
        mn = min(self.lounc)
        delta = 0.1 * (mx-mn)

        plt.axis((1848,2016,mn-delta,mx+delta))

    def plot_decadal(self,filter_width):

        nwidths = int(float(len(self.data))/float(filter_width))
        decadal = []
        decadal_years = []

        for i in range(len(self.data)-(nwidths-1)*filter_width-1, \
                       len(self.data),\
                       filter_width):
            decadal.append(np.mean(self.data[i-filter_width+1:i]))
            decadal_years.append(self.times[i])

        plt.plot(self.times,self.data)
        for i in range(0,len(decadal)):
            plt.plot([decadal_years[i]-filter_width+1,decadal_years[i]],[decadal[i],decadal[i]], color="Red")
        plt.show()


    def plot_skyscraper_diagram(self):

        plt.plot(np.zeros(50), color="White")

        for i in range(0,len(self.times)):

            y = self.times[i]
            d = self.data[i]

            color = "Silver"

            if nino_year(y) == 0:
                color = "Silver"
            elif nino_year(y) == -1:
                color = "DodgerBlue"
            elif nino_year(y) == 1:
                color = "FireBrick"

            delta = 0.45
           
            poly = Polygon(zip([y-delta,y+delta,y+delta,y-delta],\
                               [0,0,d,d]),facecolor=color,edgecolor="Black")
            plt.gca().add_patch(poly)

        plt.xlabel('Year', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)

#spoof some data points to get an appropriate legend
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='FireBrick',edgecolor="Black",label="El Nino")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Silver',edgecolor="Black",label="Neutral")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='DodgerBlue',edgecolor="Black",label="La Nina")
        plt.gca().add_patch(poly)
        plt.legend(loc='best')

        plt.plot([1949.5,2015.5],[0,0],color="Black")
        plt.axis((1949.5,2015.5,-0.32,0.62))
        plt.show()


    def plot_ranking_diagram(self):
        order = self.data[:]
        order.sort()

        plt.plot([0,0],[0,0], color="White")

        plt.xlabel('Rank', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)

#spoof some data points to get an appropriate legend
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='DarkRed',edgecolor="Black",label="2011-2014")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Red',edgecolor="Black",label="2001-2010")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='DarkOrange',edgecolor="Black",label="1991-2000")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Gold',edgecolor="Black",label="1971-1990")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Green',edgecolor="Black",label="1951-1970")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='DodgerBlue',edgecolor="Black",label="1931-1950")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Navy',edgecolor="Black",label="1911-1930")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Purple',edgecolor="Black",label="1850-1910")
        plt.gca().add_patch(poly)
        plt.legend(loc='best')

        for i in range(0,50):

            this_year = self.times[self.data.index(order[len(order)-i-1])]

            if this_year >= 1850 and this_year <=1910:
                col = "Purple"
            elif this_year >= 1911 and this_year <=1930:
                col = "Navy"
            elif this_year >= 1931 and this_year <=1950:
                col = "DodgerBlue"
            elif this_year >= 1951 and this_year <=1970:
                col = "Green"
            elif this_year >= 1971 and this_year <=1990:
                col = "Gold"
            elif this_year >= 1991 and this_year <=2000:
                col = "DarkOrange"
            elif this_year >= 2001 and this_year <=2010:
                col = "Red"
            elif this_year >= 2011 and this_year <=2014:
                col = "DarkRed"
            

#plot uncertainty range as coloured bar
            lo = self.lounc[self.data.index(order[len(order)-i-1])]
            hi = self.hiunc[self.data.index(order[len(order)-i-1])]
            mi = self.data[self.data.index(order[len(order)-i-1])]
            poly = Polygon(zip([i+1-0.5,i+1+0.5,i+1+0.5,i+1-0.5],\
                               [lo,lo,hi,hi]),facecolor=col,edgecolor="White")
            plt.gca().add_patch(poly)

#plot medians as black dash
            plt.plot([i+1-0.4,i+1+0.4],[mi,mi],color="Black")

#add year lables to each coloured bar
            plt.text(i+0.7,hi+0.05,str(int(this_year)),fontsize=10,rotation=90)

        plt.axis((0,51,-0.22,0.79))
        plt.show()



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
                
