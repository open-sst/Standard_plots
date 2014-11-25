import numpy as np
import math
from datetime import date, time, datetime
import calendar
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

    def pull_year(self,year):

        thisyr = []
        for i in range(0,len(self.data)):
            if self.years[i] == year:
                thisyr.append(self.data[i])
 
        return thisyr

    def pull_month(self,month):

        thismn = []
        for i in range(0,len(self.data)):
            if self.months[i] == month:
                thismn.append(self.data[i])
 
        return thismn


    def pull_year_to_date(self, year):

        year_to_date = []
        thisyr = self.pull_year(year)

        for i in range(0,len(thisyr)):
            year_to_date.append(np.mean(thisyr[0:i+1]))

        return year_to_date

    def pull_year_to_date_with_fill(self,year,option=0):

        year_to_date = []
        thisyr = self.pull_year(year)
        n = len(thisyr)

        if n < 12:
            for m in range(n+1,13):
                allmonths = self.pull_month(m)
                allmonths.sort()
                nmonths = len(allmonths)
                if option == 0:
                    thisyr.append(allmonths[nmonths-1])
                elif option == 1:
                    thisyr.append(allmonths[nmonths-10])
                elif option == 2:
                    thisyr.append(allmonths[nmonths-3])
                elif option == 3:
                    thisyr.append(np.mean(allmonths[nmonths-10:nmonths-1]))
                    
        for i in range(0,len(thisyr)):
            year_to_date.append(np.mean(thisyr[0:i+1]))

        return year_to_date

    def plot_scenario(self,title):

        ytd = self.pull_year_to_date_with_fill(2014)
        xax = range(1,13)
 
        for i in range(0,4):
            ytd = self.pull_year_to_date_with_fill(2014,i)
            plt.plot(xax,ytd,color="Black",linewidth=2)

        

        years_of_interest = [2010,2005,1998,2013,2003,2014]
        colours_of_int = ['FireBrick','SteelBlue','DarkGray','DarkGray','DarkGray','Gold']
        i=0
        for y in years_of_interest:
            ytd = self.pull_year_to_date(y)
            xax = range(1,len(ytd)+1)
            plt.plot(xax,ytd,color=colours_of_int[i],linewidth=4,label=str(y))
            i+=1

        
        plt.title(title)
        plt.axis((0,16,0.37,0.70))
        plt.xlabel('Month', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
        plt.legend(loc='best')
        
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
                       
    def annualise(self, tomonth = 12):
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
            if self.months[i] <= tomonth:
                y = int(self.years[i])
                annual_data[y-fyear] += self.data[i]
                annual_data_ct[y-fyear] += 1.0
                annual_years[y-fyear] = y

        for i in range(0,len(annual_data)):
            if annual_data_ct[i] > 0:
                annual_data[i] /= annual_data_ct[i]
            else:
                annual_data[i] =  -99

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
#        plt.savefig('test.eps')
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

def nino_year(year):

    #based on CPC Nino 3.4 apparently.
    result = 0 #neutral

    elninos = [1958, 1966, 1973, 1983, 1987, 1988, 1998, 2003, 2010]
    laninas = [1950, 1955, 1956, 1974, 1976, 1989, 1999, 2000, 2008, 2011]
 
    if year in laninas:
        result = -1
    elif year in elninos:
        result = 1

    return result
    

