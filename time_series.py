import numpy as np
import math
from datetime import date, time, datetime
import calendar
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from math import ceil
from scipy import linalg
from scipy.stats import mode
import random

def lowess(x, y, f=2./3., iter=3):

    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations."""

    n = len(x)
    r = int(ceil(f*n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:,None] - x[None,:]) / h), 0.0, 1.0)
    w = (1 - w**3)**3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:,i]
            b = np.array([np.sum(weights*y), np.sum(weights*y*x)])
            A = np.array([[np.sum(weights), np.sum(weights*x)],
                   [np.sum(weights*x), np.sum(weights*x*x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1]*x[i]
        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta**2)**2
    return yest

def gauss(n=11,sigma=1):
    """Create approximate Gaussian filter"""
    r = range(-int(n/2),int(n/2)+1)
    return [1 / (sigma * np.sqrt(2*np.pi)) * np.exp(-float(x)*float(x)/(2*sigma*sigma)) for x in r]

class monthly_time_series:

    """Class for dealing with monthly time series"""
    
    def __init__(self,years,months,data):
        self.years = years
        self.months = months
        self.data = data
        self.make_time_axis()

    def make_time_axis(self):
#an internal routine that turns year and month into a single time axis
        self.time_axis = self.years[:]
        for i in range(0,len(self.time_axis)):
            self.time_axis[i] += (self.months[i]-1.)/12.

    def add_month(self,year,month,data):
#add a single month to the end of the series
        self.years.append(year)
        self.months.append(month)
        self.data.append(data)
        self.make_time_axis()

    def diff(series1,series2):
#take the difference between the intersection of two series
        n1 = len(series1.data)
        n2 = len(series2.data)

        years = []
        months = []
        data = []

        for i in range(0,n1):
            for j in range(0,n2):
                if series1.years[i] == series2.years[j] and \
                   series1.months[i] == series2.months[j]:
                    years.append(series1.years[i])
                    months.append(series1.months[i])
                    data.append(series1.data[i]-series2.data[j])

        diff = monthly_time_series(years,months,data)
        return diff

    def combine_monthly_series(had,ncdc,giss):
#take an average of three data sets. This is really only going to
#mean anything if those three data sets are HadCRUT, NCDC and GISS
#global temperatures
        comb_year = []
        comb_month = []
        comb_data = []

        for i in range(0,len(ncdc.data)):
            assert ncdc.years[i] == giss.years[i]
            assert ncdc.years[i] == had.years[i+360]
            assert ncdc.months[i] == giss.months[i]
            assert ncdc.months[i] == had.months[i+360]
            comb_year.append(ncdc.years[i])
            comb_month.append(ncdc.months[i])
            comb_data.append(np.mean([ncdc.data[i],giss.data[i],had.data[i+360]]))

        combined = monthly_time_series(comb_year,comb_month,comb_data)

        return combined



    def plot_ts(self):
#do simple time series plot of monthly series
        plt.plot(self.time_axis,self.data)
        plt.plot(self.time_axis,np.zeros(len(self.time_axis)),color="Black")
        #plt.plot([2014,2014],[-1,1],color="Red")
        plt.show()

    def pull_year(self,year):
#extract data for all the months for one specified year
        thisyr = []
        for i in range(0,len(self.data)):
            if self.years[i] == year:
                thisyr.append(self.data[i])
 
        return thisyr

    def pull_month(self,month):
#extract data for all years for one specified month
        thismn = []
        for i in range(0,len(self.data)):
            if self.months[i] == month:
                thismn.append(self.data[i])
 
        return thismn


    def pull_data(self,year,month):
#get data value for a given year and month
        val = -99.
        for i in range(0,len(self.data)):
            if self.years[i] == year and self.months[i]:
                val = self.data[i]
        return val

    def print_month_rank_table(self,nin):
#given an el nino times series, nin, print the
#warmest 10 years for each month
        for m in range(1,13):
            allm = self.pull_month(m)
            allm.sort()
            n = len(allm)
            for rank in range(0,11):
                ind = self.data.index(allm[n-rank])
                nino = nin.pull_data(self.years[ind],self.months[ind])
                print ("%2i %2i %4i %7.3f %7.3f" % (m,rank,self.years[ind],self.data[ind],nino))
            print ("")

    def pull_year_to_date(self, year):
#pull a year of data and then run a cumulative average through it
#returned list first elements is Jan avg, 2nd is Jan-Feb average,
#all the way up to Jan-Dec average (if year goes to December)
        year_to_date = []
        thisyr = self.pull_year(year)

        for i in range(0,len(thisyr)):
            year_to_date.append(np.mean(thisyr[0:i+1]))

        return year_to_date

    def pull_year_to_date_with_fill(self,year,option=0):
#pull a year of data, fill the remainder of the year with "scenarios"
#and then run a cumulative average through it
#returned list first elements is Jan avg, 2nd is Jan-Feb average,
#all the way up to Jan-Dec average (if year goes to December)
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
                    thisyr.append(allmonths[nmonths-3])
                elif option == 2:
                    thisyr.append(allmonths[nmonths-10])
                elif option == 3:
                    thisyr.append(np.mean(allmonths[nmonths-10:nmonths-1]))
                elif option == 4:
                    thisyr.append(np.mean(thisyr[0:n]))
                    
        for i in range(0,len(thisyr)):
            year_to_date.append(np.mean(thisyr[0:i+1]))

        return year_to_date


    def plot_single_years(self,title):
#extract the warmest, 3rd warmest, 5th, 10th and 20th warmest months
#and plot them in grey for context
#then plot certain selected years in colour
        sample_yr_1 = []
        sample_yr_2 = []
        sample_yr_3 = []
        sample_yr_4 = []
        sample_yr_5 = []

        for m in range(1,13):
            allmonths = self.pull_month(m)
            allmonths.sort()
            n = len(allmonths)
            sample_yr_1.append(allmonths[n-1])
            sample_yr_2.append(allmonths[n-3])
            sample_yr_3.append(allmonths[n-5])
            sample_yr_4.append(allmonths[n-10])
            sample_yr_5.append(allmonths[n-20])

        xax = range(1,13)
        plt.plot(xax,sample_yr_1,color="DarkGray",linewidth=6)
        plt.plot(xax,sample_yr_2,color="DarkGray",linewidth=6)
        plt.plot(xax,sample_yr_3,color="DarkGray",linewidth=6)
        plt.plot(xax,sample_yr_4,color="DarkGray",linewidth=6)
        plt.plot(xax,sample_yr_5,color="DarkGray",linewidth=6)
       
        years_of_interest = [2010,2005,1998,2014]
        colours_of_int = ['FireBrick','SteelBlue','Pink','Gold']
        i = 0
        for y in years_of_interest:
            thisyr = self.pull_year(y)
            xax = range(1,len(thisyr)+1)
            plt.plot(xax,thisyr,color=colours_of_int[i],linewidth=4)
            i += 1

        plt.axis((0,13,0.21,0.75))
        plt.show()

    def plot_scenario(self,title):
#plot scenarios for 2014
        ytd = self.pull_year_to_date_with_fill(2014)
        xax = range(1,13)

        #loop over 5 different "scenarios" for remainder of year
        for i in range(0,5):
            ytd = self.pull_year_to_date_with_fill(2014,i)
            print title,i,ytd[11]
            plt.plot(xax,ytd,color="Black",linewidth=2)

        years_of_interest = [2010,2005,1998,2013,2003,2002,2007,2006,2009,2014]
        colours_of_int = ['FireBrick','SteelBlue','DarkGray','DarkGray','DarkGray','DarkGray','DarkGray','DarkGray','DarkGray','Gold']
        i=0
        for y in years_of_interest:
            ytd = self.pull_year_to_date(y)
            xax = range(1,len(ytd)+1)
            plt.plot(xax,ytd,color=colours_of_int[i],linewidth=4,label=str(y))
            i+=1
        
        plt.title(title)
        if title == "Cowtan and Way Hybrid" or title=="ERA":
            #these data sets have different baselines
            plt.axis((0,16,0.09,0.45))
        else:
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
        annual_lounc = []
        annual_hiunc = []
        annual_data_ct = []

        for i in range(fyear,lyear+1):
            annual_years.append(0)
            annual_data.append(0)
            annual_lounc.append(0)
            annual_hiunc.append(0)
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
                annual_lounc[i] /= annual_data_ct[i]
                annual_hiunc[i] /= annual_data_ct[i]
            else:
                annual_data[i] =  -99
                annual_lounc[i] =  -99
                annual_hiunc[i] =  -99

        return time_series(annual_years,annual_data,annual_lounc,annual_hiunc)




######################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#********************************************************************#
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
######################################################################



class time_series:

    """class to do anual time series with uncertainty ranges"""

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

    def add_offset(self,offset):
        #little routine for shifting entire series up or down by a constant offset
        for i in range(0,len(self.times)):
            self.data[i] += offset
            self.lounc[i] += offset
            self.hiunc[i] += offset


    def get_rank_of_year(self,year):
        #this tells you the nominal rank of a particular year 1 is warmest
        years = self.times[:]
        anoms = self.data[:]
        sorted_years = [years for (anoms,years) in sorted(zip(anoms,years))]
        rank = len(self.times)-sorted_years.index(year)
        return rank
    
    def draw_sample(self):
        #assuming that uncertainty ranges is 95% confidence, that the errors
        #have a gaussian distribution and that they are uncorrelated: create
        #a single "realisation" of the data set
        times = self.times[:]
        data = self.data[:]
        lounc = self.lounc[:]
        hiunc = self.hiunc[:]
        for i in range(0,len(times)):
            draw = random.gauss(0,(hiunc[i]-lounc[i])/3.92)
            data[i] += draw
            lounc[i] += draw
            hiunc[i] += draw
        return time_series(times,data,lounc,hiunc)

    def rebaseline(self,year1,year2):
        #change baseline for anomalies (should probably use the offset command here...!)
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
        print "Year,Anomaly,LowerUncertaintyBound,HigherUncertaintyBound"
        for i in range(0,len(self.data)):
            print("%4d,%7.3f,%7.3f,%7.3f" % (self.times[i],self.data[i],self.lounc[i],self.hiunc[i]))

    def print_period_avg(self,year1,year2):
        #print average of time series from year1 to year2
        index1 = self.times.index(year1)
        index2 = self.times.index(year2)

        print("")
        print "Long Term Averages"
        print ("%4i-%4i avg = %7.3f" % (year1,year2,np.mean(self.data[index1:index2+1])))


    def plot_binomial_smooth(self,f=0.3):
        #plot time series, binomially smoothed and loess smoothed series
        #f is the loess smoothing parameter; higher is smoother
        filt = gauss(21,2.264)

        n = len(self.data)

        smo = []
        smaxis = []

        for i in range(4,n-4):
            sump = 0.0
            filtsum = 0.0
            for j in range(0,21):
                if i-10+j < n and i-10+j > 0:
                    sump += filt[j]*self.data[i-10+j]
                    filtsum += filt[j]
            smaxis.append(self.times[i])
            smo.append(sump/filtsum)

        yest = lowess(np.asarray(self.times,float), np.asarray(self.data), f=f, iter=10)

        plt.plot(self.times,self.data)
        plt.plot(self.times,yest,linewidth=4)
        plt.plot(smaxis,smo,linewidth=4)
        plt.show()

    def print_running_mean(self,filter_width):
        #print running mean of series
        print self.name+" Running mean"
        for i in range(filter_width,len(self.data)):
            print ("%3d %4d %7.3f" % (i,self.times[i],np.mean(self.data[i-filter_width+1:i+1])))
            if i == len(self.data)-1:
                print self.data[i-filter_width+1:i+1]
                print self.times[i-filter_width+1:i+1]

    def plot_running_mean(self,filter_width):
        running_mean = []
        running_axis = []
        for i in range(filter_width,len(self.data)):
            running_mean.append(np.mean(self.data[i-filter_width+1:i+1]))
            running_axis.append(np.mean(self.times[i-filter_width+1:i+1]))

        plt.plot(running_axis,running_mean, linewidth=4)
#        plt.xlabel('Year', fontsize=18)
#        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
#        plt.axis((1880+filter_width-1.5,2015.5,-0.62,0.62))
#        plt.show()
    
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
            decadal.append(np.mean(self.data[i-filter_width+1:i+1]))
            decadal_years.append(self.times[i])

        plt.plot(self.times,self.data)
        for i in range(0,len(decadal)):
            plt.plot([decadal_years[i]-filter_width+1,decadal_years[i]],[decadal[i],decadal[i]], color="Red")
        plt.show()

    def modded_skyscraper_diagram(self):

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
                               [0,0,d,d]),facecolor=color,edgecolor="Black",zorder=2)
            plt.gca().add_patch(poly)

        poly = Polygon(zip([y-delta,y+delta,y+delta,y-delta],\
                           [0,0,d,d]),facecolor="White",edgecolor="Black",zorder=2)
        plt.gca().add_patch(poly)


        plt.xlabel('Year', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)

        yyy = 2017.5
        barthick = 0.5
        l = yyy-4*barthick
        r = yyy-3*barthick
        d = 0.582718018519
        poly = Polygon(zip([l,r,r,l],\
                               [0,0,d,d]),facecolor='Yellow',edgecolor="Black",zorder=2,label="Warmest on record for rest of year")
        plt.gca().add_patch(poly)
        plt.plot([1997,2020],[d,d],color="Black")
        
        l = yyy-3*barthick
        r = yyy-2*barthick
        d = 0.566804962963
        poly = Polygon(zip([l,r,r,l],\
                               [0,0,d,d]),facecolor='Orange',edgecolor="Black",zorder=2,label="3rd warmest on record")
        plt.gca().add_patch(poly)

        l = yyy-2*barthick
        r = yyy-1*barthick
        d = 0.560179222222
        poly = Polygon(zip([l,r,r,l],\
                               [0,0,d,d]),facecolor='Maroon',edgecolor="Black",zorder=2,label="Average of 10 warmest on record")
        plt.gca().add_patch(poly)

        
        l = yyy-1*barthick
        r = yyy-0*barthick
        d = 0.547015935185
        poly = Polygon(zip([l,r,r,l],\
                               [0,0,d,d]),facecolor='SaddleBrown',edgecolor="Black",zorder=2,label="10th warmest on record")
        plt.gca().add_patch(poly)
        plt.plot([1997,2020],[d,d],color="Black")


        plt.legend(loc='lower right')

        plt.plot([1997.5,2018.5],[0,0],color="Black")
        plt.axis((1997.5,2018.5,-0.2,0.6))
        plt.show()

    def plot_ts_highlight_nino(self,ninos):

        delta = 0.5

        plt.plot(self.times,self.data,color="Black")

        for i in range(0,len(self.times)):

            y = self.times[i]
            d = self.data[i]

            color = "Beige"

            if y > 1950:
                if ninos.data[ninos.times.index(y)] == 0:
                    color = "Beige"
                elif ninos.data[ninos.times.index(y)] == -1:
                    color = "DodgerBlue"
                    poly = Polygon(zip([y-delta,y+delta,y+delta,y-delta],\
                                       [-0.6,-0.6,0.6,0.6]),facecolor=color,edgecolor="None",zorder=2,alpha=0.5)
                    plt.gca().add_patch(poly)
                elif ninos.data[ninos.times.index(y)] == 1:
                    color = "Red"
                    poly = Polygon(zip([y-delta,y+delta,y+delta,y-delta],\
                                       [-0.6,-0.6,0.6,0.6]),facecolor=color,edgecolor="None",zorder=2,alpha=0.5)
                    plt.gca().add_patch(poly)

        plt.xlabel('Year', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
        plt.axis((1949.5,2015.5,-0.32,0.62))
        plt.show()

    def plot_skyscraper_diagram(self,ninos):

        plt.plot(np.zeros(50), color="White")
        for i in range(1960,2020,10):
            plt.plot([i,i],[-10,10],color="DarkGray",zorder=1)

        for i in range(0,len(self.times)):

            y = self.times[i]
            d = self.data[i]

            color = "Beige"

            if y > 1950:
                if ninos.data[ninos.times.index(y)] == 0:
                    color = "Beige"
                elif ninos.data[ninos.times.index(y)] == -1:
                    color = "DodgerBlue"
                elif ninos.data[ninos.times.index(y)] == 1:
                    color = "FireBrick"

#            if nino_year(y) == 0:
#                color = "Silver"
#            elif nino_year(y) == -1:
#                color = "DodgerBlue"
#            elif nino_year(y) == 1:
#                color = "FireBrick"

            delta = 0.45
           
            poly = Polygon(zip([y-delta,y+delta,y+delta,y-delta],\
                               [0,0,d,d]),facecolor=color,edgecolor="Black",zorder=2)
            plt.gca().add_patch(poly)

        plt.xlabel('Year', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)

#draw off-screen polygons to get an appropriate legend
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='FireBrick',edgecolor="Black",label="El Nino")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='Silver',edgecolor="Black",label="Neutral")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor='DodgerBlue',edgecolor="Black",label="La Nina")
        plt.gca().add_patch(poly)
        plt.legend(loc='best')

        plt.plot([1949.5,2015.5],[0,0],color="Black")
        plt.axis((1949.5,2015.5,-0.32,0.62))
#        plt.savefig('skyscraper.eps')
        plt.show()


    def plot_ranking_diagram(self):
#this is the original and beautiful rainbow ranking diagram
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

        highesthigh = -999.999
        lowestlow = 999.999

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

            if lo < lowestlow:
                lowestlow = lo
            if hi > highesthigh:
                highesthigh = hi
            
            poly = Polygon(zip([i+1-0.5,i+1+0.5,i+1+0.5,i+1-0.5],\
                               [lo,lo,hi,hi]),facecolor=col,edgecolor="White")
            plt.gca().add_patch(poly)

#plot medians as black dash
            plt.plot([i+1-0.4,i+1+0.4],[mi,mi],color="Black")

#add year lables to each coloured bar
            plt.text(i+0.7,hi+0.105,str(int(this_year)),fontsize=10,rotation=90)

        plt.axis((0,51,lowestlow-0.2,highesthigh+0.2))
        plt.show()



    def plot_alt_ranking_diagram(self,iters=1000):

#this plots time on x-axis and rank on the y-axis. The polygons show the range
#of ranks for each year given uncertainties
#Assumes that errors are uncorrelated from one year to the next
#iters is the number of samples to draw to determine the percentiles of the
#ranks to display (2.5%, 25%, 50%, 75% and 97.5%)
        
        y1 = int(min(self.times))
        y2 = int(max(self.times))

        nranks = len(self.times)

        
        
        for yr in range(y1,y2+1):
            ranks = []
            for i in range(0,iters):
                ranks.append(self.draw_sample().get_rank_of_year(yr))

            lolo = -1*np.percentile(ranks,2.5) + 0.5
            lo   = -1*np.percentile(ranks,25) + 0.5
            hi   = -1*np.percentile(ranks,75) - 0.5
            hihi = -1*np.percentile(ranks,97.5) - 0.5
            md   = -1*np.percentile(ranks,50)

            barthick = 0.46
    
            poly = Polygon(zip([yr-barthick,yr+barthick,yr+barthick,yr-barthick],\
                               [lolo,lolo,hihi,hihi]), facecolor="PowderBlue",\
                           edgecolor="PowderBlue")
            plt.gca(zorder=2).add_patch(poly)
            poly = Polygon(zip([yr-barthick,yr+barthick,yr+barthick,yr-barthick],\
                               [lo,lo,hi,hi]), facecolor="DodgerBlue",\
                           edgecolor="DodgerBlue")
            plt.gca(zorder=2).add_patch(poly)
            plt.plot([yr-0.5,yr+0.5],[md,md],color="Black")

            print yr,int(round(mode(ranks)[0][0]))

        plt.xlabel('Year', fontsize=18)
        plt.ylabel('-1*Rank', fontsize=18)
        plt.axis((y1-1.5, y2+1.5, -1-nranks,4))
        plt.show()


    def plot_ranking_diagram_mono(self,scheme=1,legend=1):
#original style ranking diagram with a variety of user switchable
#colour schemes
        order = self.data[:]
        order.sort()

        if scheme == 1:
            colscheme = ['#91003f','#ce1256','#e7298a','#df65b0','#c994c7','#d4b9da','#e7e1ef','#f7f4f9']
        if scheme == 2:
            colscheme = ['#8888FF','#8888FF','#8888FF','#8888FF','#8888FF','#8888FF','#8888FF','#8888FF']
        if scheme == 3:
            colscheme = ['DarkRed','Red','DarkOrange','Gold','Green','DodgerBlue','Navy','Purple']
        if scheme == 4:
            colscheme = ['#084594','#2171b5','#4292c6','#6baed6','#9ecae1','#c6dbef','#deebf7','#f7fbff']
        if scheme == 5:
            colscheme = ['#2222FF','#8888FF','#8888FF','#8888FF','#8888FF','#8888FF','#8888FF','#8888FF']
        if scheme == 6:
            colscheme = ['#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4']
        if scheme == 7:
            colscheme = ['#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e']
        plt.plot([0,0],[0,0], color="White")

        plt.xlabel('Rank', fontsize=18)
        plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)

#draw some off-screen polygons to get an appropriate legend
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[0],edgecolor="Black",label="2011-2014")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[1],edgecolor="Black",label="2001-2010")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[2],edgecolor="Black",label="1991-2000")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[3],edgecolor="Black",label="1971-1990")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[4],edgecolor="Black",label="1951-1970")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[5],edgecolor="Black",label="1931-1950")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[6],edgecolor="Black",label="1911-1930")
        plt.gca().add_patch(poly)
        poly = Polygon(zip([0,0,0,0],[0,0,0,0]),facecolor=colscheme[7],edgecolor="Black",label="1850-1910")
        plt.gca().add_patch(poly)

        if legend == 1:
            plt.legend(loc='best')

        highesthigh = -999.999
        lowestlow = 999.999

        for i in range(0,50):

            this_year = self.times[self.data.index(order[len(order)-i-1])]

            if this_year >= 1850 and this_year <=1910:
                col = colscheme[7]
            elif this_year >= 1911 and this_year <=1930:
                col = colscheme[6]
            elif this_year >= 1931 and this_year <=1950:
                col = colscheme[5]
            elif this_year >= 1951 and this_year <=1970:
                col = colscheme[4]
            elif this_year >= 1971 and this_year <=1990:
                col = colscheme[3]
            elif this_year >= 1991 and this_year <=2000:
                col = colscheme[2]
            elif this_year >= 2001 and this_year <=2010:
                col = colscheme[1]
            elif this_year >= 2011 and this_year <=2014:
                col = colscheme[0]
            

#plot uncertainty range as coloured bar
            lo = self.lounc[self.data.index(order[len(order)-i-1])]
            hi = self.hiunc[self.data.index(order[len(order)-i-1])]
            mi = self.data[self.data.index(order[len(order)-i-1])]

            if lo < lowestlow:
                lowestlow = lo
            if hi > highesthigh:
                highesthigh = hi
            
            poly = Polygon(zip([i+1-0.5,i+1+0.5,i+1+0.5,i+1-0.5],\
                               [lo,lo,hi,hi]),facecolor=col,edgecolor="White")
            plt.gca().add_patch(poly)

#plot medians as black dash
            plt.plot([i+1-0.4,i+1+0.4],[mi,mi],color="Black",linewidth=2)

#add year lables to each coloured bar
            plt.text(i+0.7,hi+0.105,str(int(this_year)),fontsize=10,rotation=90)

        plt.axis((0,51,lowestlow-0.05,highesthigh+0.2))
        plt.show()


def nino_year(year):
    #based on CPC Nino 3.4. Apparently.

    result = 0 #neutral

    elninos = [1958, 1966, 1973, 1983, 1987, 1988, 1998, 2003, 2010]
    laninas = [1950, 1955, 1956, 1974, 1976, 1989, 1999, 2000, 2008, 2011]

    if year in laninas:
        result = -1
    elif year in elninos:
        result = 1

    return result
    

