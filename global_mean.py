import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

class time_series:

    def __init__(self,times,data,lounc,hiunc):
        self.times = times
        self.data = data
        self.lounc = lounc
        self.hiunc = hiunc

    def rebaseline(self,year1,year2):
        ind1 = self.times.index(year1)
        ind2 = self.times.index(year2)
        clim = np.mean(self.data[ind1:ind2])
        for i in range(0,len(self.data)):
            self.data[i] -= clim
            self.lounc[i] -= clim
            self.hiunc[i] -= clim

    def print_ordered_ts(self,topten):
        order = self.data[:]
        order.sort()
        for i in range(len(order)-topten,len(order)):
            print len(order)-i,self.times[self.data.index(order[i])], \
                               self.data[self.data.index(order[i])], \
                               self.lounc[self.data.index(order[i])], \
                               self.hiunc[self.data.index(order[i])]

    def print_ts(self):
        for i in range(0,len(self.data)):
            print i,self.times[i],self.data[i]


    def print_running_mean(self,filter_width):
        for i in range(filter_width-1,len(self.data)):
            print i,self.times[i],np.mean(self.data[i-filter_width+1:i])

    def plot_ts(self, color):
        plt.plot(self.times, self.data, linewidth=2.0, color=color)

    def plot_ts_with_unc(self, colora, colorbk):
        plt.plot(self.times, self.data, linewidth=2.0, color=colora)
        plt.fill_between(self.times, self.lounc, self.hiunc,
                facecolor=colorbk,color=colorbk, alpha=0.5,
                label='1 sigma range')
        plt.plot(self.times, self.hiunc, linewidth=1.0, color=colora, alpha=0.5)
        plt.plot(self.times, self.lounc, linewidth=1.0, color=colora, alpha=0.5)
        plt.axis((1848,2016,-0.79,0.79))

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


    def plot_ranking_diagram(self):
        order = self.data[:]
        order.sort()

        plt.plot(np.zeros(50), color="White")

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
            

            lo = self.lounc[self.data.index(order[len(order)-i-1])]
            hi = self.hiunc[self.data.index(order[len(order)-i-1])]
            mi = self.data[self.data.index(order[len(order)-i-1])]
            
            poly = Polygon(zip([i+1-0.5,i+1+0.5,i+1+0.5,i+1-0.5],\
                               [lo,lo,hi,hi]),facecolor=col,edgecolor="White")
            plt.gca().add_patch(poly)

            plt.plot([i+1-0.4,i+1+0.4],[mi,mi],color="Black")

        plt.axis((0,51,-0.22,0.79))
        plt.show()




def combine_series(series1, series2, series3):
        
    comb_year = []
    comb_anom = []
    comb_upper_unc = []
    comb_lower_unc = []

    for year in range(1880,2015):
            
        ind1 = series1.times.index(year)
        ind2 = series2.times.index(year)
        ind3 = series3.times.index(year)

        av = np.mean([series1.data[ind1], series2.data[ind2], series3.data[ind3] ])

        comb_year.append(year)
        comb_anom.append(av)
        comb_upper_unc.append(av+series1.hiunc[ind1]-series1.data[ind1])
        comb_lower_unc.append(av+series1.lounc[ind1]-series1.data[ind1])

        combined = time_series(comb_year,
                               comb_anom,
                               comb_lower_unc,
                               comb_upper_unc)


    return combined



def read_hadcrut4():
    f = open('Data/HadCRUT.4.3.0.0.annual_ns_avg.txt', 'r')

    hadcrut_year = []
    hadcrut_anom = []
    hadcrut_upper_unc = []
    hadcrut_lower_unc = []

    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        hadcrut_year.append(float(columns[0]))
        hadcrut_anom.append(float(columns[1]))
        hadcrut_upper_unc.append(float(columns[10]))
        hadcrut_lower_unc.append(float(columns[11]))
    
    f.close()
    had_ts = time_series(hadcrut_year,
                         hadcrut_anom,
                         hadcrut_lower_unc,
                         hadcrut_upper_unc)
    return had_ts


def read_ncdc():
    f = open('Data/aravg.ann.land_ocean.90S.90N.v3.5.4.201409.asc', 'r')

    ncdc_year = []
    ncdc_anom = []
    ncdc_lounc = []
    ncdc_hiunc = []

    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        ncdc_year.append(float(columns[0]))
        ncdc_anom.append(float(columns[1]))
        ncdc_lounc.append(float(columns[1]) - 2 * np.sqrt(float(columns[2])) )
        ncdc_hiunc.append(float(columns[1]) + 2 * np.sqrt(float(columns[2])) )

    f.close()

    ncdc_ts = time_series(ncdc_year,
                          ncdc_anom,
                          ncdc_lounc,
                          ncdc_hiunc)

    return ncdc_ts

def read_giss_block(f, block_length, giss_year, giss_anom):
    line = f.readline()
    line = f.readline()
    for j in range(1,block_length+1):
        line = f.readline()
        columns = line.split()
        giss_year.append(float(columns[0]))
        giss_anom.append(float(columns[13])/100.)
    return (giss_year,giss_anom)

def read_giss():
    f = open('Data/GLB.Ts+dSST.txt','r')
    
    giss_year = []
    giss_anom = []

#read header information and discard
    for i in range(1,7):
        f.readline()

#read first block of 21 years
    giss_year, giss_anom = read_giss_block(f, 21, giss_year, giss_anom)
#read five blocks of 20 years   
    for i in range(1,6):
        giss_year, giss_anom = read_giss_block(f, 20, giss_year, giss_anom)
#final block has less than 20
    giss_year, giss_anom = read_giss_block(f, 13, giss_year, giss_anom)

#final year is incomplete so calculate from monthlies
    g = f.readline()
    columns = g.split()
    giss_year.append(float(columns[0]))
    final_year = []
    for i in range(1,13):
        if columns[i] != '****':
            final_year.append(float(columns[i])/100.)
    giss_anom.append(np.mean(final_year))

    giss_lounc = []
    giss_hiunc = []

    for i in range(0,len(giss_year)):
        if giss_year[i] >= 1880 and giss_year[i] <= 1900:
            giss_lounc.append(giss_anom[i]-0.08)
            giss_hiunc.append(giss_anom[i]+0.08)
        elif giss_year[i] >= 1901 and giss_year[i] <= 1950:
            giss_lounc.append(giss_anom[i]-0.05)
            giss_hiunc.append(giss_anom[i]+0.05)
        elif giss_year[i] >= 1951 and giss_year[i] <= 3000:
            giss_lounc.append(giss_anom[i]-0.05)
            giss_hiunc.append(giss_anom[i]+0.05)
    
    giss_ts = time_series(giss_year,
                          giss_anom,
                          giss_lounc,
                          giss_hiunc)

    return giss_ts


###########################
##
##  MAIN STRIP
##
###########################

had_ts = read_hadcrut4()


ncdc_ts = read_ncdc()
ncdc_ts.rebaseline(1961,1990)

giss_ts = read_giss()
giss_ts.rebaseline(1961,1990)

had_ts.plot_ts_with_unc('Black','AliceBlue')
giss_ts.plot_ts('Red')
ncdc_ts.plot_ts('DeepSkyBlue')

combined_ts = combine_series(had_ts, ncdc_ts, giss_ts)

plt.savefig('gmt.png', bbox_inches='tight')
plt.show()

combined_ts.plot_ts_with_unc('Black','AliceBlue')
plt.show()

ncdc_ts.plot_ranking_diagram()
giss_ts.plot_ranking_diagram()

combined_ts.plot_ranking_diagram()
combined_ts.plot_decadal(10)

combined_ts.print_ordered_ts(5)

#had_ts.print_ordered_ts(5)
#ncdc_ts.print_ordered_ts(5)
#giss_ts.print_ordered_ts(5)
