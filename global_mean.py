import numpy as np
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
    

def combine_series(series1, series2, series3):
    #average together three annual time series and return an annual time series
    #hardcoded to do this from 1880 to 2014
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


def read_hadley(filename):
#read hadley format annual data sets and make annual time series out of them
    f = open(filename, 'r')

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
        hadcrut_upper_unc.append(float(columns[11]))
        hadcrut_lower_unc.append(float(columns[10]))
    
    f.close()
    had_ts = time_series(hadcrut_year,
                         hadcrut_anom,
                         hadcrut_lower_unc,
                         hadcrut_upper_unc)

    return had_ts


def read_hadsst3(version):
    return read_hadley('Data/HadSST.'+version+'_annual_globe_ts.txt')

def read_crutem4(version):
    return read_hadley('Data/CRUTEM.'+version+'.global_n+s')

def read_hadcrut4(version):
    return read_hadley('Data/HadCRUT.'+version+'.annual_ns_avg.txt')


def read_ncdc_format(filename):
    f = open(filename, 'r')
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


def read_ncdc_sst(version):
    return read_ncdc_format('Data/aravg.ann.ocean.90S.90N.'+version+'.asc')

def read_ncdc_lsat(version):
    return read_ncdc_format('Data/aravg.ann.land.90S.90N.'+version+'.asc')

def read_ncdc(version):
    return read_ncdc_format('Data/aravg.ann.land_ocean.90S.90N.'+version+'.asc')



def read_giss_block_monthly(f, block_length, giss_year, giss_month, giss_anom):
#burn first two pointless lines of each block
    line = f.readline()
    line = f.readline()
    for j in range(1,block_length+1):
        line = f.readline()
        columns = line.split()
        for i in range(1,13):
            giss_year.append(float(columns[0]))
            giss_month.append(float(i))
            giss_anom.append(float(columns[i])/100.)
            
    return (giss_year, giss_month, giss_anom)

def read_giss_monthly():
    f = open('Data/GLB.Ts+dSST.txt','r')
    
    giss_year = []
    giss_month = []
    giss_anom = []
  
#read header information and discard
    for i in range(1,7):
        f.readline()

#read first block of 21 years
    giss_year, giss_month, giss_anom = read_giss_block_monthly(f, 21, giss_year, giss_month, giss_anom)
#read five blocks of 20 years   
    for i in range(1,6):
        giss_year, giss_month, giss_anom = read_giss_block_monthly(f, 20, giss_year, giss_month, giss_anom)
#final block has less than 20
    giss_year, giss_month, giss_anom = read_giss_block_monthly(f, 13, giss_year, giss_month, giss_anom)

#final year is incomplete so calculate from monthlies
    g = f.readline()
    columns = g.split()

    for i in range(1,13):
        if columns[i] != '****':
            giss_year.append(float(columns[0]))
            giss_month.append(float(i))
            giss_anom.append(float(columns[i])/100.)
    
    giss_ts = monthly_time_series(giss_year,
                                  giss_month,
                                  giss_anom)

    return giss_ts

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


def read_cowtan_and_way(version):
    f = open('Data/had4_krig_annual_'+version+'.txt','r')
    
    cw_year = []
    cw_anom = []
    cw_lounc = []
    cw_hiunc = []

    for line in f:
        columns = line.split()
        cw_year.append(int(columns[0]))
        cw_anom.append(float(columns[1]))
        cw_lounc.append(float(columns[1]) - 2* float(columns[2]))
        cw_hiunc.append(float(columns[1]) + 2* float(columns[2]))
        
    cw_ts = time_series(cw_year,
                          cw_anom,
                          cw_lounc,
                          cw_hiunc)

    return cw_ts


###########################
###########################
###########################
###########################
##
##      MAIN STRIP
##
###########################
###########################
###########################
###########################


ncdc_version = "v3.5.4.201409"
hadsst_version = "3.1.1.0"
crutem_version = "4.3.0.0"
hadcrut_version = "4.3.0.0"
cowtan_and_way_version = "v2_0_0"


print("GLOBAL AVERAGE TEMPERATURES")

cw_ts = read_cowtan_and_way(cowtan_and_way_version)
latest_year = np.mean([0.594, 0.372, 0.589, 0.683, 0.688, 0.590, 0.505, 0.703, 0.718])
cw_ts.add_year(2014,latest_year,latest_year-0.1,latest_year+0.1)
cw_ts.add_name("Cowtan and Way")

had_ts = read_hadcrut4(hadcrut_version)
had_ts.add_name("HadCRUT."+hadcrut_version)

ncdc_ts = read_ncdc(ncdc_version)
ncdc_ts.rebaseline(1961,1990)
ncdc_ts.add_name("MLOST")

#giss_ts = read_giss()
#giss_ts.rebaseline(1961,1990)
giss_ts_monthly = read_giss_monthly()
giss_ts_monthly.rebaseline(1961,1990)
giss_ts = giss_ts_monthly.annualise()
giss_ts.add_name("GISTEMP")

combined_ts = combine_series(had_ts, ncdc_ts, giss_ts)
combined_ts.add_name("Combined")

combined_ts.plot_ranking_diagram()

combined_ts.plot_skyscraper_diagram()

#make diagram showing all three data sets
had_ts.plot_ts_with_unc('Black','AliceBlue')
giss_ts.plot_ts('Red')
ncdc_ts.plot_ts('DeepSkyBlue')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.legend(loc='best')
plt.xticks(range(1860,2020,20), fontsize = 16)
plt.yticks(np.arange(-0.8,0.8,0.2), fontsize = 16)
#plt.savefig('gmt.png', bbox_inches='tight')
plt.show()


combined_ts.plot_decadal(10)

combined_ts.print_ordered_ts(5)

had_ts.print_ordered_ts(5)
ncdc_ts.print_ordered_ts(5)
giss_ts.print_ordered_ts(5)
cw_ts.print_ordered_ts(5)

print("")
print("Global average SST")

ncdc_sst_ts = read_ncdc_sst(ncdc_version)
ncdc_sst_ts.rebaseline(1961,1990)
ncdc_sst_ts.add_name("ERSSTv3")
ncdc_sst_ts.print_ordered_ts(5)

hadsst_ts = read_hadsst3(hadsst_version)
latest_year = np.mean([0.342, 0.314, 0.347, 0.478, 0.477, 0.563, 0.551, 0.644, 0.574, 0.529])
hadsst_ts.add_year(2014,latest_year,latest_year-0.1,latest_year+0.1)
hadsst_ts.add_name("HadSST."+hadsst_version)
hadsst_ts.print_ordered_ts(5)

hadsst_ts.plot_ts_with_unc('Black','AliceBlue')
ncdc_sst_ts.plot_ts_with_unc('DeepSkyBlue','Yellow')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.legend(loc='best')
plt.show()

print("")
print("Global average LSAT")

ncdc_lsat_ts = read_ncdc_lsat(ncdc_version)
ncdc_lsat_ts.rebaseline(1961,1990)
ncdc_lsat_ts.add_name("NCDC LSAT")
ncdc_lsat_ts.print_ordered_ts(10)

print("")

crutem4_ts = read_crutem4(crutem_version)
crutem4_ts.add_name("CRUTEM."+crutem_version)
crutem4_ts.print_ordered_ts(10)

crutem4_ts.plot_ts_with_unc('Black','AliceBlue')
ncdc_lsat_ts.plot_ts_with_unc('DeepSkyBlue','Yellow')
plt.xlabel('Year', fontsize=18)
plt.ylabel('Anomaly relative to 1961-1990 (K)', fontsize=18)
plt.legend(loc='best')
plt.show()


