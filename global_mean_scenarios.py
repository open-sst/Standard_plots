import numpy as np
#import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from time_series import *
from read_data_sets import *

def combine_monthly_series(had,ncdc,giss):

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


ncdc_version = "v3.5.4.201410"
hadsst_version = "3.1.1.0"
crutem_version = "4.3.0.0"
hadcrut_version = "4.3.0.0"
cowtan_and_way_version = "v2_0_0"

had_monthly = read_hadcrut4_monthly(hadcrut_version)

ncdc_monthly = read_ncdc_monthly(ncdc_version)
ncdc_monthly.rebaseline(1961,1990)

giss_monthly = read_giss_monthly()
giss_monthly.rebaseline(1961,1990)

combinat = combine_monthly_series(had_monthly,ncdc_monthly,giss_monthly)

cw_monthly = read_cowtan_and_way_monthly(cowtan_and_way_version)

ncdc_monthly.plot_scenario("NCDC")
giss_monthly.plot_scenario("GISS")
had_monthly.plot_scenario("HadCRUT4")
cw_monthly.plot_scenario("Cowtan and Way")
combinat.plot_scenario("Combined")

